import gymnasium as gym
import numpy as np
from gymnasium.envs.registration import register
import os
import argparse
from sbx import SAC
from utils.env_util import make_vec_env
from utils.AllocationEvalCallback import EvalCallback
# from stable_baselines3.common.vec_env import VecNormalize
from utils.vec_normalize import VecNormalize

# 注册自定义环境
# 确保在代码开始时执行，以便 make_vec_env 可以找到它
# register(
#     id='cocco-v1',
#     entry_point='v1_cocco_env:CoccoEnv',
# )
register(
    id='cocco-v2',
    entry_point='v2_cocco_env:CoccoEnvV2',
)

def run(args):
    """
    主运行函数，执行包含奖励归一化的课程学习流程。
    """
    # --- 路径和模型名称设置 ---
    log_dir = './models/cocco_sac_curriculum/'
    base_model_name = 'cocco_sac_run'

    # 自动创建新的运行文件夹，方便实验管理
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    
    if args.run_num is None:
        # 自动确定运行编号
        existing_runs = [d for d in os.listdir(log_dir) if d.startswith(base_model_name)]
        run_num = len(existing_runs)
    else:
        run_num = args.run_num
        
    model_dir = os.path.join(log_dir, f"{base_model_name}_{run_num}")
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)
    print(f"模型和日志将保存到: {model_dir}")
    
    # 定义阶段一产出的文件路径
    phase1_model_path = os.path.join(model_dir, "phase1_model.zip")
    # VecNormalize 的统计数据也必须保存，以确保阶段二奖励尺度的一致性
    vec_normalize_stats_path = os.path.join(model_dir, "vecnormalize_phase1.pkl")

    # --- 课程学习流程 ---

    # ==========================================================================
    # 阶段一: 预训练 - 环境更简单 (无随机死亡)
    # ==========================================================================
    if not args.skip_phase1:
        print("\n===== 开始阶段一: 预训练 (无随机死亡) =====")
        
        # 阶段一环境参数：关闭随机死亡
        env_kwargs_p1 = {'disable_random_death': True}
        
        # 1. 创建原始的 VecEnv 实例
        vec_env_p1_raw = make_vec_env("cocco-v2", n_envs=24, seed=42, env_kwargs=env_kwargs_p1)
        eval_env_p1_raw = make_vec_env("cocco-v2", n_envs=1, seed=123, env_kwargs=env_kwargs_p1)

        # 2. 使用 VecNormalize 封装环境以进行奖励归一化
        #    这是稳定训练的关键步骤
        vec_env_p1 = VecNormalize(vec_env_p1_raw, norm_obs=True, norm_obs_index_mask=[2,3],norm_reward=False)
        #    评估环境也必须封装，以确保奖励尺度匹配
        #    设置 training=False 防止它在评估中更新自己的统计数据
        eval_env_p1 = VecNormalize(eval_env_p1_raw,  norm_obs=True, norm_obs_index_mask=[2,3], norm_reward=False)

            # --- 核心修改 2: 增加 VecNormalize 的预热阶段 ---
        # 3. 预热 VecNormalize
        burn_in_steps = 2000
        print(f"\n开始对 VecNormalize 进行 {burn_in_steps} 步的预热...")
        vec_env_p1.reset()
        for _ in range(burn_in_steps):
            # 使用随机动作进行探索
            random_actions = np.array([vec_env_p1.action_space.sample() for _ in range(vec_env_p1.num_envs)])
            vec_env_p1.step(random_actions)
        print("VecNormalize 预热完成，统计数据已初步建立。")

        # 从环境中获取真实的折现因子
        model_beta = vec_env_p1.get_attr("beta")[0]

        # 定义VFI文件路径
        VFI_MAT_FILE = 'result/vfi_results_crra_nopps_solo.mat'


        # 3. 设置评估回调函数
        #    传入封装后的评估环境 eval_env_p1
        eval_callback_p1 = EvalCallback(
            eval_env_p1,
            best_model_save_path=os.path.join(model_dir, "best_model_phase1"),
            log_path=os.path.join(model_dir, "logs_phase1"),
            eval_freq=1000, # 评估频率可以根据需要调整
            n_eval_episodes=300,
            deterministic=True,
        )

        # 4. 创建并配置 SAC 模型
        model = SAC(
            'MlpPolicy',
            vec_env_p1, # Agent 在封装后的环境中训练
            verbose=1,
            gamma=model_beta,
            learning_rate=5e-5,
            buffer_size=500_000,
            # ent_coef=0.01, # <--- 核心修改：从 'auto' 改为固定的较小值
            batch_size=256,
            # ent_coef='auto',
            train_freq=(1, 'step'), # 每收集1步经验，就触发一次更新流程
            gradient_steps=1,       # <--- *** 在该流程中，只执行1次梯度更新 ***
            learning_starts=0,
            policy_kwargs=dict(net_arch=[32,32]),
            tensorboard_log=f"./tf_logs/cocco_sac/run_{run_num}_phase1"
        )
        
        # 5. 开始训练
        model.learn(total_timesteps=args.timesteps_phase1, callback=eval_callback_p1, progress_bar=True)
        
        # 6. 保存模型和 VecNormalize 统计数据
        model.save(phase1_model_path)
        vec_env_p1.save(vec_normalize_stats_path)
        print(f"阶段一模型已保存到: {phase1_model_path}")
        print(f"阶段一 VecNormalize 统计数据已保存到: {vec_normalize_stats_path}")

    else:
        print("\n===== 跳过阶段一，直接进入阶段二 =====")
        if not os.path.exists(phase1_model_path) or not os.path.exists(vec_normalize_stats_path):
            print(f"错误: 找不到阶段一模型或 VecNormalize 统计数据。请先完整运行阶段一。")
            return


    # ==========================================================================
    # 阶段二: 微调 - 环境更真实 (有随机死亡)
    # ==========================================================================
    print("\n===== 开始阶段二: 微调 (有随机死亡) =====")
    aaa
    # 阶段二环境参数：开启随机死亡
    env_kwargs_p2 = {'disable_random_death': False}
    
    # 1. 创建阶段二的原始 VecEnv
    vec_env_p2_raw = make_vec_env("cocco-v2", n_envs=64, seed=42, env_kwargs=env_kwargs_p2)
    eval_env_p2_raw = make_vec_env("cocco-v2", n_envs=1, seed=123, env_kwargs=env_kwargs_p2)

    # 2. **关键**: 加载阶段一的统计数据来初始化 VecNormalize
    #    这确保了奖励尺度在两个阶段之间是连续的，对稳定微调至关重要
    vec_env_p2 = VecNormalize.load(vec_normalize_stats_path, vec_env_p2_raw)
    eval_env_p2 = VecNormalize.load(vec_normalize_stats_path, eval_env_p2_raw)
    # 再次确保评估环境不更新统计数据
    eval_env_p2.training = False

    model_beta = vec_env_p2.get_attr("beta")[0]
    
    # 3. 为阶段二创建新的回调函数
    eval_callback_p2 = EvalCallback(
        eval_env_p2,
        best_model_save_path=os.path.join(model_dir, "best_model_final"),
        log_path=os.path.join(model_dir, "logs_final"),
        eval_freq=15000,
        n_eval_episodes=500,
        deterministic=True,
        render=False,
        discount_factor=model_beta
    )

    # 4. 加载阶段一的模型，并设置为使用新的封装环境
    model = SAC.load(phase1_model_path, env=vec_env_p2)
    # 可以选择为微调设置一个更小的学习率
    # model.learning_rate = 1e-5 
    
    print("已加载阶段一模型和 VecNormalize 统计数据，开始微调...")
    
    # 5. 开始微调
    model.learn(
        total_timesteps=args.timesteps_phase2,
        callback=eval_callback_p2,
        progress_bar=True,
        reset_num_timesteps=False # 关键：不重置时间步计数器，使日志连续
    )
    
    # 6. 保存最终产出
    final_model_path = os.path.join(model_dir, "final_model.zip")
    final_stats_path = os.path.join(model_dir, "vecnormalize_final.pkl")
    model.save(final_model_path)
    vec_env_p2.save(final_stats_path)
    print(f"最终模型已保存到 {final_model_path}")
    print(f"最终 VecNormalize 统计数据已保存到 {final_stats_path}")


def parse_args():
    """
    解析命令行参数。
    """
    parser = argparse.ArgumentParser("Cocco Model Training with Curriculum Learning and Reward Normalization")
    parser.add_argument("--skip-phase1", action="store_true", help="跳过阶段一，直接从已保存的phase1模型开始微调。")
    parser.add_argument("--timesteps-phase1", type=int, default=5_000_000, help="阶段一（预训练）的总时间步数。")
    parser.add_argument("--timesteps-phase2", type=int, default=500_000, help="阶段二（微调）的总时间步数。")
    
    parser.add_argument("--run-num", type=int, default=None, help="手动指定运行编号，用于继续之前的训练或覆盖。")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args)