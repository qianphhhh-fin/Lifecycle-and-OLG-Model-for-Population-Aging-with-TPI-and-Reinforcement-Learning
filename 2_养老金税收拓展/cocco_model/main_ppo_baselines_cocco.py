import gymnasium as gym
import numpy as np
from gymnasium.envs.registration import register
import os
import argparse
from sbx import PPO
from utils.env_util import make_vec_env
from utils.AllocationEvalCallback import EvalCallback
from utils.vec_normalize import VecNormalize

# 注册自定义环境
# 确保在代码开始时执行，以便 make_vec_env 可以找到它
register(
    id='cocco-v2',
    entry_point='v2_cocco_env:CoccoEnvV2',
)

def run(args):
    """
    主运行函数，执行包含奖励归一化的课程学习流程。
    使用 PPO 算法，并配置其以模拟 REINFORCE 算法。
    """
    # --- 路径和模型名称设置 ---
    log_dir = './models/cocco_ppo_reinforce/'
    base_model_name = 'cocco_ppo_run'

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
        vec_env_p1 = VecNormalize(vec_env_p1_raw, norm_obs=True, norm_obs_index_mask=[2], norm_reward=False)
        eval_env_p1 = VecNormalize(eval_env_p1_raw,  norm_obs=True, norm_obs_index_mask=[2], norm_reward=False)

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

        # 4. 设置评估回调函数
        eval_callback_p1 = EvalCallback(
            eval_env_p1,
            best_model_save_path=os.path.join(model_dir, "best_model_phase1"),
            log_path=os.path.join(model_dir, "logs_phase1"),
            eval_freq=5000,
            n_eval_episodes=300,
            deterministic=True,
        )

        # 5. 创建并配置一个“类REINFORCE”的PPO模型
        model = PPO(
            'MlpPolicy',
            vec_env_p1,
            verbose=1,
            gamma=model_beta,
            learning_rate=1e-4, # 策略梯度方法通常需要较小的学习率
            
            # --- 核心修改：让PPO表现得像REINFORCE ---
            n_steps=1024, # 确保 n_steps 远大于一个生命周期的长度
            batch_size=128,      # 使用一个较大的批量
            n_epochs=10,         # 在收集到的数据上多迭代几次
            gae_lambda=1.0,      # *** gae_lambda=1.0 意味着不使用自举，完全依赖蒙特卡洛回报G ***
            clip_range=0.2,     
            ent_coef=0.001,      # 一个很小的熵，用于防止策略过早崩溃
            vf_coef=0.5,        
            max_grad_norm=0.5,
            
            policy_kwargs=dict(net_arch=dict(pi=[64, 64], vf=[64, 64])),
            tensorboard_log=f"./tf_logs/cocco_ppo_reinforce/run_{run_num}_phase1"
        )
        
        # 6. 开始训练
        model.learn(total_timesteps=args.timesteps_phase1, callback=eval_callback_p1, progress_bar=True)
        
        # 7. 保存模型和 VecNormalize 统计数据
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
    vec_env_p2 = VecNormalize.load(vec_normalize_stats_path, vec_env_p2_raw)
    eval_env_p2 = VecNormalize.load(vec_normalize_stats_path, eval_env_p2_raw)
    # 再次确保评估环境不更新统计数据
    eval_env_p2.training = False
    
    # 3. 为阶段二创建新的回调函数
    eval_callback_p2 = EvalCallback(
        eval_env_p2,
        best_model_save_path=os.path.join(model_dir, "best_model_final"),
        log_path=os.path.join(model_dir, "logs_final"),
        eval_freq=15000,
        n_eval_episodes=500,
        deterministic=True,
    )

    # 4. 加载阶段一的模型，并设置为使用新的封装环境
    model = PPO.load(phase1_model_path, env=vec_env_p2)
    
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
    parser = argparse.ArgumentParser("Cocco Model Training with PPO (REINFORCE-like) and Curriculum Learning")
    parser.add_argument("--skip-phase1", action="store_true", help="跳过阶段一，直接从已保存的phase1模型开始微调。")
    parser.add_argument("--timesteps-phase1", type=int, default=5_000_000, help="阶段一（预训练）的总时间步数。")
    parser.add_argument("--timesteps-phase2", type=int, default=500_000, help="阶段二（微调）的总时间步数。")
    parser.add_argument("--run-num", type=int, default=None, help="手动指定运行编号，用于继续之前的训练或覆盖。")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args)