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
    

    # 自动确定运行编号
    existing_runs = [d for d in os.listdir(log_dir) if d.startswith(base_model_name)]
    run_num = len(existing_runs)

        
    model_dir = os.path.join(log_dir, f"{base_model_name}_{run_num}")
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)
    print(f"模型和日志将保存到: {model_dir}")
    

    # --- 课程学习流程 ---

    # ==========================================================================
    # 阶段一: 预训练 - 环境更简单 (无随机死亡)
    # ==========================================================================        
    # 阶段一环境参数：关闭随机死亡
    env_kwargs_p1 = {'disable_random_death': True}
    
    # 1. 创建原始的 VecEnv 实例
    vec_env_p1_raw = make_vec_env("cocco-v2", n_envs=24, seed=42, env_kwargs=env_kwargs_p1)
    eval_env_p1_raw = make_vec_env("cocco-v2", n_envs=1, seed=123, env_kwargs=env_kwargs_p1)

    # 2. 使用 VecNormalize 封装环境以进行奖励归一化
    #    MODIFIED: norm_obs_index_mask=[1] 表示不归一化观测向量的第2个元素 (年龄)
    #    MODIFIED: norm_reward=False, 与 MATLAB 行为一致，不归一化奖励
    vec_env_p1 = VecNormalize(vec_env_p1_raw, norm_obs=True, norm_obs_index_mask=[1], norm_reward=False)
    eval_env_p1 = VecNormalize(eval_env_p1_raw, norm_obs=True, norm_obs_index_mask=[1], norm_reward=False)
    
    # 3. 预热 VecNormalize
    burn_in_steps = 100
    print(f"\n开始对 VecNormalize 进行 {burn_in_steps} 步的预热...")
    vec_env_p1.reset()
    for _ in range(burn_in_steps):
        random_actions = np.array([vec_env_p1.action_space.sample() for _ in range(vec_env_p1.num_envs)])
        vec_env_p1.step(random_actions)
    print("VecNormalize 预热完成，统计数据已初步建立。")

    # 从环境中获取真实的折现因子
    model_beta = vec_env_p1.get_attr("beta")[0]

    # 3. 设置评估回调函数
    eval_callback_p1 = EvalCallback(
        eval_env_p1,
        best_model_save_path=os.path.join(model_dir, "best_model"),
        log_path=os.path.join(model_dir, "best_model"),
        eval_freq=100, 
        n_eval_episodes=1,
        deterministic=True,
    )

    # 4. 创建并配置 SAC 模型
    model = SAC(
        'MlpPolicy',
        vec_env_p1, 
        verbose=1,
        gamma=model_beta,
        learning_rate=1e-4,
        buffer_size=500_000,
        batch_size=256,
        train_freq=(1, 'step'), 
        gradient_steps=1,
        learning_starts=0,
        # ent_coef=0.1,  # <--- 在这里设置一个固定的较小值
        policy_kwargs=dict(net_arch=[128,128]), # MODIFIED: 与 MATLAB 网络结构更接近
        tensorboard_log=f"./tf_logs/cocco_sac/run_{run_num}_phase1"
    )
    
    # 5. 开始训练
    model.learn(total_timesteps=args.timesteps_phase1, callback=eval_callback_p1, progress_bar=True)
    

    # 6. 保存最终产出
    final_model_path = os.path.join(model_dir, "final_model.zip")
    final_stats_path = os.path.join(model_dir, "vecnormalize_final.pkl")
    model.save(final_model_path)
    vec_env_p1.save(final_stats_path)
    print(f"最终模型已保存到 {final_model_path}")
    print(f"最终 VecNormalize 统计数据已保存到 {final_stats_path}")



def parse_args():
    """
    解析命令行参数。
    """
    parser = argparse.ArgumentParser("Cocco Model Training with Curriculum Learning and Reward Normalization")
    parser.add_argument("--skip-phase1", action="store_true", help="跳过阶段一，直接从已保存的phase1模型开始微调。")
    parser.add_argument("--timesteps-phase1", type=int, default=5_000_000, help="阶段一（预训练）的总时间步数。")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args)