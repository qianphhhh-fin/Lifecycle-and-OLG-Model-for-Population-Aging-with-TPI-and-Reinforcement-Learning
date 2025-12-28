import gymnasium as gym
import numpy as np
from gymnasium.envs.registration import register
import os
import argparse
from sbx import PPO  # <--- 核心修改: 导入 PPO
from utils.env_util import make_vec_env
from utils.AllocationEvalCallback_aug import EvalCallback
from utils.vec_normalize import VecNormalize

# --- [V3 核心修改] ---
# 注册新的增强环境 'cocco-v3-aug'
try:
    register(
        id='cocco-v3-aug',
        entry_point='v3_cocco_env_aug:CoccoEnvV3_Augmented',
    )
except gym.error.Error:
    print("Environment 'cocco-v3-aug' is already registered. Skipping.")


def run(args):
    """
    主运行函数，使用 PPO 算法训练基于 Bisi et al. (2022) 状态增强思想的环境。
    PPO 作为一种 On-Policy 算法，在理论上更适合处理这种稀疏奖励和数据非平稳性问题。
    """
    # --- 路径和模型名称设置 ---
    log_dir = './models/cocco_ppo_augmented/' # <--- 修改路径以区分实验
    base_model_name = 'cocco_ppo_aug_run'

    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    
    existing_runs = [d for d in os.listdir(log_dir) if d.startswith(base_model_name)]
    run_num = len(existing_runs)
        
    model_dir = os.path.join(log_dir, f"{base_model_name}_{run_num}")
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)
    print(f"模型和日志将保存到: {model_dir}")

    # --- 环境设置 (与SAC版本一致) ---
    env_kwargs = {'disable_random_death': True}
    vec_env_raw = make_vec_env("cocco-v3-aug", n_envs=24, seed=42, env_kwargs=env_kwargs)
    eval_env_raw = make_vec_env("cocco-v3-aug", n_envs=100, seed=123, env_kwargs=env_kwargs)

    # --- 归一化设置 (与SAC版本一致) ---
    # 观测空间: [X_norm, P, normalized_age, v, w]
    # 屏蔽索引 2 (age) 和 4 (w) 的归一化
    vec_env = VecNormalize(vec_env_raw, norm_obs=True, norm_obs_index_mask=[2, 4], norm_reward=True)
    eval_env = VecNormalize(eval_env_raw, norm_obs=True, norm_obs_index_mask=[2, 4], norm_reward=False)
    
    # 预热 VecNormalize
    burn_in_steps = 2000
    print(f"\n开始对 VecNormalize 进行 {burn_in_steps} 步的预热...")
    vec_env.reset()
    for _ in range(burn_in_steps):
        random_actions = np.array([vec_env.action_space.sample() for _ in range(vec_env.num_envs)])
        vec_env.step(random_actions)
    print("VecNormalize 预热完成。")

    # --- 评估回调函数设置 (与SAC版本一致) ---
    # 使用我们适配过的 Callback，可以正确评估 v3 环境并报告可比的终身效用
    eval_callback = EvalCallback(
        eval_env,
        best_model_save_path=os.path.join(model_dir, "best_model"),
        log_path=os.path.join(model_dir, "best_model"),
        eval_freq=1000,
        n_eval_episodes=300,
        deterministic=True,
    )

    # --- [核心修改] PPO 模型配置 ---
    # 根据 Bisi 论文的理论，增强后的 MDP 的折现因子应为 1。
    model_gamma = 1.0

    model = PPO(
        'MlpPolicy',
        vec_env,
        verbose=1,
        gamma=model_gamma,          # <--- 关键理论参数: 必须为 1.0
        learning_rate=1e-4,         # PPO 的一个标准学习率
        n_steps=2048,               # 每次更新前，每个环境收集的步数
        batch_size=64,              # mini-batch 大小
        n_epochs=10,                # 每次更新时，在收集到的数据上迭代的次数
        gae_lambda=0.95,            # GAE 参数，用于信用分配
        clip_range=0.2,             # PPO 裁剪参数
        ent_coef=0.01,              # 熵系数，鼓励探索，对稀疏奖励问题有益
        vf_coef=0.5,                # 值函数损失的系数
        policy_kwargs={
            "net_arch": dict(pi=[256, 256], vf=[256, 256]), # PPO策略和价值网络
        },
        tensorboard_log=f"./tf_logs/cocco_ppo_augmented/run_{run_num}" # <--- 修改Tensorboard路径
    )
    
    # --- 开始训练 ---
    print("\n--- 开始在增强环境 (Augmented Environment) 上训练 PPO 模型 ---")
    print(f"关键参数: gamma={model_gamma}, 奖励结构=稀疏, 算法=PPO (On-Policy)")
    model.learn(total_timesteps=args.timesteps, callback=eval_callback, progress_bar=True)
    
    # --- 保存最终产出 ---
    final_model_path = os.path.join(model_dir, "final_model.zip")
    final_stats_path = os.path.join(model_dir, "vecnormalize_final.pkl")
    model.save(final_model_path)
    vec_env.save(final_stats_path)
    print(f"\n训练完成。")
    print(f"最终模型已保存到 {final_model_path}")
    print(f"最终 VecNormalize 统计数据已保存到 {final_stats_path}")

def parse_args():
    parser = argparse.ArgumentParser("Cocco Model Training on Augmented State Environment with PPO")
    parser.add_argument("--timesteps", type=int, default=5_000_000, help="训练总时间步数。")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    run(args)