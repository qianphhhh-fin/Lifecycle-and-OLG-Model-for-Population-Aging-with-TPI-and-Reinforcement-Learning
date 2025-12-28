import gymnasium as gym
import numpy as np
from gymnasium.envs.registration import register
import os
import argparse
from sbx import SAC
# from stable_baselines3 import SAC
from utils.env_util import make_vec_env
from utils.AllocationEvalCallback_aug import EvalCallback
from utils.vec_normalize import VecNormalize
import optax

# --- [V3 核心修改] ---
# 1. 导入新的增强环境
# 2. 注册一个新的环境 ID 'cocco-v3-aug'
try:
    register(
        id='cocco-v3-aug',
        entry_point='v3_cocco_env_aug:CoccoEnvV3_Augmented',
    )
except gym.error.Error:
    print("Environment 'cocco-v3-aug' is already registered. Skipping.")


def run(args):
    """
    主运行函数，用于训练基于 Bisi et al. (2022) 状态增强思想的环境。
    """
    # --- 路径和模型名称设置 ---
    log_dir = './models/cocco_sac_augmented/'
    base_model_name = 'cocco_sac_aug_run'

    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    
    existing_runs = [d for d in os.listdir(log_dir) if d.startswith(base_model_name)]
    run_num = len(existing_runs)
        
    model_dir = os.path.join(log_dir, f"{base_model_name}_{run_num}")
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)
    print(f"模型和日志将保存到: {model_dir}")

    # --- 环境设置 ---
    # 为了简化和加速初始训练，我们可以在一个没有随机死亡的"简单"版本上进行
    env_kwargs = {'disable_random_death': True}
    
    # --- [V3 核心修改] ---
    # 使用新注册的环境 'cocco-v3-aug'
    vec_env_raw = make_vec_env("cocco-v3-aug", n_envs=24, seed=42, env_kwargs=env_kwargs)
    eval_env_raw = make_vec_env("cocco-v3-aug", n_envs=100, seed=123, env_kwargs=env_kwargs)

    # --- 归一化设置 ---
    # **[V3 核心修改]** 调整归一化掩码以适应新的5维观测空间
    # 观测空间: [X_norm, P, normalized_age, v, w]
    # - normalized_age (索引 2) 已经归一化到 [-1, 1]，无需再次归一化。
    # - w (索引 4) 是 beta^t，范围在 [0, 1]，也无需归一化。
    # - v (索引 3) 是累计效用，其尺度会变化，需要归一化。
    # 因此，我们屏蔽索引 2 和 4。
    # 同时，由于奖励是稀疏的（只在最后一步出现），归一化奖励至关重要。
    vec_env = VecNormalize(vec_env_raw, norm_obs=True, norm_obs_index_mask=[2, 4], norm_reward=True)
    # 评估环境的奖励不进行归一化，以便在回调中计算真实的终身效用
    eval_env = VecNormalize(eval_env_raw, norm_obs=True, norm_obs_index_mask=[2, 4], norm_reward=False)
    
    # 预热 VecNormalize
    burn_in_steps = 2000
    print(f"\n开始对 VecNormalize 进行 {burn_in_steps} 步的预热...")
    vec_env.reset()
    for _ in range(burn_in_steps):
        random_actions = np.array([vec_env.action_space.sample() for _ in range(vec_env.num_envs)])
        vec_env.step(random_actions)
    print("VecNormalize 预热完成。")

    # --- 评估回调函数设置 ---
    eval_callback = EvalCallback(
        eval_env,
        best_model_save_path=os.path.join(model_dir, "best_model"),
        log_path=os.path.join(model_dir, "best_model"),
        eval_freq=1000, # 稀疏奖励环境可能需要更频繁的评估来观察进展
        n_eval_episodes=300,
        deterministic=True,
    )

    # --- SAC 模型配置 ---
    # **[V3 核心修改]** 根据 Bisi 论文的理论，增强后的 MDP 的折现因子应为 1。
    # 原始的折现因子 beta 已经被内化到了环境的状态转移和累计效用 v 的计算中。
    # 这是至关重要的一步！
    model_gamma = 1.0

    # 学习率调度器
    lr_schedule = optax.constant_schedule(1e-4)

    model = SAC(
        'MlpPolicy',
        vec_env, 
        verbose=1,
        gamma=model_gamma,  # <--- 关键修改！
        learning_rate=lr_schedule,
        buffer_size=1_000_000,
        batch_size=256,
        train_freq=(1, 'step'), 
        gradient_steps=1,
        learning_starts=10000,
        # 稀疏奖励下，熵正则化对于鼓励探索更为重要，可能需要微调
        # target_entropy='auto', 
        policy_kwargs={
            "net_arch": [256, 256],
        },
        tensorboard_log=f"./tf_logs/cocco_sac_augmented/run_{run_num}"
    )
    
    # --- 开始训练 ---
    # 稀疏奖励环境通常需要更长的训练时间才能看到明显效果
    print("\n--- 开始在增强环境 (Augmented Environment) 上训练 SAC 模型 ---")
    print(f"关键参数: gamma={model_gamma}, 奖励结构=稀疏")
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
    parser = argparse.ArgumentParser("Cocco Model Training on Augmented State Environment")
    parser.add_argument("--timesteps", type=int, default=5_000_000, help="训练总时间步数。")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    run(args)