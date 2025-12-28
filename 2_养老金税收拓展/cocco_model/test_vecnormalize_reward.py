import gymnasium as gym
from gymnasium.envs.registration import register
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
from stable_baselines3.common.env_util import make_vec_env
from stable_baselines3.common.vec_env import VecNormalize

    # 设置字体路径
plt.rcParams['font.sans-serif'] = ['SimHei']  # 使用简体中文字体
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题

# --- 注册你的自定义环境 ---
# 确保这个脚本能找到 CoccoEnvV2 类
try:
    from v2_cocco_env import CoccoEnvV2
    register(
        id='cocco-v2',
        entry_point='v2_cocco_env:CoccoEnvV2',
    )
    print("环境 'cocco-v2' 注册成功。")
except ImportError:
    print("错误: 无法从 'v2_cocco_env.py' 导入 CoccoEnvV2。")
    print("请确保 test_vecnormalize_reward.py 和 v2_cocco_env.py 在同一个目录下。")
    exit()

# test_vecnormalize_reward.py

def run_reward_test(
    num_timesteps: int = 50000, 
    n_envs: int = 8,
    burn_in_steps: int = 10000 # 新增参数：老化/预热步数
):
    """
    运行测试，收集并可视化原始奖励和归一化奖励的分布。
    新增了 burn_in_steps 参数，以获得更准确的稳定期分布。
    """
    print(f"将在 {num_timesteps + burn_in_steps} 个总时间步上运行测试...")
    print(f"其中，前 {burn_in_steps} 步为 VecNormalize 的'老化'阶段，不计入统计。")

    # --- 1. 创建环境 ---
    raw_env = make_vec_env('cocco-v2', n_envs=n_envs)
    norm_env_to_wrap = make_vec_env('cocco-v2', n_envs=n_envs)
    norm_env = VecNormalize(norm_env_to_wrap, norm_obs=True, norm_reward=True, gamma=0.95)

    # --- 2. 老化/预热阶段 ---
    # 此阶段只运行 norm_env，让其内部统计数据收敛
    norm_env.reset()
    for _ in range(burn_in_steps):
        random_actions = np.array([norm_env.action_space.sample() for _ in range(n_envs)])
        norm_env.step(random_actions)
    print("'老化'阶段完成。VecNormalize 的统计数据已初步建立。")
    
    # --- 3. 收集数据阶段 ---
    raw_rewards_collected = []
    normalized_rewards_collected = []
    
    # 初始化环境以进行公平比较
    raw_env.reset()
    # 注意：norm_env 不需要重置，我们要保留它在老化阶段学到的统计数据
    
    print("开始收集用于分析的奖励数据...")
    for step in range(num_timesteps):
        random_actions = np.array([raw_env.action_space.sample() for _ in range(n_envs)])
        
        _, raw_rewards, _, _= raw_env.step(random_actions)
        raw_rewards_collected.extend(raw_rewards.flatten())
        
        # norm_env 从它在老化阶段结束时的状态继续运行
        _, norm_rewards, _, _= norm_env.step(random_actions)
        normalized_rewards_collected.extend(norm_rewards.flatten())

        if (step + 1) % 10000 == 0:
            print(f"  已收集 {step + 1}/{num_timesteps} 步的数据...")
    
    print("数据收集完成。")
    raw_env.close()
    norm_env.close()

    # --- 4. 数据分析与可视化 (此部分不变) ---
    raw_df = pd.DataFrame(raw_rewards_collected, columns=['reward'])
    norm_df = pd.DataFrame(normalized_rewards_collected, columns=['reward'])

    print("\n--- 原始奖励 (Raw Rewards) 统计 (稳定期) ---")
    print(raw_df.describe())
    
    print("\n--- 归一化奖励 (Normalized Rewards) 统计 (稳定期) ---")
    print(norm_df.describe())

    fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=False)
    fig.suptitle("奖励分布对比 (随机策略, 已进行'老化'处理)", fontsize=16)

    sns.histplot(data=raw_df, x='reward', kde=True, ax=axes[0], bins=50)
    axes[0].set_title('原始奖励 (Raw Reward) 分布')
    axes[0].set_xlabel('奖励值')
    axes[0].set_ylabel('频数')
    axes[0].set_yscale('log')
    axes[0].grid(True)

    sns.histplot(data=norm_df, x='reward', kde=True, ax=axes[1], bins=50)
    axes[1].set_title("VecNormalize 归一化奖励分布 (稳定期)")
    axes[1].set_xlabel('归一化后的奖励值')
    axes[1].set_ylabel('频数')
    axes[1].grid(True)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    output_filename = "reward_distribution_comparison_post_burn_in.png"
    plt.savefig(output_filename)
    print(f"\n图像已保存到: {output_filename}")
    
    plt.show()

if __name__ == "__main__":
    # 运行测试，其中10000步用于老化，50000步用于统计
    run_reward_test(num_timesteps=50000, burn_in_steps=10000)