import os
import argparse
import gymnasium as gym
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.interpolate import interp1d

from sbx import SAC
# from stable_baselines3 import SAC
from utils.vec_normalize import VecNormalize
from stable_baselines3.common.vec_env import DummyVecEnv

# 注册环境，确保 gym.make 可以找到它
from gymnasium.envs.registration import register
try:
    register(id='cocco-v2', entry_point='v2_cocco_env:CoccoEnvV2')
except gym.error.Error:
    print("Environment 'cocco-v2' is already registered. Skipping.")



def evaluate_and_compare(args):
    """
    主评估与比较函数。
    """
    print("开始执行 RL vs. VFI 评估流程...")
    
    # --- 1. 加载 VFI 基准策略 ---
    print(f"  正在加载 VFI 基准策略从: {args.vfi_file}")
    mat_data = loadmat(args.vfi_file)
    vfi_results = mat_data['vfi_results'][0, 0]
    cS = mat_data['cS'][0, 0]
    C_policy = vfi_results['C_policy']
    A_policy = vfi_results['A_policy']
    x_grid = cS['gcash'].flatten()
    
    vfi_policies = {'c': [], 'a': []}
    tn = C_policy.shape[1]

    for t in range(tn):
        c_interp = interp1d(x_grid, C_policy[:, t], kind='cubic', bounds_error=True)
        a_interp = interp1d(x_grid, A_policy[:, t], kind='cubic', bounds_error=True)
        vfi_policies['c'].append(c_interp)
        vfi_policies['a'].append(a_interp)
    print("  VFI 策略加载并插值完成。")

    # --- 2. 加载 RL 智能体和环境 ---
    print(f"  正在加载 RL 智能体从: {args.rl_model}")
    rl_model = SAC.load(args.rl_model)
    
    print(f"  正在加载归一化统计数据从: {args.stats}")
    dummy_env = DummyVecEnv([lambda: gym.make("cocco-v2")])
    norm_env = VecNormalize.load(args.stats, dummy_env)
    norm_env.training = False
    print("  RL 智能体和统计数据加载完成。")

    # --- 3. 运行模拟 [MODIFIED FOR CORRECT UTILITY CALCULATION] ---
    raw_env = gym.make("cocco-v2", disable_random_death=True)
    gamma = raw_env.unwrapped.gamma
    beta = raw_env.unwrapped.beta

    def get_policy_performance(policy_type):
        all_trajectories = []
        all_lifetime_utilities = []
        
        for i in range(args.n_sim):
            obs_raw, info = raw_env.reset()
            
            trajectory_data = {'age': [], 'wealth': [], 'consumption': [], 'income': [], 'c_prop': [], 'alpha': []}
            period_utilities = []
            is_done = False
            
            while not is_done:
                age_t = raw_env.unwrapped.age
                
                if policy_type == 'RL':
                    obs_norm = norm_env.normalize_obs(obs_raw)
                    action, _ = rl_model.predict(obs_norm, deterministic=True)
                else: # VFI
                    current_x = info['normalized_coh']
                    current_t_idx = raw_env.unwrapped.age - raw_env.unwrapped.tb
                    C_prop_vfi = vfi_policies['c'][current_t_idx](current_x)
                    alpha_vfi = vfi_policies['a'][current_t_idx](current_x)
                    action = np.array([C_prop_vfi/current_x, alpha_vfi])

                obs_raw, _, is_done, _, info = raw_env.step(action)
                
                C_t = info['absolute_consumption']
                utility_t = (C_t**(1-gamma)) / (1-gamma)
                period_utilities.append(utility_t)

                trajectory_data['age'].append(age_t)
                trajectory_data['wealth'].append(info['absolute_wealth'])
                trajectory_data['income'].append(info['absolute_income'])
                trajectory_data['consumption'].append(C_t)
                trajectory_data['c_prop'].append(info['c_prop'])
                trajectory_data['alpha'].append(info['risky_share'])
            
            # 计算该条路径的终身效用
            discounts = beta ** np.arange(len(period_utilities))
            lifetime_utility = np.sum(discounts * np.array(period_utilities))
            all_lifetime_utilities.append(lifetime_utility)
            all_trajectories.append(pd.DataFrame(trajectory_data))

        # 计算生命周期剖面
        mean_profiles = pd.concat(all_trajectories).groupby('age').mean().reset_index()
        # 计算平均终身效用
        mean_utility = np.mean(all_lifetime_utilities)
        
        return mean_profiles, mean_utility

    print("\n  [RL] 正在进行模拟...")
    raw_env.reset(seed=123)
    np.random.seed(123)
    rl_profiles, rl_utility = get_policy_performance('RL')

    print("  [VFI] 正在进行模拟...")
    raw_env.reset(seed=123)
    np.random.seed(123)
    vfi_profiles, vfi_utility = get_policy_performance('VFI')

    print(f"  模拟完成。 RL Utility: {rl_utility:.4f} | VFI Utility: {vfi_utility:.4f}")


    # --- 4. 绘图 ---
    print("\n  正在生成详细的比较图表...")
    if not os.path.exists(args.save_path):
        os.makedirs(args.save_path)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    retirement_age = raw_env.unwrapped.tr - 1

    ax1 = axes[0, 0]
    ax1.plot(rl_profiles['age'], rl_profiles['c_prop'], 'r--', linewidth=2, label='c-prop (RL)')
    ax1.plot(vfi_profiles['age'], vfi_profiles['c_prop'], 'r-', linewidth=1.5, label='c-prop (VFI)')
    ax1.axvline(x=retirement_age, color='grey', linestyle='--', label='Retirement')
    ax1.set_title('c-prop')
    ax1.set_xlabel('Age'); ax1.set_ylabel('Proportion of Cash-on-Hand')
    ax1.legend(loc='best'); ax1.grid(True)
    ax1.set_ylim(-0.05, 1.05)

    ax2 = axes[0, 1]
    ax2.plot(rl_profiles['age'], rl_profiles['alpha'], 'm--', linewidth=2, label='Risky Share (RL)')
    ax2.plot(vfi_profiles['age'], vfi_profiles['alpha'], 'm-', linewidth=1.5, label='Risky Share (VFI)')
    ax2.axvline(x=retirement_age, color='grey', linestyle='--', label='Retirement')
    ax2.set_title('Investment Policy (alpha)')
    ax2.set_xlabel('Age'); ax2.set_ylabel('Share in Risky Asset')
    ax2.legend(loc='best'); ax2.grid(True)
    ax2.set_ylim(-0.05, 1.05)
    
    ax3 = axes[1, 0]
    ax3.plot(rl_profiles['age'], rl_profiles['wealth'], 'b--', linewidth=2, label='Wealth (RL)')
    ax3.plot(vfi_profiles['age'], vfi_profiles['wealth'], 'b-', linewidth=1.5, label='Wealth (VFI)')
    ax3.set_title('Mean Lifecycle Profiles (Absolute Units)')
    ax3.set_xlabel('Age'); ax3.set_ylabel('Value')
    ax3.legend(loc='upper left'); ax3.grid(True)

    axes[1, 1].axis('off')

    base_title = 'RL vs. VFI Solution Comparison (Identical Shocks)'
    utility_title = f'RL Utility: {rl_utility:.4f}  |  VFI Utility: {vfi_utility:.4f}'
    fig.suptitle(f"{base_title}\n{utility_title}", fontsize=14)
    fig.tight_layout(rect=[0, 0.03, 1, 0.94])

    output_filename = os.path.join(args.save_path, 'rl_vs_vfi_comparison_detailed.png')
    plt.savefig(output_filename)
    print(f"  评估图表已保存到: {output_filename}")

# MODIFIED: 替换 argparse 的辅助类
class Args:
    def __init__(self):
        self.rl_model = "models/cocco_sac_curriculum/cocco_sac_run_6/best_model/best_model.zip"
        self.stats = "models/cocco_sac_curriculum/cocco_sac_run_6/best_model/vecnormalize_best.pkl"
        self.vfi_file = "result/vfi_benchmark_results.mat"
        self.save_path = "results_final_comparison"
        self.n_sim = 1000

if __name__ == "__main__":
    # MODIFIED: 直接实例化 Args 类
    args = Args()
    evaluate_and_compare(args)