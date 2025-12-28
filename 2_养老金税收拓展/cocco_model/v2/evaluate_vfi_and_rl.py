import os
import argparse
import gymnasium as gym
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.interpolate import interp1d

from sbx import SAC
from utils.vec_normalize import VecNormalize
from stable_baselines3.common.vec_env import DummyVecEnv

# 注册环境，确保 gym.make 可以找到它
from gymnasium.envs.registration import register
try:
    register(id='cocco-v2', entry_point='v2_cocco_env:CoccoEnvV2')
except gym.error.Error:
    print("Environment 'cocco-v2' is already registered. Skipping.")


def run_simulation(raw_env, n_sim, policy_type, rl_model=None, norm_env=None, vfi_policies=None):
    """
    在原始环境中运行模拟的核心函数。
    """
    # MODIFIED: Access attributes via .unwrapped
    tn = raw_env.unwrapped.tn 
    all_trajectories = []

    for i in range(n_sim):
        # NOTE: raw_env.reset() returns obs from the wrapper, which is fine
        obs_raw, info = raw_env.reset()
        
        trajectory = {
            'age': [], 'wealth': [], 'consumption': [], 'income': [], 
            'c_prop': [], 'alpha': []
        }
        
        is_done = False
        t = 0
        while not is_done:
            # MODIFIED: Access attributes via .unwrapped
            age_t = raw_env.unwrapped.age
            wealth_t = raw_env.unwrapped.W
            income_t = raw_env.unwrapped.Y
            
            if policy_type == 'RL':
                obs_norm = norm_env.normalize_obs(obs_raw)
                action, _ = rl_model.predict(obs_norm, deterministic=True)
            
            elif policy_type == 'VFI':
                # MODIFIED: Access attributes via .unwrapped
                current_x = info['normalized_coh']
                current_t_idx = raw_env.unwrapped.age - raw_env.unwrapped.tb
                
                C_vfi = vfi_policies['c'][current_t_idx](current_x)
                alpha_vfi = vfi_policies['a'][current_t_idx](current_x)
                
                c_prop_vfi = C_vfi / current_x if current_x > 1e-6 else 0.0
                action = np.array([c_prop_vfi, alpha_vfi])

            # NOTE: raw_env.step() is also called on the wrapper
            obs_raw, _, is_done, _, info = raw_env.step(action)
            
            trajectory['age'].append(age_t)
            trajectory['wealth'].append(wealth_t)
            trajectory['income'].append(income_t)
            trajectory['consumption'].append(info['absolute_consumption'])
            trajectory['c_prop'].append(info['c_prop'])
            trajectory['alpha'].append(info['risky_share'])
            t += 1

        df = pd.DataFrame(trajectory)
        all_trajectories.append(df)

    combined_df = pd.concat(all_trajectories)
    mean_profiles = combined_df.groupby('age').mean().reset_index()
    
    return mean_profiles


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
        c_interp = interp1d(x_grid, C_policy[:, t], kind='cubic', bounds_error=False, fill_value="extrapolate")
        a_interp = interp1d(x_grid, A_policy[:, t], kind='cubic', bounds_error=False, fill_value="extrapolate")
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

    # --- 3. 运行模拟 ---
    raw_env = gym.make("cocco-v2", disable_random_death=True) # 使用无死亡环境进行比较

    print("\n  [RL] 正在进行模拟...")
    raw_env.reset(seed=123)
    np.random.seed(123)
    rl_profiles = run_simulation(raw_env, args.n_sim, 'RL', rl_model=rl_model, norm_env=norm_env)
    # MODIFIED: Access attributes via .unwrapped
    gamma = raw_env.unwrapped.gamma
    beta = raw_env.unwrapped.beta
    rl_utility = np.sum( (rl_profiles['consumption']**(1-gamma))/(1-gamma) * (beta**np.arange(len(rl_profiles))) )

    print("  [VFI] 正在进行模拟...")
    raw_env.reset(seed=123)
    np.random.seed(123)
    vfi_profiles = run_simulation(raw_env, args.n_sim, 'VFI', vfi_policies=vfi_policies)
    # MODIFIED: Use variables defined above
    vfi_utility = np.sum( (vfi_profiles['consumption']**(1-gamma))/(1-gamma) * (beta**np.arange(len(vfi_profiles))) )

    print(f"  模拟完成。 RL Utility: {rl_utility:.4f} | VFI Utility: {vfi_utility:.4f}")

    # --- 4. 绘图 ---
    print("\n  正在生成详细的比较图表...")
    if not os.path.exists(args.save_path):
        os.makedirs(args.save_path)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    # MODIFIED: Access attributes via .unwrapped
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
        self.rl_model = "models/cocco_sac_curriculum/cocco_sac_run_3/best_model/best_model.zip"
        self.stats = "models/cocco_sac_curriculum/cocco_sac_run_3/best_model/vecnormalize_best.pkl"
        self.vfi_file = "result/vfi_benchmark_results.mat"
        self.save_path = "results_final_comparison"
        self.n_sim = 1000

if __name__ == "__main__":
    # MODIFIED: 直接实例化 Args 类
    args = Args()
    evaluate_and_compare(args)