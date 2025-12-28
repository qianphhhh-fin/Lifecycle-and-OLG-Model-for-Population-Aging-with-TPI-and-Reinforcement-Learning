import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.interpolate import interp1d
import gymnasium as gym

# 假设 v2_cocco_env.py 在同一个目录下或在 Python 路径中
from v2_cocco_env import CoccoEnvV2

def analyze_and_plot_results(rewards, wealths, perm_incomes, norm_wealths, lifecycle_data):
    """
    对模拟结果进行统计分析和绘图。
    """
    print("\n" + "="*30)
    print("      SIMULATION RESULTS ANALYSIS")
    print("="*30)

    # --- 1. 统计量分析 ---
    def print_stats(data, name):
        data = np.array(data)
        print(f"\n--- Distribution of: {name} ---")
        print(f"  Mean:      {np.mean(data):.4f}")
        print(f"  Std Dev:   {np.std(data):.4f}")
        print(f"  Min:       {np.min(data):.4f}")
        print(f"  Median:    {np.median(data):.4f}")
        print(f"  Max:       {np.max(data):.4f}")
        print(f"  1% Pctile: {np.percentile(data, 1):.4f}")
        print(f"  99% Pctile:{np.percentile(data, 99):.4f}")
        print("-"*(len(name) + 23))

    print_stats(rewards, "Raw Reward (Utility)")
    print_stats(wealths, "Absolute Wealth (W_t)")
    print_stats(perm_incomes, "Permanent Income Level (P_t)")
    print_stats(norm_wealths, "Normalized Wealth (w_t = W_t/P_t)")

    # --- 2. 绘图 ---
    print("\nGenerating plots...")
    fig = plt.figure(figsize=(16, 10))

    # a. 奖励分布
    ax1 = plt.subplot(2, 2, 1)
    ax1.hist(rewards, bins=100, density=True, label='Reward Distribution')
    ax1.axvline(np.mean(rewards), color='r', linestyle='--', label=f'Mean: {np.mean(rewards):.2f}')
    ax1.axvline(np.median(rewards), color='g', linestyle='-', label=f'Median: {np.median(rewards):.2f}')
    ax1.set_title('Distribution of Raw Step Rewards (Utility)')
    ax1.set_xlabel('Reward Value')
    ax1.set_ylabel('Probability Density')
    ax1.set_yscale('log')
    ax1.legend()
    ax1.grid(True)

    # b. 生命周期策略
    mean_profiles = pd.concat(lifecycle_data).groupby('age').mean()
    retirement_age = 26 - 1 # tr - 1

    ax2 = plt.subplot(2, 2, 2)
    ax2.plot(mean_profiles.index, mean_profiles['c_prop'], 'r-', lw=2, label='c_prop (VFI)')
    ax2.set_title('Mean Lifecycle Consumption Policy (c_prop)')
    ax2.set_xlabel('Age'); ax2.set_ylabel('Proportion of Cash-on-Hand')
    ax2.axvline(retirement_age, color='grey', linestyle='--', label='Retirement')
    ax2.grid(True); ax2.legend(); ax2.set_ylim(0, 1.05)
    
    ax3 = plt.subplot(2, 2, 3)
    ax3.plot(mean_profiles.index, mean_profiles['alpha'], 'm-', lw=2, label='alpha (VFI)')
    ax3.set_title('Mean Lifecycle Investment Policy (alpha)')
    ax3.set_xlabel('Age'); ax3.set_ylabel('Share in Risky Asset')
    ax3.axvline(retirement_age, color='grey', linestyle='--', label='Retirement')
    ax3.grid(True); ax3.legend(); ax3.set_ylim(-0.05, 1.05)
    
    # c. 绝对财富路径
    ax4 = plt.subplot(2, 2, 4)
    ax4.plot(mean_profiles.index, mean_profiles['wealth'], 'b-', lw=2, label='Wealth (VFI)')
    ax4.set_title('Mean Lifecycle Absolute Wealth')
    ax4.set_xlabel('Age'); ax4.set_ylabel('Absolute Value')
    ax4.axvline(retirement_age, color='grey', linestyle='--', label='Retirement')
    ax4.grid(True); ax4.legend()


    plt.tight_layout()
    plt.show()


def main():
    """
    主执行函数。
    """
    # --- 1. 参数设置 ---
    VFI_MAT_FILE = 'result/vfi_benchmark_results.mat'
    N_EPISODES = 1000
    DISABLE_RANDOM_DEATH = False # 在有随机死亡的完整环境中测试

    if not os.path.exists(VFI_MAT_FILE):
        raise FileNotFoundError(f"VFI结果文件未找到: {VFI_MAT_FILE}。请先运行 main_noPPS.m 生成它。")

    # --- 2. 加载和准备 VFI 策略 ---
    print(f"Loading VFI policies from: {VFI_MAT_FILE}")
    mat_data = loadmat(VFI_MAT_FILE)
    vfi_results = mat_data['vfi_results'][0, 0]
    cS = mat_data['cS'][0, 0]
    # 策略是在归一化网格 cS['gcash'] 上定义的
    C_policy_norm = vfi_results['C_policy']
    A_policy = vfi_results['A_policy']
    norm_cash_grid = cS['gcash'].flatten()
    
    vfi_policies = {'c_norm': [], 'a': []}
    tn = C_policy_norm.shape[1]
    for t in range(tn):
        c_interp = interp1d(norm_cash_grid, C_policy_norm[:, t], kind='cubic', bounds_error=False, fill_value="extrapolate")
        a_interp = interp1d(norm_cash_grid, A_policy[:, t], kind='cubic', bounds_error=False, fill_value="extrapolate")
        vfi_policies['c_norm'].append(c_interp)
        vfi_policies['a'].append(a_interp)
    print("VFI policies loaded and interpolated.")

    # --- 3. 运行模拟 ---
    env = CoccoEnvV2(disable_random_death=DISABLE_RANDOM_DEATH)
    print(f"\nStarting simulation of {N_EPISODES} lifecycles using VFI policy...")

    # 初始化数据收集器
    all_rewards = []
    all_abs_wealths = []
    all_perm_incomes = []
    all_norm_wealths = []
    all_lifecycle_data = []

    for i_episode in range(N_EPISODES):
        if (i_episode + 1) % 20 == 0:
            print(f"  Simulating episode {i_episode + 1}/{N_EPISODES}")
            
        obs, info = env.reset()
        is_done = False
        
        trajectory = {'age': [], 'c_prop': [], 'alpha': [], 'wealth': []}

        while not is_done:
            # a. 获取当前状态 (绝对单位)
            age_t =  env.unwrapped.age
            t_idx = age_t - env.unwrapped.tb

            trajectory['wealth'].append(info['absolute_wealth'])

            # b. 将状态转换为 VFI 策略所需的归一化单位
            # VFI 策略的状态是归一化现金持有量 x_t = X_t / P_t
            norm_cash_on_hand_t = info['normalized_coh']

            # c. 从插值函数获取最优决策 (归一化消费 c_t 和 alpha_t)
            norm_c_t = vfi_policies['c_norm'][t_idx](norm_cash_on_hand_t)
            alpha_t = vfi_policies['a'][t_idx](norm_cash_on_hand_t)

            # d. 将 VFI 决策转换回环境 step 函数所需的动作 (消费倾向 c_prop)
            c_prop_t = norm_c_t/norm_cash_on_hand_t
            
            
            action = np.array([c_prop_t, alpha_t])

            # e. 在环境中执行动作
            obs, reward, is_done, _, info = env.step(action)

            c_prop_t = info['c_prop']
            alpha_t = info['risky_share']
            
            # f. 收集数据
            all_rewards.append(reward)
            all_abs_wealths.append(info['absolute_wealth'])
            all_perm_incomes.append(info['raw_permant_shock'])
            all_norm_wealths.append(info['normalized_coh'])
            
            trajectory['age'].append(age_t)
            trajectory['c_prop'].append(c_prop_t)
            trajectory['alpha'].append(alpha_t)
            

        all_lifecycle_data.append(pd.DataFrame(trajectory))
        
    print("Simulation finished.")
    
    # --- 4. 分析与可视化 ---
    analyze_and_plot_results(
        all_rewards, all_abs_wealths, all_perm_incomes, all_norm_wealths, all_lifecycle_data
    )

if __name__ == "__main__":
    main()