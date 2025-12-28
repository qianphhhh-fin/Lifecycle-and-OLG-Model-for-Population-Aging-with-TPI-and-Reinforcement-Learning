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
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

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

    # --- 3. 运行模拟 ---
    raw_env = gym.make("cocco-v2", disable_random_death=True)
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
                    action = np.array([C_prop_vfi, alpha_vfi])

                obs_raw, utility_t, is_done, _, info = raw_env.step(action)
                
                C_t = info['absolute_consumption']
                period_utilities.append(utility_t)

                trajectory_data['age'].append(age_t)
                trajectory_data['wealth'].append(info['absolute_wealth'])
                trajectory_data['income'].append(info['absolute_income'])
                trajectory_data['consumption'].append(C_t)
                trajectory_data['c_prop'].append(info['c_prop'])
                trajectory_data['alpha'].append(info['risky_share'])
            
            discounts = beta ** np.arange(len(period_utilities))
            lifetime_utility = np.sum(discounts * np.array(period_utilities))
            all_lifetime_utilities.append(lifetime_utility)
            all_trajectories.append(pd.DataFrame(trajectory_data))

        mean_profiles = pd.concat(all_trajectories).groupby('age').mean().reset_index()
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

    # --- 4. 绘图 (高质量、符合出版要求) ---
    print("\n  正在生成符合出版要求的比较图表...")
    if not os.path.exists(args.save_path):
        os.makedirs(args.save_path)

    # --- Matplotlib 全局设置 (用于学术图表) ---
    # 为确保中文字体正确显示，如果系统没有 "Times New Roman"，可能需要配置
    plt.rcParams['font.sans-serif'] = ['SimHei'] # 例如，使用黑体
    plt.rcParams['axes.unicode_minus'] = False  # 解决负号'-'显示为方块的问题
    plt.rcParams.update({
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 9,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "figure.dpi": 600,
    })

    # --- 尺寸转换 (cm to inch) ---
    fig_width_cm = 14.99
    fig_height_cm = 4.5
    fig_width_in = fig_width_cm / 2.54
    fig_height_in = fig_height_cm / 2.54

    # --- 创建图表 ---
    fig, axes = plt.subplots(1, 2, figsize=(fig_width_in, fig_height_in), constrained_layout=True)
    retirement_age = raw_env.unwrapped.tr - 1
    age_grid = rl_profiles['age']

    # --- 子图 (a): 风险资产配置比例路径 ---
    ax1 = axes[0]
    ax1.plot(vfi_profiles['age'], vfi_profiles['alpha'], color='black', linestyle='-', linewidth=1.2, label='VFI')
    ax1.plot(rl_profiles['age'], rl_profiles['alpha'], color='#0072B2', linestyle='--', linewidth=1.2, label='RL') # 使用蓝色
    ax1.axvline(x=retirement_age, color='grey', linestyle=':', linewidth=1)
    ax1.set_title('(a) 风险资产配置比例路径')
    ax1.set_xlabel('年龄')
    ax1.set_ylabel('风险资产份额')
    ax1.legend(loc='upper right', frameon=False)
    ax1.grid(True, which='major', linestyle='--', linewidth='0.5', alpha=0.7)
    ax1.set_xlim(age_grid.min(), age_grid.max())
    ax1.set_ylim(-0.05, 1.05)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # --- 子图 (b): 消费和财富路径 ---
    ax2 = axes[1]
    # 财富路径 (使用一种颜色，不同线型)
    ax2.plot(vfi_profiles['age'], vfi_profiles['wealth'], color='#D55E00', linestyle='-', linewidth=1.2, label='财富 (VFI)') # 橙色
    ax2.plot(rl_profiles['age'], rl_profiles['wealth'], color='#D55E00', linestyle='--', linewidth=1.2, label='财富 (RL)') # 橙色
    # 消费路径 (使用另一种颜色，不同线型)
    ax2.plot(vfi_profiles['age'], vfi_profiles['consumption'], color='#0072B2', linestyle='-', linewidth=1.2, label='消费 (VFI)') # 蓝色
    ax2.plot(rl_profiles['age'], rl_profiles['consumption'], color='#0072B2', linestyle='--', linewidth=1.2, label='消费 (RL)') # 蓝色

    ax2.axvline(x=retirement_age, color='grey', linestyle=':', linewidth=1)
    ax2.set_title('(b) 消费和财富路径')
    ax2.set_xlabel('年龄')
    ax2.set_ylabel('绝对值')
    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(handles, labels, loc='lower center', frameon=False, ncol=2, handlelength=2.0, columnspacing=1.0)
    ax2.grid(True, which='major', linestyle='--', linewidth='0.5', alpha=0.7)
    ax2.set_xlim(age_grid.min(), age_grid.max())
    ax2.set_ylim(bottom=0)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    # --- 保存图表 ---
    output_filename = os.path.join(args.save_path, 'rl_vs_vfi_lifecyle_comparison_paper.png')
    plt.savefig(output_filename, dpi=600, bbox_inches='tight')
    plt.close(fig) 
    print(f"  论文质量图表已保存到: {output_filename}")

# MODIFIED: 替换 argparse 的辅助类
class Args:
    def __init__(self):
        self.rl_model = "models/cocco_sac_curriculum/cocco_sac_run_11/best_model/best_model.zip"
        self.stats = "models/cocco_sac_curriculum/cocco_sac_run_11/best_model/vecnormalize_best.pkl"
        self.vfi_file = "result/vfi_benchmark_results.mat"
        self.save_path = "fig"
        self.n_sim = 1000

if __name__ == "__main__":
    # MODIFIED: 直接实例化 Args 类
    args = Args()
    evaluate_and_compare(args)