import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import scipy.io
from scipy.interpolate import interp1d
import pandas as pd
import gymnasium as gym # 用于注册环境

# 导入并注册您的Cocco环境
# 这使得RL模型加载时能找到环境的元数据
gym.register(
    id='cocco-v2',
    entry_point='v2_cocco_env:CoccoEnvV2',
)
from v2_cocco_env import CoccoEnvV2
from sbx import SAC # 假设您使用SBX库的SAC模型

def load_vfi_policies(mat_file_path):
    """加载并插值VFI策略 (与之前相同)"""
    print(f"正在从 {mat_file_path} 加载VFI结果...")
    try:
        mat_contents = scipy.io.loadmat(mat_file_path)
    except FileNotFoundError:
        print(f"错误: 找不到MATLAB结果文件 '{mat_file_path}'。")
        return None, None
    vfi_results = mat_contents['vfi_results'][0, 0]
    cS = mat_contents['cS'][0, 0]
    c_policy_grid = vfi_results['C_policy']
    a_policy_grid = vfi_results['A_policy']
    cash_grid = cS['gcash'].flatten()
    tb = int(cS['tb'].flatten()[0])
    tn = int(cS['tn'].flatten()[0])
    policies = []
    for t_idx in range(tn):
        c_interpolator = interp1d(
            cash_grid, c_policy_grid[:, t_idx], kind='linear',
            bounds_error=False, fill_value=(c_policy_grid[0, t_idx], c_policy_grid[-1, t_idx])
        )
        a_interpolator = interp1d(
            cash_grid, a_policy_grid[:, t_idx], kind='linear',
            bounds_error=False, fill_value=(a_policy_grid[0, t_idx], a_policy_grid[-1, t_idx])
        )
        policies.append({'c_policy': c_interpolator, 'a_policy': a_interpolator})
    print("VFI策略加载并插值完成。")
    return policies

def run_simulation(policy, env_params, nsim=10000):
    """
    统一模拟器，在绝对单位 (primitive units) 下运行，以匹配 v2_cocco_env。
    
    Args:
        policy: 策略对象。可以是包含插值器的VFI策略list，也可以是RL模型。
        env_params: 从CoccoEnvV2实例中提取的环境参数字典。
        nsim: 模拟的个体数量。
    """
    print(f"开始模拟策略: {policy.__class__.__name__ if not isinstance(policy, list) else 'VFI Policy'}")
    
    # --- 1. 初始化 ---
    tn = env_params['tn']
    tb = env_params['tb']
    tr = env_params['tr']
    td = env_params['td']
    work_periods = tr - tb

    # 预分配模拟数组 (tn, nsim) - 全部使用绝对单位
    sim_W_abs = np.zeros((tn, nsim))
    sim_C_abs = np.zeros((tn, nsim))
    sim_Y_abs = np.zeros((tn, nsim))
    sim_A = np.zeros((tn, nsim))
    sim_P_abs = np.ones((tn, nsim)) # 这是随机部分 P_t = exp(v_t)
    
    # --- 2. 预生成所有外生冲击 (与之前相同) ---
    print("  正在生成外生冲击路径 (使用Antithetic Variates)...")
    nsim_half = nsim // 2
    nsim_half_plus = nsim_half + (nsim % 2)
    z_shocks_half = np.random.randn(work_periods, nsim_half_plus) * env_params['smav']
    u_shocks_log_half = np.random.randn(work_periods, nsim_half_plus) * env_params['smay']
    r_shocks_half = np.random.randn(tn, nsim_half_plus) * env_params['sigr']
    z_shocks = np.hstack([z_shocks_half, -z_shocks_half])[:, :nsim]
    u_shocks_log = np.hstack([u_shocks_log_half, -u_shocks_log_half])[:, :nsim]
    r_shocks = np.hstack([r_shocks_half, -r_shocks_half])[:, :nsim]
    R = env_params['rf'] + env_params['mu'] + r_shocks
    
    # --- 3. Vectorized模拟主循环 (按时间t推进) ---
    print("  正在模拟生命周期决策 (绝对单位)...")
    ages = np.arange(tb, tb + tn)
    normalized_ages = 2 * (ages - tb) / (td - tb) - 1

    for t in tqdm(range(tn)):
        current_age = ages[t]
        
        # a. 计算当期收入 Y_t
        if current_age < tr:
            u_shock_t = np.exp(u_shocks_log[t, :]) if t < work_periods else 1.0
            det_Y_t = env_params['f_y'][t]
            Y_gross_t = det_Y_t * sim_P_abs[t, :] * u_shock_t
            sim_Y_abs[t, :] = Y_gross_t * (1 - env_params['tau_y'])
        else:
            P_retirement = sim_P_abs[work_periods, :]
            sim_Y_abs[t, :] = env_params['ret_fac'] * env_params['f_y'][tr-tb] * P_retirement
        
        X_t_abs = sim_W_abs[t, :] + sim_Y_abs[t, :]
        
        # c. 查询策略 (VFI和RL的分支点)
        if isinstance(policy, list): # --- VFI 策略 ---
            # *** 核心修正: 进行正确的单位转换 (方法2, 'hat' normalization) ***
            
            # 1. 计算总永久收入 P_total = f(t) * P_t
            # 注意: 退休后没有确定性增长 f(t)，但为避免除零，P_total 可视为 P_retirement
            if current_age < tr:
                P_total_t = env_params['f_y'][t] * sim_P_abs[t, :]
            else:
                # 退休期没有f(t), VFI策略的输入是绝对财富除以退休时的P
                P_total_t = env_params['f_y'][tr-tb] * sim_P_abs[t, :] 
            
            # 2. 将绝对现金持有量转换为VFI策略所需的'hat'单位
            # 避免除以过小的值
            P_total_t[P_total_t < 1e-10] = 1e-10
            X_t_hat = X_t_abs / P_total_t
            
            c_policy_interp = policy[t]['c_policy']
            a_policy_interp = policy[t]['a_policy']
            
            # 3. 查询策略得到 'hat'单位的消费
            c_t_hat = c_policy_interp(X_t_hat)
            
            # 4. 将 'hat'单位的消费转换回绝对消费
            sim_C_abs[t, :] = c_t_hat * P_total_t
            sim_A[t, :] = a_policy_interp(X_t_hat)

        else: # --- RL 策略 (逻辑不变) ---
            obs = np.stack([
                sim_W_abs[t, :], 
                sim_P_abs[t, :], 
                np.full(nsim, normalized_ages[t])
            ], axis=1).astype(np.float32)
            action, _ = policy.predict(obs, deterministic=True)
            c_prop = action[:, 0]
            sim_A[t, :] = action[:, 1]
            sim_C_abs[t, :] = c_prop * X_t_abs

        # d. 应用约束
        sim_A[t, :] = np.clip(sim_A[t, :], 0.0, 1.0)
        sim_C_abs[t, :] = np.clip(sim_C_abs[t, :], 1e-6, X_t_abs * 0.99999)

        # e. 计算储蓄并更新下一期状态
        if t < tn - 1:
            S_t_abs = X_t_abs - sim_C_abs[t, :]
            portfolio_return_t = sim_A[t, :] * R[t, :] + (1 - sim_A[t, :]) * env_params['rf']
            sim_W_abs[t+1, :] = S_t_abs * portfolio_return_t
            
            next_age = ages[t+1]
            if next_age < tr:
                z_shock_t_plus_1 = np.exp(z_shocks[t, :]) if t < work_periods -1 else 1.0
                sim_P_abs[t+1, :] = sim_P_abs[t, :] * z_shock_t_plus_1
            else:
                sim_P_abs[t+1, :] = sim_P_abs[t, :]

    # --- 4. 结果汇总与计算 (逻辑不变) ---
    print("  汇总结果并计算指标...")
    mean_profiles = pd.DataFrame({
        'age': ages,
        'wealth': np.mean(sim_W_abs, axis=1),
        'consumption': np.mean(sim_C_abs, axis=1),
        'income': np.mean(sim_Y_abs, axis=1),
        'alpha': np.mean(sim_A, axis=1),
    })
    
    gamma = env_params['gamma']
    beta = env_params['beta']
    if gamma == 1: utility = np.log(sim_C_abs)
    else: utility = (sim_C_abs**(1 - gamma)) / (1 - gamma)
    
    discounted_utility = np.sum(utility * (beta ** np.arange(tn))[:, np.newaxis], axis=0)
    mean_reward = np.mean(discounted_utility)
    std_reward = np.std(discounted_utility)
    
    print("\n" + "="*50)
    print(f"策略评估结果: {policy.__class__.__name__ if not isinstance(policy, list) else 'VFI Policy'}")
    print("="*50)
    print(f"  - 平均折现终身效用: {mean_reward:.4f}")
    print(f"  - 效用标准差:         {std_reward:.4f}")
    print("="*50)

    return mean_profiles

def plot_comparison(vfi_df, rl_df, env_params):
    """绘制VFI和RL策略的对比图。"""
    fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
    
    retirement_age = env_params['tr']
    ages = vfi_df['age']
    
    # 子图1: 财富、消费、收入
    axes[0].plot(ages, vfi_df['consumption'], 'r-', label='Consumption (VFI)')
    axes[0].plot(ages, rl_df['consumption'], 'r--', label='Consumption (RL)')
    axes[0].plot(ages, vfi_df['wealth'], 'b-', label='Wealth (VFI)')
    axes[0].plot(ages, rl_df['wealth'], 'b--', label='Wealth (RL)')
    axes[0].plot(ages, vfi_df['income'], 'k:', label='Income (Mean)')
    axes[0].axvline(x=retirement_age, color='grey', linestyle='--', label='Retirement')
    axes[0].set_title('Mean Lifecycle Profiles (Absolute Values)')
    axes[0].set_ylabel('Value')
    axes[0].legend()
    axes[0].grid(True)
    
    # 子图2: 风险资产配置
    axes[1].plot(ages, vfi_df['alpha'], 'm-', label='Risky Share (VFI)')
    axes[1].plot(ages, rl_df['alpha'], 'm--', label='Risky Share (RL)')
    axes[1].axvline(x=retirement_age, color='grey', linestyle='--', label='Retirement')
    axes[1].set_title('Mean Portfolio Choice')
    axes[1].set_xlabel('Age')
    axes[1].set_ylabel('Share')
    axes[1].set_ylim(0, 1.05)
    axes[1].legend()
    axes[1].grid(True)
    
    fig.suptitle("VFI vs. RL Policy Comparison", fontsize=16)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

if __name__ == "__main__":
    MAT_FILE_PATH = 'result/vfi_results_crra_nopps_solo.mat'
    # ！！！重要：请将下面的路径修改为您训练好的最终RL模型！！！
    RL_MODEL_PATH = './models/cocco_sac_curriculum/cocco_sac_run_102/best_model_phase1/best_model.zip'
    
    # 1. 创建一个环境实例以提取参数
    # 使用有随机死亡的真实环境设置
    env = CoccoEnvV2(disable_random_death=False)
    env_params = {
        'tb': env.tb, 'tr': env.tr, 'td': env.td, 'tn': env.tn,
        'beta': env.beta, 'gamma': env.gamma,
        'smay': env.smay, 'smav': env.smav,
        'rf': env.rf, 'mu': env.mu, 'sigr': env.sigr,
        'ret_fac': env.ret_fac, 'f_y': env.f_y,
        'tau_y':env.tau_y
    }

    # 2. 加载VFI策略
    vfi_policies = load_vfi_policies(MAT_FILE_PATH)
    if not vfi_policies:
        exit()

    # 3. 加载RL模型
    try:
        rl_model = SAC.load(RL_MODEL_PATH)
    except FileNotFoundError:
        print(f"错误: 找不到RL模型文件 '{RL_MODEL_PATH}'。请检查路径。")
        exit()

    # 4. 使用统一的模拟器运行评估
    vfi_profiles_df = run_simulation(vfi_policies, env_params, nsim=20000)
    rl_profiles_df = run_simulation(rl_model, env_params, nsim=20000)
    
    # 5. 绘图对比
    plot_comparison(vfi_profiles_df, rl_profiles_df, env_params)