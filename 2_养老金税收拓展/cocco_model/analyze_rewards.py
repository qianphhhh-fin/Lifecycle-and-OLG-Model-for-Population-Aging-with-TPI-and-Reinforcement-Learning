import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm
import warnings

# 忽略来自环境的无关警告
warnings.filterwarnings("ignore", category=UserWarning, module='gymnasium')

# 导入您的Cocco环境
try:
    from v2_cocco_env import CoccoEnvV2
except ImportError:
    print("错误: 无法找到 v2_cocco_env.py。请确保此脚本与您的环境文件在同一个目录下。")
    exit()

def analyze_distribution(env, policy_fn, num_steps=50000):
    """
    在环境中运行指定策略，并收集奖励值和状态变量。
    
    :param env: Gymnasium环境实例。
    :param policy_fn: 一个函数，输入环境实例的观测值，输出一个动作。
    :param num_steps: 模拟的总步数。
    :return: 一个包含奖励和状态变量的Pandas DataFrame。
    """
    data_collected = {
        'reward': [],
        'cash_on_hand': [],
        'permanent_income': [],
        'age': []
    }
    
    obs, info = env.reset()
    
    print(f"正在使用 '{policy_fn.__name__}' 策略进行模拟...")
    for _ in tqdm(range(num_steps)):
        action = policy_fn(obs) # 策略函数现在接收观测值
        obs, reward, terminated, truncated, info = env.step(action)
        
        data_collected['reward'].append(reward)
        # 从观测值中解包状态变量
        cash_on_hand, permanent_income, normalized_age = obs
        # 将年龄从[-1, 1]的归一化范围转换回真实年龄
        age = round(((normalized_age + 1) / 2) * (env.td - env.tb) + env.tb)

        data_collected['cash_on_hand'].append(cash_on_hand)
        data_collected['permanent_income'].append(permanent_income)
        data_collected['age'].append(age)
        
        if terminated or truncated:
            obs, info = env.reset()
            
    return pd.DataFrame(data_collected)

def random_policy(obs):
    """一个完全随机的策略。现在它不依赖env实例。"""
    # 假设动作空间是Box(2,)
    return np.random.rand(2).astype(np.float32)

def heuristic_policy(obs):
    """一个固定的启发式策略，模拟一个理性的、但非最优的代理。"""
    # 总是消费可用资源的20%，并将剩余储蓄的40%投入风险资产
    return np.array([0.2, 0.4], dtype=np.float32)

def print_report_and_judge(df, title):
    """
    计算统计数据，打印报告，并判断奖励和状态分布是否合理。
    """
    print("\n" + "="*60)
    print(f"分布分析报告: {title}")
    print("="*60)
    
    if df.empty:
        print("未能收集到任何数据。")
        return

    # --- 1. 奖励分析 ---
    rewards = df['reward'].values
    print("\n--- 奖励 (Reward) 分布 ---")
    mean_rew, std_rew, min_rew, max_rew = np.mean(rewards), np.std(rewards), np.min(rewards), np.max(rewards)
    p1, p5, p50, p95, p99 = np.percentile(rewards, [1, 5, 50, 95, 99])

    print(f"  - 均值: {mean_rew:.4e}, 标准差: {std_rew:.4e}")
    print(f"  - 范围: [{min_rew:.4e}, {max_rew:.4e}]")
    print(f"  - 分位数 (1%, 5%, 50%, 95%, 99%):")
    print(f"    [{p1:.3f}, {p5:.3f}, {p50:.3f}, {p95:.3f}, {p99:.3f}]")

    # 自动化判断
    good_range = [-20.0, 5.0] # 稍微放宽范围，因为CRRA效用非对称
    is_reasonable = (p1 >= good_range[0])
    
    print("\n  自动化判断 (奖励):")
    if is_reasonable:
        print("    [结论] \033[92m合理 (PASSED)\033[0m")
        print("    [分析] 奖励的下限在可控范围内，降低了梯度爆炸的风险。")
    else:
        print("    [结论] \033[91m可能不合理 (FAILED)\033[0m")
        print(f"    [分析] 存在极端负奖励 (低至 {p1:.2e})，可能导致Critic Loss爆炸或梯度消失。")
        print("    [建议] 这是CRRA效用的固有特性，是正常的。关键是确保模型能够学习避免这些区域。")
        print("           如果模型无法收敛，可以考虑临时增加一个奖励下限进行调试。")

    # --- 2. 状态变量分析 ---
    print("\n--- 状态变量 (State Variables) 分布 ---")
    for col in ['cash_on_hand', 'permanent_income', 'age']:
        var_data = df[col].values
        mean_var, std_var, min_var, max_var = np.mean(var_data), np.std(var_data), np.min(var_data), np.max(var_data)
        p5, p50, p95 = np.percentile(var_data, [5, 50, 95])
        
        print(f"\n  * {col.replace('_', ' ').title()}:")
        print(f"    - 均值: {mean_var:.3f}, 标准差: {std_var:.3f}")
        print(f"    - 范围: [{min_var:.3f}, {max_var:.3f}]")
        print(f"    - 分位数 (5%, 50%, 95%): [{p5:.3f}, {p50:.3f}, {p95:.3f}]")

    print("\n" + "="*60)
    
    return df

def plot_distributions(df, title):
    """绘制奖励和状态变量分布的直方图。"""
    cols_to_plot = ['reward', 'cash_on_hand', 'permanent_income']
    fig, axes = plt.subplots(len(cols_to_plot), 1, figsize=(12, 6 * len(cols_to_plot)))
    fig.suptitle(f'Distribution Analysis: {title}', fontsize=16)

    for i, col in enumerate(cols_to_plot):
        data = df[col].values
        ax = axes[i]
        ax.hist(data, bins=100, alpha=0.75, color='royalblue', density=True)
        ax.set_title(f'{col.replace("_", " ").title()} Distribution')
        ax.set_xlabel('Value')
        ax.set_ylabel('Density')
        ax.grid(True)
        
        # 绘制统计线
        p5, p50, p95 = np.percentile(data, [5, 50, 95])
        ax.axvline(p50, color='red', linestyle='--', label=f'Median: {p50:.2f}')
        ax.axvline(p5, color='orange', linestyle=':', label=f'5th Percentile: {p5:.2f}')
        ax.axvline(p95, color='orange', linestyle=':', label=f'95th Percentile: {p95:.2f}')
        ax.legend()
        
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()


if __name__ == "__main__":
    # **重要**: 此脚本会使用您在 v2_cocco_env.py 中定义的最新奖励函数和事件顺序
    env = CoccoEnvV2(disable_random_death=False)
    
    # 1. 使用随机策略进行分析
    random_df = analyze_distribution(env, random_policy)
    print_report_and_judge(random_df, "随机策略 (模拟训练初期)")
    
    # 2. 使用启发式策略进行分析
    heuristic_df = analyze_distribution(env, heuristic_policy)
    print_report_and_judge(heuristic_df, "启发式策略 (模拟训练后期)")
    
    # 3. 绘制直方图
    plot_distributions(random_df, "随机策略")
    plot_distributions(heuristic_df, "启发式策略")