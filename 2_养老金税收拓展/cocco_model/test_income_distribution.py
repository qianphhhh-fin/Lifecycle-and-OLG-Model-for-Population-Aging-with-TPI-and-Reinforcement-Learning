import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

# 确保 v1_cocco_env.py 在同一个目录下
try:
    from v1_cocco_env import CoccoEnv
except ImportError:
    print("错误: 无法找到 v1_cocco_env.py。请确保此脚本与您的环境文件在同一个目录下。")
    exit()

def simulate_income_process(num_episodes=10000, fixed_c_prop=0.75):
    """
    运行多轮模拟，专注于收集绝对收入 Y 的数据。

    Args:
        num_episodes (int): 要模拟的生命周期总数。
        fixed_c_prop (float): 一个固定的、用于驱动环境的消费比例。
                              这只是为了让模拟进行下去，对收入Y的产生过程没有影响。
    """
    print(f"开始模拟 {num_episodes} 个生命周期以收集收入数据...")
    
    # 初始化环境
    env = CoccoEnv()
    
    # 用于存储所有收入数据的列表
    all_incomes = []
    
    # 用于时序图分析的数据
    income_by_age = [[] for _ in range(env.tn)]

    for _ in tqdm(range(num_episodes)):
        env.reset()
        done = False
        
        while not done:
            # 使用一个固定的、简单的动作来推进环境
            # 这个动作本身不重要，我们只关心环境产生的 info['absolute_income']
            action = np.array([fixed_c_prop, 0.5]) # 消费75%，投资50%
            
            _obs, _reward, done, _truncated, info = env.step(action)
            
            # 提取并存储绝对收入
            income = info.get('absolute_income', 0)
            
            # 忽略初始重置时的0收入
            if income > 0:
                all_incomes.append(income)
                age_idx = info['age'] - env.tb
                if 0 <= age_idx < env.tn:
                    income_by_age[age_idx].append(income)

    print("模拟完成。")
    return np.array(all_incomes), income_by_age, env.tb

def analyze_and_plot(incomes, income_by_age, start_age):
    """
    分析收入数据并绘制图表。
    """
    if len(incomes) == 0:
        print("警告：没有收集到任何收入数据。")
        return

    print("\n" + "="*50)
    print("绝对收入 (Y) 分布统计分析")
    print("="*50)
    
    # 计算关键的百分位数
    min_income = np.min(incomes)
    p01 = np.percentile(incomes, 0.01)
    p1 = np.percentile(incomes, 1)
    p5 = np.percentile(incomes, 5)
    p10 = np.percentile(incomes, 10)
    median_income = np.median(incomes)
    mean_income = np.mean(incomes)
    
    print(f"  - 观测总数: {len(incomes)}")
    print(f"  - 均值:       {mean_income:.4f}")
    print(f"  - 中位数:     {median_income:.4f}")
    print("-" * 50)
    print("  尾部风险分析 (低收入端):")
    print(f"  - 最小值:     {min_income:.6f}")
    print(f"  - 0.01百分位: {p01:.6f} (万分之一的概率收入低于此值)")
    print(f"  - 1百分位:    {p1:.6f} (百分之一的概率收入低于此值)")
    print(f"  - 5百分位:    {p5:.6f}")
    print(f"  - 10百分位:   {p10:.6f}")
    print("="*50)
    
    print("\n建议: `min_consumption` 应设置在 0.01 百分位附近或更低，")
    print(f"       例如 0.01，以覆盖绝大多数极端情况。")

    # --- 绘图 ---
    fig = plt.figure(figsize=(15, 6))

    # 1. 收入分布直方图 (对数刻度)
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.hist(incomes, bins=100, color='skyblue', edgecolor='black', alpha=0.7)
    ax1.set_yscale('log')
    ax1.set_title('绝对收入 (Y) 的分布 (Y轴对数刻度)')
    ax1.set_xlabel('绝对收入')
    ax1.set_ylabel('频数 (对数)')
    ax1.grid(True, which="both", ls="--")
    # 添加关键百分位数的垂直线
    ax1.axvline(p01, color='red', linestyle='--', label=f'0.01%ile ({p01:.4f})')
    ax1.axvline(p1, color='orange', linestyle='--', label=f'1%ile ({p1:.4f})')
    ax1.legend()

    # 2. 按年龄划分的收入百分位数
    ax2 = fig.add_subplot(1, 2, 2)
    ages = np.arange(start_age, start_age + len(income_by_age))
    
    # 计算每个年龄的百分位数
    p1_by_age = [np.percentile(age_data, 1) if age_data else np.nan for age_data in income_by_age]
    median_by_age = [np.median(age_data) if age_data else np.nan for age_data in income_by_age]
    
    ax2.plot(ages, median_by_age, 'b-', label='中位数收入')
    ax2.plot(ages, p1_by_age, 'r--', label='1百分位收入 (最差1%情况)')
    ax2.set_title('收入随年龄的变化')
    ax2.set_xlabel('年龄')
    ax2.set_ylabel('绝对收入')
    ax2.legend()
    ax2.grid(True, ls="--")
    ax2.set_yscale('log')
    
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # 运行模拟
    # 增加模拟次数可以得到更精确的尾部分布估计
    all_incomes_data, income_by_age_data, start_age_data = simulate_income_process(num_episodes=20000)
    
    # 分析并可视化结果
    analyze_and_plot(all_incomes_data, income_by_age_data, start_age_data)