import numpy as np
import matplotlib.pyplot as plt

    # 设置字体路径
plt.rcParams['font.sans-serif'] = ['SimHei']  # 使用简体中文字体
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题

# ==============================================================================
# 1. 定义模型参数和决策环境
# ==============================================================================

# 经济学参数 (与你的Cocco环境一致)
GAMMA = 3.84          # 相对风险规避系数
MIN_CONSUMPTION = 0.01 # 实际的最低消费水平

# 单步决策问题设置
CASH_ON_HAND = 10.0  # 假设主体今天有10单位的资源
GRANULARITY = 10000  # 决策空间的精细度，越高越精确

# 创建一个消费决策空间 (所有可能的消费选择)
# 从最低消费到花光所有资源
consumption_choices = np.linspace(MIN_CONSUMPTION, CASH_ON_HAND, GRANULARITY)

# ==============================================================================
# 2. 定义两个奖励函数
# ==============================================================================

def calculate_crra_utility(C, gamma):
    """计算原始的CRRA效用函数 U(C)"""
    if gamma == 1:
        # 对数效用的特殊情况
        return np.log(C)
    else:
        return (C**(1 - gamma)) / (1 - gamma)

def calculate_transformed_reward(C, gamma, min_consumption):
    """
    计算经过对数变换后的奖励函数 R(C) = f(U(C))
    这里的 f(x) = log(x - x_min + 1)
    """
    # 步骤 A: 计算原始效用 U(C)
    raw_utility = calculate_crra_utility(C, gamma)
    
    # 步骤 B: 计算效用下限 U_min = U(C_min)
    min_utility = calculate_crra_utility(min_consumption, gamma)
    
    # 步骤 C: 平移效用值
    shifted_utility = raw_utility - min_utility
    
    # 步骤 D: 取对数 (这就是严格单调递增函数 f)
    transformed_reward = np.log(shifted_utility + 1.0)
    
    return transformed_reward

# ==============================================================================
# 3. 进行实验：为每个可能的消费选择计算两种奖励
# ==============================================================================

# 计算每个消费选择对应的原始效用值
utility_values = calculate_crra_utility(consumption_choices, GAMMA)

# 计算每个消费选择对应的变换后奖励值
transformed_reward_values = calculate_transformed_reward(consumption_choices, GAMMA, MIN_CONSUMPTION)

# ==============================================================================
# 4. 寻找最优决策并比较
# ==============================================================================

# 找到最大化 U(C) 的最优消费
max_utility_index = np.argmax(utility_values)
optimal_C_for_utility = consumption_choices[max_utility_index]
max_utility = utility_values[max_utility_index]

# 找到最大化 R(C) 的最优消费
max_reward_index = np.argmax(transformed_reward_values)
optimal_C_for_reward = consumption_choices[max_reward_index]
max_transformed_reward = transformed_reward_values[max_reward_index]

# --- 打印最终结论 ---
print("="*60)
print("       最优策略对奖励函数的单调变换的不变性证明")
print("="*60)
print(f"参数设置: gamma = {GAMMA}, 可用资源 = {CASH_ON_HAND}\n")

print(f"1. 基于原始 CRRA 效用 U(C):")
print(f"   - 最优消费 C* = {optimal_C_for_utility:.4f}")
print(f"   - 对应的最大效用 U(C*) = {max_utility:.4f}\n")

print(f"2. 基于对数变换后的奖励 R(C):")
print(f"   - 最优消费 C** = {optimal_C_for_reward:.4f}")
print(f"   - 对应的最大奖励 R(C**) = {max_transformed_reward:.4f}\n")

print("-" * 60)
print("结论:")
if np.isclose(optimal_C_for_utility, optimal_C_for_reward):
    print("✅ 证明成功: 两个奖励函数导出的最优消费策略完全相同。")
    print("   这表明，虽然奖励函数的形状和数值改变了，但最优决策保持不变。")
else:
    print("❌ 证明失败: 两个奖励函数导出的最优消费策略不同。")
print("="*60)
print("\n注意: 这个单步决策的结论可以通过贝尔曼方程归纳到多步动态规划问题中。")


# ==============================================================================
# 5. 可视化结果
# ==============================================================================

fig, (ax1, ax2) = plt.subplots(
    2, 1, figsize=(12, 10), sharex=True,
    gridspec_kw={'hspace': 0.3}
)
fig.suptitle("奖励函数变换对最优策略的影响", fontsize=16)

# --- 子图 1: 原始 CRRA 效用 ---
ax1.plot(consumption_choices, utility_values, label='U(C) - 原始CRRA效用', color='royalblue')
ax1.axvline(optimal_C_for_utility, color='red', linestyle='--', label=f'最优消费 C* = {optimal_C_for_utility:.4f}')
ax1.set_title("1. 原始奖励函数 U(C)", fontsize=14)
ax1.set_ylabel("效用值")
ax1.legend()
ax1.grid(True)
ax1.annotate(
    f'最大值点\nU(C*) = {max_utility:.2f}',
    xy=(optimal_C_for_utility, max_utility),
    xytext=(optimal_C_for_utility + 1, max_utility - 0.5),
    arrowprops=dict(facecolor='black', shrink=0.05),
    fontsize=12,
)


# --- 子图 2: 变换后的奖励 ---
ax2.plot(consumption_choices, transformed_reward_values, label='R(C) - 变换后的奖励', color='darkorange')
ax2.axvline(optimal_C_for_reward, color='red', linestyle='--', label=f'最优消费 C** = {optimal_C_for_reward:.4f}')
ax2.set_title("2. 变换后的奖励函数 R(C) = log(U(C) - U_min + 1)", fontsize=14)
ax2.set_xlabel("消费选择 (C)")
ax2.set_ylabel("变换后的奖励值")
ax2.legend()
ax2.grid(True)
ax2.annotate(
    f'最大值点\nR(C**) = {max_transformed_reward:.2f}',
    xy=(optimal_C_for_reward, max_transformed_reward),
    xytext=(optimal_C_for_reward + 1, max_transformed_reward - 0.5),
    arrowprops=dict(facecolor='black', shrink=0.05),
    fontsize=12,
)

plt.show()