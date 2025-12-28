import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from main_olg_v9_utils import OLG_V9_Utils
from simplified.main_olg_v9_utils_simplified import OLGEnvV9SACSimplified
import numpy as np

# 初始化参数
cS = OLG_V9_Utils.ParameterValues_HuggettStyle()
leLogGridV, leTrProbM, leProb1V = OLG_V9_Utils.EarningProcess_olgm(cS)
leGridV = np.exp(leLogGridV)

paramS_rl = {
    'leLogGridV': leLogGridV,
    'leTrProbM': leTrProbM,
    'leProb1V': leProb1V,
    'leGridV': leGridV,
    'ageEffV_new': cS.ageEffV_new
}

M_fixed = {
    'R_k_net_factor': 1.03,
    'w_gross': 2.0,
    'TR_total': 0.1,
    'b_payg_avg_retiree': 0.4,
    'tau_l': 0.15,
    'theta_payg_actual': 0.12
}

# 创建环境
env = OLGEnvV9SACSimplified(cS, paramS_rl, M_fixed, training_mode=False)

# 重置环境
obs, info = env.reset(seed=42)
print("🔍 调试储蓄计算过程")
print("="*50)
print(f"初始状态:")
print(f"  k = {env.current_k_val}")
print(f"  k_pps = {env.current_k_pps_val}")
print(f"  age_idx = {env.current_age_idx}")
print(f"  eps_idx = {env.current_eps_idx}")
print(f"  ε值 = {leGridV[env.current_eps_idx-1]:.4f}")

# 测试第一步
action = np.array([0.128, 0.279])  # 从测试输出中的第一步动作
print(f"\n步骤1 动作: PPS={action[0]:.3f}, 储蓄={action[1]:.3f}")

# 手动计算各步骤
print("\n详细计算过程:")

# 1. 计算当期收入
age_efficiency = cS.ageEffV_new[env.current_age_idx - 1]
epsilon_val = leGridV[env.current_eps_idx - 1]
gross_labor_income = M_fixed['w_gross'] * age_efficiency * epsilon_val
print(f"1. 毛劳动收入 = w * age_eff * ε = {M_fixed['w_gross']:.2f} * {age_efficiency:.4f} * {epsilon_val:.4f} = {gross_labor_income:.4f}")

# 2. PPS缴费计算
max_cpps_by_frac = gross_labor_income * cS.pps_max_contrib_frac
max_permissible_cpps = min(cS.pps_annual_contrib_limit, max_cpps_by_frac)
actual_c_pps = action[0] * max_permissible_cpps
print(f"2. PPS缴费:")
print(f"   最大允许缴费 = min({cS.pps_annual_contrib_limit}, {gross_labor_income:.4f} * {cS.pps_max_contrib_frac}) = {max_permissible_cpps:.4f}")
print(f"   实际PPS缴费 = {action[0]:.3f} * {max_permissible_cpps:.4f} = {actual_c_pps:.4f}")

# 3. 创建临时参数对象来计算收入
class TempParamSHH:
    def __init__(self, tau_l, theta_payg_actual_for_hh, pps_tax_deferral_active, ageEffV_new):
        self.tau_l = tau_l
        self.theta_payg_actual_for_hh = theta_payg_actual_for_hh
        self.pps_tax_deferral_active = pps_tax_deferral_active
        self.ageEffV_new = ageEffV_new
        self.current_pps_withdrawal = 0
        # 添加缺失的属性
        self.tau_k = 0.0  # 资本税率，简化版设为0
        self.pps_tax_rate_withdrawal = 0.15  # PPS提取税率

paramS_hh = TempParamSHH(
    M_fixed['tau_l'],
    M_fixed['theta_payg_actual'],
    cS.pps_active,
    cS.ageEffV_new
)

# 计算总收入
b_payg_this_age = 0  # 年轻时没有PAYG福利
total_income, _, _ = OLG_V9_Utils.HHIncome_Huggett(
    env.current_k_val,
    M_fixed['R_k_net_factor'],
    M_fixed['w_gross'],
    M_fixed['TR_total'],
    b_payg_this_age,
    actual_c_pps,
    env.current_age_idx - 1,
    paramS_hh,
    cS,
    epsilon_val
)

print(f"3. 总可支配收入 = {total_income:.4f}")

# 4. 储蓄和消费计算
consumption_floor_spending = cS.cFloor * (1 + cS.tau_c)
resources_for_kprime_and_c = total_income - consumption_floor_spending
print(f"4. 储蓄和消费:")
print(f"   最低消费支出 = {cS.cFloor:.3f} * (1 + {cS.tau_c:.3f}) = {consumption_floor_spending:.4f}")
print(f"   可用于储蓄和额外消费的资源 = {total_income:.4f} - {consumption_floor_spending:.4f} = {resources_for_kprime_and_c:.4f}")

if resources_for_kprime_and_c >= 0:
    k_next = action[1] * resources_for_kprime_and_c
    k_next = max(cS.kMin, min(k_next, resources_for_kprime_and_c))
    print(f"   下期资产 = {action[1]:.3f} * {resources_for_kprime_and_c:.4f} = {k_next:.4f}")
else:
    k_next = cS.kMin
    print(f"   资源不足，下期资产 = {cS.kMin}")

k_next = max(cS.kMin, min(k_next, cS.kMax))
consumption_expenditure = total_income - k_next
c_consumption = max(cS.cFloor, consumption_expenditure / (1 + cS.tau_c))

print(f"   最终下期资产 = {k_next:.4f}")
print(f"   消费支出 = {total_income:.4f} - {k_next:.4f} = {consumption_expenditure:.4f}")
print(f"   实际消费 = max({cS.cFloor:.3f}, {consumption_expenditure:.4f} / (1 + {cS.tau_c:.3f})) = {c_consumption:.4f}")

# 执行实际步骤验证
print(f"\n实际环境步骤验证:")
obs_next, reward, terminated, truncated, info = env.step(action)
print(f"环境计算结果:")
print(f"  下期k = {env.current_k_val:.4f}")
print(f"  下期k_pps = {env.current_k_pps_val:.4f}")
print(f"  消费 = {info.get('consumption', 'N/A')}")
print(f"  PPS缴费 = {info.get('c_pps', 'N/A')}")
print(f"  奖励 = {reward:.6f}") 