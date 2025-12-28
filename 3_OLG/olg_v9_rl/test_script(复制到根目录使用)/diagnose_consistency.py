# --- diagnose_consistency.py ---

"""
======================================================================
== OLG V9 - VFI与RL环境一致性诊断脚本 ==
======================================================================

目的:
验证 OLGEnvV9SAC (RL环境) 的物理过程与 VFI 优化问题的内在逻辑
是否完全一致。这对于保证RL与VFI比较的公平性至关重要。

方法:
1. 设定一个固定的宏观经济环境和具体的个体状态 (k, k_pps, age, ε)。
2. VFI视角:
   - 手动模拟VFI求解器的一步。
   - 对一个给定的决策 (k', c_pps)，计算其总价值。
3. RL视角:
   - 将VFI的决策“翻译”成RL的动作(action)。
   - 调用RL环境的 step() 函数。
4. 交叉验证:
   - 比较两个视角下计算出的关键中间变量(资源、消费、效用、下期资产)
     是否完全相等。

如果所有断言(assert)都通过，则证明两个框架在物理层面上是一致的。
"""

import numpy as np
import warnings
from main_olg_v9_utils import OLG_V9_Utils, OLGEnvV9SAC, TempParamSHH

# 抑制不必要的警告
warnings.filterwarnings('ignore', category=UserWarning)

def run_consistency_check():
    """执行VFI与RL环境的一致性检查"""
    print("\n" + "="*70)
    print("== 正在运行VFI与RL环境的一致性诊断脚本 ==")
    print("="*70)

    # --- 1. 初始化和参数设置 ---
    print("\n--- [1] 初始化模型参数和环境 ---")
    
    # 使用标准参数
    cS = OLG_V9_Utils.ParameterValues_HuggettStyle()
    
    # 模拟VFI求解时需要的参数
    leLogGridV, leTrProbM, leProb1V = OLG_V9_Utils.EarningProcess_olgm(cS)
    paramS_rl = {
        'leLogGridV': leLogGridV, 'leTrProbM': leTrProbM, 'leProb1V': leProb1V,
        'leGridV': np.exp(leLogGridV), 'ageEffV_new': cS.ageEffV_new
    }
    
    # 设定一个固定的宏观环境 (用于测试)
    M_fixed = {
        'R_k_net_factor': 1.03, 'w_gross': 2.0, 'TR_total': 0.1,
        'b_payg_avg_retiree': 0.4, 'tau_l': 0.15, 'theta_payg_actual': 0.12
    }
    print(f"固定宏观环境 (M_fixed): {M_fixed}")

    # 创建一个RL环境实例，作为我们的“物理引擎”
    # rng_M 在这里只是为了满足构造函数，实际使用的是M_fixed
    rng_M_dummy = {key: [val, val] for key, val in M_fixed.items()}
    env = OLGEnvV9SAC(cS, paramS_rl, rng_M_dummy, training_mode=False)
    
    # --- 2. 设定一个具体的测试状态和决策 ---
    print("\n--- [2] 设定测试状态和决策 ---")
    
    # a. 个体状态 (State)
    test_age_idx = 5          # 0-based, 工作期
    test_k_state = 10.0       # 当前非PPS资产
    test_k_pps_state = 5.0    # 当前PPS资产
    test_eps_idx = 3          # 1-based, 效率冲击索引
    
    # b. 个体决策 (Decision / Action)
    test_k_prime_decision = 11.0  # 决策: 下一期非PPS资产
    test_cpps_decision = 0.5    # 决策: 当期PPS缴费
    
    print(f"测试状态: Age={test_age_idx+1}, k={test_k_state}, k_pps={test_k_pps_state}, ε_idx={test_eps_idx}")
    print(f"测试决策: k'={test_k_prime_decision}, c_pps={test_cpps_decision}")
    
    # --- 3. VFI视角: 手动计算决策后果 ---
    print("\n--- [3] VFI视角: 手动计算决策后果 ---")
    
    # a. 准备调用HHIncome_Huggett所需的参数
    b_payg_vfi = 0 # 因为是工作期
    paramS_hh_vfi = TempParamSHH(M_fixed['tau_l'], M_fixed['theta_payg_actual'], cS.pps_active, cS.ageEffV_new)
    epsilon_val_vfi = paramS_rl['leGridV'][test_eps_idx-1]
    
    # b. 计算总资源
    resources_vfi, _, _ = OLG_V9_Utils.HHIncome_Huggett(
        k_now_val=test_k_state, 
        R_k_net_factor=M_fixed['R_k_net_factor'],
        w_gross=M_fixed['w_gross'],
        TR_total=M_fixed['TR_total'],
        b_payg_val=b_payg_vfi,
        c_pps_chosen=test_cpps_decision,
        a_idx=test_age_idx,
        paramS_hh=paramS_hh_vfi,
        cS=cS,
        epsilon_val=epsilon_val_vfi
    )
    
    # c. 计算消费和效用
    consumption_expenditure_vfi = resources_vfi - test_k_prime_decision
    c_vfi = max(cS.cFloor, consumption_expenditure_vfi / (1 + cS.tau_c))
    _, u_vfi = OLG_V9_Utils.CES_utility(c_vfi, cS.sigma, cS)

    # d. 计算下一期资产
    k_next_vfi = test_k_prime_decision
    pps_return_factor_vfi = 1 + ((M_fixed['R_k_net_factor'] - 1) + cS.pps_return_rate_premium)
    k_pps_next_vfi = (test_k_pps_state + test_cpps_decision) * pps_return_factor_vfi
    
    print(f"  - VFI计算: Resources = {resources_vfi:.6f}")
    print(f"  - VFI计算: Consumption (c) = {c_vfi:.6f}")
    print(f"  - VFI计算: Utility u(c) = {u_vfi:.6f}")
    print(f"  - VFI计算: k_next = {k_next_vfi:.6f}, k_pps_next = {k_pps_next_vfi:.6f}")
    
    # --- 4. RL视角: 通过环境的step函数计算后果 ---
    print("\n--- [4] RL视角: 通过环境的step()计算后果 ---")
    
    # a. 手动将环境设置为我们的测试状态
    env.current_age_idx = test_age_idx + 1 # env使用1-based age
    env.current_k_val = test_k_state
    env.current_k_pps_val = test_k_pps_state
    env.current_eps_idx = test_eps_idx
    env.current_M = M_fixed
    env.current_bV_payg.fill(0) # 确保PAYG福利为0
    
    # b. 将VFI决策“翻译”成RL的action
    #   Action 1: PPS缴费比例
    #   首先计算出在该状态下最大可能的PPS缴费
    _, max_permissible_cpps = env._calculate_pps_contribution(prop_pps_contrib=1.0)
    # 计算比例
    action_cpps_prop = test_cpps_decision / max_permissible_cpps if max_permissible_cpps > 1e-9 else 0.0

    #   Action 2: 消费比例
    #   首先计算出在给定c_pps决策下的可用资源
    resources_rl_pre_c = env._calculate_resources_after_pps(test_cpps_decision)
    #   计算消费支出
    consumption_expenditure_rl = resources_rl_pre_c - test_k_prime_decision
    #   计算消费占“可自由支配”资源的比例
    c_floor_spending = cS.cFloor * (1 + cS.tau_c)
    spendable_resources = resources_rl_pre_c - c_floor_spending
    consumption_above_floor_spending = consumption_expenditure_rl - c_floor_spending
    action_consump_prop = consumption_above_floor_spending / spendable_resources if spendable_resources > 1e-9 else 0.0
    
    # 组合成最终的action
    test_action = np.array([action_cpps_prop, action_consump_prop])
    print(f"  - 翻译成RL Action: [prop_pps, prop_consump] = [{test_action[0]:.6f}, {test_action[1]:.6f}]")

    # c. 执行一步，并从info字典中获取结果
    #    我们不关心obs, reward, terminated, truncated，只关心info
    _, _, _, _, info = env.step(test_action)
    
    # d. RL环境计算出的结果
    resources_rl = resources_rl_pre_c # 这是step函数内部计算的总资源
    c_rl = info['consumption']
    u_rl = info['pure_utility']
    k_next_rl = info['k_prime']
    k_pps_next_rl = env.current_k_pps_val # step函数执行后，环境的k_pps值已经更新为下一期的值

    print(f"  - RL环境计算: Resources = {resources_rl:.6f}")
    print(f"  - RL环境计算: Consumption (c) = {c_rl:.6f}")
    print(f"  - RL环境计算: Utility u(c) = {u_rl:.6f}")
    print(f"  - RL环境计算: k_next = {k_next_rl:.6f}, k_pps_next = {k_pps_next_rl:.6f}")

    # --- 5. 交叉验证和断言 ---
    print("\n--- [5] 交叉验证结果 ---")
    
    try:
        # 使用np.isclose进行浮点数比较，容忍极小的计算误差
        tolerance = 1e-9
        
        assert np.isclose(resources_vfi, resources_rl, atol=tolerance), "总资源计算不一致！"
        print("✅ [通过] 总资源 (Resources) 计算一致。")
        
        assert np.isclose(c_vfi, c_rl, atol=tolerance), "消费 (c) 计算不一致！"
        print("✅ [通过] 消费 (c) 计算一致。")
        
        assert np.isclose(u_vfi, u_rl, atol=tolerance), "效用 (u(c)) 计算不一致！"
        print("✅ [通过] 效用 (u(c)) 计算一致。")
        
        assert np.isclose(k_next_vfi, k_next_rl, atol=tolerance), "下一期非PPS资产 (k') 计算不一致！"
        print("✅ [通过] 下一期非PPS资产 (k') 计算一致。")

        assert np.isclose(k_pps_next_vfi, k_pps_next_rl, atol=tolerance), "下一期PPS资产 (k_pps') 计算不一致！"
        print("✅ [通过] 下一期PPS资产 (k_pps') 计算一致。")

        print("\n" + "*"*70)
        print("🎉 恭喜！所有检查点通过。VFI优化问题与RL环境的物理过程一致。")
        print("*"*70)

    except AssertionError as e:
        print("\n" + "!"*70)
        print(f"❌ [失败] 诊断失败: {e}")
        print("!"*70)
        print("请检查 `HHIncome_Huggett` 和 `OLGEnvV9SAC.step` 中的相关计算逻辑。")

if __name__ == '__main__':
    run_consistency_check()