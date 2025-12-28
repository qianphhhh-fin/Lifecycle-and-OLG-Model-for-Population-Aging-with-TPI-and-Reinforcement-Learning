import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
import numpy as np
import matlab.engine
from main_olg_v9_utils import OLG_V9_Utils, OLGEnvV9SACSimp, TempParamSHH
from compare_rl_and_vfi_matlab_simplified import RLVFIComparatorSimplified
from scipy.interpolate import RegularGridInterpolator

def get_vfi_policy_at_point(vfi_results, age_idx, eps_idx, k_now, k_pps_now):
    """辅助函数：从VFI策略中查询一个点的决策 (k', c_pps)"""
    # ... (与之前版本相同)
    kPolM, cPpsM = vfi_results['kPolM'], vfi_results['cPpsPolM_choice']
    cS = vfi_results['cS_python_obj']
    points = (cS.kGridV, cS.kppsGridV)
    age_idx, eps_idx = int(age_idx), int(eps_idx)
    k_interp = RegularGridInterpolator(points, kPolM[:,:,eps_idx,age_idx], bounds_error=False, fill_value=None)
    cpps_interp = RegularGridInterpolator(points, cPpsM[:,:,eps_idx,age_idx], bounds_error=False, fill_value=None)
    k_prime = k_interp((k_now, k_pps_now)).item()
    cpps = cpps_interp((k_now, k_pps_now)).item()
    return k_prime, cpps

def main():
    print("--- VFI 核心计算逻辑对比调试 (最终诊断版) ---")
    
    # --- 1. 设置要测试的状态点 ---
    # 使用您输出中 "手动版" 的状态作为我们唯一的、真实的输入
    TRACE_AGE = 5
    K_NOW_TEST = 0.00000000
    K_PPS_NOW_TEST = 1.65545379
    EPS_VAL_TEST = 0.50245600
    CPPS_DECISION_TEST = 0.54380292
    K_PRIME_DECISION_TEST = 0.00000000
    
    # --- 2. 获取必要的参数 ---
    comparator = RLVFIComparatorSimplified()
    vfi_results = comparator.run_matlab_vfi()
    cS_obj = vfi_results['cS_python_obj']
    paramS_dict = vfi_results['paramS_vfi_dict']
    M_fixed = comparator.M_FIXED
    
    # 找到 epsilon_val 对应的索引
    eps_idx_test = np.argmin(np.abs(paramS_dict['leGridV'] - EPS_VAL_TEST))

    print("\n" + "="*80)
    print(f"将在完全相同的输入点上进行对比计算:")
    print(f"  k={K_NOW_TEST}, k_pps={K_PPS_NOW_TEST}, age={TRACE_AGE}, eps={EPS_VAL_TEST}")
    print(f"  决策: cpps={CPPS_DECISION_TEST}, k'={K_PRIME_DECISION_TEST}")
    print("="*80)

    # --- 3. [方法A: 手动创建paramS对象] ---
    print("\n--- [方法A] 手动创建 paramS 对象并调用 ---")
    
    class TempParamS: pass
    paramS_A = TempParamS()
    paramS_A.tau_l = paramS_dict['tau_l']
    paramS_A.theta_payg_actual_for_hh = paramS_dict['theta_payg_actual_for_hh']
    paramS_A.pps_tax_deferral_active = paramS_dict['pps_tax_deferral_active']
    paramS_A.tau_k = cS_obj.tau_k
    # 工作期，withdrawal 相关的属性不重要
    paramS_A.current_pps_withdrawal = 0 
    
    b_payg = vfi_results['bV_payg_eq'][TRACE_AGE]
    res_A, _, _ = OLG_V9_Utils.HHIncome_Huggett(
        K_NOW_TEST, M_fixed['R_k_net_factor'], M_fixed['w_gross'], M_fixed['TR_total'],
        b_payg, CPPS_DECISION_TEST, TRACE_AGE, paramS_A, cS_obj, EPS_VAL_TEST
    )
    c_A = (res_A - K_PRIME_DECISION_TEST) / (1 + cS_obj.tau_c)
    
    print(f"  paramS_A.tau_l = {paramS_A.tau_l}")
    print(f"  Output: resources_after_pps={res_A:.8f}, c={c_A:.8f}")

    # --- 4. [方法B: 从Env中创建paramS对象] ---
    print("\n--- [方法B] 使用 Env 的方式创建 paramS 对象并调用 ---")

    # 这是 env._calculate_resources_after_pps 内部会做的事情
    paramS_B = TempParamSHH(
        M_fixed['tau_l'],
        M_fixed['theta_payg_actual'],
        cS_obj.pps_active,
        cS_obj.ageEffV_new
    )
    # 工作期，withdrawal 设为0
    paramS_B.current_pps_withdrawal = 0
    
    res_B, _, _ = OLG_V9_Utils.HHIncome_Huggett(
        K_NOW_TEST, M_fixed['R_k_net_factor'], M_fixed['w_gross'], M_fixed['TR_total'],
        b_payg, CPPS_DECISION_TEST, TRACE_AGE, paramS_B, cS_obj, EPS_VAL_TEST
    )
    c_B = (res_B - K_PRIME_DECISION_TEST) / (1 + cS_obj.tau_c)

    print(f"  paramS_B.tau_l = {paramS_B.tau_l}")
    print(f"  Output: resources_after_pps={res_B:.8f}, c={c_B:.8f}")
    
    print("\n🔍 对比完成。")
    del comparator

if __name__ == "__main__":
    main()