"""
OLG V9 模型主求解器

该脚本支持三种求解方式：
1. Python VFI (--solver vfi): 使用Python实现的值函数迭代
2. MATLAB VFI (--solver matlab_vfi): 调用main_olg_v8_utils.m中的MATLAB VFI函数
3. SAC强化学习 (--solver sac): 使用预训练的SAC模型

使用示例：
- Python VFI: python main_olg_v9.py --solver vfi
- MATLAB VFI: python main_olg_v9.py --solver matlab_vfi  
- SAC求解: python main_olg_v9.py --solver sac --sac_model_path path/to/model.zip

注意：使用MATLAB VFI需要安装MATLAB引擎 for Python
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import torch as th
from typing import Dict, Any, Tuple, Callable
import time
import logging
import os

# 添加MATLAB引擎支持
try:
    import matlab.engine
    MATLAB_AVAILABLE = True
    print("✅ MATLAB引擎可用")
except ImportError:
    MATLAB_AVAILABLE = False
    print("⚠️ MATLAB引擎不可用，将只使用Python VFI")

# 尝试导入SBX，如果失败则使用SB3
try:
    from sbx import SAC as SBX_SAC
    SBX_AVAILABLE = True
    print("✅ SBX (Stable Baselines Jax) 可用")
except ImportError:
    SBX_AVAILABLE = False
    print("⚠️ SBX 不可用，将使用 SB3")

# 导入SB3作为备选
from stable_baselines3 import SAC as SB3_SAC

from main_olg_v8_utils import OLGV9Utils

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# OnnxablePolicy for loading SB3 models with specific architectures
class OnnxablePolicy(th.nn.Module):
    def __init__(self, actor: th.nn.Module):
        super().__init__()
        self.actor = actor

    def forward(self, observation: th.Tensor) -> th.Tensor:
        # The deterministic action is the mean of the distribution
        return self.actor(observation)

def matlab_vfi_solver_wrapper() -> Callable:
    """
    创建一个使用MATLAB VFI求解器的包装函数
    该包装函数调用main_olg_v8_utils.m中的HHSolution_VFI_Huggett函数
    """
    if not MATLAB_AVAILABLE:
        raise ImportError("MATLAB引擎不可用，无法使用MATLAB VFI求解器")
    
    # 启动MATLAB引擎
    logging.info("启动MATLAB引擎...")
    eng = matlab.engine.start_matlab()
    
    # 添加当前目录到MATLAB路径，以便找到main_olg_v8_utils.m
    eng.addpath(os.getcwd(), nargout=0)
    
    logging.info("MATLAB引擎已启动并配置完成")

    def solve_hh_problem_matlab_vfi(R_k_net_factor: float, MPL_gross: float, TR_total: float, 
                                   bV_new: np.ndarray, paramS: Dict[str, Any], cS: Dict[str, Any]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        使用MATLAB VFI求解器求解家庭问题
        调用main_olg_v8_utils.HHSolution_VFI_Huggett函数
        """
        try:
            logging.info("调用MATLAB VFI求解器...")
            
            # 将Python参数转换为MATLAB格式
            matlab_R_k_net_factor = float(R_k_net_factor)
            matlab_MPL_gross = float(MPL_gross)
            matlab_TR_total = float(TR_total)
            matlab_bV_new = matlab.double(bV_new.tolist())
            
            # 辅助函数：转换Python值到MATLAB格式
            def convert_to_matlab(value):
                if isinstance(value, np.ndarray):
                    if value.ndim == 1:
                        return matlab.double(value.tolist())
                    elif value.ndim == 2:
                        return matlab.double(value.tolist())
                    elif value.ndim == 3:
                        # 对于3D数组，转换为cell数组
                        return matlab.double(value.tolist())
                    elif value.ndim == 4:
                        # 对于4D数组，转换为cell数组
                        return matlab.double(value.tolist())
                    else:
                        return matlab.double(value.tolist())
                elif isinstance(value, (int, float, np.integer, np.floating)):
                    return float(value)
                elif isinstance(value, bool):
                    return bool(value)
                elif isinstance(value, str):
                    return value
                elif isinstance(value, list):
                    return matlab.double(value) if all(isinstance(x, (int, float)) for x in value) else value
                else:
                    return value
            
            # 构建MATLAB参数结构体
            matlab_paramS = {}
            for key, value in paramS.items():
                try:
                    matlab_paramS[key] = convert_to_matlab(value)
                except Exception as e:
                    logging.warning(f"无法转换paramS[{key}]: {e}, 使用原值")
                    matlab_paramS[key] = value
            
            matlab_cS = {}
            for key, value in cS.items():
                try:
                    matlab_cS[key] = convert_to_matlab(value)
                except Exception as e:
                    logging.warning(f"无法转换cS[{key}]: {e}, 使用原值")
                    matlab_cS[key] = value
            
            # 调用MATLAB函数
            logging.info("执行MATLAB VFI计算...")
            logging.info(f"参数: R_k_net_factor={matlab_R_k_net_factor:.4f}, MPL_gross={matlab_MPL_gross:.4f}")
            logging.info(f"参数: TR_total={matlab_TR_total:.4f}, bV_new长度={len(bV_new)}")
            
            # 检查MATLAB类是否可用
            try:
                # 测试MATLAB类方法是否存在
                result = eng.main_olg_v8_utils.HHSolution_VFI_Huggett(
                    matlab_R_k_net_factor,
                    matlab_MPL_gross, 
                    matlab_TR_total,
                    matlab_bV_new,
                    matlab_paramS,
                    matlab_cS,
                    nargout=4
                )
            except Exception as matlab_error:
                logging.error(f"MATLAB函数调用失败: {matlab_error}")
                # 尝试直接调用函数而不是类方法
                try:
                    logging.info("尝试直接调用MATLAB函数...")
                    result = eng.HHSolution_VFI_Huggett(
                        matlab_R_k_net_factor,
                        matlab_MPL_gross, 
                        matlab_TR_total,
                        matlab_bV_new,
                        matlab_paramS,
                        matlab_cS,
                        nargout=4
                    )
                except Exception as direct_error:
                    logging.error(f"直接调用MATLAB函数也失败: {direct_error}")
                    raise matlab_error
            
            # 将MATLAB结果转换回Python numpy数组
            cPolM_q = np.array(result[0])
            kPolM = np.array(result[1]) 
            cPpsPolM_choice = np.array(result[2])
            valM = np.array(result[3])
            
            logging.info("MATLAB VFI求解完成")
            
            return cPolM_q, kPolM, cPpsPolM_choice, valM
            
        except Exception as e:
            logging.error(f"MATLAB VFI求解器调用失败: {str(e)}")
            logging.info("回退到Python VFI求解器...")
            # 回退到Python VFI求解器
            return OLGV9Utils.hh_solution_vfi_huggett(R_k_net_factor, MPL_gross, TR_total, bV_new, paramS, cS)

    return solve_hh_problem_matlab_vfi

def sac_solver_wrapper(model: Any) -> Callable:
    """
    Creates a wrapper function that uses a pre-trained SAC model to solve the household problem.
    This wrapper matches the signature of the VFI solver.
    
    Args:
        model: The loaded SAC model object.
    """

    def solve_hh_problem_sac(R_k_net_factor: float, MPL_gross: float, TR_total: float, 
                             bV_new: np.ndarray, paramS: Dict[str, Any], cS: Dict[str, Any]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, None]:
        """
        The actual solver function that uses the loaded SAC model.
        It uses the cS and paramS passed from the GE solver.
        """
        nk, nkpps, nw, aD_new = cS['nk'], cS['nkpps'], cS['nw'], cS['aD_new']
        k_grid, kpps_grid = cS['kGridV'], cS['kppsGridV']
        
        cPolM = np.zeros((nk, nkpps, nw, aD_new))
        kPolM = np.zeros((nk, nkpps, nw, aD_new))
        cPpsPolM = np.zeros((nk, nkpps, nw, aD_new))

        # Iterate over the state space to get policy for each state
        for a in range(aD_new):
            for i_e in range(nw):
                for ik_pps in range(nkpps):
                    for ik in range(nk):
                        # Construct observation based on how the RL agent was trained
                        k_val = k_grid[ik]
                        k_pps_val = kpps_grid[ik_pps]
                        epsilon_val = paramS['leGridV'][i_e]

                        # --- Construct observation vector for the RL agent ---
                        # This vector MUST EXACTLY match the structure and normalization
                        # used in the training environment (OLGEnvV8SAC).
                        
                        # Get all relevant macro variables from the GE loop
                        b_payg_val = bV_new[a]
                        tau_l_val = paramS['tau_l']
                        theta_payg_val = paramS['theta_payg_actual_for_hh']

                        # Normalize macro variables. Ranges are inferred from the training script.
                        # The assumed formula is: (value - midpoint) / half_range
                        norm_R = (R_k_net_factor - 1.03) / 0.02 
                        norm_MPL = (MPL_gross - 2.0) / 0.5
                        norm_TR = (TR_total - 0.1) / 0.1 if TR_total > 1e-9 else 0.0
                        norm_b = (b_payg_val - 0.45) / 0.35 if a >= cS['aR_new'] else -1.0 # Use -1 for non-retirees
                        norm_tau_l = (tau_l_val - 0.15) / 0.1
                        norm_theta = (theta_payg_val - 0.125) / 0.075

                        # Micro-state variables
                        norm_k = k_val / cS['k_max']
                        norm_k_pps = k_pps_val / cS['k_pps_max']
                        norm_age = a / cS['aD_new']
                        norm_shock = (paramS['leLogGridV'][i_e] - cS['le_mu']) / 3.0

                        obs_vec = np.array([
                            norm_k,
                            norm_k_pps,
                            norm_age,
                            norm_shock,
                            norm_R,
                            norm_MPL,
                            norm_TR,
                            norm_b,
                            norm_tau_l,
                            norm_theta
                        ])
                        
                        # Safeguard check: if agent expects different input size, truncate.
                        # A better solution is to ensure the training environment and this
                        # context are perfectly aligned.
                        if model.observation_space.shape[0] != len(obs_vec):
                           logging.warning(
                               f"Observation size mismatch! Agent expects {model.observation_space.shape[0]}, "
                               f"but GE context provides {len(obs_vec)}. Truncating observation."
                           )
                           obs_vec = obs_vec[:model.observation_space.shape[0]]

                        obs_tensor = th.as_tensor(obs_vec, dtype=th.float32).unsqueeze(0)
                        
                        action, _ = model.predict(obs_tensor, deterministic=True)
                        
                        # Rescale actions from [-1, 1] to the actual policy space
                        k_prime_norm, c_pps_norm = action[0]
                        
                        k_prime = (k_prime_norm + 1) / 2 * (cS['k_max'] - cS['k_min']) + cS['k_min']
                        
                        # Max c_pps depends on income, using prices from the GE solver
                        _, _, labor_inc_gross = OLGV9Utils.hh_income_huggett(
                            k_val, R_k_net_factor, MPL_gross, TR_total, bV_new[a], 0, a, paramS, cS, epsilon_val
                        )
                        max_pps_contrib = 0
                        if a < cS['aR_new']:
                           max_pps_contrib = min(cS['pps_annual_contrib_limit'], 
                                                 cS['pps_max_contrib_frac'] * labor_inc_gross)
                        
                        c_pps = (c_pps_norm + 1) / 2 * max_pps_contrib

                        kPolM[ik, ik_pps, i_e, a] = k_prime
                        cPpsPolM[ik, ik_pps, i_e, a] = c_pps

                        # Consumption is derived from the budget constraint
                        y_disp, _, _ = OLGV9Utils.hh_income_huggett(
                            k_val, R_k_net_factor, MPL_gross, TR_total, bV_new[a], c_pps, a, paramS, cS, epsilon_val
                        )
                        c = (y_disp - k_prime) / (1 + cS['tau_c'])
                        cPolM[ik, ik_pps, i_e, a] = max(c, 1e-6)

        return cPolM, kPolM, cPpsPolM, None # Return None for value function

    return solve_hh_problem_sac

def main():
    parser = argparse.ArgumentParser(description="Solve the OLG V9 model using VFI or a trained SAC agent.")
    parser.add_argument('--solver', type=str, choices=['vfi', 'matlab_vfi', 'sac'], default='vfi',
                        help="The method to solve the household's problem ('vfi', 'matlab_vfi', or 'sac').")
    parser.add_argument('--sac_model_path', type=str, default=os.path.join('py', 'best_model_sbx', 'best_model.zip'),
                        help="Path to the pre-trained SAC model file (e.g., 'best_model.zip').")
    parser.add_argument('--use_sbx', type=bool, default=True,
                        help="Whether to use SBX (if available) or SB3 for SAC models.")
    args = parser.parse_args()
    
    # --- Header ---
    logging.info('=== OLG 模型 V9 (Python, 内生PPS缴费, 固定 Rho_prime_payg) ===')
    logging.info(f'    求解器 (Solver): {args.solver.upper()}')
    logging.info('    (Rho_prime_payg 固定, TR_gov=0, tau_l 内生, theta_payg 有上限)')
    logging.info('    (VFI 状态: k, k_pps, eps; PPS缴费: 内生选择，有比例和绝对上限)')
    
    # --- 1. Initialize parameters ---
    logging.info("\n--- 1. 初始化模型参数 ---")
    cS = OLGV9Utils.set_parameters('default')
    paramS = {} # Initialize empty paramS dictionary

    # Set derived parameters (beta, delta, grids, etc.)
    paramS = OLGV9Utils.parameter_values_huggett_style(cS, paramS)
    
    logging.info(f"模型已使用 '{args.solver}' 求解器进行初始化。")
    logging.info(f"贴现因子 (beta, 模型周期): {paramS['beta']:.4f}, 折旧率 (delta, 模型周期): {paramS['delta']:.4f}")
    logging.info(f'>>> V9: 固定 PAYG 替代率 (rho_prime_payg_fixed): {cS["rho_prime_payg_fixed"]:.3f}')
    logging.info(f"参数已加载。nk={cS['nk']}, nkpps={cS['nkpps']}, nw={cS['nw']}, nPpsChoiceGrid={cS['n_pps_choice_grid_points']}。")
    logging.info(f"年度年龄范围: {cS['age1_orig']}-{cS['ageLast_orig']}。模型年龄组数: {cS['aD_new']}。")
    logging.info(f"固定税率: tau_k={cS['tau_k']:.2f}, tau_c={cS['tau_c']:.2f}。G/Y={cS['gov_exp_frac_Y']:.2f}, B/Y={cS['gov_debt_frac_Y']:.2f}。")
    logging.info(f"PAYG 税率上限 (theta_payg_max): {cS['theta_payg_max']:.3f}")
    logging.info(f"所得税率 tau_l 范围: [{cS['tau_l_min']:.3f}, {cS['tau_l_max']:.3f}], 总劳动税上限: {cS['max_total_labor_tax']:.3f}")
    logging.info(f"PPS 年度缴费上限 (绝对值): {cS['pps_annual_contrib_limit']:.2f}, 比例上限 (法定): {cS['pps_max_contrib_frac']:.2f}")

    # --- 2. Population dynamics ---
    logging.info("\n--- 2. 模拟人口动态 ---")
    popS = OLGV9Utils.init_population(cS)
    popS = OLGV9Utils.population_dynamics(popS, cS)
    Z_ss_group_mass, Z_ss_annual, _, _ = OLGV9Utils.detect_steady_state_population(popS, cS)
    
    # Normalize the group mass to get the distribution for ageMassV
    if np.sum(Z_ss_group_mass) > 1e-9:
        Z_ss_group_dist = Z_ss_group_mass / np.sum(Z_ss_group_mass)
    else:
        Z_ss_group_dist = np.zeros_like(Z_ss_group_mass)

    paramS['Z_ss_counts'] = Z_ss_group_dist * cS['nSim'] # Use distribution for counts
    paramS['ageMassV'] = Z_ss_group_dist # ageMassV must be a distribution
    paramS['Z_ss_norm_annual'] = Z_ss_annual
    paramS['mass_workers_group'] = np.sum(paramS['ageMassV'][:cS['aR_new']])
    
    if len(popS['totalPop_history']) > 1 and popS['totalPop_history'][-2] > 1e-9:
        pop_growth_factor_per_group = popS['totalPop_history'][-1] / popS['totalPop_history'][-2]
        paramS['popGrowthForDebt'] = pop_growth_factor_per_group**(1/cS['yearStep']) - 1
    else:
        paramS['popGrowthForDebt'] = cS['popGrowth_orig']
    
    logging.info(f"人口参数计算完毕。年化稳态人口增长率 (用于债务): {paramS['popGrowthForDebt']:.4f}")
    logging.info(f"稳态工人占比 (基于年龄组人口): {paramS['mass_workers_group']:.4f}")

    # --- 3. Pre-calculate labor supply ---
    logging.info("\n--- 3. 预计算劳动 ---")
    paramS['leLogGridV'], paramS['leTrProbM'], paramS['leProb1V'] = OLGV9Utils.earning_process_olgm(cS)
    paramS['leGridV'] = np.exp(paramS['leLogGridV'])
    paramS['ageEffV_new'] = cS['ageEffV_new']
    eIdxM = OLGV9Utils.labor_endow_simulation_olgm(cS, paramS)
    L_per_capita, _ = OLGV9Utils.labor_supply_huggett(eIdxM, cS, paramS, paramS['ageMassV'])
    paramS['L_per_capita'] = L_per_capita
    logging.info(f"总劳动供给 (L, 效率单位, 总体人均): {L_per_capita:.4f}")
    
    logging.info('劳动效率状态网格值：')
    for i, val in enumerate(paramS['leGridV']):
        logging.info(f'状态 {i+1}: {val:.4f}')
    logging.info('\n初始分布概率：')
    for i, val in enumerate(paramS['leProb1V']):
        logging.info(f'状态 {i+1}: {val:.4f}')

    # --- 4. Solve for the general equilibrium ---
    logging.info(f"\n--- 4. 求解一般均衡 (固定 rho_prime_payg_fixed={cS['rho_prime_payg_fixed']:.3f}) ---")
    logging.info(f"  当前格点参数: n_k={cS['nk']}, n_kpps={cS['nkpps']}, n_w={cS['nw']}, n_PpsChoiceGrid={cS['n_pps_choice_grid_points']}")
    logging.info(f"调用均衡求解器 (Solver: {args.solver.upper()})...")
    
    # Select the household solver function based on arguments
    if args.solver == 'vfi':
        household_solver = OLGV9Utils.hh_solution_vfi_huggett
        logging.info("使用Python VFI求解器")
    elif args.solver == 'matlab_vfi':
        if not MATLAB_AVAILABLE:
            logging.error("MATLAB引擎不可用，无法使用MATLAB VFI求解器")
            logging.info("回退到Python VFI求解器")
            household_solver = OLGV9Utils.hh_solution_vfi_huggett
        else:
            household_solver = matlab_vfi_solver_wrapper()
            logging.info("使用MATLAB VFI求解器 (main_olg_v8_utils.m)")
    elif args.solver == 'sac':
        if not args.sac_model_path or not os.path.exists(args.sac_model_path):
            raise FileNotFoundError(f"SAC model not found at path: {args.sac_model_path}")
        
        # 根据可用性和用户选择确定使用哪个SAC实现
        use_sbx = args.use_sbx and SBX_AVAILABLE
        SAC_class = SBX_SAC if use_sbx else SB3_SAC
        backend_name = "SBX" if use_sbx else "SB3"
        
        logging.info(f"Loading SAC model from: {args.sac_model_path}")
        logging.info(f"Using {backend_name} SAC implementation")
        
        sac_model = SAC_class.load(args.sac_model_path)
        household_solver = sac_solver_wrapper(sac_model)
        logging.info("使用SAC强化学习求解器")
    else:
        raise NotImplementedError(f"Solver '{args.solver}' not recognized.")

    K_global_guess = 3.0 # Initial guess for aggregate capital
    solve_start_time = time.time()
    
    # Call the equilibrium solver with MATLAB-compatible signature
    result_dict = OLGV9Utils.solve_k_tau_l_for_rho_prime(
        cS['rho_prime_payg_fixed'], K_global_guess, household_solver, cS, paramS, eIdxM
    )
    
    K_eq = result_dict['K_sol']
    tau_l_eq = result_dict['tau_l_sol']
    gbc_residual_eq = result_dict['gbc_res_final']
    eq_found = result_dict['converged']
    final_eq_details = result_dict['details']
    
    solve_time = time.time() - solve_start_time
    logging.info(f'均衡求解完成。耗时: {solve_time:.2f} 秒。')

    if not eq_found:
        logging.error(f"未能为固定的 rho_prime_payg_fixed = {cS['rho_prime_payg_fixed']:.3f} 找到均衡解。")
        return

    logging.info(f'  均衡求解器返回状态: eq_found = {eq_found}')
    logging.info(f'  均衡结果: K_eq = {K_eq:.4f}, tau_l_eq = {tau_l_eq:.4f}, GBC 残差 = {gbc_residual_eq:.3e}')
    
    theta_payg_optimal_calc = final_eq_details.get('theta_payg_required_before_cap', np.nan)
    logging.info(f'  (理论所需PAYG税率，未考虑上限前: {theta_payg_optimal_calc:.4f})')
    
    if abs(gbc_residual_eq) > cS.get('gbc_tol_for_internal_loop', 1e-5) * 10:
        logging.warning(f'最终均衡的GBC残差 ({gbc_residual_eq:.2e}) 较大。可能均衡解的质量不高。')

    # --- 5. Analyze and plot final equilibrium results ---
    logging.info(f"\n--- 5. 最终均衡结果与绘图 (rho_prime_payg_fixed={cS['rho_prime_payg_fixed']:.4f}, tau_l_eq={tau_l_eq:.4f}, TR_gov=0) ---")

    # Re-calculate final state for detailed reporting
    paramS_eq = paramS.copy()
    paramS_eq['tau_l'] = tau_l_eq
    paramS_eq['theta_payg_actual_for_hh'] = final_eq_details.get('theta_payg', 0.0)
    paramS_eq['pps_tax_deferral_active'] = cS['pps_active']

    R_mkt_gross_factor_eq_final = final_eq_details.get('R_mkt_gross_factor', 1.04)
    MPL_gross_eq_final = final_eq_details.get('MPL_gross', 2.0)
    r_mkt_gross_eq_final = R_mkt_gross_factor_eq_final - 1
    r_k_net_hh_eq_final = r_mkt_gross_eq_final * (1 - cS['tau_k'])
    R_k_net_factor_hh_eq_final = 1 + r_k_net_hh_eq_final
    
    avg_worker_gross_wage_eq_final = (MPL_gross_eq_final * paramS.get('L_per_capita', 0)) / paramS.get('mass_workers_group', 1)
    b_payg_eq_final = cS['rho_prime_payg_fixed'] * avg_worker_gross_wage_eq_final
    bV_eq_new_final = np.zeros(cS['aD_new'])
    if cS['aR_new'] < cS['aD_new']:
        bV_eq_new_final[cS['aR_new']:] = b_payg_eq_final
    
    T_bequest_eq_final = final_eq_details.get('T_bequest_aggr', 0.0)
    TR_total_eq_final_for_vfi = T_bequest_eq_final # Assuming TR_gov=0

    logging.info(f"最终 VFI 调用参数: MPL_gross={MPL_gross_eq_final:.4f}, tau_l={paramS_eq['tau_l']:.4f}, theta_payg_actual={paramS_eq['theta_payg_actual_for_hh']:.4f}, TR_total={TR_total_eq_final_for_vfi:.4f} (T_bequest)")

    # 根据所使用的求解器选择最终的策略函数计算方法
    if args.solver == 'matlab_vfi' and MATLAB_AVAILABLE:
        logging.info("调用最终的 MATLAB HHSolution_VFI_Huggett...")
        cPolM_plot, kPolM_plot, cPpsPolM_plot, _ = household_solver(
            R_k_net_factor_hh_eq_final, MPL_gross_eq_final, TR_total_eq_final_for_vfi, bV_eq_new_final, paramS_eq, cS
        )
    else:
        logging.info("调用最终的 Python HHSolution_VFI_Huggett (V9)...")
        cPolM_plot, kPolM_plot, cPpsPolM_plot, _ = OLGV9Utils.hh_solution_vfi_huggett(
            R_k_net_factor_hh_eq_final, MPL_gross_eq_final, TR_total_eq_final_for_vfi, bV_eq_new_final, paramS_eq, cS
        )

    logging.info("模拟最终均衡的分布 (HHSimulation_olgm V9)...")
    kHistM_eq, kPpsHistM_eq, cHistM_eq = OLGV9Utils.hh_simulation_olgm(
        kPolM_plot, cPpsPolM_plot, cPolM_plot, eIdxM,
        R_k_net_factor_hh_eq_final, MPL_gross_eq_final,
        TR_total_eq_final_for_vfi, bV_eq_new_final, paramS_eq, cS
    )

    K_nonpps_eq_agg = np.mean(kHistM_eq, axis=0) @ paramS['Z_ss_norm_annual']
    K_pps_eq_agg = np.mean(kPpsHistM_eq, axis=0) @ paramS['Z_ss_norm_annual'] if cS['pps_active'] else 0
    Actual_K_eq_final = K_nonpps_eq_agg + K_pps_eq_agg
    C_eq_final = np.mean(cHistM_eq, axis=0) @ paramS['Z_ss_norm_annual']
    Y_eq_final = cS['A'] * (Actual_K_eq_final ** cS['alpha']) * (paramS['L_per_capita'] ** (1 - cS['alpha']))
    G_eq_final = cS['gov_exp_frac_Y'] * Y_eq_final
    B_eq_final = cS['gov_debt_frac_Y'] * Y_eq_final

    logging.info('\n--- V9 最终均衡汇总 ---')
    logging.info(f'K_eq (来自均衡求解器): {K_eq:.4f}, K_eq (来自最终模拟): {Actual_K_eq_final:.4f}')
    if abs(K_eq - Actual_K_eq_final) > 2e-2 and K_eq > 1e-9:
        logging.warning(f'K_eq from solver and K from final simulation differ significantly by {abs(K_eq - Actual_K_eq_final):.3e}.')
    
    logging.info(f'均衡总生产性资本 (K*): {Actual_K_eq_final:.4f} (总体人均)')
    logging.info(f'  其中: 非PPS资本 K_non_pps: {K_nonpps_eq_agg:.4f}, PPS资本 K_pps: {K_pps_eq_agg:.4f}')
    logging.info(f"均衡总劳动 (L, 效率单位, 总体人均): {paramS['L_per_capita']:.4f}")
    logging.info(f'均衡总产出 (Y*): {Y_eq_final:.4f}')
    logging.info(f'均衡市场毛回报率因子 (R_mkt_gross*): {R_mkt_gross_factor_eq_final:.4f} (对应 r_mkt_gross*={r_mkt_gross_eq_final:.4f})')
    logging.info(f'  家庭税后资本净回报率因子 (R_k_net_hh*): {R_k_net_factor_hh_eq_final:.4f} (对应 r_k_net_hh*={r_k_net_hh_eq_final:.4f})')
    logging.info(f'均衡市场总工资率 (MPL_gross*): {MPL_gross_eq_final:.4f}')
    logging.info(f"目标PAYG替代率 (rho_prime_payg_fixed*): {cS['rho_prime_payg_fixed']:.4f} (b_payg / avg_worker_gross_wage)")
    logging.info(f"均衡内生实际PAYG税率 (theta_payg_eq*, 上限 {cS['theta_payg_max']:.3f}): {paramS_eq['theta_payg_actual_for_hh']:.4f}")
    if 'theta_payg_required_before_cap' in final_eq_details:
        logging.info(f"  (理论所需PAYG税率，未考虑上限前: {final_eq_details['theta_payg_required_before_cap']:.4f})")
    logging.info(f'均衡内生"所得"税率 (tau_l_eq*): {tau_l_eq:.4f}')
    logging.info(f"  固定资本所得税率 (tau_k): {cS['tau_k']:.2f}, 固定消费税率 (tau_c): {cS['tau_c']:.2f}")

    w_net_hh_approx_display = MPL_gross_eq_final * (1 - paramS_eq['theta_payg_actual_for_hh'] - paramS_eq['tau_l'])
    logging.info(f'均衡近似家庭平均净工资率 (w_gross * (1-theta_act-tau_l)): {w_net_hh_approx_display:.4f} (仅为说明)')

    logging.info(f'均衡PAYG福利 (b_payg*, 每位退休者): {b_payg_eq_final:.4f}')
    logging.info(f'均衡总净转移支付 (TR_total*, 总体人均, TR_gov=0): {TR_total_eq_final_for_vfi:.4f}')
    logging.info(f'  其中意外遗赠 (T_bequest*): {T_bequest_eq_final:.4f}')
    logging.info('  其中政府直接转移 (TR_gov*): 0.0000 (按设定)')
    logging.info(f'均衡政府消费 (G*): {G_eq_final:.4f} (G/Y* = {G_eq_final/Y_eq_final if Y_eq_final > 0 else 0:.3f})')
    logging.info(f'均衡政府债务 (B*): {B_eq_final:.4f} (B/Y* = {B_eq_final/Y_eq_final if Y_eq_final > 0 else 0:.3f})')
    logging.info(f'均衡 K*/Y* 比率: {Actual_K_eq_final / Y_eq_final if Y_eq_final > 0 else 0:.4f}')
    logging.info(f'均衡 C*/Y* 比率: {C_eq_final / Y_eq_final if Y_eq_final > 0 else 0:.4f}')

    achieved_replacement_rate_final = b_payg_eq_final / avg_worker_gross_wage_eq_final if avg_worker_gross_wage_eq_final > 1e-9 else 0
    logging.info(f'实际达成替代率 (b_payg / avg_worker_gross_wage): {achieved_replacement_rate_final:.4f} (应接近 rho_prime_payg_fixed*)')
    if abs(achieved_replacement_rate_final - cS['rho_prime_payg_fixed']) > 1e-3:
        logging.warning(f"最终达成的替代率与目标替代率差异较大 (差异: {abs(achieved_replacement_rate_final - cS['rho_prime_payg_fixed']):.3e})。")

    final_gbc_residual_check = OLGV9Utils.check_gbc_residual(Actual_K_eq_final, C_eq_final, Y_eq_final, 
        G_eq_final, B_eq_final, MPL_gross_eq_final, r_mkt_gross_eq_final, 
        paramS_eq['theta_payg_actual_for_hh'], tau_l_eq, b_payg_eq_final, 
        T_bequest_eq_final, 0, cS, paramS_eq)
    logging.info(f'最终GBC(一般预算)检查 @ 均衡状态: Residual = {final_gbc_residual_check:.4e}')
    if abs(final_gbc_residual_check) > 1e-2:
        logging.warning(f'最终GBC残差 ({final_gbc_residual_check:.3e}) 较大。')

    # Plotting
    logging.info("\n绘制最终均衡的策略函数...")
    plot_a_idx = cS['aR_new'] // 2
    plot_ie_idx = cS['nw'] // 2
    plot_ikpps_indices = np.round(np.linspace(0, cS['nkpps'] - 1, min(3, cS['nkpps']))).astype(int)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f"Equilibrium Policy Functions (Age Group {plot_a_idx}, Effic State {plot_ie_idx})")

    # Savings Policy k'
    ax = axes[0]
    for ikpps in plot_ikpps_indices:
        ax.plot(cS['kGridV'], kPolM_plot[:, ikpps, plot_ie_idx, plot_a_idx], label=f"k_pps={cS['kppsGridV'][ikpps]:.2f}")
    ax.plot(cS['kGridV'], cS['kGridV'], 'k--', label='k\'=k')
    ax.set_title("Savings Policy k'(k | k_pps)")
    ax.set_xlabel("Current Assets (k)")
    ax.set_ylabel("Next Period Assets (k')")
    ax.legend()
    ax.grid(True)

    # PPS Contribution Policy
    ax = axes[1]
    plot_a_idx_contrib = min(5, cS['aR_new']-1)
    for ikpps in plot_ikpps_indices:
        ax.plot(cS['kGridV'], cPpsPolM_plot[:, ikpps, plot_ie_idx, plot_a_idx_contrib], label=f"k_pps={cS['kppsGridV'][ikpps]:.2f}")
    ax.set_title(f"PPS Contribution (Age Group {plot_a_idx_contrib})")
    ax.set_xlabel("Current Assets (k)")
    ax.set_ylabel("PPS Contribution (c_pps)")
    ax.legend()
    ax.grid(True)

    # Consumption Policy
    ax = axes[2]
    for ikpps in plot_ikpps_indices:
        ax.plot(cS['kGridV'], cPolM_plot[:, ikpps, plot_ie_idx, plot_a_idx], label=f"k_pps={cS['kppsGridV'][ikpps]:.2f}")
    ax.set_title("Consumption Policy c(k | k_pps)")
    ax.set_xlabel("Current Assets (k)")
    ax.set_ylabel("Consumption (c)")
    ax.legend()
    ax.grid(True)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    # save figure
    plt.savefig('main_olg_v9_results.png', dpi=300)
    logging.info("Policy function plots saved to 'main_olg_v9_results.png'")
    plt.show()


if __name__ == '__main__':
    main() 