# --- START OF FILE main_olg_v8.py ---

# OLG 模型 V8 (内生PPS缴费决策, PPS所得税递延, VFI w k_pps状态):
# 目标: 求解给定 PAYG 替代率 (rho_prime_payg_fixed) 下的均衡
# PPS缴费: 个体优化选择PPS缴费额，但受收入比例上限和年度绝对上限约束。
# 其他特性同Baseline:
#   - PPS 缴费可从所得税前扣除 (所得税率为 tau_l)。
#   - tau_l 内生调整以平衡政府一般预算 (TR_gov = 0)。
#   - PAYG 税率 (theta_payg) 内生决定，但有上限 cS.theta_payg_max。
#   - VFI 状态变量仍然包含 k_pps (PPS资产)。
# 🤖 新增功能: 支持使用RL模型替代VFI进行内层求解

import numpy as np
import matplotlib.pyplot as plt
import time
import os
# 解决OpenMP库冲突问题
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
import pickle
from main_olg_v9_utils import OLG_V9_Utils, ModelParameters

def main():
    plt.close('all')
    print('=== OLG 模型 V9 (内生PPS缴费, 固定 Rho_prime_payg, VFI/RL w k_pps) ===')
    print('    (Rho_prime_payg 固定, TR_gov=0, tau_l 内生, theta_payg 有上限)')
    print('    (VFI 状态: k, k_pps, eps; PPS缴费: 内生选择，有比例和绝对上限)')
    print('🤖 新增功能: 支持使用RL模型替代VFI进行内层求解')

    # --- 🤖 选择求解方法 ---
    print('\n🤖 选择求解方法:')
    print('  1. VFI (传统值函数迭代)')
    print('  2. RL (使用训练好的强化学习模型)')
    
    use_rl_solver = False
    rl_model = None
    rl_config = None
    
    # 自动检测最佳RL模型路径
    best_model_path = './py/best_model_sbx_full/best_model.zip'
    if not os.path.exists(best_model_path):
        best_model_path = './py/final_sac_agent_olg_sbx_full.zip'

    if os.path.exists(best_model_path):
        print(f'  🎯 检测到可用的RL模型: {best_model_path}')
        try:
            user_choice = input('  是否使用RL模型进行求解? (y/n, 默认n): ').strip().lower()
            if user_choice in ['y', 'yes', '是']:
                use_rl_solver = True
                print(f'  ✅ 选择使用RL模型: {best_model_path}')
                
                print('  🔄 加载SBX RL模型和配置...')
                try:
                    from sbx import SAC as SBX_SAC
                    rl_model = SBX_SAC.load(best_model_path)
                    
                    config_path = best_model_path.replace('.zip', '_config.pkl')
                    with open(config_path, 'rb') as f:
                        rl_config = pickle.load(f)
                    
                    print('  ✅ RL模型和配置加载成功。')
                except Exception as e:
                    print(f'  ❌ RL模型加载失败: {e}')
                    print('  🔄 回退到VFI求解')
                    use_rl_solver = False
            else:
                print('  ✅ 选择使用传统VFI求解')
        except (EOFError, KeyboardInterrupt):
            print('\n  ✅ 默认使用传统VFI求解')
    else:
        print('  ⚠️ 未检测到可用的RL模型，将使用传统VFI求解')

    # --- 1. 初始化参数 ---
    print('\n--- 1. 初始化参数 ---')
    cS = OLG_V9_Utils.ParameterValues_HuggettStyle()
    paramS = ModelParameters()

    cS.rho_prime_payg_fixed = 0.2
    print(f'>>> V9: 固定 PAYG 替代率 (rho_prime_payg_fixed): {cS.rho_prime_payg_fixed:.3f}')
    
    # --- 2. 模拟人口动态 ---
    print('\n--- 2. 模拟人口动态 ---')
    popS = OLG_V9_Utils.initPopulation(cS)
    popS = OLG_V9_Utils.populationDynamics(popS, cS)
    Z_ss, _, _, _ = OLG_V9_Utils.detectSteadyStatePopulation(popS, cS)
    paramS.Z_ss_counts = Z_ss
    Z_ss_total = np.sum(Z_ss)
    paramS.ageMassV = Z_ss / Z_ss_total if Z_ss_total > 1e-9 else np.zeros(cS.aD_new)
    paramS.mass_workers_group = np.sum(paramS.ageMassV[:cS.aR_new])
    
    Z_ss_norm_annual = np.zeros(cS.aD_orig)
    if Z_ss_total > 1e-9:
        for a_new, indices in enumerate(cS.physAgeMap):
            if indices:
                Z_ss_norm_annual[indices] = paramS.ageMassV[a_new] / len(indices)
    paramS.Z_ss_norm_annual = Z_ss_norm_annual / np.sum(Z_ss_norm_annual)
    paramS.popGrowthForDebt = (popS.totalPop[-1] / popS.totalPop[-2])**(1/cS.yearStep) - 1 if len(popS.totalPop) > 1 and popS.totalPop[-2]>0 else 0.01

    # --- 3. 预计算劳动 ---
    print('\n--- 3. 预计算劳动 ---')
    paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V = OLG_V9_Utils.EarningProcess_olgm(cS)
    paramS.leGridV = np.exp(paramS.leLogGridV.flatten())
    paramS.ageEffV_new = cS.ageEffV_new
    eIdxM = OLG_V9_Utils.LaborEndowSimulation_olgm_AgeGroup(cS, paramS)
    _, L_per_capita = OLG_V9_Utils.LaborSupply_Huggett(eIdxM, cS, paramS, paramS.ageMassV)
    paramS.L_per_capita = max(L_per_capita, 1e-6)

    # --- 4. 求解一般均衡 ---
    print(f'\n--- 4. 求解一般均衡 (固定 rho_prime_payg_fixed={cS.rho_prime_payg_fixed:.3f}) ---')
    K_global_guess = 2.0
    
    solve_start_time = time.time()
    if use_rl_solver and rl_model is not None:
        print('🤖 调用RL模型均衡求解器 solve_K_tau_l_for_rho_prime_rl...')
        # 将必要参数放入rl_config, 如果它们不存在的话
        if 'cS' not in rl_config: rl_config['cS'] = cS
        if 'paramS_for_rl' not in rl_config: 
            rl_config['paramS_for_rl'] = {'leGridV': paramS.leGridV}
        
        K_eq, tau_l_eq, gbc_residual_eq, eq_found, final_eq_solution_details = \
            OLG_V9_Utils.solve_K_tau_l_for_rho_prime_rl(
                cS.rho_prime_payg_fixed, K_global_guess, cS, paramS, eIdxM, rl_model, rl_config
            )
        solver_method = 'RL_model'
    else:
        print('📊 调用传统VFI均衡求解器 solve_K_tau_l_for_rho_prime_vfi...')
        K_eq, tau_l_eq, gbc_residual_eq, eq_found, final_eq_solution_details = \
            OLG_V9_Utils.solve_K_tau_l_for_rho_prime_vfi(
                cS.rho_prime_payg_fixed, K_global_guess, cS, paramS, eIdxM
            )
        solver_method = 'VFI'
    solve_time = time.time() - solve_start_time

     # (此部分代码从 solver_method = ... 之后开始)
    
    print(f'均衡求解完成。耗时: {solve_time:.2f} 秒。')
    print(f'  求解方法: {solver_method}')
    print(f'  均衡结果: K_eq = {K_eq:.4f}, tau_l_eq = {tau_l_eq:.4f}, GBC 残差 = {gbc_residual_eq:.3e}')

    if not eq_found or np.isnan(K_eq) or np.isnan(tau_l_eq):
        raise ValueError(f'未能为固定的 rho_prime_payg_fixed = {cS.rho_prime_payg_fixed:.3f} 找到均衡解。')
    if abs(gbc_residual_eq) > cS.gbc_tol_for_internal_loop * 10:
        print(f'警告: 最终均衡的GBC残差 ({gbc_residual_eq:.2e}) 较大。')
    
    # --- 5. 分析和绘制最终均衡结果 ---
    print(f'\n--- 5. 最终均衡结果与绘图 (rho_prime_payg_fixed={cS.rho_prime_payg_fixed:.3f}, tau_l_eq={tau_l_eq:.4f}, TR_gov=0) ---')
    
    # a. 准备最终模拟所需的参数
    paramS_eq = paramS
    paramS_eq.tau_l = tau_l_eq
    paramS_eq.theta_payg_actual_for_hh = final_eq_solution_details.get('theta_payg', 0.0)
    paramS_eq.pps_tax_deferral_active = cS.pps_active
    
    R_mkt_gross_factor_eq, MPL_gross_eq = OLG_V9_Utils.HHPrices_Huggett(K_eq, paramS.L_per_capita, cS)
    r_mkt_gross_eq = R_mkt_gross_factor_eq - 1
    r_k_net_hh_eq = r_mkt_gross_eq * (1 - cS.tau_k)
    R_k_net_factor_hh_eq = 1 + r_k_net_hh_eq

    b_payg_eq = final_eq_solution_details.get('b_payg', 0.0)
    bV_eq = np.zeros(cS.aD_new)
    bV_eq[cS.aR_new:] = b_payg_eq

    TR_total_eq = 0.0 # 与基准RL设定对齐，无意外遗赠

    print(f'最终 {solver_method} 调用参数: MPL_gross={MPL_gross_eq:.4f}, tau_l={paramS_eq.tau_l:.4f}, theta_payg_actual={paramS_eq.theta_payg_actual_for_hh:.4f}, TR_total={TR_total_eq:.4f}')
    
    # b. 再次运行模拟以获取最终的生命周期路径
    if use_rl_solver:
        print('🤖 使用RL模型进行最终均衡分析...')
        kHistM_eq, kPpsHistM_eq, cHistM_eq, cppsHistM_eq = OLG_V9_Utils.HHSimulation_olgm_rl(
            rl_model, rl_config, eIdxM, R_k_net_factor_hh_eq, MPL_gross_eq, TR_total_eq, bV_eq, paramS_eq, cS
        )
    else:
        print('📊 调用最终的 HHSolution_VFI_Huggett (V9)...')
        cPolM_eq, kPolM_eq, cPpsPolM_choice_eq, _ = OLG_V9_Utils.HHSolution_VFI_Huggett(
            R_k_net_factor_hh_eq, MPL_gross_eq, TR_total_eq, bV_eq, paramS_eq, cS
        )

        print('模拟最终均衡的分布 (HHSimulation_olgm V9)...')
        kHistM_eq, kPpsHistM_eq, cHistM_eq, cppsHistM_eq = OLG_V9_Utils.HHSimulation_olgm(
            kPolM_eq, cPpsPolM_choice_eq, cPolM_eq, eIdxM,
            R_k_net_factor_hh_eq, MPL_gross_eq, TR_total_eq, bV_eq, paramS_eq, cS
        )

    # c. 计算最终均衡的宏观总量
    K_nonpps_eq_agg = np.dot(np.mean(kHistM_eq, axis=0), paramS.ageMassV)
    K_pps_eq_agg = 0
    if cS.pps_active and cS.pps_in_K and kPpsHistM_eq.size > 0:
        K_pps_eq_agg = np.dot(np.mean(kPpsHistM_eq, axis=0), paramS.ageMassV)
    
    Actual_K_eq_final = K_nonpps_eq_agg + K_pps_eq_agg
    C_eq_final = np.dot(np.mean(cHistM_eq, axis=0), paramS.ageMassV)
    Y_eq_final = cS.A * (Actual_K_eq_final**cS.alpha) * (paramS.L_per_capita**(1-cS.alpha))
    G_eq_final = cS.gov_exp_frac_Y * Y_eq_final
    B_eq_final = cS.gov_debt_frac_Y * Y_eq_final
    
    # d. 打印详细的均衡结果
    print(f'\n--- V9 最终均衡汇总 ({solver_method}) ---')
    print(f'K_eq (来自均衡求解器): {K_eq:.4f}, K_eq (来自最终模拟): {Actual_K_eq_final:.4f}')
    if abs(K_eq - Actual_K_eq_final) > 2e-2 and K_eq > 1e-9:
        print(f'警告: K_eq from solver and K from final simulation differ by {abs(K_eq - Actual_K_eq_final):.3e}.')

    print(f'均衡总生产性资本 (K*): {Actual_K_eq_final:.4f} (非PPS: {K_nonpps_eq_agg:.4f}, PPS: {K_pps_eq_agg:.4f})')
    print(f'均衡总劳动 (L*): {paramS.L_per_capita:.4f}')
    print(f'均衡总产出 (Y*): {Y_eq_final:.4f}')
    print(f'均衡市场毛回报率因子 (R_mkt_gross*): {R_mkt_gross_factor_eq:.4f} (对应 r_mkt_gross*={r_mkt_gross_eq:.4f})')
    print(f'  家庭税后资本净回报率因子 (R_k_net_hh*): {R_k_net_factor_hh_eq:.4f} (对应 r_k_net_hh*={r_k_net_hh_eq:.4f})')
    print(f'均衡市场总工资率 (MPL_gross*): {MPL_gross_eq:.4f}')
    print(f'目标PAYG替代率 (rho_prime_payg_fixed*): {cS.rho_prime_payg_fixed:.4f}')
    print(f'均衡内生实际PAYG税率 (theta_payg_eq*, 上限 {cS.theta_payg_max:.3f}): {paramS_eq.theta_payg_actual_for_hh:.4f}')
    if 'theta_payg_required_before_cap' in final_eq_solution_details:
        print(f'  (理论所需PAYG税率，未考虑上限前: {final_eq_solution_details["theta_payg_required_before_cap"]:.4f})')
    print(f'均衡内生"所得"税率 (tau_l_eq*): {tau_l_eq:.4f}')
    
    # e. 绘图 - 只有VFI方法才有策略函数可以绘制
    if not use_rl_solver:
        print('\n绘制最终均衡的策略函数...')
        plot_a_idx = min(round(cS.aR_new / 2), cS.aD_new - 1)
        plot_ie_idx = round(cS.nw / 2)
        
        plot_nkpps_to_show = min(3, cS.nkpps)
        if cS.nkpps > 0:
            plot_ikpps_indices = np.round(np.linspace(0, cS.nkpps - 1, plot_nkpps_to_show)).astype(int)
        else:
            plot_ikpps_indices = []

        figure_title_suffix = f'年龄组 {plot_a_idx} (约{int(cS.age1_orig + plot_a_idx * cS.yearStep)}岁), 效率状态 {plot_ie_idx}'
        
        if cS.nk > 1 and len(plot_ikpps_indices) > 0:
            # 非PPS储蓄策略 k'(k | k_pps)
            fig1, ax1 = plt.subplots(num=f'VFI: 非PPS储蓄策略 k\'(k | k_pps): {figure_title_suffix}')
            colors = plt.cm.viridis(np.linspace(0, 1, plot_nkpps_to_show))
            for i, ikpps in enumerate(plot_ikpps_indices):
                k_prime_slice = np.squeeze(kPolM_eq[:, ikpps, plot_ie_idx, plot_a_idx])
                ax1.plot(cS.kGridV, k_prime_slice, label=f'k_pps={cS.kppsGridV[ikpps]:.2f}', color=colors[i])
            ax1.plot(cS.kGridV, cS.kGridV, 'k--', label='k\'=k (45度线)')
            ax1.set_xlabel('当前非PPS资产 k'); ax1.set_ylabel('下一期非PPS资产 k\'')
            ax1.set_title(f'VFI: 非PPS储蓄策略 k\'(k | k_pps)\n{figure_title_suffix}')
            ax1.legend(loc='best'); ax1.grid(True)

            # PPS缴费策略 c_pps(k | k_pps)
            if cS.pps_active:
                plot_a_idx_pps = min(4, cS.aR_new - 1)
                figure_title_suffix_pps = f'年龄组 {plot_a_idx_pps} (约{int(cS.age1_orig + plot_a_idx_pps * cS.yearStep)}岁), 效率状态 {plot_ie_idx}'
                fig2, ax2 = plt.subplots(num=f'VFI: PPS缴费策略 c_pps(k | k_pps): {figure_title_suffix_pps}')
                for i, ikpps in enumerate(plot_ikpps_indices):
                    cpps_slice = np.squeeze(cPpsPolM_choice_eq[:, ikpps, plot_ie_idx, plot_a_idx_pps])
                    ax2.plot(cS.kGridV, cpps_slice, label=f'k_pps={cS.kppsGridV[ikpps]:.2f}', color=colors[i])
                ax2.set_xlabel('当前非PPS资产 k'); ax2.set_ylabel('PPS缴费 c_pps')
                ax2.set_title(f'VFI: PPS缴费策略\n{figure_title_suffix_pps}')
                ax2.legend(loc='best'); ax2.grid(True)

            # 消费策略 c(k | k_pps)
            fig3, ax3 = plt.subplots(num=f'VFI: 消费策略 c(k | k_pps): {figure_title_suffix}')
            for i, ikpps in enumerate(plot_ikpps_indices):
                c_slice = np.squeeze(cPolM_eq[:, ikpps, plot_ie_idx, plot_a_idx])
                ax3.plot(cS.kGridV, c_slice, label=f'k_pps={cS.kppsGridV[ikpps]:.2f}', color=colors[i])
            ax3.set_xlabel('当前非PPS资产 k'); ax3.set_ylabel('消费量 c')
            ax3.set_title(f'VFI: 消费策略 c(k | k_pps)\n{figure_title_suffix}')
            ax3.legend(loc='best'); ax3.grid(True)
    else:
        print('\n🤖 RL模型求解完成，无策略函数可视化（RL模型为黑盒）。')
        print('  💡 如需分析策略函数，请使用VFI求解方法。')
    
    print(f'\n--- V9 OLG 模型 (内生PPS缴费, 固定 Rho_prime_payg, {solver_method}) 分析完成 ---')
    if not use_rl_solver:
        plt.show()


if __name__ == "__main__":
    main()