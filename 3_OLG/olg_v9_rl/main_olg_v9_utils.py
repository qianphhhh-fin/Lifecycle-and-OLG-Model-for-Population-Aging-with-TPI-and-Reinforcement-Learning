# --- START OF FILE main_olg_v9_utils.py ---

# =====================================================================
# === OLG 模型 V8 工具函数库: 内生PPS缴费决策的数值实现 (Python版) ===
# =====================================================================
#
# 理论基础：
# 本工具函数库实现了v8.tex理论模型中描述的所有核心算法和数值方法
# 包括值函数迭代(VFI)、家庭问题求解、一般均衡算法等关键组件
#
# 主要功能模块：
# 1. 参数设定函数 - 校准模型所有参数
# 2. 人口动态函数 - 模拟人口老龄化过程
# 3. 劳动禀赋过程 - Tauchen方法离散化AR(1)过程
# 4. 家庭效用函数 - CRRA效用和边际效用计算
# 5. VFI核心算法 - 内生PPS缴费的值函数迭代求解
# 6. 宏观经济函数 - 要素价格和市场出清条件
# 7. 一般均衡求解器 - 迭代求解K和tau_l的均衡值
#
# V8模型核心创新：
# - HHSolutionByAge_VFI_Huggett_v8: 实现内生PPS缴费选择的VFI
# - 对每个状态(k,k_pps,ε,a)，个体选择最优的(c,k',c_pps)组合
# - PPS缴费受双重约束：收入比例上限和年度绝对上限
# - 通过连续优化(scipy.optimize.minimize)实现PPS缴费的优化选择
# =====================================================================

import numpy as np
import scipy.stats as stats
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import minimize, fminbound
import warnings
import time
import matplotlib.pyplot as plt
from typing import Dict, Any, Tuple, List, Optional, Union # <--- 添加或确保这一行存在
# Suppress warnings for cleaner output, similar to MATLAB's warning('off', ...)
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=FutureWarning)

class ModelParameters:
    """A simple class to hold model parameters, mimicking a MATLAB struct."""
    def __init__(self):
        pass

class OLG_V9_Utils:

    @staticmethod
    def ParameterValues_HuggettStyle(**kwargs):
        """
        ParameterValues_HuggettStyle - 设置OLG模型V8的所有参数
        对应v8.tex第4节"参数设定和参数校准"
        """
        print('V8: 开始设置参数...')
        cS = ModelParameters()

        # --- 人口结构基础参数（对应v8.tex人口动态设定） ---
        cS.age1_orig = 20              # 模型起始年龄（岁）
        cS.ageLast_orig = 98           # 模型终止年龄（岁）
        cS.ageRetire_orig = 65         # 退休年龄（岁）
        cS.popGrowth_orig = 0.012      # 原始人口增长率
        cS.aD_orig = cS.ageLast_orig - cS.age1_orig + 1        # 年度年龄数量 (20-98岁，共79年)
        cS.aR_idx_orig = cS.ageRetire_orig - cS.age1_orig + 1  # 🔧 修复：退休年度年龄索引 (1-based，与MATLAB一致)
        cS.aW_orig = cS.aR_idx_orig - 1                       # 🔧 修复：工作年数 (不包含退休年龄，与MATLAB一致)
        cS.physAgeV_orig = np.arange(cS.age1_orig, cS.ageLast_orig + 1)   # 年度年龄向量

        # --- 年度死亡率数据（基于中国生命表） ---
        cS.d_orig = np.array([0.00159,0.00169,0.00174,0.00172,0.00165,0.00156,0.00149,0.00145,0.00145,0.00149,0.00156,0.00163,0.00171,0.00181,0.00193,0.00207,0.00225,0.00246,0.00270,0.00299,0.00332,0.00368,0.00409,0.00454,0.00504,0.00558,0.00617,0.00686,0.00766,0.00865,0.00955,0.01058,0.01162,0.01264,0.01368,0.01475,0.01593,0.01730,0.01891,0.02074,0.02271,0.02476,0.02690,0.02912,0.03143,0.03389,0.03652,0.03930,0.04225,0.04538,0.04871,0.05230,0.05623,0.06060,0.06542,0.07066,0.07636,0.08271,0.08986,0.09788,0.10732,0.11799,0.12895,0.13920,0.14861,0.16039,0.17303,0.18665,0.20194,0.21877,0.23601,0.25289,0.26973,0.28612,0.30128,0.31416,0.32915,0.34450,0.36018])
        if len(cS.d_orig) != cS.aD_orig:
            raise ValueError('年度死亡率数据 d_orig 长度与年龄跨度不匹配')
        cS.s_orig = 1 - cS.d_orig      # 年度存活率

        # --- 年龄组聚合参数（将年度年龄聚合为5年期年龄组） ---
        cS.yearStep = 5                # 每个年龄组跨度（年）
        cS.aD_new = int(np.ceil(cS.aD_orig / cS.yearStep))    # 🔧 修复：年龄组数量 (与MATLAB一致，使用ceil)
        cS.aR_new = int(np.ceil(cS.aW_orig / cS.yearStep))    # 工作年龄组数量

        # 建立年度年龄到年龄组的映射关系 (0-based)
        cS.physAgeMap = [[] for _ in range(cS.aD_new)]
        for a in range(cS.aD_new):
            startIdx = a * cS.yearStep
            endIdx = min((a + 1) * cS.yearStep, cS.aD_orig)
            cS.physAgeMap[a] = list(range(startIdx, endIdx))

        # 计算各年龄组代表性年龄
        cS.physAgeV_new = np.zeros(cS.aD_new)
        for a in range(cS.aD_new):
            cS.physAgeV_new[a] = cS.physAgeV_orig[cS.physAgeMap[a][0]]

        # --- [核心修改] 计算年龄组间转移存活率 ---
        # 与MATLAB的简洁逻辑完全对齐
        cS.s_1yr_transitionV = np.zeros(cS.aD_new)
        for a in range(cS.aD_new - 1):
            lastYearIdxInGroup = cS.physAgeMap[a][-1]
            # s_orig的长度是 aD_orig，其有效索引是 0 到 aD_orig-1
            # lastYearIdxInGroup 是正确的0-based年度索引
            cS.s_1yr_transitionV[a] = cS.s_orig[lastYearIdxInGroup]
        # 最后一个年龄组的转移存活率默认为0，这在初始化时已经完成，无需额外操作。

        # --- 初始人口分布（基于2023年中国人口结构） ---
        cS.initial_pop = np.array([76.2, 86.4, 113.8, 98.6, 86.6, 102.7, 112.0, 99.0, 64.0, 66.9, 44.1, 25.4, 14.9, 6.8, 1.7, 0.2])
        if len(cS.initial_pop) != cS.aD_new:
            warnings.warn('initial_pop长度与年龄组数不匹配，已重设为均匀分布。')
            cS.initial_pop = np.ones(cS.aD_new) * (100 / cS.aD_new)

        # --- 年龄组间存活率（用于人口动态模拟） ---
        # 对应v8.tex公式中的β_{surv,a-1,t-1}参数
        # 🔧 修复：与MATLAB版本保持一致，确保长度为aD_new-1
        beta_surv_pop = [0.998, 0.996, 0.994, 0.992, 0.988, 0.984, 0.980, 0.976, 
                        0.970, 0.960, 0.945, 0.920, 0.880, 0.800, 0.680]  # 15个元素，对应16个年龄组
        if len(beta_surv_pop) != cS.aD_new - 1:
            raise ValueError(f'年龄组间存活率 beta_surv_pop 的长度对于 {cS.aD_new} 个年龄组不正确。应为 {cS.aD_new - 1}。')        
        cS.survivalProbV_popdyn = np.array(beta_surv_pop + [0])  # 最后一个年龄组存活率为0

        # --- 人口动态收敛参数 ---
        cS.bgp_tolerance = 0.001       # 稳态收敛容忍度
        cS.bgp_window = 5              # 稳态检测窗口期
        cS.max_periods = 50            # 最大模拟期数

        # --- 家庭偏好参数 ---
        cS.sigma = 1.5            # 相对风险厌恶系数γ
        cS.beta = 0.97            # 主观贴现因子β_disc
        cS.cFloor = 0.0001        # 最低消费约束
        cS.nSim = 1000            # 蒙特卡洛模拟个体数

        # --- 生产技术参数 ---
        cS.A = 0.895944       # 全要素生产率
        cS.alpha = 0.36           # 资本产出弹性
        cS.ddk = 0.06           # 资本折旧率δ

        # --- 政府财政参数 ---
        cS.tau_k = 0.20                # 资本所得税率
        cS.tau_c = 0.10                # 消费税率
        cS.gov_exp_frac_Y = 0.15       # 政府支出占GDP比例
        cS.gov_debt_frac_Y = 0.60      # 政府债务占GDP比例

        # --- 劳动效率冲击过程参数 ---
        cS.leSigma1 = 0.38**0.5         # 初期效率分布标准差
        cS.leShockStd = 0.045**0.5      # 效率冲击标准差
        cS.lePersistence = 0.96        # AR(1)持续性参数
        cS.leWidth = 4                 # Tauchen方法的标准差倍数
        cS.nw = 5                      # 效率状态网格点数

        # --- 资产网格参数 ---
        cS.tgKY = 3                    # 目标资本产出比
        cS.tgWage = (1-cS.alpha)*cS.A*((cS.tgKY/cS.A)**(cS.alpha/(1-cS.alpha)))
        cS.nk = 40                # 非PPS资产状态网格点数
        cS.nkpps = 40                 # PPS资产状态网格点数
        cS.nkprime = 40             # 非PPS资产决策网格点数
        cS.npps = 40            # PPS资产决策网格点数

            # 2. 使用 kwargs 覆盖默认值
        if 'nk' in kwargs: cS.nk = kwargs['nk']
        if 'nkpps' in kwargs: cS.nkpps = kwargs['nkpps']
        if 'nkprime' in kwargs: cS.nkprime = kwargs['nkprime']
        if 'npps' in kwargs: cS.npps = kwargs['npps']



        cS.kMin = 0                    # 非PPS资产下界
        cS.kMax = 15 * cS.tgWage       # 非PPS资产上界
        power = 1.5                    # 网格密度参数
        kGridV = cS.kMin + (cS.kMax - cS.kMin) * (np.linspace(0, 1, cS.nk)**power)
        if cS.nk > 0:
            kGridV[0] = cS.kMin
        cS.kGridV = kGridV

        # --- 年龄效率剖面 ---
        ageEffV_orig_temp = np.zeros(100)
        ageEffV_orig_temp[19:72] = np.concatenate([
            np.linspace(0.3, 1.5, 36 - 20 + 1),
            1.5 * np.ones(47 - 37 + 1),
            np.linspace(1.5, 0.2, 65 - 48 + 1),
            np.linspace(0.18, 0, 72 - 66 + 1)
        ])
        cS.ageEffV_orig = ageEffV_orig_temp[cS.age1_orig -1 : cS.ageLast_orig]
        if len(cS.ageEffV_orig) != cS.aD_orig:
            raise ValueError('ageEffV_orig 年度年龄效率剖面长度不匹配')
        
        cS.ageEffV_new = np.zeros(cS.aD_new)
        for a in range(cS.aD_new):
            cS.ageEffV_new[a] = np.mean(cS.ageEffV_orig[cS.physAgeMap[a]])

        # === V8模型核心：PPS制度参数设计 ===
        cS.use_continuous_optimization = True

        # --- PPS制度基础参数 ---
        cS.pps_active = True
        cS.pps_tax_rate_withdrawal = 0.03
        cS.pps_return_rate_premium = 0.00
        cS.pps_withdrawal_rate = 0.15
        cS.pps_in_K = True
        cS.pps_bequeathable = True

        # --- V8模型关键创新：PPS缴费约束参数 ---
        cS.pps_contrib_limit = 9999
        cS.pps_max_contrib_frac = 1
        cS.pps_contribution_age_max_idx = cS.aR_idx_orig - 1
        cS.pps_withdrawal_age_min_idx = cS.aR_idx_orig


        # --- PPS资产网格 ---
        
        cS.kppsMin = 0
        cS.kppsMax = cS.kMax / 2
        if cS.nkpps > 0:
            cS.kppsMax = max(cS.kppsMax, 1e-3)
        power_kpps = 1.5
        if cS.nkpps > 1:
            kppsGridV_temp = cS.kppsMin + (cS.kppsMax - cS.kppsMin) * (np.linspace(0, 1, cS.nkpps)**power_kpps)
            kppsGridV_temp[0] = cS.kppsMin
        elif cS.nkpps == 1:
            kppsGridV_temp = np.array([cS.kppsMin])
        else:
            kppsGridV_temp = np.array([])
        cS.kppsGridV = kppsGridV_temp

        # --- 一般均衡求解参数 ---
        cS.max_iter_K_tau_l = 100
        cS.tol_K_tau_l = 1e-4
        cS.damp_K_v5 = 0.1
        cS.damp_tau_l_v5 = 0.1
        cS.gbc_tol_for_internal_loop = 1e-3

        # --- 收敛检测参数 ---
        cS.max_stagnation_iters = 10
        cS.min_norm_improvement_frac = 1e-3
        cS.max_tau_l_boundary_strikes = 5

        # --- PAYG税率约束参数 ---
        cS.tau_l_init_guess = 0.1509
        cS.tau_l_min = 0.00
        cS.tau_l_max = 0.3
        cS.max_total_labor_tax = 1.0
        cS.theta_payg_max = 1.0

        # --- 参考性PPS缴费时间表（V8中不再直接使用） ---
        cS.pps_fixed_contrib_schedule_frac = np.zeros(cS.aD_new)
        num_working_age_groups = cS.aR_new
        if num_working_age_groups > 0:
            if num_working_age_groups == 1:
                cS.pps_fixed_contrib_schedule_frac[0] = 0.05
            elif num_working_age_groups > 1:
                mid_point1 = int(np.ceil(num_working_age_groups / 3))
                mid_point2 = int(np.ceil(2 * num_working_age_groups / 3))
                if mid_point1 > 0:
                    cS.pps_fixed_contrib_schedule_frac[0:mid_point1] = np.linspace(0.02, 0.06, mid_point1)
                if mid_point2 > mid_point1:
                    cS.pps_fixed_contrib_schedule_frac[mid_point1:mid_point2] = np.linspace(0.06, 0.10, mid_point2 - mid_point1)
                if num_working_age_groups > mid_point2:
                    cS.pps_fixed_contrib_schedule_frac[mid_point2:num_working_age_groups] = np.linspace(0.10, 0.04, num_working_age_groups - mid_point2)
        if cS.aR_new < cS.aD_new:
            cS.pps_fixed_contrib_schedule_frac[cS.aR_new:] = 0

        # --- 数值优化参数 ---
        # Note: In Python, options are passed directly to the optimizer function.
        # cS.fminbnd_TolX = 1e-6
        # cS.fminbnd_Display = 'none'


        
        print('V8: 完整参数已设置完毕。')
        print(f'    PPS约束：收入比例上限={cS.pps_max_contrib_frac*100:.1f}%, 年度绝对上限={cS.pps_contrib_limit:.2f}')
        print(f'    PPS可遗赠 (cS.pps_bequeathable): {cS.pps_bequeathable}')
        return cS

# 在 main_olg_v9_utils.py 的 OLG_V9_Utils class 内部

    @staticmethod
    def generateGrids(cS):
        """
        [新增, 与MATLAB对齐] 根据当前的网格参数设置，重新生成资产网格。
        """
        # 重新生成非PPS资产网格 (kGridV)
        power_k = 1.5
        if cS.nk > 1:
            kGridV_temp = cS.kMin + (cS.kMax - cS.kMin) * (np.linspace(0, 1, cS.nk)**power_k)
            kGridV_temp[0] = cS.kMin
        elif cS.nk == 1:
            kGridV_temp = np.array([cS.kMin])
        else:
            kGridV_temp = np.array([])
        cS.kGridV = kGridV_temp.reshape(-1, 1) # 确保是列向量

        # 重新生成PPS资产网格 (kppsGridV)
        power_kpps = 1.5
        if cS.nkpps > 1:
            kppsGridV_temp = cS.kppsMin + (cS.kppsMax - cS.kppsMin) * (np.linspace(0, 1, cS.nkpps)**power_kpps)
            kppsGridV_temp[0] = cS.kppsMin
        elif cS.nkpps == 1:
            kppsGridV_temp = np.array([cS.kppsMin])
        else:
            kppsGridV_temp = np.array([])
        cS.kppsGridV = kppsGridV_temp.reshape(-1, 1) # 确保是列向量
        
        # print(f"Python端网格已重新生成：nk={cS.nk}, nkpps={cS.nkpps}")
        return cS

    @staticmethod
    def initPopulation(cS):
        """
        initPopulation - 初始化人口结构
        """
        popS = ModelParameters()
        initial_total = np.sum(cS.initial_pop)

        # 🔧 修复：cS.initial_pop是年龄组数据，应与cS.aD_new比较，不是cS.aD_orig
        if initial_total > 0 and len(cS.initial_pop) == cS.aD_new:
            popS.Z = (cS.initial_pop / initial_total * 100).reshape(-1, 1)
        else:
            warnings.warn('初始人口数据不匹配或总和为零。将设置为均匀的初始年龄组人口分布。')
            popS.Z = np.full((cS.aD_new, 1), 100 / cS.aD_new)

        popS.totalPop = np.array([np.sum(popS.Z[:, 0])])

        if popS.totalPop[0] > 1e-9:
            popS.ageDist = popS.Z[:, 0] / popS.totalPop[0]
        else:
            popS.ageDist = np.zeros(cS.aD_new)
        popS.initialAgeDist = popS.ageDist.reshape(-1, 1)
        print(f'初始年龄组人口已设置。总人口={popS.totalPop[0]:.2f} (代表百分比基数)。')
        return popS

    @staticmethod
    def populationDynamics(popS, cS):
        """
        populationDynamics - 模拟人口动态演进
        """
        max_periods_sim = cS.max_periods
        Z_history = np.zeros((cS.aD_new, max_periods_sim + 1))
        totalPop_history = np.zeros(max_periods_sim + 1)
        ageDist_history = np.zeros((cS.aD_new, max_periods_sim + 1))


        # 已经是年龄组数据
        Z_history[:, 0] = popS.Z.flatten()
        totalPop_history[0] = popS.totalPop[0]
        ageDist_history[:, 0] = popS.ageDist.flatten()
        
        print(f'人口动态模拟开始 (年龄组, 最大期数 = {max_periods_sim})...')
        bgp_reached_flag = False
        actual_periods_run = max_periods_sim

        for t in range(max_periods_sim):
            if t % 10 == 0 or t == 0:
                print(f'  模拟人口期数 {t+1} (年龄组)')

            Z_current_period = Z_history[:, t]
            Z_next_period = np.zeros(cS.aD_new)

            if t < 5:
                time_varying_growth_rate = -0.01 - 0.003 * (t+1)
            else:
                time_varying_growth_rate = -0.03 - 0.004 * min(t - 5, 10)
            
            Z_next_period[0] = Z_current_period[0] * (1 + time_varying_growth_rate)
            Z_next_period[0] = max(0, Z_next_period[0])

            for a in range(1, cS.aD_new):
                survival_prob_group = cS.survivalProbV_popdyn[a - 1]
                Z_next_period[a] = Z_current_period[a - 1] * survival_prob_group
                Z_next_period[a] = max(0, Z_next_period[a])

            Z_history[:, t + 1] = Z_next_period
            totalPop_history[t + 1] = np.sum(Z_next_period)
            if totalPop_history[t + 1] > 1e-9:
                ageDist_history[:, t + 1] = Z_next_period / totalPop_history[t + 1]
            else:
                ageDist_history[:, t + 1] = 0
                totalPop_history[t + 1] = 0

            current_check_period_idx = t + 1
            if current_check_period_idx >= cS.bgp_window:
                stable = True
                for w_idx in range(cS.bgp_window):
                    hist_idx1 = current_check_period_idx - w_idx
                    hist_idx2 = current_check_period_idx - w_idx - 1
                    change = np.linalg.norm(ageDist_history[:, hist_idx1] - ageDist_history[:, hist_idx2])
                    if change >= cS.bgp_tolerance:
                        stable = False
                        break
                if stable:
                    print(f'\n人口稳态 (年龄组) 在模拟期数 {t+1} (对应历史数据索引 {current_check_period_idx+1}) 达到。')
                    bgp_reached_flag = True
                    actual_periods_run = t + 1
                    break
        
        final_period_idx_to_store = actual_periods_run + 1
        popS.Z = Z_history[:, :final_period_idx_to_store]
        popS.totalPop = totalPop_history[:final_period_idx_to_store]
        popS.ageDist = ageDist_history[:, :final_period_idx_to_store]
        
        depRatio_history = np.zeros(actual_periods_run)
        for th_loop in range(actual_periods_run):
            Z_t_for_depratio = Z_history[:, th_loop + 1]
            working_pop = np.sum(Z_t_for_depratio[:cS.aR_new])
            retired_pop = np.sum(Z_t_for_depratio[cS.aR_new:])
            if working_pop > 1e-9:
                depRatio_history[th_loop] = retired_pop / working_pop
            else:
                depRatio_history[th_loop] = np.inf
        popS.dependencyRatio = depRatio_history
        
        print(f'人口动态模拟完成。运行期数: {actual_periods_run}。达到BGP: {bgp_reached_flag}')
        if not bgp_reached_flag:
            print(f'警告: 人口稳态 (年龄组) 未在 {max_periods_sim} 期内达到。')
        
        return popS

    @staticmethod
    def detectSteadyStatePopulation(popS, cS):
        actual_periods_in_data = popS.Z.shape[1]
        bgp_reached = False
        bgp_period = actual_periods_in_data - 1 # 0-based

        if actual_periods_in_data < cS.bgp_window + 1:
            print(f'人口模拟期数过短 ({actual_periods_in_data} 数据点)，无法进行稳态检查 (窗口期 = {cS.bgp_window})。')
        else:
            print(f'检查人口稳态 (年龄组, 最近 {cS.bgp_window} 期)...')
            for t_check_end_idx in range(actual_periods_in_data - 1, cS.bgp_window - 1, -1):
                stable = True
                for w_idx in range(cS.bgp_window):
                    idx1 = t_check_end_idx - w_idx
                    idx2 = t_check_end_idx - w_idx - 1
                    change = np.linalg.norm(popS.ageDist[:, idx1] - popS.ageDist[:, idx2])
                    if change >= cS.bgp_tolerance:
                        stable = False
                        break
                if stable:
                    bgp_reached = True
                    bgp_period = t_check_end_idx -1
                    print(f'人口稳态 (年龄组) 从模拟期数 {bgp_period+1} (数据索引 {t_check_end_idx+1}) 开始检测到 (稳定窗口结束于此)。')
                    break
            if not bgp_reached:
                print('未检测到人口稳态 (年龄组)。将使用最终期数据。')
                bgp_period = actual_periods_in_data - 2

        ss_data_index = min(bgp_period + 1, popS.Z.shape[1]-1)
        Z_ss = popS.Z[:, ss_data_index]

        Z_ss_norm = np.zeros(cS.aD_new)
        if np.sum(Z_ss) > 1e-9:
            Z_ss_norm = Z_ss / np.sum(Z_ss)
        
        dependency_ratio_ss = np.nan
        if hasattr(popS, 'dependencyRatio') and len(popS.dependencyRatio) > 0:
            valid_dep_ratio_index = min(max(0, bgp_period), len(popS.dependencyRatio) - 1)
            dependency_ratio_ss = popS.dependencyRatio[valid_dep_ratio_index]
        if np.isnan(dependency_ratio_ss):
            working_pop_ss = np.sum(Z_ss[:cS.aR_new])
            retired_pop_ss = np.sum(Z_ss[cS.aR_new:])
            if working_pop_ss > 1e-9:
                dependency_ratio_ss = retired_pop_ss / working_pop_ss
            else:
                dependency_ratio_ss = np.inf
            if not hasattr(popS, 'dependencyRatio') or len(popS.dependencyRatio) == 0:
                warnings.warn('抚养比历史未找到或过短，已基于Z_ss重新计算。')
        
        plt.figure(num='V8: 初始 vs 稳态/最终 年龄组人口分布')
        group_indices = np.arange(cS.aD_new)
        width = 0.4
        
        # 🔧 修复：现在initialAgeDist已经是年龄组数据，直接使用
        if hasattr(popS, 'initialAgeDist') and popS.initialAgeDist is not None:
            initial_dist_flat = popS.initialAgeDist.flatten()
            if len(initial_dist_flat) == cS.aD_new:
                # 已经是年龄组数据，直接使用
                plt.bar(group_indices - width/2, initial_dist_flat * 100, width, label='初始年龄组分布')
            else:
                warnings.warn(f'初始年龄分布长度 ({len(initial_dist_flat)}) 与年龄组数 ({cS.aD_new}) 不匹配。')
        
        plt.bar(group_indices + width/2, Z_ss_norm * 100, width, label=f'稳态年龄组分布 (模拟期 {bgp_period+1})', color='r')
        plt.xlabel(f'年龄组索引 (0 至 {cS.aD_new-1})')
        plt.ylabel('占总人口百分比 (%)')
        plt.title(f'V8: 初始 vs 稳态/最终 年龄组人口分布 (稳态代表模拟期 t={bgp_period+1})')
        plt.legend(loc='best')
        plt.xticks(group_indices)
        plt.grid(True)
        plt.draw()
        plt.pause(0.01) # To ensure the plot is drawn
        print('已绘制初始与稳态/最终年龄组人口分布图。')
        
        return Z_ss, dependency_ratio_ss, bgp_reached, bgp_period

    @staticmethod
    def tauchen(N_states, persistence_rho, shock_sigma, mean_val_mu, num_std_dev_width):
        std_y_unconditional = np.sqrt(shock_sigma**2 / (1 - persistence_rho**2))
        
        y_max_boundary = num_std_dev_width * std_y_unconditional
        y_grid_centered = np.linspace(-y_max_boundary, y_max_boundary, N_states)

        step_size_d = 0
        if N_states > 1:
            step_size_d = y_grid_centered[1] - y_grid_centered[0]

        trProbM_calc = np.zeros((N_states, N_states))
        if N_states == 1:
            trProbM_calc[0,0] = 1.0
        else:
            for iRow in range(N_states):
                mean_next_y_conditional = persistence_rho * y_grid_centered[iRow]
                # First column
                trProbM_calc[iRow, 0] = stats.norm.cdf((y_grid_centered[0] - mean_next_y_conditional + step_size_d/2) / shock_sigma)
                # Last column
                trProbM_calc[iRow, N_states-1] = 1 - stats.norm.cdf((y_grid_centered[N_states-1] - mean_next_y_conditional - step_size_d/2) / shock_sigma)
                # Middle columns
                for iCol in range(1, N_states - 1):
                    upper_bound_cdf = stats.norm.cdf((y_grid_centered[iCol] - mean_next_y_conditional + step_size_d/2) / shock_sigma)
                    lower_bound_cdf = stats.norm.cdf((y_grid_centered[iCol] - mean_next_y_conditional - step_size_d/2) / shock_sigma)
                    trProbM_calc[iRow, iCol] = upper_bound_cdf - lower_bound_cdf
        
        row_sums_check = np.sum(trProbM_calc, axis=1)
        row_sums_check[row_sums_check <= 1e-9] = 1
        trProbM_out = trProbM_calc / row_sums_check[:, np.newaxis]

        unconditional_mean_shift = mean_val_mu / (1 - persistence_rho)
        y_grid_out = y_grid_centered + unconditional_mean_shift
        return y_grid_out.flatten(), trProbM_out

    @staticmethod
    def EarningProcess_olgm(cS):
        """
        EarningProcess_olgm - 生成劳动效率冲击的离散Markov链
        """
        print('劳动禀赋过程参数已生成 (Tauchen)。')
        lePersistence = 0.90
        leShockStd = 0.15
        Tauchen_q = 2.0
        
        leLogGridV_raw, leTrProbM = OLG_V9_Utils.tauchen(cS.nw, lePersistence, leShockStd, 0, Tauchen_q)
        leLogGridV = leLogGridV_raw - np.mean(leLogGridV_raw)

        leGridV_test = np.exp(leLogGridV)
        efficiency_ratio = leGridV_test[-1] / leGridV_test[0]

        max_acceptable_ratio = 5.0
        if efficiency_ratio > max_acceptable_ratio:
            compression_factor = np.log(max_acceptable_ratio) / np.log(efficiency_ratio)
            leLogGridV = leLogGridV * compression_factor
            print(f'效率分布压缩因子: {compression_factor:.3f}')

        # 计算平稳分布
        try:
            eigenvals, eigenvecs = np.linalg.eig(leTrProbM.T)
            idx = np.argmin(np.abs(eigenvals - 1))
            leProb1V = np.real(eigenvecs[:, idx])
            leProb1V = leProb1V / np.sum(leProb1V)
            leProb1V = np.abs(leProb1V)
        except:
            leProb1V = np.ones(cS.nw) / cS.nw
            for _ in range(1000):
                leProb1V_new = leTrProbM.T @ leProb1V
                if np.linalg.norm(leProb1V_new - leProb1V) < 1e-10:
                    break
                leProb1V = leProb1V_new
        
        if np.sum(leProb1V) > 1e-9:
            leProb1V /= np.sum(leProb1V)
        else:
            warnings.warn('EarningProcess_olgm: 平稳分布概率和过小，已重置为均匀分布。')
            leProb1V = np.ones(cS.nw) / cS.nw

        leGridV_final = np.exp(leLogGridV)
        efficiency_ratio_final = leGridV_final[-1] / leGridV_final[0]
        mean_efficiency = np.dot(leGridV_final, leProb1V)

        print(f'AR(1)参数: ρ={lePersistence:.3f}, σ={leShockStd:.3f}, q={Tauchen_q:.1f}')
        print(f'效率网格范围: [{leGridV_final[0]:.4f}, {leGridV_final[-1]:.4f}], 比值={efficiency_ratio_final:.2f}')
        print(f'平均效率: {mean_efficiency:.4f}')
        return leLogGridV, leTrProbM, leProb1V


    @staticmethod
    def MarkovChainSimulation_AgeGroup(num_simulations, cS, initial_prob_dist_p0V, transition_matrix_trProbM):
        """
        MarkovChainSimulation_AgeGroup - 直接模拟年龄组的马尔可夫链
        
        Args:
            num_simulations: 模拟个体数量
            cS: 包含年龄组信息的参数对象
            initial_prob_dist_p0V: 初始状态分布
            transition_matrix_trProbM: 转移概率矩阵
            
        Returns:
            eIdxM_group_out: 年龄组效率冲击索引矩阵 (nSim × aD_new)
        """
        # 设置随机数种子确保可重复性
        np.random.seed(433)
        
        num_states = len(initial_prob_dist_p0V)
        num_age_groups = cS.aD_new
        
        # 参数验证
        if transition_matrix_trProbM.shape[0] != num_states or transition_matrix_trProbM.shape[1] != num_states:
            raise ValueError('MarkovChainSimulation_AgeGroup: 转移矩阵维度与初始分布长度不匹配。')
        
        # 归一化概率分布
        if abs(np.sum(initial_prob_dist_p0V) - 1) > 1e-5:
            initial_prob_dist_p0V = initial_prob_dist_p0V / np.sum(initial_prob_dist_p0V)
        
        row_sums = np.sum(transition_matrix_trProbM, axis=1)
        if np.any(np.abs(row_sums - 1) > 1e-5):
            row_sums[row_sums <= 1e-9] = 1
            transition_matrix_trProbM = transition_matrix_trProbM / row_sums[:, np.newaxis]
        
        # 生成年龄组随机数
        random_numbers_group = np.random.rand(num_simulations, num_age_groups)
        
        # 计算累积概率
        cumulative_initial_prob_cP0 = np.cumsum(initial_prob_dist_p0V)
        cumulative_transition_prob_cPT = np.cumsum(transition_matrix_trProbM, axis=1)
        if num_states > 0:
            cumulative_transition_prob_cPT[:, -1] = 1.0
        
        # 初始化结果矩阵
        eIdxM_group_out = np.zeros((num_simulations, num_age_groups), dtype=np.uint16)
        
        # 第一个年龄组：使用初始分布
        if num_simulations > 0 and num_age_groups > 0 and num_states > 0:
            initial_indices = np.searchsorted(cumulative_initial_prob_cP0, random_numbers_group[:, 0])
            # 🔧 确保索引在有效范围内
            initial_indices = np.clip(initial_indices, 0, num_states - 1)
            eIdxM_group_out[:, 0] = initial_indices
        
        # 后续年龄组：使用转移概率
        for a_group in range(1, num_age_groups):
            current_state_indices = eIdxM_group_out[:, a_group - 1]
            
            # 验证状态索引有效性
            valid_indices = (current_state_indices >= 0) & (current_state_indices < num_states)
            if not np.all(valid_indices):
                print(f'警告: MarkovChainSimulation_AgeGroup: 在年龄组 {a_group} 检测到无效状态索引。已重置为状态0。')
                current_state_indices[~valid_indices] = 0
                eIdxM_group_out[:, a_group - 1] = current_state_indices
            
            # 计算下一期状态
            cPt_for_next_state = cumulative_transition_prob_cPT[current_state_indices, :]
            for i_sim in range(num_simulations):
                next_state = np.searchsorted(
                    cPt_for_next_state[i_sim], random_numbers_group[i_sim, a_group]
                )
                # 🔧 确保索引在有效范围内
                next_state = np.clip(next_state, 0, num_states - 1)
                eIdxM_group_out[i_sim, a_group] = next_state
        
        print(f'年龄组马尔可夫链模拟完成 ({num_simulations} 个体, {num_age_groups} 年龄组)。')
        return eIdxM_group_out

    @staticmethod
    def LaborEndowSimulation_olgm_AgeGroup(cS, paramS):
        """
        LaborEndowSimulation_olgm_AgeGroup - 直接生成年龄组的劳动禀赋路径
        
        Args:
            cS: 包含年龄组参数的对象
            paramS: 包含马尔可夫链参数的对象
            
        Returns:
            eIdxM_group: 年龄组效率冲击索引矩阵 (nSim × aD_new)
            
        功能：直接为年龄组生成效率冲击序列，避免年度到年龄组的转换
        """
        eIdxM_group = OLG_V9_Utils.MarkovChainSimulation_AgeGroup(cS.nSim, cS,
                                                                 paramS.leProb1V, paramS.leTrProbM)
        print(f'年龄组劳动禀赋路径已生成 ({cS.nSim} 个体, {cS.aD_new} 年龄组)。')
        return eIdxM_group

    @staticmethod
    def LaborSupply_Huggett(eIdxM_input, cS, paramS, Z_ss_norm_group):
        """
        LaborSupply_Huggett - 计算劳动供给（适配年龄组模拟）
        
        Args:
            eIdxM_input: 效率冲击矩阵，可以是年度 (n_sim x aD_orig) 或年龄组 (n_sim x aD_new)
            cS: 模型参数
            paramS: 参数结构体
            Z_ss_norm_group: 年龄组人口分布
        """
        nSim = eIdxM_input.shape[0]
        

            # 年龄组数据，直接使用
        # print('LaborSupply_Huggett: 直接处理年龄组效率冲击数据。')
        eIdxM_group = eIdxM_input


        HHlaborM_group = np.zeros((nSim, cS.aD_new))
        leGridV_col_local = paramS.leGridV.flatten()
        
        # 直接在年龄组层面计算劳动供给
        for a_group in range(cS.aD_new):
            if a_group < cS.aR_new:  # 工作年龄组
                current_eIdx = eIdxM_group[:, a_group]
                # 🔧 添加边界检查，确保索引不超出leGridV范围
                current_eIdx = np.clip(current_eIdx, 0, len(leGridV_col_local) - 1)
                labor_eff_for_valid = leGridV_col_local[current_eIdx]
                HHlaborM_group[:, a_group] = cS.ageEffV_new[a_group] * labor_eff_for_valid

        # 计算总体人均有效劳动供给
        L_total_eff_pc_sum = 0
        if cS.aR_new > 0:
            mean_labor_per_working_group = np.mean(HHlaborM_group[:, :cS.aR_new], axis=0)
            L_total_eff_pc_sum = np.dot(mean_labor_per_working_group, Z_ss_norm_group[:cS.aR_new])
            
        L_total_eff_pc = max(0, L_total_eff_pc_sum)
        # print(f'家庭劳动供给已计算（年龄组版本）。总体人均有效劳动供给 (L_eff_pc) = {L_total_eff_pc:.4f}')
        return HHlaborM_group, L_total_eff_pc

    @staticmethod
    def HHPrices_Huggett(K_productive, L_total_eff, cS):
        """
        HHPrices_Huggett - 根据边际生产力计算要素价格
        """
        if K_productive <= 0:
            K_productive = 1e-6
        if L_total_eff <= 0:
            L_total_eff = 1e-6
            
        Y_gross = cS.A * (K_productive**cS.alpha) * (L_total_eff**(1-cS.alpha))
        MPK_gross_val = cS.alpha * Y_gross / K_productive
        MPL_gross = (1-cS.alpha) * Y_gross / L_total_eff
        R_market_gross_factor = 1 + MPK_gross_val - cS.ddk
        R_market_gross_factor = max(1.0 + 1e-6, R_market_gross_factor)
        return R_market_gross_factor, MPL_gross

    @staticmethod
    def CES_utility(cM_quantity, sigma_crra, cS_common):
        """
        CES_utility - 计算CRRA效用函数和边际效用
        """
        # Convert to numpy array to handle both scalar and array inputs
        cM_quantity = np.asarray(cM_quantity)
        scalar_input = cM_quantity.ndim == 0
        
        # Ensure we have at least 1D array for consistent processing
        if scalar_input:
            cM_quantity = cM_quantity.reshape(1)
        
        c_adjusted_quantity = np.maximum(cS_common.cFloor, cM_quantity)
        
        is_valid_consumption = (cM_quantity >= cS_common.cFloor)
        utilM = np.full(cM_quantity.shape, -np.inf)
        muM = np.full(cM_quantity.shape, np.inf)

        if abs(sigma_crra - 1) < 1e-6: # Log utility
            utilM[is_valid_consumption] = np.log(c_adjusted_quantity[is_valid_consumption])
            muM[is_valid_consumption] = 1.0 / c_adjusted_quantity[is_valid_consumption]
        else: # CRRA utility
            utilM[is_valid_consumption] = (c_adjusted_quantity[is_valid_consumption]**(1-sigma_crra)) / (1-sigma_crra)
            muM[is_valid_consumption] = c_adjusted_quantity[is_valid_consumption]**(-sigma_crra)

        # Penalty for consumption below floor
        if np.any(~is_valid_consumption):
             utilM[~is_valid_consumption] = -1e10 - (cS_common.cFloor - cM_quantity[~is_valid_consumption])*1e10
        
        # Convert back to scalar if input was scalar
        if scalar_input:
            muM = muM[0]
            utilM = utilM[0]
            
        return muM, utilM

    @staticmethod
    def HHIncome_Huggett(k_now_val, R_k_net_factor, w_gross, TR_total, b_payg_val, 
                         c_pps_chosen, a_idx, paramS_hh, cS, epsilon_val):
        """
        [最终修正版] 计算家庭收入，与MATLAB的`HHIncome_Huggett`物理过程完全一致。

        核心修正:
        - R_k_net_factor 被正确地解释为税后回报因子, 即 (1 + r_net)。
        - 资本收入 = k_now_val * (R_k_net_factor - 1)，不再错误地重复计算资本税。
        - 确保工作期和退休期的非资本收入计算逻辑清晰且正确。

        Args:
            k_now_val (float): 当前非PPS资产存量。
            R_k_net_factor (float): 市场税后资本回报因子 (1 + r_net)。
            w_gross (float): 市场工资率。
            TR_total (float): 政府 lump-sum 转移支付。
            b_payg_val (float): 当前年龄组的PAYG养老金。
            c_pps_chosen (float): 当前选择的PPS缴费金额。
            a_idx (int): 当前年龄组索引 (0-based)。
            paramS_hh (object): 包含税率等家庭特定参数的对象。
            cS (object): 包含模型公共参数的对象。
            epsilon_val (float): 当前的个人劳动效率冲击。

        Returns:
            tuple: (
                resources_for_c_and_k_prime,      # 可用于消费和储蓄的总资源
                actual_pps_contribution_expenditure, # 实际的PPS支出
                capital_income_net_of_tax           # 税后资本收入
            )
        """
        
        # --- 1. 计算非资本收入和PPS支出 ---
        
        # a_idx是0-based Python索引, cS.aR_new是1-based长度
        if a_idx < cS.aR_new:  # 工作年龄组
            # a. 计算税后劳动收入
            age_efficiency = cS.ageEffV_new[a_idx]
            labor_income_gross = w_gross * age_efficiency * epsilon_val
            
            # 如果PPS缴费可以税收递延，则从应税收入中扣除
            pps_deduction = 0.0
            if hasattr(paramS_hh, 'pps_tax_deferral_active') and paramS_hh.pps_tax_deferral_active:
                pps_deduction = c_pps_chosen
                
            labor_income_taxable = labor_income_gross - pps_deduction
            labor_income_tax = paramS_hh.tau_l * max(0, labor_income_taxable)
            labor_income_net_of_tax = labor_income_gross - labor_income_tax
            
            # b. 记录实际的PPS支出
            actual_pps_contribution_expenditure = c_pps_chosen
            
            # c. 工作期非资本总收入
            non_capital_income = labor_income_net_of_tax + TR_total + b_payg_val
            
        else:  # 退休年龄组
            # a. 退休后无PPS缴费
            actual_pps_contribution_expenditure = 0.0
            
            # b. 计算税后PPS提取收入
            pps_withdrawal_net_of_tax = 0.0
            # paramS_hh.current_pps_withdrawal 是由环境在step开始时计算并传入的当期应提取总额
            if hasattr(paramS_hh, 'current_pps_withdrawal') and paramS_hh.current_pps_withdrawal > 0:
                pps_withdrawal_gross = paramS_hh.current_pps_withdrawal
                pps_withdrawal_tax = cS.pps_tax_rate_withdrawal * pps_withdrawal_gross
                pps_withdrawal_net_of_tax = pps_withdrawal_gross - pps_withdrawal_tax
            
            # c. 退休期非资本总收入
            non_capital_income = TR_total + b_payg_val + pps_withdrawal_net_of_tax
        
        # --- 2. [核心修正] 计算税后资本收入 ---
        # R_k_net_factor 是 (1 + r_net)，所以 k * (R_k_net_factor - 1) 直接就是税后资本收入。
        # 这个简单的表达式避免了所有关于如何计算 r_gross, r_net, tau_k 的混淆。
        capital_income_net_of_tax = k_now_val * (R_k_net_factor - 1)
        
        # --- 3. 计算可用于消费和新储蓄的总资源 ---
        # 总资源 = 期初财富 + 所有税后收入 - 当期PPS支出
        resources_for_c_and_k_prime = (k_now_val + 
                                       capital_income_net_of_tax + 
                                       non_capital_income - 
                                       actual_pps_contribution_expenditure)
        
        return resources_for_c_and_k_prime, actual_pps_contribution_expenditure, capital_income_net_of_tax

    @staticmethod
    def HHSolution_VFI_Huggett(R_k_net_factor_vfi, w_gross_vfi, TR_total_vfi, bV_payg_vfi, paramS_vfi, cS):
        """
        HHSolution_VFI_Huggett - V8模型的值函数迭代主求解器
        """
        cPolM_q = np.zeros((cS.nk, cS.nkpps, cS.nw, cS.aD_new))
        kPolM = np.zeros((cS.nk, cS.nkpps, cS.nw, cS.aD_new))
        cPpsPolM_choice = np.zeros((cS.nk, cS.nkpps, cS.nw, cS.aD_new))
        valM = np.full((cS.nk, cS.nkpps, cS.nw, cS.aD_new), -np.inf)
        
        # Backward induction
        for a_idx in range(cS.aD_new - 1, -1, -1):
            vPrime_kkppse_next = None
            if a_idx < cS.aD_new - 1:
                vPrime_kkppse_next = valM[:, :, :, a_idx + 1]
            
            eps_grid_for_vfi = paramS_vfi.leGridV
            
            # This is the Python version of the fmincon-based solver
            cPolM_q[:,:,:,a_idx], kPolM[:,:,:,a_idx], cPpsPolM_choice[:,:,:,a_idx], valM[:,:,:,a_idx] = \
                OLG_V9_Utils.HHSolutionByAge_VFI_Huggett_v8_fmincon(
                    a_idx, vPrime_kkppse_next, R_k_net_factor_vfi, w_gross_vfi, TR_total_vfi, 
                    bV_payg_vfi[a_idx], paramS_vfi, cS, eps_grid_for_vfi
                )
        return cPolM_q, kPolM, cPpsPolM_choice, valM

    @staticmethod
    def HHSolutionByAge_VFI_Huggett_v8_fmincon(a_idx, vPrime_kkppse_next, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, paramS_age, cS, epsilon_grid):
        """
        HHSolutionByAge_VFI_Huggett_v8_fmincon - 基于fmincon的连续优化VFI算法
        Python equivalent using scipy.optimize.minimize with 'SLSQP'
        """
        cPol_age_q = np.zeros((cS.nk, cS.nkpps, cS.nw))
        kPol_age = np.zeros((cS.nk, cS.nkpps, cS.nw))
        cPpsPol_age_choice = np.zeros((cS.nk, cS.nkpps, cS.nw))
        val_age = np.full((cS.nk, cS.nkpps, cS.nw), -np.inf)

        # Last period solution
        if a_idx == cS.aD_new - 1:
            k_grid_nd, kpps_grid_nd, eps_grid_nd = np.meshgrid(cS.kGridV, cS.kppsGridV, epsilon_grid, indexing='ij')
            
            resources_batch = np.zeros_like(k_grid_nd)
            for ik, ikpps, ie in np.ndindex(k_grid_nd.shape):
                resources_batch[ik, ikpps, ie], _, _ = OLG_V9_Utils.HHIncome_Huggett(
                    k_grid_nd[ik, ikpps, ie], R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val,
                    0, a_idx, paramS_age, cS, eps_grid_nd[ik, ikpps, ie]
                )

            total_resources = resources_batch
            if cS.pps_active:
                total_resources += kpps_grid_nd * (1 - cS.pps_tax_rate_withdrawal)
            
            cPol_age_q = np.maximum(cS.cFloor, total_resources / (1 + cS.tau_c))
            kPol_age[:] = cS.kMin
            cPpsPol_age_choice[:] = 0
            _, val_age = OLG_V9_Utils.CES_utility(cPol_age_q, cS.sigma, cS)
            return cPol_age_q, kPol_age, cPpsPol_age_choice, val_age

        # Other periods
        EV_matrix = np.zeros((cS.nk, cS.nkpps, cS.nw))
        for ie_current in range(cS.nw):
            transition_probs = paramS_age.leTrProbM[ie_current, :]
            EV_slice = np.sum(vPrime_kkppse_next * transition_probs.reshape(1, 1, -1), axis=2)
            EV_matrix[:, :, ie_current] = EV_slice

        EV_interpolants = []
        for ie_current in range(cS.nw):
            if cS.nk > 1 and cS.nkpps > 1:
                # 最接近MATLAB griddedInterpolant('linear','linear')的设置：
                # - method='linear': 双线性插值，对应MATLAB第一个'linear'
                # - bounds_error=False: 允许超出边界，对应MATLAB默认行为
                # - fill_value=None: 使用边界值外推，接近MATLAB第二个'linear'的效果
                interp = RegularGridInterpolator((cS.kGridV, cS.kppsGridV), EV_matrix[:, :, ie_current], 
                                                 method='linear', bounds_error=False, fill_value=None)
            elif cS.nk > 1:
                # 一维情况：线性插值和边界外推
                interp = RegularGridInterpolator((cS.kGridV,), EV_matrix[:, 0, ie_current],
                                                 method='linear', bounds_error=False, fill_value=None)
            elif cS.nkpps > 1:
                # 一维情况：线性插值和边界外推
                interp = RegularGridInterpolator((cS.kppsGridV,), EV_matrix[0, :, ie_current],
                                                 method='linear', bounds_error=False, fill_value=None)
            else:
                interp = lambda k, kpps: EV_matrix[0, 0, ie_current]
            EV_interpolants.append(interp)

        # Loop over all states
        for ik, ikpps, ie in np.ndindex(cS.nk, cS.nkpps, cS.nw):
            k_state = cS.kGridV[ik]
            k_pps_state = cS.kppsGridV[ikpps]
            epsilon_state = epsilon_grid[ie]

            # 🔧 修复：使用年龄组逻辑替代年度年龄逻辑
            is_pps_eligible = (a_idx < cS.aR_new and cS.pps_active)
            
            max_permissible_cpps = 0
            if is_pps_eligible:
                age_efficiency = cS.ageEffV_new[a_idx]
                current_gross_labor_income = w_gross_age * age_efficiency * epsilon_state
                if current_gross_labor_income > 1e-6:
                    max_cpps_by_frac = current_gross_labor_income * cS.pps_max_contrib_frac
                    max_permissible_cpps = min(cS.pps_contrib_limit, max_cpps_by_frac)
                    max_permissible_cpps = max(0, max_permissible_cpps)
            
            # Define objective for scipy.optimize.minimize
            obj_func = lambda x_prop: OLG_V9_Utils.fmincon_objective_helper_proportional(
                x_prop, k_state, k_pps_state, epsilon_state, a_idx, ie,
                R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val,
                paramS_age, cS, EV_interpolants, max_permissible_cpps
            )

            bounds = [(0, 1), (0, 1)]
            x0_pps_prop = 0.5
            if max_permissible_cpps < 1e-9:
                x0_pps_prop = 0
                bounds[0] = (0, 0)
            x0 = [x0_pps_prop, 0.5]

            res = minimize(obj_func, x0, method='SLSQP', bounds=bounds, 
                           options={'ftol': 1e-7, 'maxiter': 500})
            
            if res.success:
                pps_prop_opt, k_prime_prop_opt = res.x
                optimal_value_val = -res.fun

                optimal_cpps_val = pps_prop_opt * max_permissible_cpps
                
                resources_after_pps_opt, _, _ = OLG_V9_Utils.HHIncome_Huggett(
                    k_state, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val,
                    optimal_cpps_val, a_idx, paramS_age, cS, epsilon_state)

                consumption_floor_spending_opt = cS.cFloor * (1 + cS.tau_c)
                resources_for_kprime_c_above_floor_opt = resources_after_pps_opt - consumption_floor_spending_opt
                
                if resources_for_kprime_c_above_floor_opt >= 0:
                    optimal_k_prime_val = k_prime_prop_opt * resources_for_kprime_c_above_floor_opt
                    optimal_k_prime_val = max(cS.kMin, min(optimal_k_prime_val, resources_for_kprime_c_above_floor_opt))
                else:
                    optimal_k_prime_val = cS.kMin
                
                optimal_k_prime_val = max(cS.kMin, min(optimal_k_prime_val, cS.kMax))
                
                consumption_expenditure_opt = resources_after_pps_opt - optimal_k_prime_val
                optimal_c_val = max(cS.cFloor, consumption_expenditure_opt / (1 + cS.tau_c))
            else:
                # Fallback to a simpler discrete search if optimizer fails
                optimal_c_val, optimal_k_prime_val, optimal_cpps_val, optimal_value_val = \
                    OLG_V9_Utils.fallback_discrete_solution(
                        k_state, k_pps_state, epsilon_state, a_idx, ie,
                        R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val,
                        paramS_age, cS, EV_interpolants, max_permissible_cpps)
            
            val_age[ik, ikpps, ie] = optimal_value_val
            cPol_age_q[ik, ikpps, ie] = optimal_c_val
            kPol_age[ik, ikpps, ie] = optimal_k_prime_val
            cPpsPol_age_choice[ik, ikpps, ie] = optimal_cpps_val
            
        return cPol_age_q, kPol_age, cPpsPol_age_choice, val_age

    @staticmethod
    def fmincon_objective_helper_proportional(x_prop, k_state, k_pps_state, epsilon_state, a_idx, ie, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, paramS_age, cS, EV_interpolants, max_permissible_cpps):
        pps_proportion, k_prime_proportion = x_prop

        actual_c_pps = pps_proportion * max_permissible_cpps
        
        resources_after_pps, _, _ = OLG_V9_Utils.HHIncome_Huggett(
            k_state, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val,
            actual_c_pps, a_idx, paramS_age, cS, epsilon_state)

        consumption_floor_spending = cS.cFloor * (1 + cS.tau_c)
        resources_for_kprime_and_c_above_floor = resources_after_pps - consumption_floor_spending

        if resources_for_kprime_and_c_above_floor >= 0:
            actual_k_prime = k_prime_proportion * resources_for_kprime_and_c_above_floor
            actual_k_prime = max(cS.kMin, min(actual_k_prime, resources_for_kprime_and_c_above_floor))
        else:
            actual_k_prime = cS.kMin
        
        actual_k_prime = max(cS.kMin, min(actual_k_prime, cS.kMax))
        
        consumption_expenditure = resources_after_pps - actual_k_prime
        current_c = max(cS.cFloor, consumption_expenditure / (1 + cS.tau_c))
        
        _, current_utility = OLG_V9_Utils.CES_utility(current_c, cS.sigma, cS)
        if not np.isfinite(current_utility):
            return 1e12 + abs(current_c - cS.cFloor) * 1e10

        pps_withdrawal = 0
        # 🔧 修复：使用年龄组逻辑替代年度年龄逻辑
        if a_idx >= cS.aR_new and cS.pps_active:
            pps_withdrawal = k_pps_state * cS.pps_withdrawal_rate

        pps_return_factor = 1 + ((R_k_net_factor_age - 1) + cS.pps_return_rate_premium)
        k_pps_prime = (k_pps_state + actual_c_pps - pps_withdrawal) * pps_return_factor
        k_pps_prime = max(cS.kppsMin, min(cS.kppsMax, k_pps_prime))

        expected_future_value = -np.inf
        if a_idx < cS.aD_new - 1:
            # 调用方式要匹配MATLAB的griddedInterpolant(query_point)行为
            if cS.nk > 1 and cS.nkpps > 1:
                # 二维插值：点格式需要是[k, k_pps]，RegularGridInterpolator返回标量
                query_point = np.array([actual_k_prime, k_pps_prime])
                result = EV_interpolants[ie](query_point)
                expected_future_value = result.item() if hasattr(result, 'item') else result
            elif cS.nk > 1:
                # 一维插值：只有k维度
                result = EV_interpolants[ie](np.array([actual_k_prime]))
                expected_future_value = result.item() if hasattr(result, 'item') else result
            elif cS.nkpps > 1:
                # 一维插值：只有k_pps维度  
                result = EV_interpolants[ie](np.array([k_pps_prime]))
                expected_future_value = result.item() if hasattr(result, 'item') else result
            else: # nk=1, nkpps=1
                # 常数情况
                expected_future_value = EV_interpolants[ie](actual_k_prime, k_pps_prime)

        if not np.isfinite(expected_future_value):
            expected_future_value = -1e11

        s_transition = cS.s_1yr_transitionV[a_idx]
        total_value = current_utility + cS.beta * s_transition * expected_future_value

        return -total_value if np.isfinite(total_value) else 1e12

    @staticmethod
    def fallback_discrete_solution(k_state, k_pps_state, epsilon_state, a_idx, ie, 
                                  R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, 
                                  paramS_age, cS, EV_interpolants, max_permissible_cpps):
        """
        fallback_discrete_solution - 当连续优化失败时的离散网格搜索后备方案
        """
        # 简单的网格搜索
        best_value = -np.inf
        best_c = cS.cFloor
        best_k_prime = cS.kMin
        best_cpps = 0
        
        # 在可行域内搜索
        n_grid_points = 5  # 简化网格
        cpps_grid = np.linspace(0, max_permissible_cpps, n_grid_points)
        
        for cpps_test in cpps_grid:
            resources_after_pps, _, _ = OLG_V9_Utils.HHIncome_Huggett(
                k_state, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val,
                cpps_test, a_idx, paramS_age, cS, epsilon_state)
            
            resources_for_k_and_c = max(0, resources_after_pps - cS.cFloor * (1 + cS.tau_c))
            k_prime_grid = np.linspace(cS.kMin, min(cS.kMax, resources_for_k_and_c), n_grid_points)
            
            for k_prime_test in k_prime_grid:
                consumption_expenditure = resources_after_pps - k_prime_test
                c_test = max(cS.cFloor, consumption_expenditure / (1 + cS.tau_c))
                
                _, utility = OLG_V9_Utils.CES_utility(c_test, cS.sigma, cS)
                
                if np.isfinite(utility) and utility > best_value:
                    best_value = utility
                    best_c = c_test
                    best_k_prime = k_prime_test
                    best_cpps = cpps_test
        
        return best_c, best_k_prime, best_cpps, best_value

    @staticmethod
    def HHSimulation_olgm(kPolM_4D_input, cPpsPolM_choice_4D_input, cPolM_consump_q_4D_input, eIdxM_input, R_k_net_factor_hh_sim, w_gross_sim_price, TR_total_sim_transfer, bV_payg_sim_benefit, paramS_sim_household, cS_common_sim):
        """
        HHSimulation_olgm - 基于最优策略模拟家庭生命周期路径（年龄组版本）
        
        重要更新：改为按年龄组模拟，而不是年度年龄模拟
        - 模拟维度：从aD_orig（年度）改为aD_new（年龄组）
        - 效率冲击：自动检测并处理年度或年龄组数据
        - 年龄判断：直接使用年龄组索引，无需复杂映射
        - 存活概率：使用年龄组间转移存活率
        
        Args:
            eIdxM_input: 效率冲击矩阵，可以是年度 (nSim × aD_orig) 或年龄组 (nSim × aD_new)
        """
        
        nSim_sim = eIdxM_input.shape[0]
        

        # 年龄组数据，直接使用
        print('HHSimulation_olgm: 直接使用年龄组效率冲击数据。')
        eIdxM_group = eIdxM_input

        # 初始化结果存储（年龄组尺度）
        kHistM_out = np.zeros((nSim_sim, cS_common_sim.aD_new))
        kPpsHistM_out = np.zeros((nSim_sim, cS_common_sim.aD_new))
        cHistM_out = np.zeros((nSim_sim, cS_common_sim.aD_new))
        cppsHistM_out = np.zeros((nSim_sim, cS_common_sim.aD_new))  # PPS缴费路径

        leGridV_col_sim = paramS_sim_household.leGridV.flatten()

        # 创建插值器
        kPolInterp_sim = [[None for _ in range(cS_common_sim.aD_new)] for _ in range(cS_common_sim.nw)]
        cPpsPolInterp_choice_sim = [[None for _ in range(cS_common_sim.aD_new)] for _ in range(cS_common_sim.nw)]
        cPolqInterp_sim = [[None for _ in range(cS_common_sim.aD_new)] for _ in range(cS_common_sim.nw)]

        for ia_interp in range(cS_common_sim.aD_new):
            for ie_interp in range(cS_common_sim.nw):
                kPol_slice = np.squeeze(kPolM_4D_input[:, :, ie_interp, ia_interp])
                cPpsPol_slice = np.squeeze(cPpsPolM_choice_4D_input[:, :, ie_interp, ia_interp])
                cPolq_slice = np.squeeze(cPolM_consump_q_4D_input[:, :, ie_interp, ia_interp])

                if cS_common_sim.nk > 1 and cS_common_sim.nkpps > 1:
                    # 🔧 使用线性插值避免网格点不足问题（nkpps=3 < 4，不支持立方插值）
                    kPolInterp_sim[ie_interp][ia_interp] = RegularGridInterpolator(
                        (cS_common_sim.kGridV, cS_common_sim.kppsGridV), kPol_slice, 
                        method='linear', bounds_error=False, fill_value=None)
                    cPpsPolInterp_choice_sim[ie_interp][ia_interp] = RegularGridInterpolator(
                        (cS_common_sim.kGridV, cS_common_sim.kppsGridV), cPpsPol_slice, 
                        method='linear', bounds_error=False, fill_value=None)
                    cPolqInterp_sim[ie_interp][ia_interp] = RegularGridInterpolator(
                        (cS_common_sim.kGridV, cS_common_sim.kppsGridV), cPolq_slice, 
                        method='linear', bounds_error=False, fill_value=None)
                elif cS_common_sim.nk > 1 and cS_common_sim.nkpps == 1:
                    # 🔧 使用线性插值（nk=5足够支持立方，但为一致性使用线性）
                    kPolInterp_sim[ie_interp][ia_interp] = RegularGridInterpolator(
                        (cS_common_sim.kGridV,), kPol_slice[:, 0], 
                        method='linear', bounds_error=False, fill_value=None)
                    cPpsPolInterp_choice_sim[ie_interp][ia_interp] = RegularGridInterpolator(
                        (cS_common_sim.kGridV,), cPpsPol_slice[:, 0], 
                        method='linear', bounds_error=False, fill_value=None)
                    cPolqInterp_sim[ie_interp][ia_interp] = RegularGridInterpolator(
                        (cS_common_sim.kGridV,), cPolq_slice[:, 0], 
                        method='linear', bounds_error=False, fill_value=None)
                elif cS_common_sim.nk == 1 and cS_common_sim.nkpps > 1:
                    # 🔧 使用线性插值（nkpps=3 < 4，不支持立方插值）
                    kPolInterp_sim[ie_interp][ia_interp] = RegularGridInterpolator(
                        (cS_common_sim.kppsGridV,), kPol_slice[0, :], 
                        method='linear', bounds_error=False, fill_value=None)
                    cPpsPolInterp_choice_sim[ie_interp][ia_interp] = RegularGridInterpolator(
                        (cS_common_sim.kppsGridV,), cPpsPol_slice[0, :], 
                        method='linear', bounds_error=False, fill_value=None)
                    cPolqInterp_sim[ie_interp][ia_interp] = RegularGridInterpolator(
                        (cS_common_sim.kppsGridV,), cPolq_slice[0, :], 
                        method='linear', bounds_error=False, fill_value=None)
                elif cS_common_sim.nk == 1 and cS_common_sim.nkpps == 1:
                    kVal = kPol_slice[0, 0] if kPol_slice.ndim == 2 else kPol_slice.item()
                    cPpsVal = cPpsPol_slice[0, 0] if cPpsPol_slice.ndim == 2 else cPpsPol_slice.item()
                    cPolqVal = cPolq_slice[0, 0] if cPolq_slice.ndim == 2 else cPolq_slice.item()
                    kPolInterp_sim[ie_interp][ia_interp] = lambda x, y=None: kVal
                    cPpsPolInterp_choice_sim[ie_interp][ia_interp] = lambda x, y=None: cPpsVal
                    cPolqInterp_sim[ie_interp][ia_interp] = lambda x, y=None: cPolqVal
                else:
                    raise ValueError('HHSimulation_olgm: nk 或 nkpps 为零，无法创建插值器。')

        pps_return_factor_sim = 1 + ((R_k_net_factor_hh_sim - 1) + cS_common_sim.pps_return_rate_premium)

        # 主循环：年龄组模拟
        for a_group in range(cS_common_sim.aD_new):  # 循环年龄组
            # 🔧 修复：获取当期资产状态（从上期的策略决策中获得）
            if a_group == 0:
                # 第一期：初始资产状态
                kNowV_group = np.full(nSim_sim, cS_common_sim.kMin)
                kPpsNowV_group = np.full(nSim_sim, cS_common_sim.kppsMin)
            else:
                # 后续期：从上期模拟结果中获取当期资产
                kNowV_group = kHistM_out[:, a_group]  # 这是上期决定的本期初始资产
                kPpsNowV_group = kPpsHistM_out[:, a_group]  # 这是上期PPS的本期初始值
            
            # 年龄组判断（简化）
            is_working_age_group = (a_group < cS_common_sim.aR_new)
            is_pps_withdrawal_eligible_group = (not is_working_age_group and cS_common_sim.pps_active)

            pps_withdrawal_pretax_group = np.zeros(nSim_sim)
            if is_pps_withdrawal_eligible_group:
                pps_withdrawal_pretax_group = kPpsNowV_group * cS_common_sim.pps_withdrawal_rate
                
                # 设置PPS提取参数
                if not hasattr(paramS_sim_household, 'pps_tax_rate_withdrawal'):
                    paramS_sim_household.pps_tax_rate_withdrawal = 0.15
                paramS_sim_household.current_pps_withdrawal = np.mean(pps_withdrawal_pretax_group)
            else:
                paramS_sim_household.current_pps_withdrawal = 0
            
            actual_cpps_group = np.zeros(nSim_sim)
            k_prime_group = np.zeros(nSim_sim)  # 🔧 新增：下期非PPS资产选择

            # 按效率状态循环进行插值
            for ie_sim_idx in range(cS_common_sim.nw):
                simIdxV_for_this_e = np.where(eIdxM_group[:, a_group] == ie_sim_idx)[0]
                if len(simIdxV_for_this_e) == 0:
                    continue

                # 网格约束处理
                kNow_clamped = np.clip(kNowV_group[simIdxV_for_this_e], 
                                     cS_common_sim.kGridV[0], cS_common_sim.kGridV[-1])
                kPpsNow_clamped = np.clip(kPpsNowV_group[simIdxV_for_this_e], 
                                        cS_common_sim.kppsGridV[0], cS_common_sim.kppsGridV[-1])

                # 🔧 修复：插值获取策略（k'策略，PPS缴费策略，消费策略）
                if cS_common_sim.nk > 1 and cS_common_sim.nkpps > 1:
                    points = np.vstack((kNow_clamped, kPpsNow_clamped)).T
                    k_prime_group[simIdxV_for_this_e] = kPolInterp_sim[ie_sim_idx][a_group](points)
                    actual_cpps_group[simIdxV_for_this_e] = cPpsPolInterp_choice_sim[ie_sim_idx][a_group](points)
                    cHistM_out[simIdxV_for_this_e, a_group] = cPolqInterp_sim[ie_sim_idx][a_group](points)
                elif cS_common_sim.nk > 1 and cS_common_sim.nkpps == 1:
                    k_prime_group[simIdxV_for_this_e] = kPolInterp_sim[ie_sim_idx][a_group](kNow_clamped)
                    actual_cpps_group[simIdxV_for_this_e] = cPpsPolInterp_choice_sim[ie_sim_idx][a_group](kNow_clamped)
                    cHistM_out[simIdxV_for_this_e, a_group] = cPolqInterp_sim[ie_sim_idx][a_group](kNow_clamped)
                elif cS_common_sim.nk == 1 and cS_common_sim.nkpps > 1:
                    k_prime_group[simIdxV_for_this_e] = kPolInterp_sim[ie_sim_idx][a_group](kPpsNow_clamped)
                    actual_cpps_group[simIdxV_for_this_e] = cPpsPolInterp_choice_sim[ie_sim_idx][a_group](kPpsNow_clamped)
                    cHistM_out[simIdxV_for_this_e, a_group] = cPolqInterp_sim[ie_sim_idx][a_group](kPpsNow_clamped)
                else:  # nk=1, nkpps=1
                    if len(simIdxV_for_this_e) > 0:
                        k_prime_group[simIdxV_for_this_e] = kPolInterp_sim[ie_sim_idx][a_group](kNow_clamped[0])
                        actual_cpps_group[simIdxV_for_this_e] = cPpsPolInterp_choice_sim[ie_sim_idx][a_group](kNow_clamped[0])
                        cHistM_out[simIdxV_for_this_e, a_group] = cPolqInterp_sim[ie_sim_idx][a_group](kNow_clamped[0])

            # 🔧 修复：记录当期资产和决策
            if a_group == 0:
                # 第一期：记录初始资产
                kHistM_out[:, a_group] = kNowV_group  # 当期初始资产
                kPpsHistM_out[:, a_group] = kPpsNowV_group  # 当期初始PPS资产
            # 对于后续期，资产已在上期更新时设置好，无需重复设置

            # PPS逻辑处理
            if cS_common_sim.pps_active:
                can_contribute_pps_group = is_working_age_group
                if not can_contribute_pps_group:
                    actual_cpps_group[:] = 0
                else:
                    actual_cpps_group = np.maximum(0, actual_cpps_group)
                
                # 🔧 修复：PPS资产演化到下期
                if a_group < cS_common_sim.aD_new - 1:  # 非最后一期
                    kPpsHistM_out[:, a_group + 1] = ((kPpsNowV_group + actual_cpps_group - 
                                                    pps_withdrawal_pretax_group) * pps_return_factor_sim)
                    kPpsHistM_out[:, a_group + 1] = np.clip(kPpsHistM_out[:, a_group + 1], 
                                                           cS_common_sim.kppsMin, cS_common_sim.kppsMax)
            else:
                if a_group < cS_common_sim.aD_new - 1:
                    kPpsHistM_out[:, a_group + 1] = 0

            # 🔧 修复：非PPS资产演化到下期
            if a_group < cS_common_sim.aD_new - 1:
                kHistM_out[:, a_group + 1] = np.clip(k_prime_group, 
                                                   cS_common_sim.kMin, cS_common_sim.kMax)
            
            # 确保消费非负
            cHistM_out[:, a_group] = np.maximum(cS_common_sim.cFloor, cHistM_out[:, a_group])
            cppsHistM_out[:, a_group] = actual_cpps_group

        return kHistM_out, kPpsHistM_out, cHistM_out, cppsHistM_out


# --- 在 main_olg_v9_utils.py 中，替换 HHSimulation_olgm_rl 函数 ---

# 在 main_olg_v9_utils.py 中，替换此函数

    @staticmethod
    def HHSimulation_olgm_rl(rl_model: Any, rl_config: Dict, eIdxM_input: np.ndarray, 
                            R_k_net_factor_hh_sim: float, w_gross_sim_price: float, 
                            TR_total_sim_transfer: float, bV_payg_sim_benefit: np.ndarray, 
                            paramS_sim_household: Any, cS_common_sim: Any) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        [v3 - 终极对齐版] 使用RL模型进行家庭模拟，并实现与VFI完全一致的最后一期清算逻辑。
        """
        n_sim = eIdxM_input.shape[0]
        cS_obj = cS_common_sim
        paramS_rl = rl_config['paramS_for_rl']
        sim_env = OLGEnvV9SAC(cS=cS_obj, paramS_rl=paramS_rl, rng_M=rl_config['rng_M'], training_mode=False)
        
        b_payg_avg_retiree = np.mean(bV_payg_sim_benefit[bV_payg_sim_benefit > 0]) if np.any(bV_payg_sim_benefit > 0) else 0.0
        M_current = {
            'R_k_net_factor': R_k_net_factor_hh_sim, 'w_gross': w_gross_sim_price,
            'TR_total': TR_total_sim_transfer, 'b_payg_avg_for_obs': b_payg_avg_retiree,
            'tau_l': paramS_sim_household.tau_l, 'theta_payg_actual': paramS_sim_household.theta_payg_actual_for_hh
        }
        
        k_paths, kpps_paths, c_paths, cpps_paths = [], [], [], []

        for i_sim in range(n_sim):
            sim_env.current_age_idx = 1
            sim_env.current_k_val = cS_obj.kMin
            sim_env.current_k_pps_val = cS_obj.kppsMin
            
            sim_env.current_M = M_current
            sim_env.current_bV_payg = bV_payg_sim_benefit
            
            k_path, kpps_path, c_path, cpps_path = [], [], [], []
            
            for age_idx in range(cS_obj.aD_new):
                sim_env.current_age_idx = age_idx + 1
                sim_env.current_eps_idx = eIdxM_input[i_sim, age_idx]

                k_path.append(sim_env.current_k_val)
                kpps_path.append(sim_env.current_k_pps_val)
                
                is_last_period = (age_idx == cS_obj.aD_new - 1)
                
                if not is_last_period:
                    obs = sim_env._get_observation()
                    action, _ = rl_model.predict(obs, deterministic=True)
                    _, _, terminated, _, info = sim_env.step(action)
                    
                    c_path.append(info['consumption'])
                    cpps_path.append(info['c_pps'])
                    
                    if terminated: break
                else:
                    # --- [核心修复] 最后一期特殊清算逻辑，与VFI完全对齐 ---
                    
                    # 1. 在最后一期，PPS缴费必须为0
                    c_pps_last = 0.0
                    cpps_path.append(c_pps_last)
                    
                    # 2. 计算除了PPS清算财富之外的所有资源
                    #    我们直接调用HHIncome_Huggett，就像VFI求解器做的那样
                    paramS_hh_last = TempParamSHH(M_current['tau_l'], M_current['theta_payg_actual'], cS_obj.pps_active, cS_obj.ageEffV_new)
                    paramS_hh_last.current_pps_withdrawal = 0 # 在这里不考虑按比例提取
                    
                    epsilon_val_last = paramS_rl['leGridV'][sim_env.current_eps_idx - 1]
                    b_payg_last = sim_env.current_bV_payg[age_idx]

                    base_resources, _, _ = OLG_V9_Utils.HHIncome_Huggett(
                        sim_env.current_k_val, M_current['R_k_net_factor'], M_current['w_gross'],
                        M_current['TR_total'], b_payg_last, c_pps_last, age_idx, 
                        paramS_hh_last, cS_obj, epsilon_val_last
                    )
                    
                    # 3. 计算税后清算的PPS财富
                    liquidated_pps_wealth = 0.0
                    if cS_obj.pps_active:
                        liquidated_pps_wealth = sim_env.current_k_pps_val * (1 - cS_obj.pps_tax_rate_withdrawal)
                    
                    # 4. 计算总财富
                    total_wealth_to_spend = base_resources + liquidated_pps_wealth
                    
                    # 5. 计算最终消费
                    final_c = max(cS_obj.cFloor, total_wealth_to_spend / (1 + cS_obj.tau_c))
                    c_path.append(final_c)
                    
                    break # 结束生命周期

            k_paths.append(k_path)
            kpps_paths.append(kpps_path)
            c_paths.append(c_path)
            cpps_paths.append(cpps_path)
                
        return (np.array(k_paths), np.array(kpps_paths),
                np.array(c_paths), np.array(cpps_paths))

# --- 在 OLG_V9_Utils 类的静态方法中添加以下新函数 ---

# --- 在 OLG_V9_Utils 类的静态方法中，替换/添加以下函数 ---

    @staticmethod
    def HHSimulation_olgm_rl_simplified(rl_model: Any, cS_obj: Any, paramS_rl: Dict, 
                                        M_fixed: Dict, n_sim: int, 
                                        eIdxM_group: np.ndarray) -> Dict: # [修正] 更改函数签名
        """
        [v3 - 最终统一随机性版] 在一个固定的宏观环境下，模拟简化版RL智能体的生命周期路径。
        这个版本强制使用外部提供的效率路径 (eIdxM_group)，以确保与VFI的比较绝对公平。

        Args:
            rl_model: 训练好的简化版RL模型。
            cS_obj: 完整的模型参数对象。
            paramS_rl: RL相关的派生参数。
            M_fixed: 固定的宏观环境字典。
            n_sim: 模拟的个体数量。
            eIdxM_group: [新增] 预先生成的效率冲击路径 (n_sim x aD_new)，期望是1-based索引。

        Returns:
            一个包含所有模拟路径的字典 (k_path, kpps_path, c_path, cpps_path)。
        """
        print("🤖 运行简化版RL模拟 (强制统一效率路径)...")
        
        # 1. 创建与训练时完全一致的简化版评估环境
        eval_env = OLGEnvV9SACSimp(
            cS=cS_obj,
            paramS_rl=paramS_rl,
            M_fixed=M_fixed,
            training_mode=False
        )
        
        # [修正] 2. 循环模拟 n_sim 个生命周期，但不再依赖env.reset()的随机性
        k_paths, kpps_paths, c_paths, cpps_paths = [], [], [], []
        
        for i_sim in range(n_sim):
            # [修正] a. 手动管理生命周期状态，不调用 env.reset()
            eval_env.current_age_idx = 1
            eval_env.current_k_val = cS_obj.kMin
            eval_env.current_k_pps_val = cS_obj.kppsMin
            # 累积生存概率在env内部管理，每次step都会更新，但每次新生命周期开始前需要重置
            eval_env.cumulative_survival_prob = 1.0
            
            k_path, kpps_path, c_path, cpps_path = [], [], [], []
            
            for age_idx in range(cS_obj.aD_new): # 0-based
                # [修正] b. 强制设定当前状态，特别是使用“剧本”中的效率冲击
                eval_env.current_age_idx = age_idx + 1
                # k, k_pps 的值是从上一步循环中继承过来的，是正确的
                eval_env.current_eps_idx = eIdxM_group[i_sim, age_idx] # 强制使用剧本

                k_path.append(eval_env.current_k_val)
                kpps_path.append(eval_env.current_k_pps_val)
                
                # [修正] c. 从强制设定的状态中获取观测
                obs = eval_env._get_observation()
                
                # d. 模型根据当前状态预测动作
                action, _ = rl_model.predict(obs, deterministic=True)
                
                # [修正] e. 执行一步。env.step现在将在我们强制设定的状态(k,k_pps,age,eps)下执行
                # 注意：step内部会根据转移概率更新到下一个eps，但这没关系，因为在下一个
                # 循环开始时，我们会再次用剧本强制覆盖它。
                _, _, terminated, _, info = eval_env.step(action)
                
                # f. 记录决策结果
                c_path.append(info['consumption'])
                cpps_path.append(info['c_pps'])

                if terminated:
                    break # 如果生命周期提前结束（虽然在这个模型里不会）
            
            k_paths.append(k_path)
            kpps_paths.append(kpps_path)
            c_paths.append(c_path)
            cpps_paths.append(cpps_path)

        print("✅ 模拟完成。")
        return {
            "k_path_rl": np.array(k_paths),
            "kpps_path_rl": np.array(kpps_paths),
            "c_path_rl": np.array(c_paths),
            "cpps_path_rl": np.array(cpps_paths),
        }

# --- 在 OLG_V9_Utils 类的静态方法中，替换为以下最终版本 ---

 # 在 main_olg_v9_utils.py 中，替换此函数

# 在 main_olg_v9_utils.py 中，替换此函数

    @staticmethod
    def HHSimulation_olgm_vfi_simplified(vfi_results: Dict, n_sim: int, 
                                        eIdxM_group: np.ndarray) -> Dict: 
        """
        [v7 - 最终对齐版] 通过插值消费策略，消除“双重离散化”误差。
        
        此版本与MATLAB的HHSimulation_olgm在逻辑上完全对齐，可以公平地评估VFI策略。
        它接收完整的策略矩阵(k', c_pps, c)，并直接插值执行，而不是重新计算消费。
        """
        print("🎯 运行VFI模拟 (最终对齐版，插值消费)...")
        
        from scipy.interpolate import RegularGridInterpolator
        
        # 1. 提取参数和完整的策略矩阵
        cS_obj = vfi_results['cS_python_obj']
        M_fixed = vfi_results['M_test']
        kPolM = vfi_results['kPolM']
        cPpsPolM_choice = vfi_results['cPpsPolM_choice']
        cPolM = vfi_results.get('cPolM') # 使用.get()以保持向后兼容性

        if cPolM is None:
            raise ValueError("错误: vfi_results 中缺少消费策略矩阵 'cPolM'。请确保Python调用端已返回它。")

        # 2. [核心] 为所有三个决策变量创建插值器
        k_prime_interps = [[None for _ in range(cS_obj.aD_new)] for _ in range(cS_obj.nw)]
        c_pps_interps = [[None for _ in range(cS_obj.aD_new)] for _ in range(cS_obj.nw)]
        c_interps = [[None for _ in range(cS_obj.aD_new)] for _ in range(cS_obj.nw)]

        for ia in range(cS_obj.aD_new):
            for ie in range(cS_obj.nw):
                k_prime_slice = np.squeeze(kPolM[:, :, ie, ia])
                c_pps_slice = np.squeeze(cPpsPolM_choice[:, :, ie, ia])
                c_slice = np.squeeze(cPolM[:, :, ie, ia])

                # 动态构建插值点和检查数据形状
                points_list, data_shape = [], []
                if cS_obj.nk > 1:
                    points_list.append(np.array(cS_obj.kGridV).flatten())
                    data_shape.append(cS_obj.nk)
                if cS_obj.nkpps > 1:
                    points_list.append(np.array(cS_obj.kppsGridV).flatten())
                    data_shape.append(cS_obj.nkpps)
                
                points = tuple(points_list)
                expected_shape = tuple(data_shape)

                # 确保数据形状正确
                if k_prime_slice.shape != expected_shape: k_prime_slice = k_prime_slice.reshape(expected_shape)
                if c_pps_slice.shape != expected_shape: c_pps_slice = c_pps_slice.reshape(expected_shape)
                if c_slice.shape != expected_shape: c_slice = c_slice.reshape(expected_shape)

                # 创建插值器
                if not points: # 0D Case (nk=1, nkpps=1)
                    k_prime_interps[ie][ia] = lambda p: k_prime_slice.item()
                    c_pps_interps[ie][ia] = lambda p: c_pps_slice.item()
                    c_interps[ie][ia] = lambda p: c_slice.item()
                else: # 1D or 2D Case
                    k_prime_interps[ie][ia] = RegularGridInterpolator(points, k_prime_slice, method='linear', bounds_error=False, fill_value=None)
                    c_pps_interps[ie][ia] = RegularGridInterpolator(points, c_pps_slice, method='linear', bounds_error=False, fill_value=None)
                    c_interps[ie][ia] = RegularGridInterpolator(points, c_slice, method='linear', bounds_error=False, fill_value=None)

        # 3. 创建一个临时的简化版环境，仅用于状态演化（可选，但保持结构清晰）
        sim_env = OLGEnvV9SACSimp(
            cS=cS_obj,
            paramS_rl=vfi_results['paramS_vfi_dict'],
            M_fixed=M_fixed,
            training_mode=False
        )
        
        # 4. 初始化路径存储
        k_paths, kpps_paths, c_paths, cpps_paths = [], [], [], []
        
        # 5. 模拟生命周期
        for i_sim in range(n_sim):
            sim_env.current_age_idx = 1
            sim_env.current_k_val = cS_obj.kMin
            sim_env.current_k_pps_val = cS_obj.kppsMin
            k_path, kpps_path, c_path, cpps_path = [], [], [], []
            
            for age_idx in range(cS_obj.aD_new):
                sim_env.current_age_idx = age_idx + 1
                sim_env.current_eps_idx = eIdxM_group[i_sim, age_idx]
                
                k_path.append(sim_env.current_k_val)
                kpps_path.append(sim_env.current_k_pps_val)
                
                eps_vfi_idx = sim_env.current_eps_idx - 1
                
                # 准备查询点，并进行钳位以避免外插
                query_point_list = []
                if cS_obj.nk > 1:
                    query_point_list.append(np.clip(sim_env.current_k_val, cS_obj.kGridV[0], cS_obj.kGridV[-1]))
                if cS_obj.nkpps > 1:
                    query_point_list.append(np.clip(sim_env.current_k_pps_val, cS_obj.kppsGridV[0], cS_obj.kppsGridV[-1]))
                if not query_point_list: query_point_list.append(0)
                query_point = tuple(query_point_list)

                # [核心修复] 直接插值获取所有决策变量
                k_prime_vfi = k_prime_interps[eps_vfi_idx][age_idx](query_point).item()
                c_pps_vfi = c_pps_interps[eps_vfi_idx][age_idx](query_point).item()
                current_c = c_interps[eps_vfi_idx][age_idx](query_point).item()
                current_c = max(cS_obj.cFloor, current_c) # 确保满足最低消费

                # 记录决策和结果
                c_path.append(current_c)
                cpps_path.append(c_pps_vfi)
                
                # 手动更新下一期的资产状态 (物理过程)
                if age_idx < cS_obj.aD_new - 1:
                    # 更新非PPS资产
                    sim_env.current_k_val = max(cS_obj.kMin, min(cS_obj.kMax, k_prime_vfi))
                    
                    # 更新PPS资产
                    if cS_obj.pps_active:
                        pps_return_factor = 1 + ((sim_env.current_M['R_k_net_factor'] - 1) + cS_obj.pps_return_rate_premium)
                        pps_withdrawal = 0
                        if age_idx >= cS_obj.aR_new:
                            pps_withdrawal = sim_env.current_k_pps_val * cS_obj.pps_withdrawal_rate
                        next_k_pps = (sim_env.current_k_pps_val + c_pps_vfi - pps_withdrawal) * pps_return_factor
                        sim_env.current_k_pps_val = max(cS_obj.kppsMin, min(cS_obj.kppsMax, next_k_pps))
                    else:
                        sim_env.current_k_pps_val = cS_obj.kppsMin
            
            k_paths.append(k_path)
            kpps_paths.append(kpps_path)
            c_paths.append(c_path)
            cpps_paths.append(cpps_path)
                
        print("✅ VFI模拟完成。")
        
        return {
            "k_path_vfi": np.array(k_paths), "kpps_path_vfi": np.array(kpps_paths),
            "c_path_vfi": np.array(c_paths), "cpps_path_vfi": np.array(cpps_paths),
        }



    @staticmethod
    def solve_K_tau_l_for_rho_prime_vfi(rho_prime_payg_target_input, K_init_guess_input, cS_global, paramS_global_in, eIdxM_global_sim_paths):
        K_current_guess = K_init_guess_input
        tau_l_current_guess = cS_global.tau_l_init_guess
        
        # Pre-calculate some values
        mass_retirees_global = np.sum(paramS_global_in.ageMassV[cS_global.aR_new:])
        if paramS_global_in.mass_workers_group > 1e-9:
            theta_payg_required_calc = rho_prime_payg_target_input * (mass_retirees_global / paramS_global_in.mass_workers_group)
        else:
            theta_payg_required_calc = np.inf if rho_prime_payg_target_input > 1e-9 else 0
        theta_payg_required_calc = max(0, theta_payg_required_calc)
        
        solution_details_out = {'theta_payg_required_before_cap': theta_payg_required_calc}

        if theta_payg_required_calc > cS_global.theta_payg_max + 1e-5:
            if not getattr(paramS_global_in, 'suppress_initial_theta_print', False):
                print(f'  solve_K_tau_l (V8): rho_prime_target={rho_prime_payg_target_input:.4f} 导致理论theta_req={theta_payg_required_calc:.4f} > theta_max={cS_global.theta_payg_max:.3f}. 直接标记为不可行。')
            return K_init_guess_input, tau_l_current_guess, np.inf, False, solution_details_out

        stagnation_counter_ktl = 0
        prev_devNorm_ktl = np.inf
        tau_l_boundary_strike_count_ktl = 0

        if not getattr(paramS_global_in, 'suppress_inner_print_header', False):
            print(f'  solve_K_tau_l_for_rho_prime_vfi (V8): rho_prime_target={rho_prime_payg_target_input:.4f} (理论theta_req={theta_payg_required_calc:.4f}), K_init={K_current_guess:.2f}, tau_l_init={tau_l_current_guess:.3f}')
            print('  IterKTL | K_guess  | tau_l_gs | MPL_g    | theta_act| K_tot_mod| K_pps_mod| GBC_res  | K_dev    | tau_l_dev| Norm     | Improv   | Strikes  | Time (s) |')
            print('  ' + '-'*123)

        for iter_ktl_idx in range(cS_global.max_iter_K_tau_l):
            iter_timer_start = time.time()
            
            R_mkt_gross_factor, MPL_gross = OLG_V9_Utils.HHPrices_Huggett(K_current_guess, paramS_global_in.L_per_capita, cS_global)
            r_mkt_gross = R_mkt_gross_factor - 1
            
            avg_worker_gross_wage = (MPL_gross * paramS_global_in.L_per_capita) / paramS_global_in.mass_workers_group if paramS_global_in.mass_workers_group > 1e-9 else 0
            b_payg = max(0, rho_prime_payg_target_input * avg_worker_gross_wage)

            theta_payg_actual = min(theta_payg_required_calc, cS_global.theta_payg_max)
            if (theta_payg_actual + tau_l_current_guess) > cS_global.max_total_labor_tax:
                theta_payg_actual = max(0, cS_global.max_total_labor_tax - tau_l_current_guess)
            
            r_k_net_hh = r_mkt_gross * (1 - cS_global.tau_k)
            R_k_net_hh_factor = 1 + r_k_net_hh
            
            bV_payg_vec = np.zeros(cS_global.aD_new)
            bV_payg_vec[cS_global.aR_new:] = b_payg
            
            paramS_for_vfi = paramS_global_in
            paramS_for_vfi.tau_l = tau_l_current_guess
            paramS_for_vfi.theta_payg_actual_for_hh = theta_payg_actual
            paramS_for_vfi.pps_tax_deferral_active = cS_global.pps_active
            
            # Inner loop for bequest consistency
            T_bequest_model_iter = 0.01 * MPL_gross
            for _ in range(5): # Max 5 iterations for bequest
                cPolM, kPolM, cPpsPolM_choice, _ = OLG_V9_Utils.HHSolution_VFI_Huggett(
                    R_k_net_hh_factor, MPL_gross, T_bequest_model_iter, bV_payg_vec, paramS_for_vfi, cS_global
                )
                kHistM, kPpsHistM, _, _ = OLG_V9_Utils.HHSimulation_olgm(
                    kPolM, cPpsPolM_choice, cPolM, eIdxM_global_sim_paths, R_k_net_hh_factor, MPL_gross, T_bequest_model_iter, bV_payg_vec, paramS_for_vfi, cS_global
                )
                
                ageDeathMass_annual = paramS_global_in.Z_ss_norm_annual * cS_global.d_orig
                
                TotalNonPPSBequests_pc = np.dot(np.mean(kHistM, axis=0), ageDeathMass_annual)
                TotalPPSBequests_pc = 0
                if cS_global.pps_active and cS_global.pps_bequeathable and kPpsHistM.size > 0:
                    TotalPPSBequests_pc = np.dot(np.mean(kPpsHistM, axis=0), ageDeathMass_annual)
                
                TotalBequests_pc = TotalNonPPSBequests_pc + TotalPPSBequests_pc
                T_bequest_new = max(0, TotalBequests_pc / (1 + paramS_global_in.popGrowthForDebt))
                
                if abs(T_bequest_new - T_bequest_model_iter) < 1e-3 * (MPL_gross + 1e-9):
                    T_bequest_model_iter = T_bequest_new
                    break
                T_bequest_model_iter = 0.5 * T_bequest_model_iter + 0.5 * T_bequest_new
            
            # Final simulation with consistent bequests
            kHistM, kPpsHistM, cHistM, _ = OLG_V9_Utils.HHSimulation_olgm(
                kPolM, cPpsPolM_choice, cPolM, eIdxM_global_sim_paths, R_k_net_hh_factor, MPL_gross, T_bequest_model_iter, bV_payg_vec, paramS_for_vfi, cS_global
            )
            
            # 使用年龄组人口权重 paramS_global_in.ageMassV
            K_model_nonpps_sim = np.dot(np.mean(kHistM, axis=0), paramS_global_in.ageMassV)
            K_model_pps_sim = 0
            if cS_global.pps_active and cS_global.pps_in_K and kPpsHistM.size > 0:
                K_model_pps_sim = np.dot(np.mean(kPpsHistM, axis=0), paramS_global_in.ageMassV)
            C_model = np.dot(np.mean(cHistM, axis=0), paramS_global_in.ageMassV)
           
            K_model_from_sim = max(1e-6, K_model_nonpps_sim + K_model_pps_sim)

            
            Y_for_gbc = cS_global.A * (K_current_guess**cS_global.alpha) * (paramS_global_in.L_per_capita**(1-cS_global.alpha))
            G_target = cS_global.gov_exp_frac_Y * Y_for_gbc
            B_target = cS_global.gov_debt_frac_Y * Y_for_gbc

            gbc_residual = OLG_V9_Utils.check_gbc_residual(
                K_current_guess, C_model, Y_for_gbc, G_target, B_target, MPL_gross, r_mkt_gross,
                theta_payg_actual, tau_l_current_guess, b_payg, T_bequest_model_iter, 0, cS_global, paramS_global_in
            )

            K_dev = K_current_guess - K_model_from_sim
            tau_l_dev_raw = -gbc_residual / (MPL_gross * paramS_global_in.L_per_capita + 1e-9)
            current_devNorm = np.sqrt(K_dev**2 + gbc_residual**2)
            norm_improvement = prev_devNorm_ktl - current_devNorm
            elapsed_time = time.time() - iter_timer_start
            
            print(f'  {iter_ktl_idx+1:7d} | {K_current_guess:8.4f} | {tau_l_current_guess:8.4f} | {MPL_gross:8.4f} | {theta_payg_actual:8.4f} | {K_model_from_sim:8.4f} | {K_model_pps_sim:8.4f} | {gbc_residual:8.2e} | {K_dev:8.4f} | {tau_l_dev_raw:8.4f} | {current_devNorm:8.2e} | {norm_improvement:.1e} | {tau_l_boundary_strike_count_ktl:7d} | {elapsed_time:8.2f} |')

            payg_fully_funded_check = (theta_payg_actual >= theta_payg_required_calc - 1e-5)
            
            if current_devNorm < cS_global.tol_K_tau_l and abs(gbc_residual) < cS_global.gbc_tol_for_internal_loop and payg_fully_funded_check:
                print(f'  solve_K_tau_l (V8): K和tau_l成功收敛 (rho_prime_target={rho_prime_payg_target_input:.4f}, 实际theta_act={theta_payg_actual:.4f}).')
                solution_details_out.update({'R_mkt_gross_factor': R_mkt_gross_factor, 'MPL_gross': MPL_gross, 'theta_payg': theta_payg_actual, 'b_payg': b_payg, 'T_bequest_Model': T_bequest_model_iter, 'C_model': C_model, 'K_model_pps': K_model_pps_sim, 'K_model_non_pps': K_model_nonpps_sim})
                return K_model_from_sim, tau_l_current_guess, gbc_residual, True, solution_details_out

            # Update guesses
            K_current_guess = max(1e-3, K_current_guess - cS_global.damp_K_v5 * K_dev)
            tau_l_next_unconstrained = tau_l_current_guess + cS_global.damp_tau_l_v5 * tau_l_dev_raw
            tau_l_next_constrained = np.clip(tau_l_next_unconstrained, cS_global.tau_l_min, cS_global.tau_l_max)
            
            is_at_boundary = (abs(tau_l_next_constrained - cS_global.tau_l_max) < 1e-7 and tau_l_next_unconstrained >= cS_global.tau_l_max - 1e-7) or \
                             (abs(tau_l_next_constrained - cS_global.tau_l_min) < 1e-7 and tau_l_next_unconstrained <= cS_global.tau_l_min + 1e-7)

            if is_at_boundary and abs(gbc_residual) > cS_global.gbc_tol_for_internal_loop:
                tau_l_boundary_strike_count_ktl += 1
            else:
                tau_l_boundary_strike_count_ktl = 0
            
            tau_l_current_guess = tau_l_next_constrained
            
            if tau_l_boundary_strike_count_ktl >= cS_global.max_tau_l_boundary_strikes:
                print(f'  警告 (V8): tau_l 在边界 ({tau_l_current_guess:.4f}) 持续 {tau_l_boundary_strike_count_ktl} 次迭代，且GBC ({gbc_residual:.2e}) 未平衡。中止。')
                break
            
            if iter_ktl_idx > 0 and norm_improvement < (cS_global.min_norm_improvement_frac * prev_devNorm_ktl):
                stagnation_counter_ktl += 1
            else:
                stagnation_counter_ktl = 0
            
            prev_devNorm_ktl = current_devNorm
            
            if stagnation_counter_ktl >= cS_global.max_stagnation_iters:
                print(f'  警告 (V8): 在 {iter_ktl_idx+1} 次迭代后检测到范数停滞。中止。')
                break
        
        # If loop finishes without converging
        print(f'  警告 (V8): K和tau_l迭代达到最大次数 ({cS_global.max_iter_K_tau_l}) 或在该次数内未达可行解。')
        solution_details_out.update({'R_mkt_gross_factor': R_mkt_gross_factor, 'MPL_gross': MPL_gross, 'theta_payg': theta_payg_actual, 'b_payg': b_payg, 'T_bequest_Model': T_bequest_model_iter, 'C_model': C_model, 'K_model_pps': K_model_pps_sim, 'K_model_non_pps': K_model_nonpps_sim})
        return K_model_from_sim, tau_l_current_guess, gbc_residual, False, solution_details_out
       
    @staticmethod
    def solve_K_tau_l_for_rho_prime_rl(rho_prime_payg_target_input, K_init_guess_input, cS_global, paramS_global_in, eIdxM_global_sim_paths, rl_model, rl_config):
        """
        🤖 [最终修正版] 使用RL模型求解一般均衡（替代VFI）
        
        这个函数在结构上与VFI版本类似，但在内层循环中使用RL模型替代
        VFI策略函数来模拟家庭决策，从而找到市场的均衡 K 和 tau_l。
        
        修正点:
        - 修正了宏观总量聚合时的人口权重，使用 ageMassV (年龄组权重)。
        - 简化了不必要的意外遗赠循环，与基准RL模型设定对齐。
        
        Args:
            rho_prime_payg_target_input: 目标PAYG替代率
            K_init_guess_input: 初始资本存量猜测
            cS_global: 全局参数结构体
            paramS_global_in: 全局参数结构体（包含人口分布等）
            eIdxM_global_sim_paths: 全局效率冲击路径（剧本）
            rl_model: 训练好的SBX SAC模型
            rl_config: RL模型训练时的配置字典
            
        Returns:
            K_eq, tau_l_eq, gbc_residual_eq, eq_found, final_eq_solution_details
        """
        K_current_guess = K_init_guess_input
        tau_l_current_guess = cS_global.tau_l_init_guess
        
        # --- 初始检查和设置 ---
        mass_retirees_global = np.sum(paramS_global_in.ageMassV[cS_global.aR_new:])
        if paramS_global_in.mass_workers_group > 1e-9:
            theta_payg_required_calc = rho_prime_payg_target_input * (mass_retirees_global / paramS_global_in.mass_workers_group)
        else:
            theta_payg_required_calc = np.inf if rho_prime_payg_target_input > 1e-9 else 0
        theta_payg_required_calc = max(0, theta_payg_required_calc)
        
        solution_details_out = {'theta_payg_required_before_cap': theta_payg_required_calc}

        if theta_payg_required_calc > cS_global.theta_payg_max + 1e-5:
            if not getattr(paramS_global_in, 'suppress_initial_theta_print', False):
                print(f'  solve_K_tau_l (RL): rho_prime_target={rho_prime_payg_target_input:.4f} 导致理论theta_req={theta_payg_required_calc:.4f} > theta_max={cS_global.theta_payg_max:.3f}. 不可行。')
            return K_init_guess_input, tau_l_current_guess, np.inf, False, solution_details_out

        # --- 初始化迭代变量 ---
        stagnation_counter_ktl = 0
        prev_devNorm_ktl = np.inf
        tau_l_boundary_strike_count_ktl = 0

        if not getattr(paramS_global_in, 'suppress_inner_print_header', False):
            print(f'  solve_K_tau_l_for_rho_prime_rl: rho_prime_target={rho_prime_payg_target_input:.4f} (理论theta_req={theta_payg_required_calc:.4f})')
            print('  IterKTL | K_guess  | tau_l_gs | MPL_g    | theta_act| K_tot_mod| K_pps_mod| GBC_res  | K_dev    | tau_l_dev| Norm     | Improv   | Strikes  | Time (s) |')
            print('  ' + '-'*123)

        # --- 主迭代循环 ---
        for iter_ktl_idx in range(cS_global.max_iter_K_tau_l):
            iter_timer_start = time.time()
            
            # a. 计算要素价格和PAYG福利
            R_mkt_gross_factor, MPL_gross = OLG_V9_Utils.HHPrices_Huggett(K_current_guess, paramS_global_in.L_per_capita, cS_global)
            r_mkt_gross = R_mkt_gross_factor - 1
            
            avg_worker_gross_wage = (MPL_gross * paramS_global_in.L_per_capita) / paramS_global_in.mass_workers_group if paramS_global_in.mass_workers_group > 1e-9 else 0
            b_payg = max(0, rho_prime_payg_target_input * avg_worker_gross_wage)

            theta_payg_actual = min(theta_payg_required_calc, cS_global.theta_payg_max)
            if (theta_payg_actual + tau_l_current_guess) > cS_global.max_total_labor_tax:
                theta_payg_actual = max(0, cS_global.max_total_labor_tax - tau_l_current_guess)
            
            r_k_net_hh = r_mkt_gross * (1 - cS_global.tau_k)
            R_k_net_factor_hh = 1 + r_k_net_hh
            
            bV_payg_vec = np.zeros(cS_global.aD_new)
            bV_payg_vec[cS_global.aR_new:] = b_payg
            
            # b. 准备RL模拟所需的参数
            paramS_for_rl_sim = type('ParameterStruct', (), {})()
            paramS_for_rl_sim.tau_l = tau_l_current_guess
            paramS_for_rl_sim.theta_payg_actual_for_hh = theta_payg_actual

            # c. 关键步骤: 调用RL模拟器替代VFI求解器
            #    TR_total 设为0，因为在此均衡设定中无意外遗赠
            kHistM, kPpsHistM, cHistM, _ = OLG_V9_Utils.HHSimulation_olgm_rl(
                rl_model, rl_config, eIdxM_global_sim_paths,
                R_k_net_factor_hh, MPL_gross, 0.0,
                bV_payg_vec, paramS_for_rl_sim, cS_global
            )
            
            # d. [核心修正] 计算宏观总量，使用与微观数据维度一致的年龄组人口权重
            K_model_nonpps_sim = np.dot(np.mean(kHistM, axis=0), paramS_global_in.ageMassV)
            K_model_pps_sim = 0
            if cS_global.pps_active and cS_global.pps_in_K and kPpsHistM.size > 0:
                K_model_pps_sim = np.dot(np.mean(kPpsHistM, axis=0), paramS_global_in.ageMassV)
            
            K_model_from_sim = max(1e-6, K_model_nonpps_sim + K_model_pps_sim)
            C_model = np.dot(np.mean(cHistM, axis=0), paramS_global_in.ageMassV)
            
            # e. 检查政府预算约束 (GBC)
            Y_for_gbc = cS_global.A * (K_current_guess**cS_global.alpha) * (paramS_global_in.L_per_capita**(1-cS_global.alpha))
            G_target = cS_global.gov_exp_frac_Y * Y_for_gbc
            B_target = cS_global.gov_debt_frac_Y * Y_for_gbc

            paramS_for_gbc = type('ParameterStruct', (), {})()
            paramS_for_gbc.L_per_capita = paramS_global_in.L_per_capita
            paramS_for_gbc.popGrowthForDebt = paramS_global_in.popGrowthForDebt
            gbc_residual = OLG_V9_Utils.check_gbc_residual(
                K_current_guess, C_model, Y_for_gbc, G_target, B_target, MPL_gross, r_mkt_gross,
                theta_payg_actual, tau_l_current_guess, b_payg, 0, 0, cS_global, paramS_for_gbc
            )

            # f. 计算偏差并准备下一次迭代
            K_dev = K_current_guess - K_model_from_sim
            tau_l_dev_raw = -gbc_residual / (MPL_gross * paramS_global_in.L_per_capita + 1e-9)
            current_devNorm = np.sqrt(K_dev**2 + gbc_residual**2)
            norm_improvement = prev_devNorm_ktl - current_devNorm
            elapsed_time = time.time() - iter_timer_start
            
            print(f'  {iter_ktl_idx+1:7d} | {K_current_guess:8.4f} | {tau_l_current_guess:8.4f} | {MPL_gross:8.4f} | {theta_payg_actual:8.4f} | {K_model_from_sim:8.4f} | {K_model_pps_sim:8.4f} | {gbc_residual:8.2e} | {K_dev:8.4f} | {tau_l_dev_raw:8.4f} | {current_devNorm:8.2e} | {norm_improvement:.1e} | {tau_l_boundary_strike_count_ktl:7d} | {elapsed_time:8.2f} |')

            # g. 收敛检查
            payg_fully_funded_check = (theta_payg_actual >= theta_payg_required_calc - 1e-5)
            if current_devNorm < cS_global.tol_K_tau_l and abs(gbc_residual) < cS_global.gbc_tol_for_internal_loop and payg_fully_funded_check:
                print(f'  solve_K_tau_l (RL): K和tau_l成功收敛。')
                solution_details_out.update({'R_mkt_gross_factor': R_mkt_gross_factor, 'MPL_gross': MPL_gross, 'theta_payg': theta_payg_actual, 'b_payg': b_payg, 'T_bequest_Model': 0.0, 'C_model': C_model, 'K_model_pps': K_model_pps_sim, 'K_model_non_pps': K_model_nonpps_sim})
                return K_model_from_sim, tau_l_current_guess, gbc_residual, True, solution_details_out

            # h. 更新猜测值
            K_current_guess = max(1e-3, K_current_guess - cS_global.damp_K_v5 * K_dev)
            tau_l_next_unconstrained = tau_l_current_guess + cS_global.damp_tau_l_v5 * tau_l_dev_raw
            tau_l_current_guess = np.clip(tau_l_next_unconstrained, cS_global.tau_l_min, cS_global.tau_l_max)
            
            # i. 检查停滞和边界撞击
            is_at_boundary = (abs(tau_l_current_guess - cS_global.tau_l_max) < 1e-7 and tau_l_next_unconstrained >= cS_global.tau_l_max - 1e-7) or \
                             (abs(tau_l_current_guess - cS_global.tau_l_min) < 1e-7 and tau_l_next_unconstrained <= cS_global.tau_l_min + 1e-7)

            if is_at_boundary and abs(gbc_residual) > cS_global.gbc_tol_for_internal_loop:
                tau_l_boundary_strike_count_ktl += 1
            else: tau_l_boundary_strike_count_ktl = 0
            
            if tau_l_boundary_strike_count_ktl >= cS_global.max_tau_l_boundary_strikes:
                print(f'  警告 (RL): tau_l 在边界持续撞击，且GBC未平衡。中止。')
                break
            
            if iter_ktl_idx > 0 and norm_improvement < (cS_global.min_norm_improvement_frac * prev_devNorm_ktl):
                stagnation_counter_ktl += 1
            else: stagnation_counter_ktl = 0
            
            prev_devNorm_ktl = current_devNorm
            
            if stagnation_counter_ktl >= cS_global.max_stagnation_iters:
                print(f'  警告 (RL): 检测到范数停滞。中止。')
                break
        
        # --- 循环结束但未收敛 ---
        print(f'  警告 (RL): K和tau_l迭代达到最大次数或未达可行解。')
        solution_details_out.update({'R_mkt_gross_factor': R_mkt_gross_factor, 'MPL_gross': MPL_gross, 'theta_payg': theta_payg_actual, 'b_payg': b_payg, 'T_bequest_Model': 0.0, 'C_model': C_model, 'K_model_pps': K_model_pps_sim, 'K_model_non_pps': K_model_nonpps_sim})
        return K_model_from_sim, tau_l_current_guess, gbc_residual, False, solution_details_out

    @staticmethod
    def check_gbc_residual(K_val_market_input, C_val_model_input, Y_val_market_input, G_val_target_input, B_val_target_input, MPL_gross_val_input, r_mkt_gross_val_input, theta_payg_val_actual_input, tau_l_val_input, b_payg_val_per_retiree_input, T_bequest_val_pc_input, TR_gov_val_pc_input, cS_check, paramS_loc_check):
        L_per_capita_local_check = paramS_loc_check.L_per_capita
        
        LaborTaxRev_general_part_calc = tau_l_val_input * MPL_gross_val_input * L_per_capita_local_check
        CapitalTaxRev_calc = r_mkt_gross_val_input * K_val_market_input * cS_check.tau_k
        ConsumptionTaxRev_calc = C_val_model_input * cS_check.tau_c
        GeneralRevenue_calc = LaborTaxRev_general_part_calc + CapitalTaxRev_calc + ConsumptionTaxRev_calc
        
        GovConsumption_calc = G_val_target_input
        r_b_for_debt_service_calc = r_mkt_gross_val_input
        DebtService_calc = (r_b_for_debt_service_calc - paramS_loc_check.popGrowthForDebt) * B_val_target_input
        GovDirectTransfers_calc = TR_gov_val_pc_input
        
        GeneralOutlays_calc = GovConsumption_calc + DebtService_calc + GovDirectTransfers_calc
        
        gbc_residual_out = GeneralRevenue_calc - GeneralOutlays_calc
        return gbc_residual_out

    @staticmethod
    def diagnose_rl_vfi_consistency(cS, verbose=True):
        """
        🔍 诊断RL和VFI模拟的逻辑一致性
        
        检查关键的年龄判断、PPS缴费提取逻辑、折现因子等是否在RL和VFI之间保持一致
        
        Args:
            cS: 模型参数
            verbose: 是否输出详细信息
            
        Returns:
            dict: 诊断结果
        """
        if verbose:
            print("🔍 RL vs VFI 逻辑一致性诊断")
            print("=" * 50)
        
        issues = []
        
        # 1. 检查年龄相关参数
        if verbose:
            print(f"📊 年龄参数:")
            print(f"  年度年龄范围: {cS.age1_orig} - {cS.ageLast_orig} 岁")
            print(f"  年度退休年龄: {cS.ageRetire_orig} 岁 (索引 {cS.aR_idx_orig})")
            print(f"  年龄组数量: {cS.aD_new} (每组 {cS.yearStep} 年)")
            print(f"  组别退休年龄组: {cS.aR_new}")
        
        # 2. 检查PPS相关参数
        if hasattr(cS, 'pps_active') and cS.pps_active:
            if verbose:
                print(f"📊 PPS参数:")
                print(f"  PPS激活: {cS.pps_active}")
                print(f"  PPS缴费最大年龄索引: {cS.pps_contribution_age_max_idx}")
                print(f"  PPS提取最小年龄索引: {cS.pps_withdrawal_age_min_idx}")
                print(f"  PPS提取率: {cS.pps_withdrawal_rate:.3f}")
                print(f"  PPS提取税率: {cS.pps_tax_rate_withdrawal:.3f}")
                print(f"  PPS收益率溢价: {cS.pps_return_rate_premium:.3f}")
        
        # 3. 检查年龄映射一致性
        age_mapping_consistent = True
        for a_group in range(cS.aD_new):
            if a_group < len(cS.physAgeMap):
                annual_ages = cS.physAgeMap[a_group]
                if annual_ages:
                    first_annual_age = annual_ages[0]
                    # 检查退休判断一致性
                    group_is_working = (a_group < cS.aR_new)
                    annual_is_working = (first_annual_age + 1 < cS.aR_idx_orig)  # 🔧 修复：年度年龄判断 (first_annual_age是0-based，aR_idx_orig是1-based)
                    if group_is_working != annual_is_working:
                        age_mapping_consistent = False
                        issues.append(f"年龄组 {a_group} 与年度年龄 {first_annual_age} 的工作状态判断不一致")
        
        if verbose:
            print(f"📊 年龄映射一致性: {'✅ 一致' if age_mapping_consistent else '❌ 不一致'}")
        
        # 4. 检查折现因子
        beta = getattr(cS, 'beta', 0.97)
        if verbose:
            print(f"📊 折现因子: β = {beta:.4f}")
        
        # 5. 检查效用函数参数
        sigma = getattr(cS, 'sigma', 1.5)
        if verbose:
            print(f"📊 风险厌恶系数: σ = {sigma:.3f}")
        
        # 6. 总结
        if verbose:
            print("=" * 50)
            if len(issues) == 0:
                print("✅ 未发现逻辑不一致问题")
            else:
                print(f"❌ 发现 {len(issues)} 个潜在问题:")
                for i, issue in enumerate(issues, 1):
                    print(f"  {i}. {issue}")
        
        return {
            'age_mapping_consistent': age_mapping_consistent,
            'issues': issues,
            'beta': beta,
            'sigma': sigma,
            'aR_idx_orig': cS.aR_idx_orig,
            'aR_new': cS.aR_new,
            'pps_active': getattr(cS, 'pps_active', False)
        }

    @staticmethod
    def GeneralEquilibrium(kHistM, kPpsHistM, cHistM, cppsHistM, eIdxM_input, cS, paramS, Z_ss_norm, rng_M, verbose=False):
        """
        GeneralEquilibrium - 一般均衡计算（适配年龄组模拟）
        
        Args:
            kHistM, kPpsHistM, cHistM, cppsHistM: 生命周期路径矩阵
            eIdxM_input: 效率冲击矩阵，可以是年度或年龄组数据
            其他参数保持不变
        """
        


        print('处理年龄组模拟结果...')
        # 年龄组数据，直接使用
        kHistM_group = kHistM
        kPpsHistM_group = kPpsHistM
        cHistM_group = cHistM
        cppsHistM_group = cppsHistM


        # 计算劳动供给（现在使用年龄组版本）
        HHlaborM_group, L_total_eff_pc = OLG_V9_Utils.LaborSupply_Huggett(
            eIdxM_input, cS, paramS, Z_ss_norm
        )

        # 计算平均资产
        K_avg_HH = np.mean(kHistM_group[:, :-1])  # 排除最后一个年龄组（死亡）
        print(f'💰 平均非PPS资产 (K_avg): {K_avg_HH:.4f}')

        # 计算平均PPS资产
        K_pps_avg_HH = np.mean(kPpsHistM_group[:, :-1])
        print(f'🏦 平均PPS资产 (K_pps_avg): {K_pps_avg_HH:.4f}')

        # 计算平均消费
        C_avg_HH = np.mean(cHistM_group)
        print(f'🛒 平均消费 (C_avg): {C_avg_HH:.4f}')

        # 计算平均PPS缴费（仅工作年龄组）
        if cS.aR_new > 0:
            C_pps_avg_HH = np.mean(cppsHistM_group[:, :cS.aR_new])
        else:
            C_pps_avg_HH = 0
        print(f'💳 平均PPS缴费 (C_pps_avg): {C_pps_avg_HH:.4f}')

        # 计算总体经济指标
        K_total = K_avg_HH  # 总资本存量
        L_total = L_total_eff_pc  # 总有效劳动供给
        
        # 生产函数
        alpha = cS.alpha
        delta = cS.delta
        Y_total = K_total**alpha * L_total**(1-alpha)  # 总产出
        
        # 要素价格
        R_gross = alpha * Y_total / K_total if K_total > 0 else 1 + cS.r_ss
        w_gross = (1-alpha) * Y_total / L_total if L_total > 0 else cS.w_ss
        R_net = R_gross - delta  # 净利率

        if verbose:
            print(f'📊 一般均衡结果（年龄组版本）:')
            print(f'   总产出 (Y): {Y_total:.4f}')
            print(f'   总资本 (K): {K_total:.4f}')
            print(f'   总劳动 (L): {L_total:.4f}')
            print(f'   毛利率 (R_gross): {R_gross:.4f}')
            print(f'   净利率 (R_net): {R_net:.4f}')
            print(f'   工资率 (w): {w_gross:.4f}')

        # 构建结果字典
        ge_results = {
            'K_avg': K_avg_HH,
            'K_pps_avg': K_pps_avg_HH,
            'C_avg': C_avg_HH,
            'C_pps_avg': C_pps_avg_HH,
            'L_eff_pc': L_total_eff_pc,
            'Y_total': Y_total,
            'R_gross': R_gross,
            'R_net': R_net,
            'w_gross': w_gross,
            'kHistM_group': kHistM_group,
            'kPpsHistM_group': kPpsHistM_group,
            'cHistM_group': cHistM_group,
            'cppsHistM_group': cppsHistM_group,
            'HHlaborM_group': HHlaborM_group
        }

        return ge_results




# ==============================================================================
# === OLGEnvV9SAC with Exponential Utility for Reward Shaping ===
# ==============================================================================

import numpy as np
from typing import Dict, Any, Tuple, List, Optional, Union
import gymnasium as gym
from gymnasium import spaces


# --- 在 main_olg_v9_utils.py 中，替换 OLGEnvV9SAC 类 ---

class OLGEnvV9SAC(gym.Env):
    """
    [重构版] OLG 完整环境，用于 SAC Agent 训练。
    
    此版本处理动态的宏观经济变量，是进行一般均衡分析的基础。
    观测空间包含个体状态和宏观状态。
    """
    metadata = {'render.modes': ['human']}

    def __init__(self, cS: Any, paramS_rl: Dict[str, Any], rng_M: Dict[str, Any], 
                 training_mode: bool = True, 
                 reward_shaping_scheme: str = 'exponential',
                 reward_alpha: float = 0.1):
        """
        初始化完整版环境

        Args:
            cS: OLG模型参数 (作为对象)
            paramS_rl: RL相关的派生参数
            rng_M: 宏观变量的采样范围
            training_mode: bool, 是否为训练模式 (True: 加权奖励, False: 纯效用奖励)
            reward_shaping_scheme: str, 使用的奖励方案。
            reward_alpha: float, 指数化效用中的alpha参数。
        """
        if gym is None or spaces is None:
            raise ImportError("需要安装gymnasium库：pip install gymnasium")

        super().__init__()

        self.cS = cS
        self.paramS_rl = paramS_rl
        self.rng_M = rng_M
        self.training_mode = training_mode
        self.reward_shaping_scheme = reward_shaping_scheme
        self.reward_alpha = reward_alpha
        self.cumulative_survival_prob = 1.0

        # [核心] 定义完整的观测和动作空间
        self.obs_dim_full = 10  # k, k_pps, age_idx, eps_idx, M_vars(6)
        self.observation_space = spaces.Box(
            low=0.0, high=1.0, shape=(self.obs_dim_full,), dtype=np.float32
        )

        act_dim = 2  # [prop_pps_contrib, prop_consump] -> 已统一为消费比例
        self.action_space = spaces.Box(
            low=0.0, high=1.0, shape=(act_dim,), dtype=np.float32
        )

        # 初始化当前状态
        self.current_M = {}
        self.current_bV_payg = np.zeros(self.cS.aD_new)
        self.current_age_idx = 1
        self.current_k_val = self.cS.kMin
        self.current_k_pps_val = self.cS.kppsMin
        self.current_eps_idx = 1

        # 初始化完整的归一化参数
        self._init_normalization_params()
        
        # 打印初始化信息
        mode = "训练模式" if training_mode else "评估模式"
        # print(f"🏋️ OLG完整环境(OLGEnvV9SAC)初始化 - {mode}")

    def _init_normalization_params(self):
        """初始化完整的观测归一化参数"""
        self.obs_norm_min = np.array([
            self.cS.kMin, self.cS.kppsMin, 1.0, 1.0,
            self.rng_M['R_k_net_factor'][0], self.rng_M['w_gross'][0],
            self.rng_M['TR_total'][0], self.rng_M['b_payg_avg_retiree'][0],
            self.rng_M['tau_l'][0], self.rng_M['theta_payg_actual'][0]
        ])
        self.obs_norm_max = np.array([
            self.cS.kMax, self.cS.kppsMax, float(self.cS.aD_new), float(self.cS.nw),
            self.rng_M['R_k_net_factor'][1], self.rng_M['w_gross'][1],
            self.rng_M['TR_total'][1], self.rng_M['b_payg_avg_retiree'][1],
            self.rng_M['tau_l'][1], self.rng_M['theta_payg_actual'][1]
        ])
        self.obs_norm_range = self.obs_norm_max - self.obs_norm_min
        self.obs_norm_range[self.obs_norm_range < 1e-6] = 1.0
    
    def _get_observation(self) -> np.ndarray:
        """获取当前观测（并归一化）"""
        raw_obs_vec = np.array([
            self.current_k_val, self.current_k_pps_val,
            float(self.current_age_idx), float(self.current_eps_idx),
            self.current_M['R_k_net_factor'], self.current_M['w_gross'],
            self.current_M['TR_total'], self.current_M.get('b_payg_avg_for_obs', 0.0), # 使用.get以防万一
            self.current_M['tau_l'], self.current_M['theta_payg_actual']
        ])
        
        # 归一化
        obs = (raw_obs_vec - self.obs_norm_min) / self.obs_norm_range
        obs = np.clip(obs, 0, 1)
        return obs.astype(np.float32)

    def _sample_macro_vars(self):
        """在每轮开始时采样宏观变量"""
        self.current_M['R_k_net_factor'] = self.np_random.uniform(*self.rng_M['R_k_net_factor'])
        self.current_M['w_gross'] = self.np_random.uniform(*self.rng_M['w_gross'])
        self.current_M['TR_total'] = self.np_random.uniform(*self.rng_M['TR_total'])
        b_payg_avg = self.np_random.uniform(*self.rng_M['b_payg_avg_retiree'])
        self.current_M['tau_l'] = self.np_random.uniform(*self.rng_M['tau_l'])
        self.current_M['theta_payg_actual'] = self.np_random.uniform(*self.rng_M['theta_payg_actual'])
        
        # 更新PAYG福利向量和用于观测的平均值
        self.current_bV_payg.fill(0)
        if self.cS.aR_new < self.cS.aD_new:
            self.current_bV_payg[self.cS.aR_new:] = b_payg_avg
        self.current_M['b_payg_avg_for_obs'] = b_payg_avg
        
    def reset(self, seed: Optional[int] = None, options: Optional[Dict] = None) -> Tuple[np.ndarray, Dict]:
        """ 重置环境 """
        super().reset(seed=seed)

        # 采样宏观经济状态 M
        self._sample_macro_vars()

        # 初始化个体状态
        self.current_age_idx = 1
        self.current_k_val = self.cS.kMin
        self.current_k_pps_val = self.cS.kppsMin
        self.cumulative_survival_prob = 1.0
        
        # 初始效率冲击
        self.current_eps_idx = self.np_random.choice(
            self.cS.nw, p=self.paramS_rl['leProb1V']
        ) + 1

        observation = self._get_observation()
        info = {'age_idx': self.current_age_idx}
        return observation, info

    def step(self, action: np.ndarray) -> Tuple[np.ndarray, float, bool, bool, Dict]:
        """ 执行一步 """
        # [核心] 这里的逻辑现在是基于“消费比例”的动作空间
        prop_pps_contrib, prop_consump = np.clip(action, 0, 1)

        # 1. 计算决策变量
        actual_c_pps, _ = self._calculate_pps_contribution(prop_pps_contrib)
        resources_after_pps = self._calculate_resources_after_pps(actual_c_pps)
        actual_k_prime, current_c = self._calculate_consumption_and_savings(
            resources_after_pps, prop_consump
        )
        
        # 2. 计算纯效用 u(c)
        _, pure_utility = OLG_V9_Utils.CES_utility(current_c, self.cS.sigma, self.cS)
        if not np.isfinite(pure_utility):
            pure_utility = -1000.0

        # 3. 获取存活概率
        vfi_age_idx = self.current_age_idx - 1
        survival_prob = self.cS.s_1yr_transitionV[vfi_age_idx] if 0 <= vfi_age_idx < self.cS.aD_new else 0.0

        # 4. 根据模式计算奖励
        if self.training_mode:
            if self.reward_shaping_scheme == 'exponential':
                reward = np.exp(self.reward_alpha * pure_utility) * self.cumulative_survival_prob
            else:
                reward = pure_utility * self.cumulative_survival_prob
        else: 
            reward = pure_utility

        # 5. 更新状态
        terminated = self._update_state(actual_k_prime, actual_c_pps)
        
        if not terminated:
            self.cumulative_survival_prob *= survival_prob

        observation = self._get_observation()

        info = {
            "pure_utility": pure_utility,
            "consumption": current_c,
            "k_prime": actual_k_prime,
            "c_pps": actual_c_pps,
        }

        truncated = False
        return observation, reward, terminated, truncated, info

    # --- 以下是共享的物理过程辅助函数 ---
    
    def _calculate_pps_contribution(self, prop_pps_contrib: float) -> Tuple[float, float]:
        # ... (此函数无需修改, 已在之前版本中提供) ...
        # ... 保持与 OLGEnvV9SACSimp 的父类版本一致 ...
        actual_c_pps = 0.0
        max_permissible_cpps = 0.0
        current_epsilon_val = self.paramS_rl['leGridV'][self.current_eps_idx - 1]
        is_working_age_group = (self.current_age_idx - 1 < self.cS.aR_new)
        is_pps_contribution_eligible = (is_working_age_group and self.cS.pps_active)
        if is_pps_contribution_eligible:
            age_efficiency = self.cS.ageEffV_new[self.current_age_idx - 1]
            current_gross_labor_income = self.current_M['w_gross'] * age_efficiency * current_epsilon_val
            if current_gross_labor_income > 1e-6:
                max_cpps_by_frac = current_gross_labor_income * self.cS.pps_max_contrib_frac
                max_permissible_cpps = min(self.cS.pps_contrib_limit, max_cpps_by_frac)
                max_permissible_cpps = max(0, max_permissible_cpps)
            actual_c_pps = prop_pps_contrib * max_permissible_cpps
            actual_c_pps = max(0, min(actual_c_pps, max_permissible_cpps))
        return actual_c_pps, max_permissible_cpps


    def _calculate_resources_after_pps(self, actual_c_pps: float) -> float:
        # ... (此函数无需修改, 已在之前版本中提供) ...
        # ... 保持与 OLGEnvV9SACSimp 的父类版本一致 ...
        paramS_hh_step = TempParamSHH(self.current_M['tau_l'], self.current_M['theta_payg_actual'], self.cS.pps_active, self.cS.ageEffV_new)
        is_working_age_group = (self.current_age_idx - 1 < self.cS.aR_new)
        is_pps_withdrawal_eligible = (not is_working_age_group and self.cS.pps_active)
        if is_pps_withdrawal_eligible:
            paramS_hh_step.current_pps_withdrawal = self.current_k_pps_val * self.cS.pps_withdrawal_rate
        else:
            paramS_hh_step.current_pps_withdrawal = 0
        b_payg_this_age = self.current_bV_payg[self.current_age_idx - 1]
        current_epsilon_val = self.paramS_rl['leGridV'][self.current_eps_idx - 1]
        resources_after_pps, _, _ = OLG_V9_Utils.HHIncome_Huggett(self.current_k_val, self.current_M['R_k_net_factor'], self.current_M['w_gross'], self.current_M['TR_total'], b_payg_this_age, actual_c_pps, self.current_age_idx - 1, paramS_hh_step, self.cS, current_epsilon_val)
        return resources_after_pps

    def _calculate_consumption_and_savings(self, resources_after_pps: float, prop_consump: float) -> Tuple[float, float]:
        # [核心] 基于消费比例的决策逻辑
        consumption_floor_spending = self.cS.cFloor * (1 + self.cS.tau_c)
        if resources_after_pps <= consumption_floor_spending:
            total_consumption_spending = resources_after_pps
            actual_k_prime = self.cS.kMin
        else:
            spendable_resources = resources_after_pps - consumption_floor_spending
            consumption_above_floor_spending = prop_consump * spendable_resources
            total_consumption_spending = consumption_floor_spending + consumption_above_floor_spending
            actual_k_prime = resources_after_pps - total_consumption_spending
        current_c = total_consumption_spending / (1 + self.cS.tau_c)
        actual_k_prime = max(self.cS.kMin, min(actual_k_prime, self.cS.kMax))
        current_c = max(self.cS.cFloor, current_c)
        return actual_k_prime, current_c

    def _update_state(self, actual_k_prime: float, actual_c_pps: float) -> bool:
        # ... (此函数无需修改, 已在之前版本中提供) ...
        # ... 保持与 OLGEnvV9SACSimp 的父类版本一致 ...
        self.current_k_val = actual_k_prime
        if self.cS.pps_active:
            is_working_age_group = (self.current_age_idx - 1 < self.cS.aR_new)
            is_pps_withdrawal_eligible = (not is_working_age_group and self.cS.pps_active)
            pps_withdrawal = 0
            if is_pps_withdrawal_eligible: pps_withdrawal = self.current_k_pps_val * self.cS.pps_withdrawal_rate
            pps_return_factor = 1 + ((self.current_M['R_k_net_factor'] - 1) + self.cS.pps_return_rate_premium)
            k_pps_next_unclamped = (self.current_k_pps_val + actual_c_pps - pps_withdrawal) * pps_return_factor
            self.current_k_pps_val = max(self.cS.kppsMin, min(self.cS.kppsMax, k_pps_next_unclamped))
        else: self.current_k_pps_val = self.cS.kppsMin
        terminated = False
        if self.current_age_idx < self.cS.aD_new:
            trans_probs = self.paramS_rl['leTrProbM'][self.current_eps_idx - 1, :]
            self.current_eps_idx = self.np_random.choice(len(trans_probs), p=trans_probs) + 1
            self.current_age_idx += 1
        else: terminated = True
        return terminated

# ==============================================================================
# === OLGEnvV9SACSimp: Simplified OLG Environment with Fixed Macro-variables ===
# ==============================================================================

# ==============================================================================
# === OLGEnvV9SACSimp: Simplified Env with Consumption-based Action Space ===
# ==============================================================================

import numpy as np
from typing import Dict, Any, Tuple, Optional
import gymnasium as gym
from gymnasium import spaces

# 导入基类和工具函数

class OLGEnvV9SACSimp(OLGEnvV9SAC):
    """
    [修正版] 一个简化的 OLG 环境，具有固定的宏观变量。
    
    此版本继承自 OLGEnvV9SAC，并通过以下方式进行简化：
    1. 在构造时接收固定的宏观变量 M_fixed。
    2. 重写观测空间为4维 [k, k_pps, age, ε]。
    3. 重写 _get_observation 和归一化方法以匹配4维观测。
    4. 重写 reset 方法，使其不再采样宏观变量。
    """
    
    def __init__(self, cS: Any, paramS_rl: Dict[str, Any], 
                 M_fixed: Dict[str, float],
                 **kwargs):
        """
        初始化简化版环境
        
        Args:
            cS: 模型参数对象
            paramS_rl: RL相关参数
            M_fixed: 固定的宏观变量字典
            **kwargs: 传递给父类构造函数的其他参数
        """
        # [核心修正] 1. 使用固定的 M_fixed 构建一个“假的” rng_M
        #    这个假的 rng_M 仅用于满足父类构造函数对范围的需求，其值不重要。
        rng_M_dummy = {key: [val, val] for key, val in M_fixed.items()}
        required_keys = ['R_k_net_factor', 'w_gross', 'TR_total', 'b_payg_avg_retiree', 'tau_l', 'theta_payg_actual']
        for key in required_keys:
            if key not in rng_M_dummy:
                rng_M_dummy[key] = [0.0, 1.0] # 提供一个默认范围

        # [核心修正] 2. 显式调用父类的构造函数
        #    这将初始化所有父类的方法和属性，包括 set_macro_parameters
        super().__init__(cS, paramS_rl, rng_M_dummy, **kwargs)

        # 3. 存储固定的宏观参数并立即设置
        self.M_fixed = M_fixed
        self._set_fixed_macro_parameters() # 使用一个内部方法来设置

        # 4. [重写] 定义简化的观测空间
        obs_dim = 4 # [k, k_pps, age_idx, eps_idx]
        self.observation_space = spaces.Box(
            low=0.0, high=1.0, shape=(obs_dim,), dtype=np.float32
        )
        
        # 5. [重写] 初始化简化的归一化参数
        self._init_simplified_normalization_params()
        
        # print("✅ 简化版OLG环境(OLGEnvV9SACSimp)初始化完成")

    def _set_fixed_macro_parameters(self):
        """一个内部方法，用于设置和更新固定的宏观参数"""
        self.current_M = self.M_fixed.copy()
        b_payg_avg = self.M_fixed.get('b_payg_avg_retiree', 0.0)
        self.current_bV_payg.fill(0)
        if self.cS.aR_new < self.cS.aD_new:
            self.current_bV_payg[self.cS.aR_new:] = b_payg_avg
        # 确保用于观测的键存在
        self.current_M['b_payg_avg_for_obs'] = b_payg_avg


    def _init_simplified_normalization_params(self):
        """[重写] 初始化简化的观测归一化参数"""
        self.obs_norm_min_simp = np.array([
            self.cS.kMin, self.cS.kppsMin, 1.0, 1.0
        ])
        self.obs_norm_max_simp = np.array([
            self.cS.kMax, self.cS.kppsMax, float(self.cS.aD_new), float(self.cS.nw)
        ])
        self.obs_norm_range_simp = self.obs_norm_max_simp - self.obs_norm_min_simp
        self.obs_norm_range_simp[self.obs_norm_range_simp < 1e-6] = 1.0

    def reset(self, seed: Optional[int] = None, options: Optional[Dict] = None) -> Tuple[np.ndarray, Dict]:
        """
        [重写] 重置简化版环境。
        关键区别：不再采样宏观变量，而是使用固定的 M_fixed。
        """
        super(OLGEnvV9SAC, self).reset(seed=seed) # 调用祖父类的reset来处理seed

        # 确保宏观参数是固定的
        self._set_fixed_macro_parameters()

        # 初始化个体状态
        self.current_age_idx = 1
        self.current_k_val = self.cS.kMin
        self.current_k_pps_val = self.cS.kppsMin
        self.cumulative_survival_prob = 1.0
        
        # 初始效率冲击
        self.current_eps_idx = self.np_random.choice(
            self.cS.nw, p=self.paramS_rl['leProb1V']
        ) + 1
        
        observation = self._get_observation()
        info = {'age_idx': self.current_age_idx}
        return observation, info

    def _get_observation(self) -> np.ndarray:
        """[重写] 获取简化的当前观测，并使用简化的归一化"""
        raw_obs_vec = np.array([
            self.current_k_val, self.current_k_pps_val,
            float(self.current_age_idx), float(self.current_eps_idx)
        ])
        
        # 使用简化的归一化参数
        obs = (raw_obs_vec - self.obs_norm_min_simp) / self.obs_norm_range_simp
        obs = np.clip(obs, 0, 1)
        return obs.astype(np.float32)

    # step 方法不需要重写，因为它依赖的 _calculate... 和 _update_state 方法
    # 都使用 self.current_M，而这个 M 在本类中被固定了。

# 临时参数结构体类（用于HHIncome_Huggett调用），与VFI保持一致
class TempParamSHH:
    def __init__(self, tau_l, theta_payg_actual_for_hh, pps_tax_deferral_active, ageEffV_new):
        self.tau_l = tau_l
        self.theta_payg_actual_for_hh = theta_payg_actual_for_hh
        self.pps_tax_deferral_active = pps_tax_deferral_active
        self.ageEffV_new = ageEffV_new
        self.tau_k = 0.2
        self.pps_tax_rate_withdrawal = 0.15
        self.current_pps_withdrawal = 0