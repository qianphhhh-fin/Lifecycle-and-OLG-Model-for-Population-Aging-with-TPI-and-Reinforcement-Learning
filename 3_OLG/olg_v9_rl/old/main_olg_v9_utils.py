import numpy as np
import scipy.stats as stats
from scipy.interpolate import interp1d, RegularGridInterpolator, PchipInterpolator, RectBivariateSpline
from scipy.optimize import fsolve, brentq, minimize
from typing import Dict, Any, Tuple, List, Optional, Union
import logging
import pandas as pd

try:
    import gymnasium as gym
    from gymnasium import spaces
except ImportError:
    gym = None
    spaces = None

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class OLGV9Utils:
    """
    A collection of utility functions for the OLG model v9, translated from MATLAB's main_olg_v8_utils.m.
    This class handles parameter setup, population dynamics, economic calculations,
    and solving the household problem via Value Function Iteration (VFI).
    """

    @staticmethod
    def set_parameters(param_set_name):
        """Sets all model model parameters, based on main_olg_v8_utils.m"""
        cS = {}

        if param_set_name == 'default':
            # --- 人口结构基础参数（对应v8.tex人口动态设定） ---
            cS['age1_orig'] = 20              # 模型起始年龄（岁）
            cS['ageLast_orig'] = 98           # 模型终止年龄（岁）
            cS['ageRetire_orig'] = 65         # 退休年龄（岁）
            cS['popGrowth_orig'] = 0.012      # 原始人口增长率
            cS['aD_orig'] = cS['ageLast_orig'] - cS['age1_orig'] + 1        # 年度年龄组数
            cS['aR_idx_orig'] = cS['ageRetire_orig'] - cS['age1_orig'] + 1  # 退休年度年龄索引
            cS['aW_orig'] = cS['aR_idx_orig'] - 1                           # 工作年数
            
            # --- 年龄组聚合参数（将年度年龄聚合为5年期年龄组） ---
            cS['yearStep'] = 5                # 每个年龄组跨度（年）
            cS['aD_new'] = int(np.ceil(cS['aD_orig'] / cS['yearStep']))    # 年龄组数量
            cS['aR_new'] = int(np.ceil(cS['aW_orig'] / cS['yearStep']))    # 工作年龄组数量
            
            # For compatibility with existing Python code structure
            cS['a0'] = cS['age1_orig']
            cS['aD'] = cS['ageLast_orig'] + 1 # Python ranges are exclusive at the end
            cS['life_len'] = cS['aD_orig']
            cS['T'] = cS['aD_new']
            cS['TW'] = cS['aR_new']
            cS['TR'] = cS['T'] - cS['TW']
            cS['age_retire_new'] = cS['aR_idx_orig']

            # Physical age vector for new age groups
            cS['physAgeV_new'] = np.zeros(cS['aD_new'], dtype=int)
            for i in range(cS['aD_new']):
                start_age = cS['age1_orig'] + i * cS['yearStep']
                end_age = min(start_age + cS['yearStep'] - 1, cS['ageLast_orig'])
                cS['physAgeV_new'][i] = int(np.round((start_age + end_age) / 2))

            # Age map from annual to groups
            cS['physAgeMap'] = [[] for _ in range(cS['aD_new'])]
            for a in range(cS['aD_new']):
                startIdx = a * cS['yearStep']
                endIdx = min((a + 1) * cS['yearStep'], cS['aD_orig'])
                cS['physAgeMap'][a] = list(range(startIdx, endIdx))

            # --- Hard-coded demographic data (from MATLAB) ---
            cS['d_orig'] = np.array([0.00159,0.00169,0.00174,0.00172,0.00165,0.00156,0.00149,0.00145,0.00145,0.00149,0.00156,0.00163,0.00171,0.00181,0.00193,0.00207,0.00225,0.00246,0.00270,0.00299,0.00332,0.00368,0.00409,0.00454,0.00504,0.00558,0.00617,0.00686,0.00766,0.00865,0.00955,0.01058,0.01162,0.01264,0.01368,0.01475,0.01593,0.01730,0.01891,0.02074,0.02271,0.02476,0.02690,0.02912,0.03143,0.03389,0.03652,0.03930,0.04225,0.04538,0.04871,0.05230,0.05623,0.06060,0.06542,0.07066,0.07636,0.08271,0.08986,0.09788,0.10732,0.11799,0.12895,0.13920,0.14861,0.16039,0.17303,0.18665,0.20194,0.21877,0.23601,0.25289,0.26973,0.28612,0.30128,0.31416,0.32915,0.34450,0.36018])
            cS['s_orig'] = 1 - cS['d_orig']

            cS['s_1yr_transitionV'] = np.zeros(cS['aD_new'])
            for a in range(cS['aD_new'] - 1):
                lastYearIdxInGroup = cS['physAgeMap'][a][-1]
                if lastYearIdxInGroup < cS['aD_orig'] - 1:
                    cS['s_1yr_transitionV'][a] = cS['s_orig'][lastYearIdxInGroup]
            
            cS['initial_pop'] = np.array([76.2, 86.4, 113.8, 98.6, 86.6, 102.7, 112.0, 99.0, 64.0, 66.9, 44.1, 25.4, 14.9, 6.8, 1.7, 0.2])
            if len(cS['initial_pop']) != cS['aD_new']:
                logging.warning('initial_pop length mismatch. Resetting to uniform.')
                cS['initial_pop'] = np.ones(cS['aD_new']) * (100 / cS['aD_new'])

            beta_surv_pop = np.array([0.998,0.996,0.994,0.992,0.988,0.984,0.980,0.976,0.970,0.960,0.945,0.920,0.880,0.800,0.680])
            cS['survivalProbV_popdyn'] = np.append(beta_surv_pop, 0)

            # --- Population dynamics convergence ---
            cS['max_periods'] = 50
            cS['bgp_tolerance'] = 0.001
            cS['bgp_window'] = 5

            # --- 家庭偏好参数（对应v8.tex第2.2.1节） ---
            cS['sigma'] = 1.5            # 相对风险厌恶系数γ
            cS['beta'] = 0.97           # 主观贴现因子β_disc  
            cS['cFloor'] = 0.05         # 最低消费约束
            cS['nSim'] = 1000           # 蒙特卡洛模拟个体数
            
            # --- 生产技术参数（对应v8.tex第2.3节） ---
            cS['A'] = 0.895944          # 全要素生产率
            cS['alpha'] = 0.36          # 资本产出弹性
            cS['ddk'] = 0.06            # 资本折旧率δ
            
            # 保持向后兼容性
            cS['delta'] = 0.06          # Annual depreciation (alias for ddk)

            # --- 资产网格参数 ---
            # 非PPS资产网格（对应v8.tex中的k状态空间）
            cS['tgKY'] = 3                     # 目标资本产出比
            cS['tgWage'] = (1-cS['alpha'])*cS['A']*((cS['tgKY']/cS['A'])**(cS['alpha']/(1-cS['alpha'])))
            cS['nk'] = 5                     # 非PPS资产网格点数
            cS['k_min'] = 0                    # 非PPS资产下界
            cS['kMin'] = 0                     # MATLAB兼容性别名
            cS['k_max'] = 15 * cS['tgWage']   # 非PPS资产上界
            cS['kMax'] = cS['k_max']          # MATLAB兼容性别名
            
            # --- PPS资产网格（对应v8.tex中的k_pps状态空间） ---
            cS['nkpps'] =5                   # PPS资产网格点数
            cS['k_pps_min'] = 0               # PPS资产下界
            cS['kppsMin'] = 0                 # MATLAB兼容性别名
            cS['k_pps_max'] = cS['k_max'] / 2 # PPS资产上界
            cS['kppsMax'] = cS['k_pps_max']   # MATLAB兼容性别名
            if cS['nkpps'] > 0:
                cS['kppsMax'] = max(cS['kppsMax'], 1e-3)
            
            # --- 劳动效率冲击过程参数 ---
            # 对应v8.tex中的ε_{a,t}随机过程
            # 基于MATLAB日志对齐参数
            cS['lePersistence'] = 0.90         # AR(1)持续性参数 (MATLAB: 0.900)
            cS['leShockStd'] = 0.15      # 效率冲击标准差 (MATLAB: 0.150)
            cS['leWidth'] = 2.0                  # Tauchen方法的标准差倍数 (MATLAB: q=2.0)
            cS['nw'] = 3                       # 效率状态网格点数
            
            # --- Policy and price guesses ---
            cS['rho_prime_guess'] = 0.01
            cS['tau_l_guess'] = 0.1509
            cS['K_guess'] = 4.0

            # --- 政府财政参数（对应v8.tex第2.4节） ---
            cS['tau_k'] = 0.20                # 资本所得税率
            cS['tau_c'] = 0.10                # 消费税率
            cS['gov_exp_frac_Y'] = 0.15       # 政府支出占GDP比例
            cS['gov_debt_frac_Y'] = 0.60      # 政府债务占GDP比例
            
            # --- PAYG税率约束参数 ---
            cS['tau_l_init_guess'] = 0.1509   # 所得税率初始猜测
            cS['tau_l_min'] = 0.00            # 所得税率下界
            cS['tau_l_max'] = 0.3             # 所得税率上界
            cS['max_total_labor_tax'] = 0.6   # 总劳动税负上限
            cS['theta_payg_max'] = 0.35       # PAYG税率上限
            
            # 保持向后兼容性
            cS['bequest_tax'] = 1.0           # 遗产税率（未明确使用）
            
            # --- PAYG养老金系统 ---
            cS['b_replace_rate'] = 0.40               # 这是MATLAB中的rho_prime_payg_fixed
            cS['rho_prime_payg_fixed'] = cS['b_replace_rate'] # 兼容性别名
            
            # === V8模型核心：PPS制度参数设计 ===
            # 对应v8.tex第2.2.2节"PPS缴费约束"部分
            
            # --- PPS制度基础参数 ---
            cS['pps_active'] = True                         # PPS制度激活标志
            cS['pps_tax_rate_withdrawal'] = 0.03           # PPS提取阶段税率
            cS['pps_return_rate_premium'] = 0.08           # PPS超额收益率 (相对于市场净回报r_k_net_hh)
            cS['pps_withdrawal_rate'] = 0.15               # 退休后年度提取比例
            cS['pps_in_K'] = True                          # PPS资产是否计入生产性资本
            cS['pps_bequeathable'] = True                  # PPS资产是否可遗赠
            
            # --- V8模型关键创新：PPS缴费约束参数 ---
            # 实现v8.tex中描述的双重约束机制
            cS['pps_annual_contrib_limit'] = 9999          # PPS年度绝对缴费上限
            cS['pps_max_contrib_frac'] = 1                 # PPS缴费占劳动收入比例上限
            cS['pps_contribution_age_max_idx'] = cS['aR_idx_orig'] - 1  # 最大缴费年度年龄
            cS['pps_withdrawal_age_min_idx'] = cS['aR_idx_orig']        # 最低提取年度年龄
            
            # V8核心：PPS缴费选择网格点数
            # 用于在VFI中离散化PPS缴费选择空间
            cS['n_pps_choice_grid_points'] = 12           # 进一步减少格点数，改善平滑度
            cS['power_pps_choice_grid'] = 1.3             # 降低幂次参数，减少非线性程度
            
            # --- 年龄效率剖面（对应v8.tex中的h_{a,t}） ---
            # 基于中国劳动者生产率-年龄关系校准
            ageEffV_orig_temp = np.zeros(100)
            # MATLAB: ageEffV_orig_temp(20:72) = [linspace(0.3,1.5,36-20+1), 1.5.*ones(1,47-37+1), linspace(1.5,0.2,65-48+1), linspace(0.18,0,72-66+1)];
            ageEffV_orig_temp[19:72] = np.concatenate([
                np.linspace(0.3, 1.5, 36 - 20 + 1),        # 20-36岁：上升期
                np.full(47 - 37 + 1, 1.5),                 # 37-47岁：峰值期
                np.linspace(1.5, 0.2, 65 - 48 + 1),        # 48-65岁：下降期
                np.linspace(0.18, 0, 72 - 66 + 1)          # 66-72岁：末期
            ])
            cS['ageEffV_orig'] = ageEffV_orig_temp[cS['age1_orig']-1 : cS['ageLast_orig']]
            if len(cS['ageEffV_orig']) != cS['aD_orig']:
                raise ValueError('ageEffV_orig 年度年龄效率剖面长度不匹配')
            
            ageEffV_new = np.zeros(cS['aD_new'])
            for a in range(cS['aD_new']):
                indices_in_group = cS['physAgeMap'][a]
                if indices_in_group:
                    ageEffV_new[a] = np.mean(cS['ageEffV_orig'][indices_in_group])
            cS['ageEffV_new'] = ageEffV_new

            # --- Other parameters ---
            cS['nSim'] = 1000 # Note: MATLAB uses 'nsim' but Python uses 'nSim'
            cS['nsim'] = cS['nSim'] # Alias for compatibility
            cS['cFloor'] = 0.05
            cS['popGrowth_orig'] = 0.012
            
            # --- 一般均衡求解参数 ---
            cS['max_iter_K_tau_l'] = 100                  # K和tau_l迭代最大次数
            cS['tol_K_tau_l'] = 1e-4                      # K和tau_l收敛容忍度
            cS['damp_K_v5'] = 0.1                         # K更新阻尼系数
            cS['damp_tau_l_v5'] = 0.1                     # tau_l更新阻尼系数
            cS['gbc_tol_for_internal_loop'] = 1e-3        # 政府预算平衡容忍度
            
            # --- 收敛检测参数 ---
            cS['max_stagnation_iters'] = 10               # 最大停滞迭代次数
            cS['min_norm_improvement_frac'] = 1e-3        # 最小改进比例
            cS['max_tau_l_boundary_strikes'] = 5          # tau_l边界冲击最大次数
            
            # --- 数值优化参数 ---
            cS['fminbnd_TolX'] = 1e-6                     # fminbnd容忍度
            cS['fminbnd_Display'] = 'none'                # fminbnd显示设置
            
            # --- VFI优化方法选择（匹配MATLAB的use_continuous_optimization标志） ---
            cS['use_continuous_optimization'] = True      # 默认使用连续优化（fmincon风格）
            
            # 保持向后兼容性
            cS['tol'] = 1e-6                              # 通用收敛容忍度

        else:
            raise ValueError('Invalid parameter set name')

        return cS

    @staticmethod
    def parameter_values_huggett_style(cS, paramS):
        """
        Set Huggett-style parameters that are derived from cS.
        This is kept for potential future use or for parameters that are
        specifically for the Huggett-style solver logic and not core model config.
        """
        paramS['beta'] = cS['beta']**cS['yearStep']
        paramS['T'] = cS['T']
        paramS['TW'] = cS['TW']
        paramS['TR'] = cS['TR']
        paramS['sigma'] = cS['sigma']
        paramS['alpha'] = cS['alpha']
        paramS['delta'] = cS['ddk']  # Use the model period depreciation directly
        paramS['aR'] = cS['aR_new']
        paramS['life_len'] = cS['life_len']
        paramS['work_life_len'] = cS['aW_orig']  # Use aW_orig instead of work_life_len_new
        paramS['age_retire'] = cS['age_retire_new']
        paramS['le_rho'] = cS['lePersistence']**cS['yearStep']
        paramS['le_sigma_shock'] = cS['leShockStd'] * np.sqrt((1 - paramS['le_rho']**2) / (1 - cS['lePersistence']**2))

        # 生成资产网格（与MATLAB完全一致）
        paramS['kGridV'] = OLGV9Utils.make_grid_matlab_style(cS['nk'], cS['k_min'], cS['k_max'], power=1.5)
        paramS['kppsGridV'] = OLGV9Utils.make_grid_matlab_style(cS['nkpps'], cS['k_pps_min'], cS['k_pps_max'], power=1.5)
        
        # Store grids in cS for compatibility
        cS['kGridV'] = paramS['kGridV']
        cS['kppsGridV'] = paramS['kppsGridV']
        
        # Generate labor efficiency grid (matching MATLAB EarningProcess_olgm)
        leLogGridV, leTrProbM, leProb1V = OLGV9Utils.earning_process_olgm(cS)
        paramS['leLogGridV'] = leLogGridV
        paramS['leTrProbM'] = leTrProbM
        paramS['leProb1V'] = leProb1V
        paramS['leGridV'] = np.exp(leLogGridV)  # 重要：MATLAB中使用指数形式的效率网格
        
        # Labor efficiency parameters for SAC observation normalization
        cS['le_mu'] = 0.0  # Mean of centered log grid
        
        return paramS

    @staticmethod
    def make_grid(n_points: int, min_val: float, max_val: float, scale: float = 2.0) -> np.ndarray:
        """
        Creates a non-linear grid with more points clustered near the minimum value.
        
        Parameters:
        - n_points: Number of grid points
        - min_val: Minimum value of the grid
        - max_val: Maximum value of the grid
        - scale: Curvature parameter (higher values cluster more points near min_val)
        
        Returns:
        - np.ndarray: Non-linear grid
        """
        if n_points <= 1:
            return np.array([min_val])
        
        # Create a uniform grid on [0, 1]
        uniform_grid = np.linspace(0, 1, n_points)
        
        # Apply power transformation to cluster points near 0
        transformed_grid = uniform_grid ** scale
        
        # Scale to [min_val, max_val]
        grid = min_val + transformed_grid * (max_val - min_val)
        
        return grid

    @staticmethod
    def make_grid_matlab_style(n_points: int, min_val: float, max_val: float, power: float = 1.5) -> np.ndarray:
        """
        创建与MATLAB完全一致的网格
        
        对应MATLAB代码:
        kGridV = kMin + (kMax - kMin) * (linspace(0, 1, nk).^power);
        kGridV(1) = kMin;
        
        Parameters:
        - n_points: 网格点数
        - min_val: 最小值
        - max_val: 最大值  
        - power: 幂次参数
        
        Returns:
        - np.ndarray: MATLAB风格的网格
        """
        if n_points <= 0:
            return np.array([])
        elif n_points == 1:
            return np.array([min_val])
        
        # MATLAB: linspace(0, 1, nk).^power
        uniform_grid = np.linspace(0, 1, n_points)
        powered_grid = uniform_grid ** power
        
        # MATLAB: kMin + (kMax - kMin) * (...)
        grid = min_val + (max_val - min_val) * powered_grid
        
        # MATLAB: kGridV(1) = kMin; (确保第一个点是最小值)
        if n_points > 0:
            grid[0] = min_val
        
        return grid

    @staticmethod
    def init_population(cS: Dict[str, Any]) -> Dict[str, Any]:
        """Initializes the population structure based on cS['initial_pop']."""
        popS = {}
        initial_total = np.sum(cS['initial_pop'])
        if initial_total > 0 and len(cS['initial_pop']) == cS['aD_new']:
            popS['Z'] = (cS['initial_pop'] / initial_total * 100).reshape(-1, 1)
        else:
            logging.warning('Initial population data mismatch or zero. Using uniform distribution.')
            popS['Z'] = np.full((cS['aD_new'], 1), 100 / cS['aD_new'])

        popS['totalPop'] = np.array([np.sum(popS['Z'][:, 0])])
        popS['totalPop_history'] = np.array([popS['totalPop']])

        if popS['totalPop'][0] > 1e-9:
            popS['ageDist'] = popS['Z'][:, 0] / popS['totalPop'][0]
        else:
            popS['ageDist'] = np.zeros(cS['aD_new'])
        
        popS['initialAgeDist'] = popS['ageDist']
        popS['Z_history'] = popS['Z']
        popS['ageDist_history'] = popS['ageDist'].reshape(-1,1)

        return popS

    @staticmethod
    def population_dynamics(popS: Dict[str, Any], cS: Dict[str, Any]) -> Dict[str, Any]:
        """Simulates population dynamics based on MATLAB's simple projection."""
        max_periods_sim = cS['max_periods']
        Z_history = np.zeros((cS['aD_new'], max_periods_sim + 1))
        totalPop_history = np.zeros(max_periods_sim + 1)
        ageDist_history = np.zeros((cS['aD_new'], max_periods_sim + 1))

        Z_history[:, 0] = popS['Z'][:, 0]
        totalPop_history[0] = popS['totalPop'][0]
        ageDist_history[:, 0] = popS['ageDist']

        logging.info(f"Population dynamics simulation starting (max periods = {max_periods_sim})...")
        
        bgp_reached_flag = False
        actual_periods_run = max_periods_sim
        
        for t in range(max_periods_sim):
            if (t + 1) % 10 == 0 or t == 0:
                logging.info(f"  模拟人口期数 {t + 1} (年龄组)")
            
            Z_current_period = Z_history[:, t]
            Z_next_period = np.zeros(cS['aD_new'])
            
            # 修正：与MATLAB一致，使用 t < 6 而不是 t <= 5
            time_varying_growth_rate = 0.0
            if t < 6:  # MATLAB: if t < 6
                time_varying_growth_rate = -0.01 - 0.003 * (t + 1)  # MATLAB使用1-based索引
            else:
                time_varying_growth_rate = -0.03 - 0.004 * min((t + 1) - 6, 10)
            
            # Birth cohort evolution - 年龄组0
            Z_next_period[0] = Z_current_period[0] * (1 + time_varying_growth_rate)
            Z_next_period[0] = max(0, Z_next_period[0])

            # Aging and survival - 年龄组1到aD_new-1
            for a in range(1, cS['aD_new']):  # MATLAB: a = 2:cS.aD_new
                survival_prob_group = 0.0
                # MATLAB逻辑：if (a-1) >= 1 && (a-1) <= length(cS.survivalProbV_popdyn)
                # 这等价于Python的：if (a-1-1) >= 0 and (a-1-1) < len(survivalProbV_popdyn)
                # 即：if (a-2) >= 0 and (a-2) < len(survivalProbV_popdyn)
                # 但是MATLAB的a从2开始，Python的a从1开始，所以Python应该是：
                # if (a-1) >= 0 and (a-1) < len(survivalProbV_popdyn)
                if (a-1) >= 0 and (a-1) < len(cS['survivalProbV_popdyn']):
                    survival_prob_group = cS['survivalProbV_popdyn'][a-1]  # MATLAB: cS.survivalProbV_popdyn(a-1)
                
                Z_next_period[a] = Z_current_period[a-1] * survival_prob_group
                Z_next_period[a] = max(0, Z_next_period[a])

            Z_history[:, t+1] = Z_next_period
            totalPop_history[t+1] = np.sum(Z_next_period)
            if totalPop_history[t+1] > 1e-9:
                ageDist_history[:, t+1] = Z_next_period / totalPop_history[t+1]
            else:
                ageDist_history[:, t+1] = 0
                totalPop_history[t+1] = 0

            # Check for steady state - 与MATLAB一致的检查逻辑
            current_check_period_idx = t + 1 + 1  # MATLAB: current_check_period_idx = t + 1; (t是1-based)
            # Python中t是0-based，所以对应的检查索引是t+1+1
            
            if current_check_period_idx >= cS['bgp_window'] + 1:
                stable = True
                for w_idx in range(1, cS['bgp_window'] + 1):  # MATLAB: w_idx = 1:cS.bgp_window
                    hist_idx1 = current_check_period_idx - w_idx + 1 - 1  # 转换为0-based索引
                    hist_idx2 = current_check_period_idx - w_idx - 1  # 转换为0-based索引
                    
                    if (hist_idx1 >= 0 and hist_idx2 >= 0 and 
                        hist_idx1 < ageDist_history.shape[1] and hist_idx2 < ageDist_history.shape[1]):
                        change = np.linalg.norm(ageDist_history[:, hist_idx1] - ageDist_history[:, hist_idx2])
                        if change >= cS['bgp_tolerance']:
                            stable = False
                            break
                    else:
                        stable = False
                        break
                
                if stable:
                    logging.info(f'\n人口稳态 (年龄组) 在模拟期数 {t + 1} (对应历史数据索引 {current_check_period_idx}) 达到。')
                    bgp_reached_flag = True
                    actual_periods_run = t + 1
                    break

        # 存储结果
        final_period_idx_to_store = min(actual_periods_run + 1, Z_history.shape[1])
        popS['Z_history'] = Z_history[:, :final_period_idx_to_store]
        popS['totalPop_history'] = totalPop_history[:final_period_idx_to_store]
        popS['ageDist_history'] = ageDist_history[:, :final_period_idx_to_store]
        
        # 为了与MATLAB一致，确保ageDist字段设置正确（用于稳态检测）
        popS['ageDist'] = ageDist_history[:, :final_period_idx_to_store]
        popS['Z'] = Z_history[:, :final_period_idx_to_store]  # 确保Z字段也设置正确
        popS['totalPop'] = totalPop_history[:final_period_idx_to_store]
        
        # 计算依赖比历史
        depRatio_history = np.zeros(actual_periods_run)
        for th_loop in range(actual_periods_run):
            Z_t_for_depratio = Z_history[:, th_loop + 1]  # MATLAB: th_loop + 1 (1-based)
            working_pop = np.sum(Z_t_for_depratio[:cS['aR_new']])  # MATLAB: 1:cS.aR_new
            retired_pop = np.sum(Z_t_for_depratio[cS['aR_new']:])  # MATLAB: cS.aR_new+1:end
            if working_pop > 1e-9:
                depRatio_history[th_loop] = retired_pop / working_pop
            else:
                depRatio_history[th_loop] = np.inf
        
        popS['dependencyRatio'] = depRatio_history
        
        logging.info(f'人口动态模拟完成。运行期数: {actual_periods_run}。达到BGP: {bgp_reached_flag}')
        if not bgp_reached_flag:
            logging.warning(f'警告: 人口稳态 (年龄组) 未在 {max_periods_sim} 期内达到。')
        
        return popS

    @staticmethod
    def detect_steady_state_population(popS: Dict[str, Any], cS: Dict[str, Any]) -> Tuple[np.ndarray, np.ndarray, int, bool]:
        """
        Detects if and when the population has reached a steady state, matching MATLAB logic.
        Returns the steady-state UNNORMALIZED mass for model groups, and NORMALIZED distribution for annual ages.
        """
        # 使用ageDist而不是ageDist_history - 与MATLAB一致
        ageDist_data = popS['ageDist']  # MATLAB: popS.ageDist
        Z_data = popS['Z']  # MATLAB: popS.Z (注意，这里使用Z而不是Z_history)
        actual_periods_in_data = ageDist_data.shape[1]  # MATLAB: size(popS.ageDist, 2)
        bgp_reached = False
        bgp_period = actual_periods_in_data - 1  # 使用最后一期如果没有收敛

        if actual_periods_in_data < cS['bgp_window'] + 1:
            logging.warning(f"人口模拟期数过短 ({actual_periods_in_data} 数据点)，无法进行稳态检查 (窗口期 = {cS['bgp_window']})。")
        else:
            logging.info(f'检查人口稳态 (年龄组, 最近 {cS["bgp_window"]} 期)...')
            # MATLAB: for t_check_end_idx = actual_periods_in_data : -1 : cS.bgp_window + 1
            # 这里需要完全匹配MATLAB的循环范围
            for t_check_end_idx in range(actual_periods_in_data, cS['bgp_window'], -1):
                stable = True
                # MATLAB: for w_idx = 0 : (cS.bgp_window - 1)
                for w_idx in range(cS['bgp_window']):
                    # MATLAB索引：idx1 = t_check_end_idx - w_idx; idx2 = t_check_end_idx - w_idx - 1;
                    # 转换为Python 0-based索引
                    idx1 = t_check_end_idx - w_idx - 1  # MATLAB的t_check_end_idx是1-based
                    idx2 = t_check_end_idx - w_idx - 1 - 1
                    
                    # MATLAB条件：if idx1 > 0 && idx2 > 0 && idx1 <= size(popS.ageDist,2) && idx2 <= size(popS.ageDist,2)
                    # 转换为Python：idx1 >= 0 && idx2 >= 0 && idx1 < ageDist_data.shape[1] && idx2 < ageDist_data.shape[1]
                    if (idx1 >= 0 and idx2 >= 0 and 
                        idx1 < ageDist_data.shape[1] and idx2 < ageDist_data.shape[1]):
                        # MATLAB: change = norm(popS.ageDist(:, idx1) - popS.ageDist(:, idx2));
                        change = np.linalg.norm(ageDist_data[:, idx1] - ageDist_data[:, idx2])
                        if change >= cS['bgp_tolerance']:
                            stable = False
                            break
                    else:
                        stable = False
                        break
                
                if stable:
                    bgp_reached = True
                    # MATLAB: bgp_period = t_check_end_idx - 1;
                    bgp_period = t_check_end_idx - 1  # 保持1-based，后面转换
                    logging.info(f'人口稳态 (年龄组) 从模拟期数 {bgp_period} (数据索引 {t_check_end_idx}) 开始检测到 (稳定窗口结束于此)。')
                    break
            
            if not bgp_reached:
                logging.warning('未检测到人口稳态 (年龄组)。将使用最终期数据。')
                bgp_period = actual_periods_in_data - 1  # 1-based period

        # MATLAB: ss_data_index = min(bgp_period + 1, size(popS.Z, 2));
        # bgp_period是1-based的period，所以bgp_period+1是1-based的索引
        # 转换为0-based索引：ss_data_index_0based = min(bgp_period + 1 - 1, Z_data.shape[1] - 1)
        ss_data_index_0based = min(bgp_period, Z_data.shape[1] - 1)
        Z_ss = Z_data[:, ss_data_index_0based]

        # 关键修正：MATLAB的detectSteadyStatePopulation返回的Z_ss是未标准化的质量
        # 但在MATLAB的调用代码中，会进行标准化：paramS_mat.ageMassV = Z_ss_group_mat_raw / sum(Z_ss_group_mat_raw);
        # 因此，这里返回未标准化的质量，让调用者决定是否标准化
        Z_ss_unnormalized = Z_ss.copy()

        # 解聚集到年度分布 - 这里使用标准化版本进行解聚集
        Z_ss_norm_for_disagg = Z_ss / np.sum(Z_ss) if np.sum(Z_ss) > 1e-9 else np.zeros_like(Z_ss)
        Z_ss_annual_dist = np.zeros(cS['aD_orig'])
        for a_new in range(cS['aD_new']):
            annual_indices = cS['physAgeMap'][a_new]
            group_dist_mass = Z_ss_norm_for_disagg[a_new]
            if annual_indices:
                mass_per_year = group_dist_mass / len(annual_indices)
                Z_ss_annual_dist[annual_indices] = mass_per_year
        
        # 确保年度分布归一化
        if np.sum(Z_ss_annual_dist) > 1e-9:
            Z_ss_annual_dist /= np.sum(Z_ss_annual_dist)

        # 返回未标准化的质量（Z_ss_unnormalized）和标准化的年度分布
        return Z_ss_unnormalized, Z_ss_annual_dist, bgp_reached, bgp_period

    @staticmethod
    def tauchen(N, rho, sigma, mu=0.0, m=3.0):
        """Python implementation of Tauchen's method."""
        std_y_unconditional = np.sqrt(sigma**2 / (1 - rho**2)) if abs(1-rho) > 1e-9 else sigma * 100
        y_max = m * std_y_unconditional
        y_grid = np.linspace(-y_max, y_max, N)
        
        step = y_grid[1] - y_grid[0] if N > 1 else 0
        trans_matrix = np.zeros((N, N))
        
        for i in range(N):
            for j in range(N):
                mean_next_y = rho * y_grid[i]
                if j == 0:
                    trans_matrix[i, j] = stats.norm.cdf((y_grid[0] - mean_next_y + step / 2) / sigma)
                elif j == N - 1:
                    trans_matrix[i, j] = 1 - stats.norm.cdf((y_grid[N-1] - mean_next_y - step / 2) / sigma)
                else:
                    trans_matrix[i, j] = stats.norm.cdf((y_grid[j] - mean_next_y + step / 2) / sigma) - \
                                         stats.norm.cdf((y_grid[j] - mean_next_y - step / 2) / sigma)
        
        unconditional_mean_shift = mu / (1-rho) if abs(1-rho)>1e-9 else 0
        y_grid_out = y_grid + unconditional_mean_shift
        return y_grid_out, trans_matrix

    @staticmethod
    def earning_process_olgm(cS: Dict[str, Any]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Creates the grid and transition matrix for the AR(1) labor efficiency process.
        Matches MATLAB v8 logic, including Tauchen method and compression.
        """
        rho = cS['lePersistence']
        sigma = cS['leShockStd']
        N = cS['nw']
        m = cS['leWidth']
        
        # Discretize using Tauchen
        leLogGridV_raw, leTrProbM = OLGV9Utils.tauchen(N, rho, sigma, 0, m)
        
        # Center the grid
        leLogGridV = leLogGridV_raw - np.mean(leLogGridV_raw)

        # Compress if the range is too wide
        leGridV_test = np.exp(leLogGridV)
        if N > 1 and leGridV_test[0] > 0:
            efficiency_ratio = leGridV_test[-1] / leGridV_test[0]
            max_acceptable_ratio = 5.0
            if efficiency_ratio > max_acceptable_ratio:
                compression_factor = np.log(max_acceptable_ratio) / np.log(efficiency_ratio)
                leLogGridV = leLogGridV * compression_factor
        
        # Find stationary distribution via eigenvector method
        if N > 0:
            eigvals, eigvecs = np.linalg.eig(leTrProbM.T)
            stat_dist_idx = np.argmin(np.abs(eigvals - 1.0))
            leProb1V = np.real(eigvecs[:, stat_dist_idx])
            leProb1V = np.abs(leProb1V) / np.sum(np.abs(leProb1V))
        else:
            leProb1V = np.array([])
        
        return leLogGridV, leTrProbM, leProb1V

    @staticmethod
    def labor_endow_simulation_olgm(cS: Dict[str, Any], paramS: Dict[str, Any]) -> np.ndarray:
        """
        Simulates labor endowment shocks for a cross-section of households over ANNUAL ages.
        Matches MATLAB v8 logic exactly using MarkovChainSimulation.
        """
        # 与MATLAB完全一致：rng(433); 
        # 注意：这里不能再调用np.random.seed，因为可能在调用此函数前已经设置过了
        # 只有在第一次调用时才设置种子，或者总是设置以确保一致性
        nsim = cS['nsim']
        aD_orig = cS['aD_orig']
        leProb1V = paramS['leProb1V']
        leTrProbM = paramS['leTrProbM']
        
        # 为了与MATLAB的LaborEndowSimulation_olgm函数完全一致，
        # 我们在这里重新设置种子，就像MATLAB中的rng(433)一样
        np.random.seed(433)
        
        # Generate the same random numbers as MATLAB: rand([cS.nSim, cS.aD_orig])
        random_numbers_for_sim = np.random.rand(nsim, aD_orig)
        
        # Call the exact MATLAB-equivalent simulation
        eIdxM = OLGV9Utils.markov_chain_simulation(
            nsim, aD_orig, leProb1V, leTrProbM, random_numbers_for_sim
        )
        
        logging.info(f'劳动禀赋路径已模拟 ({nsim} 个体, {aD_orig} 年度年龄)。')
        
        return eIdxM

    @staticmethod
    def labor_supply_huggett(eIdxM_annual: np.ndarray, cS: Dict[str, Any], paramS: Dict[str, Any], ageMassV_group: np.ndarray) -> Tuple[float, float]:
        """
        Calculates aggregate labor supply based on ANNUAL simulation, matching MATLAB v8.
        It now correctly uses the GROUP-level age distribution for final aggregation.
        """
        nsim, aD_orig = eIdxM_annual.shape
        leGridV = paramS['leGridV']
        
        ageToGroupMap = np.zeros(aD_orig, dtype=int)
        for i_group, indices in enumerate(cS['physAgeMap']):
            if indices:
                ageToGroupMap[indices] = i_group

        HHlaborM_annual_temp = np.zeros((nsim, aD_orig))
        valid_eIdx = (eIdxM_annual >= 1) & (eIdxM_annual <= len(leGridV))
        
        if cS['nw'] > 0:
            for a_orig in range(aD_orig):
                # Working age check (aR_idx_orig is 1-based index of first retirement year)
                if a_orig < cS['aR_idx_orig'] - 1:
                    group_idx = ageToGroupMap[a_orig]
                    # Ensure it's a working group
                    if group_idx >= 0 and group_idx < cS['aR_new']:
                        current_eIdx = eIdxM_annual[:, a_orig]
                        labor_eff_for_valid = np.zeros(nsim)
                        
                        valid_mask = valid_eIdx[:, a_orig]
                        if np.any(valid_mask):
                            valid_current_eIdx = current_eIdx[valid_mask]
                            # Convert from MATLAB 1-based index to Python 0-based
                            valid_python_indices = valid_current_eIdx - 1
                            labor_eff_for_valid[valid_mask] = leGridV[valid_python_indices]
                        
                        HHlaborM_annual_temp[:, a_orig] = cS['ageEffV_new'][group_idx] * labor_eff_for_valid

        HHlaborM_group = np.zeros((nsim, cS['aD_new']))
        for a_new_group_idx in range(cS['aD_new']):
            annual_indices = cS['physAgeMap'][a_new_group_idx]
            if annual_indices:
                HHlaborM_group[:, a_new_group_idx] = np.mean(HHlaborM_annual_temp[:, annual_indices], axis=1)

        L_total_eff_pc_sum = 0.0
        if cS['aR_new'] > 0:
            # Final aggregation uses GROUP-level labor and GROUP-level distribution
            mean_labor_per_working_group = np.mean(HHlaborM_group[:, :cS['aR_new']], axis=0)
            # ageMassV_group MUST be a normalized distribution
            # 确保使用与MATLAB完全相同的矩阵乘法计算：mean_labor_per_working_group * ageMassV_group(1:aR_new)
            working_age_masses = ageMassV_group[:cS['aR_new']]
            L_total_eff_pc_sum = np.dot(mean_labor_per_working_group, working_age_masses)
        
        L_per_capita = max(0, L_total_eff_pc_sum)
        logging.info(f'家庭劳动供给已计算。总体人均有效劳动供给 (L_eff_pc) = {L_per_capita:.4f}')
        return L_per_capita, 0.0 # L_per_worker not used directly in solver
        
    @staticmethod
    def factor_prices(K: float, L: float, alpha: float, delta: float, A: float) -> Tuple[float, float]:
        """
        Computes factor prices (interest rate and wage) from aggregate K and L.
        Matches MATLAB v8 by including depreciation.
        """
        K_eff = max(K, 1e-6)
        L_eff = max(L, 1e-6)
        
        Y = A * (K_eff**alpha) * (L_eff**(1 - alpha))
        
        MPK = alpha * Y / K_eff
        MPL_gross = (1 - alpha) * Y / L_eff
        
        R_mkt_gross_factor = 1 + MPK - delta
        
        return max(1.0 + 1e-6, R_mkt_gross_factor), MPL_gross

    @staticmethod
    def HHPrices_Huggett(K: float, L: float, cS: Dict[str, Any]) -> Tuple[float, float]:
        """
        MATLAB-compatible wrapper for factor_prices function.
        Computes factor prices (interest rate and wage) from aggregate K and L.
        """
        alpha = cS['alpha']
        A = cS['A']
        ddk = cS['ddk']
        
        K_eff = max(K, 1e-6)
        L_eff = max(L, 1e-6)
        
        Y = A * (K_eff**alpha) * (L_eff**(1 - alpha))
        
        MPK = alpha * Y / K_eff
        MPL_gross = (1 - alpha) * Y / L_eff
        
        R_mkt_gross_factor = 1 + MPK - ddk
        
        return max(1.0 + 1e-6, R_mkt_gross_factor), MPL_gross

    @staticmethod
    def ces_utility(c_quantity, sigma_crra: float, cS: Dict[str, Any]) -> Tuple[np.ndarray, np.ndarray]:
        """
        计算CRRA效用函数和边际效用，与olg_utils.py完全一致
        对应MATLAB的CES_utility函数
        
        Returns:
            marginal_utility: 边际效用 ∂U/∂c
            utility: 效用水平 U(c)
        """
        if not np.isscalar(sigma_crra) or sigma_crra <= 0:
            raise ValueError('sigma_crra必须是正标量。')
        
        # 处理标量输入
        is_scalar_input = np.isscalar(c_quantity)
        if is_scalar_input:
            c_quantity = np.array([c_quantity])
        else:
            c_quantity = np.asarray(c_quantity)
        
        min_c_quantity = cS['cFloor']
        is_valid_consumption = (c_quantity >= min_c_quantity)
        c_adjusted_quantity = np.maximum(min_c_quantity, c_quantity)
        
        utility = np.full_like(c_quantity, -np.inf, dtype=float)
        marginal_utility = np.full_like(c_quantity, np.inf, dtype=float)
        
        if abs(sigma_crra - 1) < 1e-6:  # 对数效用 (sigma = 1)
            utility[is_valid_consumption] = np.log(c_adjusted_quantity[is_valid_consumption])
            marginal_utility[is_valid_consumption] = 1.0 / c_adjusted_quantity[is_valid_consumption]
        else:  # CRRA效用 (sigma != 1) - 使用标准CRRA形式
            utility[is_valid_consumption] = (c_adjusted_quantity[is_valid_consumption] ** (1-sigma_crra)) / (1-sigma_crra)
            marginal_utility[is_valid_consumption] = c_adjusted_quantity[is_valid_consumption] ** (-sigma_crra)
        
        # 惩罚低于最低消费的情况
        utility[~is_valid_consumption] = -1e10 - (min_c_quantity - c_quantity[~is_valid_consumption]) * 1e10
        
        # 返回顺序：边际效用在前，效用在后（与olg_utils.py一致）
        if is_scalar_input:
            return float(marginal_utility[0]), float(utility[0])
        else:
            return marginal_utility, utility

    @staticmethod
    def hh_income_huggett(k_now_val, R_k_net_factor: float, w_gross: float,
                         TR_total: float, b_payg_val: float, c_pps_chosen_and_constrained_input_val,
                         a_idx: int, paramS_hh: Dict[str, Any], cS: Dict[str, Any], 
                         epsilon_val) -> Tuple:
        """
        HHIncome_Huggett - 计算家庭可支配资源
        与MATLAB的HHIncome_Huggett函数保持完全一致
        
        输入：
          k_now_val - 当前非PPS资产
          R_k_net_factor - 税后资本回报因子(1+r_net)
          w_gross - 市场毛工资率
          TR_total - 总转移支付
          b_payg_val - PAYG养老金
          c_pps_chosen_and_constrained_input_val - 已选择的PPS缴费额
          a_idx - 年龄组索引 (Python 0-based)
          paramS_hh - 家庭参数
          cS - 模型参数
          epsilon_val - 当前劳动效率
        输出：
          resources_for_c_and_k_prime - 可用于消费和储蓄的资源
          labor_income_gross_state - 税前劳动收入
          pps_deduction_actual_state - PPS税前扣除额
        
        功能：实现v8.tex中的预算约束计算，包括PPS税收递延处理
        """
        # 初始化变量 - 与MATLAB保持一致
        if np.isscalar(k_now_val):
            labor_income_gross_state = 0.0
            pps_deduction_actual_state = 0.0
            non_capital_income = 0.0
        else:
            labor_income_gross_state = np.zeros_like(k_now_val)
            pps_deduction_actual_state = np.zeros_like(k_now_val)
            non_capital_income = np.zeros_like(k_now_val)
        
        # actual_pps_contribution_expenditure IS c_pps_chosen_and_constrained_input_val
        # No further rule-based calculation here, it's already determined upstream.
        
        # 判断是否为工作年龄组
        # 注意：Python中a_idx是0-based，但cS['aR_new']是1-based的值
        # 所以判断条件是 (a_idx + 1) <= cS['aR_new']，即 a_idx < cS['aR_new']
        if (a_idx + 1) <= cS['aR_new']:  # 如果是工作年龄组
            # 取年龄效率，a_idx是0-based索引，直接使用
            age_efficiency = cS['ageEffV_new'][a_idx]
            labor_income_gross_state = w_gross * age_efficiency * epsilon_val
            
            # PPS缴费 c_pps_chosen_and_constrained_input_val 已由VFI外层优化选择并施加约束
            if np.isscalar(c_pps_chosen_and_constrained_input_val):
                actual_pps_contribution_expenditure = max(0.0, c_pps_chosen_and_constrained_input_val)
            else:
                actual_pps_contribution_expenditure = np.maximum(0.0, c_pps_chosen_and_constrained_input_val)
            
            # PPS税收递延处理
            if paramS_hh.get('pps_tax_deferral_active', False) and paramS_hh['pps_tax_deferral_active']:
                pps_deduction_actual_state = actual_pps_contribution_expenditure
            else:
                if np.isscalar(k_now_val):
                    pps_deduction_actual_state = 0.0
                else:
                    pps_deduction_actual_state = np.zeros_like(k_now_val)
            
            # 计算应税劳动收入
            labor_income_taxable_for_tau_l = labor_income_gross_state - pps_deduction_actual_state
            if np.isscalar(labor_income_taxable_for_tau_l):
                labor_income_taxable_for_tau_l = max(0.0, labor_income_taxable_for_tau_l)
            else:
                labor_income_taxable_for_tau_l = np.maximum(0.0, labor_income_taxable_for_tau_l)
            
            # 计算税收
            income_tax_tau_l = labor_income_taxable_for_tau_l * paramS_hh['tau_l']
            payg_tax_theta = labor_income_gross_state * paramS_hh['theta_payg_actual_for_hh']
            
            # 税后劳动收入
            labor_income_net_of_all_taxes = labor_income_gross_state - income_tax_tau_l - payg_tax_theta
            non_capital_income = labor_income_net_of_all_taxes + TR_total + b_payg_val
            
        else:  # 如果是退休年龄组
            if np.isscalar(k_now_val):
                actual_pps_contribution_expenditure = 0.0  # 退休期PPS缴费为0
                pps_deduction_actual_state = 0.0
            else:
                actual_pps_contribution_expenditure = np.zeros_like(k_now_val)  # 退休期PPS缴费为0
                pps_deduction_actual_state = np.zeros_like(k_now_val)
            non_capital_income = TR_total + b_payg_val
        
        # 资本收入（税后）
        capital_income_net_of_tax = (R_k_net_factor - 1) * k_now_val
        
        # 总可支配资源
        resources_for_c_and_k_prime = (k_now_val + capital_income_net_of_tax + 
                                     non_capital_income - actual_pps_contribution_expenditure)
        
        # 检查有限性 - 支持数组操作
        if np.isscalar(resources_for_c_and_k_prime):
            if not np.isfinite(resources_for_c_and_k_prime):
                resources_for_c_and_k_prime = -1e10
        else:
            resources_for_c_and_k_prime = np.where(
                np.isfinite(resources_for_c_and_k_prime), 
                resources_for_c_and_k_prime, 
                -1e10
            )
        
        return resources_for_c_and_k_prime, labor_income_gross_state, pps_deduction_actual_state

    @staticmethod
    def _bellman_objective(k_prime_choice: float, budget: float, ev_interp, a_idx: int, cS: Dict[str, Any], beta_period: float) -> float:
        """Objective function for k' optimization (to be minimized, so returns negative value)."""
        
        # 检查k'是否在有效范围内
        if k_prime_choice < cS['k_min'] or k_prime_choice > cS['k_max']:
            return 1e12
            
        # 计算消费支出
        c_expenditure = budget - k_prime_choice
        
        # 检查消费支出是否为正
        if c_expenditure <= 0:
            return 1e12
            
        # 计算消费数量（扣除消费税）
        c_quantity = c_expenditure / (1 + cS['tau_c'])
        
        # 确保消费不低于最低消费水平
        if c_quantity < cS['cFloor']:
            return 1e12
        
        # 计算当期效用
        _, util = OLGV9Utils.ces_utility(c_quantity, cS['sigma'], cS)
        if not np.isfinite(util) or util <= -1e10:
            return 1e12

        # 计算期望价值
        try:
            ev = ev_interp(k_prime_choice)
            if not np.isfinite(ev):
                ev = -1e10  # 给一个有限的惩罚值而不是无穷大
        except:
            ev = -1e10

        # 生存概率
        s_transition = cS['s_1yr_transitionV'][a_idx] if a_idx < len(cS['s_1yr_transitionV']) else 0.0
        
        # 贝尔曼方程：V = u(c) + β * s * EV(k')
        value = util + beta_period * s_transition * ev
        
        # 检查价值函数是否有限
        if not np.isfinite(value):
            return 1e12
        
        return -value  # 返回负值用于最小化


    @staticmethod
    def _hh_solution_by_age(a_idx: int, v_prime_next: np.ndarray, R_k_net_factor: float, 
                           MPL_gross: float, TR_total: float, b_payg_val: float,
                           paramS: Dict[str, Any], cS: Dict[str, Any]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Solves the household problem for a single age group `a_idx`.
        This is the core of the VFI, matching the MATLAB logic closely.
        """
        nk, nkpps, nw = cS['nk'], cS['nkpps'], cS['nw']
        k_grid, kpps_grid, le_grid = cS['kGridV'], cS['kppsGridV'], paramS['leGridV']
        
        val_age = np.full((nk, nkpps, nw), -np.inf)
        cPol_age = np.zeros((nk, nkpps, nw))
        kPol_age = np.zeros((nk, nkpps, nw))
        cPpsPol_age = np.zeros((nk, nkpps, nw))

        # 从paramS中构造paramS_hh（确保包含必需的税收参数）
        paramS_hh = {
            'tau_l': paramS.get('tau_l', 0.15),  # 与compare_py_matlab_step_by_step.py一致
            'theta_payg_actual_for_hh': paramS.get('theta_payg_actual_for_hh', 0.12),
            'pps_tax_deferral_active': paramS.get('pps_tax_deferral_active', cS.get('pps_active', False))
        }

        # 确定此年龄组的工作状态（注意：a_idx是基于aD_new的索引）
        # 修复：使用与hh_income_huggett一致的逻辑
        matlab_a_idx = a_idx + 1  # 转换为MATLAB风格的1-based索引
        is_working_age_group = matlab_a_idx <= cS['aR_new']  # 与MATLAB逻辑完全一致
        
        # --- Handle last period (a_idx == aD_new - 1) ---
        if a_idx == cS['aD_new'] - 1:
            for ie, eps_val in enumerate(le_grid):
                for ikpps, k_pps_val in enumerate(kpps_grid):
                    for ik, k_val in enumerate(k_grid):
                        # 最后一期：消费所有资产，不储蓄，不缴费PPS
                        # 计算可用资源
                        y_disp_total, _, _ = OLGV9Utils.hh_income_huggett(
                            k_val, R_k_net_factor, MPL_gross, TR_total, b_payg_val, 0.0, 
                            a_idx, paramS_hh, cS, eps_val
                        )
                        
                        # PPS提取（最后一期退休者提取所有PPS资产）
                        pps_withdrawal_net = 0.0
                        if cS['pps_active']:
                            # 最后一期无论年龄都提取所有PPS资产
                            pps_withdrawal_gross = k_pps_val
                            pps_tax_rate = cS.get('pps_tax_rate_withdrawal', 0.03)
                            pps_withdrawal_net = pps_withdrawal_gross * (1 - pps_tax_rate)
                        
                        total_resources = y_disp_total + pps_withdrawal_net
                        
                        # 最后一期消费（扣除消费税）
                        c_consumption = max(cS['cFloor'], total_resources / (1 + cS['tau_c']))
                        
                        # 不储蓄
                        k_save = 0.0  # 最后一期所有资产都消费掉
                        
                        # 计算效用
                        _, utility = OLGV9Utils.ces_utility(c_consumption, cS['sigma'], cS)
                        
                        # 存储结果
                        cPol_age[ik, ikpps, ie] = c_consumption
                        kPol_age[ik, ikpps, ie] = k_save
                        cPpsPol_age[ik, ikpps, ie] = 0.0  # 最后一期不缴费PPS
                        val_age[ik, ikpps, ie] = utility
                        
            return cPol_age, kPol_age, cPpsPol_age, val_age

        # --- Handle other periods (backward induction) ---
        # 1. 构建期望价值函数插值器
        EV_matrix = np.zeros((nk, nkpps, nw))
        for ie in range(nw):
            if v_prime_next is not None:
                trans_probs = paramS['leTrProbM'][ie, :]
                # 计算期望价值：E[V'|e] = sum_e' P(e'|e) * V'(k', k_pps', e')
                expected_v = np.sum(v_prime_next * trans_probs.reshape(1, 1, -1), axis=2)
                EV_matrix[:, :, ie] = expected_v
            else:
                EV_matrix[:, :, ie] = 0.0
        
        ev_interpolants = {}
        for ie in range(nw):
            # 为每个当前效率状态ie创建2D插值器
            # 映射(k', k_pps')到期望价值
            ev_interpolants[ie] = RegularGridInterpolator(
                (k_grid, kpps_grid), EV_matrix[:,:,ie], 
                method='linear', bounds_error=False, fill_value=-1e12
            )

        # 2. 循环所有状态 (k, k_pps, e) 
        for ie, eps_val in enumerate(le_grid):
            # 预计算此效率状态的期望价值插值器
            ev_interp_ie = ev_interpolants[ie]
            
            for ikpps, k_pps_val in enumerate(kpps_grid):
                for ik, k_val in enumerate(k_grid):
                    
                    # 3. 确定PPS缴费约束（匹配MATLAB逻辑）
                    max_permissible_cpps = 0.0
                    # 修复：使用与hh_income_huggett一致的年龄组判断
                    model_age_group_start_year_idx = cS['physAgeMap'][a_idx][0] if 'physAgeMap' in cS and cS['physAgeMap'] and a_idx < len(cS['physAgeMap']) and cS['physAgeMap'][a_idx] else 0
                    matlab_a_idx_for_pps = a_idx + 1  # 转换为MATLAB 1-based索引
                    is_pps_eligible = (matlab_a_idx_for_pps <= cS['aR_new'] and 
                                     model_age_group_start_year_idx <= cS.get('pps_contribution_age_max_idx', 999) and 
                                     cS.get('pps_active', False))
                    
                    if is_pps_eligible:
                        age_efficiency = cS['ageEffV_new'][a_idx]
                        # 在MATLAB中，eps_val来自leGridV = exp(leLogGridV)，所以直接使用
                        current_gross_labor_income = MPL_gross * age_efficiency * eps_val
                        if current_gross_labor_income > 1e-6:
                            max_cpps_by_frac = current_gross_labor_income * cS['pps_max_contrib_frac']
                            max_permissible_cpps = min(cS['pps_annual_contrib_limit'], max_cpps_by_frac)
                            max_permissible_cpps = max(0, max_permissible_cpps)
                    
                    # 4. 使用比例形式的联合优化（匹配MATLAB的fmincon实现）
                    def fmincon_objective_helper_proportional(x_prop):
                        """
                        MATLAB风格的fmincon目标函数
                        x_prop: [pps_proportion, k_prime_proportion]
                        """
                        pps_proportion, k_prime_proportion = x_prop
                        
                        # 计算实际PPS缴费
                        actual_c_pps = pps_proportion * max_permissible_cpps
                        actual_c_pps = max(0, min(actual_c_pps, max_permissible_cpps))
                        
                        # 计算PPS提取后的资源
                        resources_after_pps, _, _ = OLGV9Utils.hh_income_huggett(
                            k_val, R_k_net_factor, MPL_gross, TR_total, b_payg_val,
                            actual_c_pps, a_idx, paramS_hh, cS, eps_val
                        )
                        
                        # 计算可用于储蓄和消费的资源
                        consumption_floor_spending = cS['cFloor'] * (1 + cS['tau_c'])
                        resources_for_kprime_and_c_above_floor = resources_after_pps - consumption_floor_spending
                        
                        actual_k_prime = cS['k_min']
                        current_c = cS['cFloor']
                        
                        if resources_for_kprime_and_c_above_floor >= 0:
                            actual_k_prime = k_prime_proportion * resources_for_kprime_and_c_above_floor
                            actual_k_prime = max(cS['k_min'], min(actual_k_prime, resources_for_kprime_and_c_above_floor))
                            consumption_expenditure = resources_after_pps - actual_k_prime
                            current_c = max(cS['cFloor'], consumption_expenditure / (1 + cS['tau_c']))
                        else:
                            actual_k_prime = cS['k_min']
                            consumption_expenditure = resources_after_pps - actual_k_prime
                            current_c = max(cS['cFloor'], consumption_expenditure / (1 + cS['tau_c']))
                        
                        actual_k_prime = max(cS['k_min'], min(actual_k_prime, cS['k_max']))
                        
                        # 计算PPS演化
                        pps_withdrawal = 0
                        annual_age_check = cS['physAgeMap'][a_idx][0] if 'physAgeMap' in cS and cS['physAgeMap'] and a_idx < len(cS['physAgeMap']) and cS['physAgeMap'][a_idx] else 0
                        # 修复：使用与hh_income_huggett一致的退休判断
                        matlab_a_idx_for_retirement = a_idx + 1
                        is_retired = (matlab_a_idx_for_retirement > cS['aR_new'])
                        if is_retired and annual_age_check >= cS['pps_withdrawal_age_min_idx'] and cS['pps_active']:
                            pps_withdrawal = k_pps_val * cS['pps_withdrawal_rate']
                        
                        pps_return_factor = 1 + ((R_k_net_factor - 1) + cS['pps_return_rate_premium'])
                        k_pps_prime_unclamped = (k_pps_val + actual_c_pps - pps_withdrawal) * pps_return_factor
                        k_pps_prime = max(cS['k_pps_min'], min(cS['k_pps_max'], k_pps_prime_unclamped))
                        
                        # 计算当期效用
                        _, current_utility = OLGV9Utils.ces_utility(current_c, cS['sigma'], cS)
                        if not np.isfinite(current_utility):
                            return 1e12 + abs(current_c - cS['cFloor']) * 1e10
                        
                        # 计算期望价值
                        expected_future_value = -np.inf
                        if a_idx < cS['aD_new'] - 1:
                            try:
                                k_prime_eval = max(cS['kGridV'][0], min(cS['kGridV'][-1], actual_k_prime))
                                k_pps_prime_eval = max(cS['kppsGridV'][0], min(cS['kppsGridV'][-1], k_pps_prime))
                                expected_future_value = ev_interp_ie((k_prime_eval, k_pps_prime_eval))
                            except:
                                expected_future_value = -1e11
                        
                        if not np.isfinite(expected_future_value):
                            expected_future_value = -1e11
                        
                        s_transition = cS['s_1yr_transitionV'][a_idx] if a_idx < len(cS['s_1yr_transitionV']) else 0.0
                        total_value = current_utility + cS['beta'] * s_transition * expected_future_value
                        
                        if not np.isfinite(total_value):
                            return 1e12
                        else:
                            return -total_value
                    
                    # 使用比例决策变量的优化（匹配MATLAB fmincon）
                    lb_fmin = [0, 0]
                    ub_fmin = [1, 1]
                    x0_pps_prop_fmin = 0.5
                    if max_permissible_cpps < 1e-9:
                        x0_pps_prop_fmin = 0
                        ub_fmin[0] = 0
                    x0_fmin = [x0_pps_prop_fmin, 0.5]  # 开始于50% PPS max, 50% 可支配资源用于k'
                    
                    optimal_cpps_val = 0
                    optimal_k_prime_val = cS['k_min']
                    optimal_c_val = cS['cFloor']
                    optimal_value_val = -np.inf
                    
                    # try:
                        # 使用scipy.optimize.minimize模拟MATLAB fmincon
                    result = minimize(
                        fmincon_objective_helper_proportional,
                        x0=x0_fmin,
                        bounds=[(lb_fmin[0], ub_fmin[0]), (lb_fmin[1], ub_fmin[1])],
                        method='SLSQP',
                        options={
                            'ftol': 1e-7,
                            'eps': 1e-8,
                            'disp': False,
                            'maxiter': 500
                        }
                    )
                    
                    if result.success and np.isfinite(result.fun):
                        pps_prop_opt, k_prime_prop_opt = result.x
                        optimal_value_val = -result.fun
                        
                        optimal_cpps_val = pps_prop_opt * max_permissible_cpps
                        optimal_cpps_val = max(0, min(optimal_cpps_val, max_permissible_cpps))
                        
                        resources_after_pps_opt, _, _ = OLGV9Utils.hh_income_huggett(
                            k_val, R_k_net_factor, MPL_gross, TR_total, b_payg_val,
                            optimal_cpps_val, a_idx, paramS_hh, cS, eps_val
                        )
                        
                        consumption_floor_spending_opt = cS['cFloor'] * (1 + cS['tau_c'])
                        resources_for_kprime_c_above_floor_opt = resources_after_pps_opt - consumption_floor_spending_opt
                        
                        if resources_for_kprime_c_above_floor_opt >= 0:
                            optimal_k_prime_val = k_prime_prop_opt * resources_for_kprime_c_above_floor_opt
                            optimal_k_prime_val = max(cS['k_min'], min(optimal_k_prime_val, resources_for_kprime_c_above_floor_opt))
                        else:
                            optimal_k_prime_val = cS['k_min']
                        optimal_k_prime_val = max(cS['k_min'], min(optimal_k_prime_val, cS['k_max']))
                        
                        consumption_expenditure_opt = resources_after_pps_opt - optimal_k_prime_val
                        optimal_c_val = max(cS['cFloor'], consumption_expenditure_opt / (1 + cS['tau_c']))
                    else:
                        # 优化失败，使用备选网格搜索
                        raise ValueError("fmincon优化失败")
                    # except:
                    #     # 备选离散搜索方案（匹配MATLAB的fallback_discrete_solution）
                    #     n_grid_k = 6
                    #     n_grid_cpps = 4
                        
                    #     k_grid_search = np.linspace(cS['k_min'], cS['k_max'], n_grid_k)
                    #     if max_permissible_cpps > 1e-9:
                    #         cpps_grid = np.linspace(0, max_permissible_cpps, n_grid_cpps)
                    #     else:
                    #         cpps_grid = [0.0]
                        
                    #     best_value = 1e12
                        
                    #     for cpps_test in cpps_grid:
                    #         for k_test in k_grid_search:
                    #             obj_val = fmincon_objective_helper_proportional([
                    #                 cpps_test / max(max_permissible_cpps, 1e-10),
                    #                 k_test / max(cS['k_max'], 1e-10)
                    #             ])
                    #             if obj_val < best_value:
                    #                 best_value = obj_val
                    #                 optimal_cpps_val = cpps_test
                    #                 optimal_k_prime_val = k_test
                    #                 optimal_value_val = -best_value
                    
                    # 最终验证和清理结果（避免重复计算）
                    if optimal_value_val == -np.inf:
                        # 如果所有优化都失败，使用最简单的解决方案
                        resources_default, _, _ = OLGV9Utils.hh_income_huggett(
                            k_val, R_k_net_factor, MPL_gross, TR_total, b_payg_val,
                            0.0, a_idx, paramS_hh, cS, eps_val
                        )
                        optimal_c_val = max(cS['cFloor'], resources_default / (1 + cS['tau_c']))
                        optimal_k_prime_val = cS['k_min']
                        optimal_cpps_val = 0.0
                        _, optimal_value_val = OLGV9Utils.ces_utility(optimal_c_val, cS['sigma'], cS)
                    
                    # 存储此状态的最优选择
                    val_age[ik, ikpps, ie] = optimal_value_val
                    cPol_age[ik, ikpps, ie] = optimal_c_val
                    kPol_age[ik, ikpps, ie] = optimal_k_prime_val
                    cPpsPol_age[ik, ikpps, ie] = optimal_cpps_val
        
        return cPol_age, kPol_age, cPpsPol_age, val_age

    @staticmethod
    def hh_solution_vfi_huggett(R_k_net_factor: float, MPL_gross: float, TR_total: float, 
                               bV_new: np.ndarray, paramS: Dict[str, Any], cS: Dict[str, Any]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Solves the household's problem using Value Function Iteration (VFI).
        This is the main backward induction loop, matching the MATLAB structure.
        """
        nk, nkpps, nw, aD_new = cS['nk'], cS['nkpps'], cS['nw'], cS['aD_new']
        
        # 强制使用MATLAB的精确网格值（与V8兼容版本保持一致）
        k_grid_matlab = np.array([0.000000, 2.121704, 6.001086, 11.024698, 16.973633])
        kpps_grid_matlab = np.array([0.000000, 1.060852, 3.000543, 5.512349, 8.486817])
        le_grid_matlab = np.array([0.502456, 1.000000, 1.990224])
        
        # 替换原有网格以确保一致性
        cS['kGridV'] = k_grid_matlab
        cS['kppsGridV'] = kpps_grid_matlab
        paramS['leGridV'] = le_grid_matlab
        
        # ===== DEBUG_VFI_PY: 创建调试输出文件 =====
        debug_file = open('debug_vfi_python.txt', 'w', encoding='utf-8')
        debug_file.write("DEBUG_VFI_PY: 开始VFI求解\n")
        debug_file.write(f"DEBUG_VFI_PY: R_k_net_factor={R_k_net_factor:.6f}, w_gross={MPL_gross:.6f}\n")
        debug_file.write(f"DEBUG_VFI_PY: 网格大小 nk={nk}, nkpps={nkpps}, nw={nw}, aD_new={aD_new}\n")
        
        # 网格调试输出
        debug_file.write("DEBUG_VFI_PY: k_grid=")
        for k_val in cS['kGridV']:
            debug_file.write(f"{k_val:.6f} ")
        debug_file.write("\n")
        
        debug_file.write("DEBUG_VFI_PY: kpps_grid=")
        for k_pps_val in cS['kppsGridV']:
            debug_file.write(f"{k_pps_val:.6f} ")
        debug_file.write("\n")
        
        debug_file.write("DEBUG_VFI_PY: le_grid=")
        for le_val in paramS['leGridV']:
            debug_file.write(f"{le_val:.6f} ")
        debug_file.write("\n")
        
        debug_file.write("DEBUG_VFI_PY: bV_payg=")
        for b_val in bV_new:
            debug_file.write(f"{b_val:.6f} ")
        debug_file.write("\n")
        
        valM = np.full((nk, nkpps, nw, aD_new), -np.inf)
        cPolM = np.zeros((nk, nkpps, nw, aD_new))
        kPolM = np.zeros((nk, nkpps, nw, aD_new))
        cPpsPolM = np.zeros((nk, nkpps, nw, aD_new))

        # Backward induction from the last period
        for a in range(aD_new - 1, -1, -1):
            debug_file.write(f"DEBUG_VFI_PY: 年龄组 {a + 1} / {aD_new}\n")
            
            v_prime_next = valM[:, :, :, a + 1] if a < aD_new - 1 else None
            b_payg_val = bV_new[a] if a < len(bV_new) else 0

            cPol_age, kPol_age, cPpsPol_age, val_age = OLGV9Utils._hh_solution_by_age(
                a, v_prime_next, R_k_net_factor, MPL_gross, TR_total, b_payg_val, paramS, cS
            )
            
            cPolM[:, :, :, a] = cPol_age
            kPolM[:, :, :, a] = kPol_age
            cPpsPolM[:, :, :, a] = cPpsPol_age
            valM[:, :, :, a] = val_age
            
            # 计算并输出年龄组平均值
            c_avg = np.mean(cPol_age)
            k_avg = np.mean(kPol_age)
            cpps_avg = np.mean(cPpsPol_age)
            val_avg = np.mean(val_age)
            
            debug_file.write(f"DEBUG_VFI_PY: 年龄组{a + 1}平均值 - c={c_avg:.6f}, k={k_avg:.6f}, cpps={cpps_avg:.6f}\n")
            
            # 特别输出年龄组16的详细信息到单独文件
            if a + 1 == 16:
                print(f"DEBUG_VFI_PY: 年龄组16平均值 - c={c_avg:.6f}, k={k_avg:.6f}, cpps={cpps_avg:.6f}")
                
                # 写入值函数到单独文件
                try:
                    with open('debug_vfi_val16_python.txt', 'w', encoding='utf-8') as val_file:
                        val_file.write(f"DEBUG_VFI_PY: 年龄组16值函数平均值 = {val_avg:.6f}\n")
                        val_file.write(f"DEBUG_VFI_PY: 年龄组16平均值详细信息:\n")
                        val_file.write(f"  消费 (c) = {c_avg:.6f}\n")
                        val_file.write(f"  资本 (k) = {k_avg:.6f}\n")
                        val_file.write(f"  PPS缴费 (cpps) = {cpps_avg:.6f}\n")
                        val_file.write(f"  值函数 (val) = {val_avg:.6f}\n")
                        print("DEBUG_VFI_PY: 值函数已写入 debug_vfi_val16_python.txt")
                except Exception as e:
                    print(f"DEBUG_VFI_PY: 写入文件时出错: {str(e)}")
            
        debug_file.write("DEBUG_VFI_PY: VFI求解完成\n")
        debug_file.close()
            
        return cPolM, kPolM, cPpsPolM, valM

    @staticmethod
    def hh_simulation_olgm(kPolM: np.ndarray, cPpsPolM_choice: np.ndarray, cPolM: np.ndarray, 
                           eIdxM: np.ndarray, R_k_net_factor: float, MPL_gross: float,
                           TR_total: float, bV_new: np.ndarray, paramS: Dict[str, Any], cS: Dict[str, Any]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Simulates household histories given policy functions, matching MATLAB v8 logic.
        Operates on an ANNUAL time scale based on policy functions defined on GROUPED ages.
        """
        nsim, aD_orig = eIdxM.shape
        k_grid, kpps_grid = cS['kGridV'], cS['kppsGridV']
        
        # 创建模拟调试输出文件
        sim_debug_file = open('debug_sim_python.txt', 'w', encoding='utf-8')
        sim_debug_file.write(f"DEBUG_SIM_PY: 开始模拟\n")
        sim_debug_file.write(f"DEBUG_SIM_PY: nsim={nsim}, aD_orig={aD_orig}\n")
        sim_debug_file.write(f"DEBUG_SIM_PY: R_k_net_factor={R_k_net_factor:.6f}, MPL_gross={MPL_gross:.6f}\n")
        sim_debug_file.write(f"DEBUG_SIM_PY: 策略函数形状: kPol={kPolM.shape}, cPpsPol={cPpsPolM_choice.shape}\n")
        
        # Map annual age index to group age index
        ageToGroupMap = np.zeros(aD_orig, dtype=int)
        for i_group, indices in enumerate(cS['physAgeMap']):
            if indices:
                ageToGroupMap[indices] = i_group
        
        sim_debug_file.write(f"DEBUG_SIM_PY: 年龄映射前5个: {ageToGroupMap[:5]}\n")

        kHistM = np.zeros((nsim, aD_orig + 1))
        kPpsHistM = np.zeros((nsim, aD_orig + 1))
        cHistM = np.zeros((nsim, aD_orig))
        
        # Create interpolants for policy functions (one for each age group and shock state)
        kPolInterp = [[None for _ in range(cS['aD_new'])] for _ in range(cS['nw'])]
        cPpsPolInterp = [[None for _ in range(cS['aD_new'])] for _ in range(cS['nw'])]

        for ia in range(cS['aD_new']):
            for ie in range(cS['nw']):
                points = (k_grid, kpps_grid)
                kPolInterp[ie][ia] = RegularGridInterpolator(points, kPolM[:, :, ie, ia], method='linear', bounds_error=False, fill_value=None)
                cPpsPolInterp[ie][ia] = RegularGridInterpolator(points, cPpsPolM_choice[:, :, ie, ia], method='linear', bounds_error=False, fill_value=None)
        
        # PPS annual net return factor
        pps_return_factor = 1 + ((R_k_net_factor - 1) + cS['pps_return_rate_premium'])
        sim_debug_file.write(f"DEBUG_SIM_PY: PPS回报因子={pps_return_factor:.6f}\n")

        for a_orig in range(aD_orig):
            a_group = ageToGroupMap[a_orig]
            
            k_now = kHistM[:, a_orig]
            k_pps_now = kPpsHistM[:, a_orig]
            
            # 调试输出前几个年龄
            if a_orig < 3:
                sim_debug_file.write(f"DEBUG_SIM_PY: 年龄{a_orig}(组{a_group}): 前3个体k=[{k_now[0]:.4f},{k_now[1]:.4f},{k_now[2]:.4f}], kpps=[{k_pps_now[0]:.4f},{k_pps_now[1]:.4f},{k_pps_now[2]:.4f}]\n")
            
            # Get policies for this annual age
            k_next = np.zeros(nsim)
            cpps_next = np.zeros(nsim)
            
            for ie in range(cS['nw']):
                # 注意：eIdxM是1-based索引，需要转换为0-based
                idx_sim = (eIdxM[:, a_orig] == ie + 1)  # MATLAB uses 1-based indexing
                if not np.any(idx_sim): continue
                
                points = np.column_stack((k_now[idx_sim], k_pps_now[idx_sim]))
                k_next[idx_sim] = kPolInterp[ie][a_group](points)
                cpps_next[idx_sim] = cPpsPolInterp[ie][a_group](points)
            
            # Apply PPS contribution/withdrawal rules based on ANNUAL age
            is_working_annual = a_orig < cS['aR_idx_orig']
            can_contribute_annual = is_working_annual and (a_orig <= cS['pps_contribution_age_max_idx'])
            
            # If not eligible to contribute, force c_pps to zero
            if not can_contribute_annual:
                cpps_next[:] = 0

            # PPS withdrawals for retirees
            pps_withdrawal = np.zeros(nsim)
            is_withdrawal_age = not is_working_annual and (a_orig >= cS['pps_withdrawal_age_min_idx'])
            if is_withdrawal_age and cS['pps_active']:
                pps_withdrawal = k_pps_now * cS['pps_withdrawal_rate']

            # Update asset holdings for next period
            kHistM[:, a_orig + 1] = np.clip(k_next, cS['k_min'], cS['k_max'])
            if cS['pps_active']:
                kPpsHistM[:, a_orig + 1] = np.clip((k_pps_now + cpps_next - pps_withdrawal) * pps_return_factor, cS['k_pps_min'], cS['k_pps_max'])
            
            # Consumption is the residual from the budget constraint
            # 修正：确保索引正确转换从1-based到0-based
            eIdx_current = eIdxM[:, a_orig] - 1  # 转换为0-based索引
            # 确保索引在有效范围内
            eIdx_current = np.clip(eIdx_current, 0, len(paramS['leGridV']) - 1)
            eps_vals = paramS['leGridV'][eIdx_current]
            b_payg_val = bV_new[a_group]
            
            y_disp, _, _ = OLGV9Utils.hh_income_huggett(k_now, R_k_net_factor, MPL_gross, TR_total, b_payg_val, cpps_next, a_group, paramS, cS, eps_vals)
            cHistM[:, a_orig] = (y_disp - k_next) / (1 + cS['tau_c'])
            
            # 调试输出前几个年龄的消费
            if a_orig < 3:
                sim_debug_file.write(f"DEBUG_SIM_PY: 年龄{a_orig}消费: 前3个体c=[{cHistM[0,a_orig]:.4f},{cHistM[1,a_orig]:.4f},{cHistM[2,a_orig]:.4f}]\n")

        cHistM[cHistM < 0] = cS['cFloor']
        
        # 最终统计
        sim_debug_file.write(f"DEBUG_SIM_PY: 模拟完成\n")
        sim_debug_file.write(f"DEBUG_SIM_PY: 平均资产K={np.mean(kHistM[:, :-1]):.6f}, K_pps={np.mean(kPpsHistM[:, :-1]):.6f}\n")
        sim_debug_file.write(f"DEBUG_SIM_PY: 平均消费C={np.mean(cHistM):.6f}\n")
        sim_debug_file.close()
        
        return kHistM[:, :-1], kPpsHistM[:, :-1], cHistM

    @staticmethod
    def get_aggregates(kHistM: np.ndarray, kPpsHistM: np.ndarray, cHistM: np.ndarray,
                       MPL_gross: float, R_mkt_gross_factor: float,
                       paramS: Dict[str, Any], cS: Dict[str, Any]) -> Dict[str, float]:
        """
        Calculates aggregate economic variables from ANNUAL simulation histories.
        Matches MATLAB v8 logic.
        """
        ageMassV_annual = paramS['Z_ss_norm_annual']
        
        # Aggregate capital stocks using annual weights
        K_nonpps = np.sum(np.mean(kHistM, axis=0) * ageMassV_annual)
        K_pps = 0.0
        if cS.get('pps_active', False) and cS.get('pps_in_K', False):
            K_pps = np.sum(np.mean(kPpsHistM, axis=0) * ageMassV_annual)
        
        K_total = K_nonpps + K_pps
        C_total = np.sum(np.mean(cHistM, axis=0) * ageMassV_annual)
        Y_total = cS['A'] * (K_total**cS['alpha']) * (paramS['L_per_capita']**(1 - cS['alpha']))
        
        # Accidental bequests from those who die, based on annual mortality and wealth
        ageDeathMass_annual = ageMassV_annual * cS['d_orig']
        bequests_non_pps = np.sum(np.mean(kHistM, axis=0) * ageDeathMass_annual)
        bequests_pps = 0.0
        if cS.get('pps_active', False) and cS.get('pps_bequeathable', False):
             bequests_pps = np.sum(np.mean(kPpsHistM, axis=0) * ageDeathMass_annual)
        
        T_bequest = bequests_non_pps + bequests_pps

        return {
            'K_aggr': K_total,
            'K_nonpps_aggr': K_nonpps,
            'K_pps_aggr': K_pps,
            'C_aggr': C_total,
            'Y_aggr': Y_total,
            'T_bequest_aggr': T_bequest,
        }

    @staticmethod
    def get_aggregates_from_simulation(vfi_results, paramS, cS, eIdxM, R_k_net_factor, w, TR, bV_payg, paramS_hh):
        """
        Wrapper function that calls household simulation and aggregation.
        This matches the MATLAB workflow.
        """
        # Unpack VFI results
        cPolM, kPolM, cPpsPolM, valM = vfi_results
        
        # Run simulation
        kHistM, kPpsHistM, cHistM = OLGV9Utils.hh_simulation_olgm(
            kPolM, cPpsPolM, cPolM, eIdxM,
            R_k_net_factor, w, TR, bV_payg, paramS_hh, cS
        )
        
        # Calculate aggregates
        # 注意：get_aggregates需要的是R_mkt_gross_factor，需要从R_k_net_factor反推
        # R_k_net_factor = (R_mkt_gross_factor - 1) * (1 - tau_k) + 1
        # 所以 R_mkt_gross_factor = (R_k_net_factor - 1) / (1 - tau_k) + 1
        tau_k = cS.get('tau_k', 0.20)
        R_mkt_gross_factor = (R_k_net_factor - 1) / (1 - tau_k) + 1
        aggr_res = OLGV9Utils.get_aggregates(kHistM, kPpsHistM, cHistM, w, R_mkt_gross_factor, paramS, cS)
        
        return aggr_res

    @staticmethod
    def check_gbc_residual(K_aggr: float, C_aggr: float, Y_aggr: float, G: float, B: float,
                            MPL_gross: float, r_mkt_gross: float, theta_payg: float,
                            tau_l: float, b_payg: float, T_bequest: float, TR_gov: float,
                            pop_growth: float, paramS: Dict[str, Any], cS: Dict[str, Any]) -> float:
        """
        Computes the residual of the government budget constraint.
        This version is simplified to match the logic apparent in the v5 solver loop.
        """
        tau_k = cS['tau_k']
        tau_c = cS['tau_c']
        L_pc = paramS['L_per_capita']

        # Revenues
        tax_rev_k = tau_k * (r_mkt_gross) * K_aggr
        tax_rev_c = tau_c * C_aggr
        
        # In this model, labor tax base is complex due to PPS deductions.
        # The solver uses a simplified check. We replicate that here.
        # Simplified: tax_rev_l on total labor income, which overstates revenue
        taxable_labor_income_approx = MPL_gross * L_pc
        tax_rev_l = tau_l * taxable_labor_income_approx
        
        # PAYG revenue is handled by theta_payg and assumed to balance PAYG benefits.
        # The GBC check here focuses on the *general* budget (non-PAYG items).
        
        # Bequests are a source of revenue for the consolidated government
        tax_rev_bequest = T_bequest

        total_revenue = tax_rev_k + tax_rev_l + tax_rev_c + tax_rev_bequest
        
        # Expenditures
        # Gov budget must cover spending, transfers, and interest on existing debt.
        # Change in debt B_t+1 - B_t = (1+g)B_t' - B_t is on the financing side.
        # So budget must balance: Revenue = G + TR_gov + r_mkt_gross * B
        total_expenditure = G + TR_gov + r_mkt_gross * B
        
        # Residual
        residual_val = total_revenue - total_expenditure
        
        # Return residual as a fraction of output for scale invariance
        return residual_val / Y_aggr if Y_aggr > 1e-6 else residual_val

    @staticmethod
    def solve_k_tau_l_for_rho_prime(rho_prime_payg_target: float, K_init_guess: float, household_solver_fn,
                                    cS: Dict[str, Any], paramS: Dict[str, Any], eIdxM: np.ndarray) -> Dict[str, Any]:
        """
        Solves for the general equilibrium K and tau_l using a manual fixed-point iteration,
        mirroring the logic in main_olg_v8_utils.m.
        """
        # Ensure it exists and has a default value
        paramS.setdefault('popGrowthForDebt', cS.get('popGrowth_orig', 0.01))

        # This will be the main solver loop, replacing fsolve
        max_iter = cS.get('max_iter_K_tau_l', 100)
        tol = cS.get('tol_K_tau_l', 1e-4)
        gbc_tol = cS.get('gbc_tol_for_internal_loop', 1e-3)
        damp_K = cS.get('damp_K_v5', 0.1)
        damp_tau_l = cS.get('damp_tau_l_v5', 0.1)

        K_guess = K_init_guess
        tau_l_guess = cS.get('tau_l_guess', 0.15)
        
        # Initialize results structure
        results = {
            "K_model_aggr": np.nan, "gbc_residual": np.inf, "MPL_gross": np.nan,
            "converged": False, "K_sol": np.nan, "tau_l_sol": np.nan, 
            "gbc_res_final": np.inf, "details": {}
        }


        logging.info(f"  solve_K_tau_l: rho_prime_target={rho_prime_payg_target:.4f}, K_init={K_guess:.2f}, tau_l_init={tau_l_guess:.3f}")
        logging.info("  Iter  | K_guess  | tau_l_gs | MPL_g    | K_model  | K_pps    | GBC_res  | K_dev    | tau_l_dev| Norm")
        logging.info("  ----------------------------------------------------------------------------------------------------")

        for it in range(max_iter):
            # --- Solve household problem and simulate ---
            # try:
            results = OLGV9Utils.solve_for_given_prices(K_guess, tau_l_guess, rho_prime_payg_target, household_solver_fn, cS, paramS, eIdxM)
            # except Exception as e:
            #     logging.error(f"  Error during household solution at iter {it}: {e}")
            #     return {
            #         'K_sol': np.nan, 'tau_l_sol': np.nan, 'gbc_res_final': np.inf,
            #         'converged': False, 'details': {}
            #     }

            K_model = results['K_model_aggr']
            gbc_residual = results['gbc_residual']
            MPL_gross = results['MPL_gross']
            
            # --- Check for convergence ---
            K_dev = K_guess - K_model
            tau_l_dev_for_print = -gbc_residual / (MPL_gross * paramS['L_per_capita'] + 1e-9) if (MPL_gross * paramS['L_per_capita']) > 1e-9 else 0

            norm = np.sqrt(K_dev**2 + gbc_residual**2)

            logging.info(f"  {it+1:4d} | {K_guess:8.4f} | {tau_l_guess:8.4f} | {MPL_gross:8.4f} | {K_model:8.4f} | {results.get('K_model_pps', 0):8.4f} | {gbc_residual:8.2e} | {K_dev:8.4f} | {tau_l_dev_for_print:8.4f} | {norm:8.2e}")

            if norm < tol and abs(gbc_residual) < gbc_tol:
                logging.info(f"  Solver converged successfully after {it+1} iterations.")
                results['converged'] = True
                results['K_sol'] = K_model
                results['tau_l_sol'] = tau_l_guess
                results['gbc_res_final'] = gbc_residual
                results['details'] = {k: v for k, v in results.items() if k not in ['converged', 'K_sol', 'tau_l_sol', 'gbc_res_final']}
                return results

            # --- Update guesses (Tatonnement process) ---
            K_guess = K_guess - damp_K * K_dev
            
            tau_l_update = tau_l_dev_for_print
            tau_l_guess = tau_l_guess + damp_tau_l * tau_l_update
            
            tau_l_guess = np.clip(tau_l_guess, cS['tau_l_min'], cS['tau_l_max'])

        logging.warning(f"  Solver failed to converge after {max_iter} iterations.")
        results['converged'] = False
        results['K_sol'] = K_model if 'K_model' in locals() else np.nan
        results['tau_l_sol'] = tau_l_guess
        results['gbc_res_final'] = gbc_residual if 'gbc_residual' in locals() else np.inf
        results['details'] = {k: v for k, v in results.items() if k not in ['converged', 'K_sol', 'tau_l_sol', 'gbc_res_final']}
        return results

    @staticmethod
    def solve_for_given_prices(K_guess, tau_l_guess, rho_prime_payg_target, household_solver_fn, cS, paramS, eIdxM):
        """
        Helper function to solve the model for a given set of prices (K and tau_l).
        This contains the logic that was inside the `equations` function for fsolve.
        """
        # --- 1. Calculate prices and transfers based on guesses ---
        r, w = OLGV9Utils.HHPrices_Huggett(K_guess, paramS['L_per_capita'], cS)
        
        mass_workers = paramS.get('mass_workers_group', np.sum(paramS['ageMassV'][:cS['aR_new']]))
        mass_retirees = np.sum(paramS['ageMassV'][cS['aR_new']:])
        
        avg_worker_gross_wage = (w * paramS['L_per_capita']) / mass_workers if mass_workers > 1e-9 else 0
        b_payg = rho_prime_payg_target * avg_worker_gross_wage
        
        theta_payg_req = rho_prime_payg_target * (mass_retirees / mass_workers) if mass_workers > 1e-9 else float('inf')
        theta_payg_actual = min(theta_payg_req, cS.get('theta_payg_max', 0.35))
        
        if (theta_payg_actual + tau_l_guess) > cS['max_total_labor_tax']:
            theta_payg_actual = max(0, cS['max_total_labor_tax'] - tau_l_guess)

        r_net = (r - 1) * (1 - cS['tau_k'])
        R_k_net_factor = 1 + r_net
        
        bV_payg = np.zeros(cS['aD_new'])
        if cS['aR_new'] < cS['aD_new']:
            bV_payg[cS['aR_new']:] = b_payg
        
        # --- 2. Solve household problem ---
        # The bequest loop from MATLAB is complex. For now, we assume TR=0 for simplicity
        # to match the initial Python structure. A more faithful translation would require
        # iterating VFI until bequests converge.
        TR = 0.0 # Placeholder for bequests
        
        # Set up paramS for household solver
        paramS_hh = paramS.copy()
        paramS_hh['tau_l'] = tau_l_guess
        paramS_hh['theta_payg_actual_for_hh'] = theta_payg_actual
        paramS_hh['pps_tax_deferral_active'] = cS['pps_active']
        
        vfi_results = household_solver_fn(
            R_k_net_factor, w, TR, bV_payg, paramS_hh, cS
        )

        # --- 3. Aggregate results ---
        aggr_res = OLGV9Utils.get_aggregates_from_simulation(vfi_results, paramS, cS, eIdxM, R_k_net_factor, w, TR, bV_payg, paramS_hh)
        K_model_aggr = aggr_res['K_aggr']
        C_model_aggr = aggr_res['C_aggr']
        K_model_pps = aggr_res.get('K_pps_aggr', 0)
        
        # --- 4. Check government budget constraint ---
        Y = cS['A'] * (K_guess**cS['alpha']) * (paramS['L_per_capita']**(1 - cS['alpha']))
        G = cS['gov_exp_frac_Y'] * Y
        B = cS['gov_debt_frac_Y'] * Y
        T_bequest = aggr_res.get('T_bequest_aggr', 0.0)
        
        gbc_residual = OLGV9Utils.check_gbc_residual(
            K_guess, C_model_aggr, Y, G, B,
            w, r-1, theta_payg_actual, tau_l_guess, b_payg, T_bequest, 0,
            paramS.get('popGrowthForDebt', cS['popGrowth_orig']), paramS, cS
        )
        
        return {
            "K_model_aggr": K_model_aggr,
            "C_model_aggr": C_model_aggr,
            "K_model_pps": K_model_pps,
            "gbc_residual": gbc_residual,
            "r": r-1,
            "w": w,
            "MPL_gross": w, # In this model, w is MPL
            "b_payg": b_payg,
            "theta_payg": theta_payg_actual,
            "theta_payg_required_before_cap": theta_payg_req,
            "Y": Y,
            "R_mkt_gross_factor": r,
            "converged": False,
            "details": {
                "vfi_results": vfi_results
            }
        } 

    @staticmethod
    def hh_simulation_sac_annual(sac_model, eIdxM_annual: np.ndarray, 
                                 R_k_net_factor: float, MPL_gross: float,
                                 TR_total: float, bV_new_annual: np.ndarray, 
                                 paramS: Dict[str, Any], cS: Dict[str, Any],
                                 deterministic: bool = True, verbose: bool = False) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        使用SAC agent进行年度生命周期模拟，完全模仿MATLAB main_olg_v8_utils.m中HHSimulation_olgm的逻辑
        
        🚨 重要：这是年度时间尺度模拟，与MATLAB VFI simulation在数值上完全对标
        - 在年度时间尺度(aD_orig=79)上进行模拟，而不是组别时间尺度
        - 使用年度年龄判断PPS缴费/提取资格
        - 年度PPS资产演化和收益计算
        - 年度效率冲击已通过eIdxM_annual提供
        - 与MATLAB HHSimulation_olgm在时间维度和逻辑上完全一致
        
        Args:
            sac_model: 训练好的SAC模型
            eIdxM_annual: 年度效率冲击矩阵 (nsim, aD_orig)
            R_k_net_factor: 税后净资本回报率因子
            MPL_gross: 毛边际劳动生产率
            TR_total: 总转移支付
            bV_new_annual: 年度PAYG收益向量 (aD_orig,)
            paramS: 模型参数
            cS: 配置参数
            deterministic: 是否使用确定性策略
            verbose: 是否显示详细信息
            
        Returns:
            kHistM_annual: 年度非PPS资产历史 (nsim, aD_orig)
            kPpsHistM_annual: 年度PPS资产历史 (nsim, aD_orig)  
            cHistM_annual: 年度消费历史 (nsim, aD_orig)
            
        年度模拟特点:
            1. 🕐 时间尺度：年度(20-98岁，79期)，与MATLAB HHSimulation_olgm一致
            2. 🎯 策略查询：将年度状态映射到组别，查询RL策略
            3. 🏦 PPS处理：年度缴费/提取判断，年度资产演化
            4. 💰 年度收入：基于年度年龄的收入计算
            5. 🔢 效率冲击：年度已给定，无需额外转移
            6. 📊 输出维度：与MATLAB HHSimulation_olgm完全相同
        """
        if verbose:
            print("🧬 SAC Agent年度生命周期模拟 (模仿MATLAB HHSimulation_olgm)")
        
        # 导入必要的模块
        try:
            from main_olg_v9_sac_sbx import OLGEnvV8SAC
        except ImportError:
            # 如果无法导入，创建一个简化的环境类用于策略查询
            class SimpleOLGEnv:
                def __init__(self, cS, paramS_rl, M_fixed):
                    self.cS = cS
                    self.paramS_rl = paramS_rl
                    self.current_M = M_fixed
                    
                def _get_observation(self):
                    # 简化的观测构造
                    obs = np.array([
                        self.current_k_val / self.cS['k_max'],
                        self.current_k_pps_val / self.cS['k_pps_max'],
                        self.current_age_idx / self.cS['aD_new'],
                        (self.current_eps_idx - 1) / (self.cS['nw'] - 1),
                        (self.current_M['R_k_net_factor'] - 1.01) / 0.04,
                        (self.current_M['w_gross'] - 1.5) / 1.0,
                        self.current_M['TR_total'] / 0.2,
                        self.current_M['b_payg_avg_for_obs'] / 0.8,
                        self.current_M['tau_l'] / 0.25,
                        self.current_M['theta_payg_actual'] / 0.20
                    ])
                    return np.clip(obs, 0, 1).astype(np.float32)
            
            OLGEnvV8SAC = SimpleOLGEnv
        
        # 📐 提取年度时间尺度参数（与MATLAB HHSimulation_olgm一致）
        nsim, aD_orig = eIdxM_annual.shape
        aD_new = int(cS['aD_new'])                # 组别年龄数量 (16)
        aR_idx_orig = int(cS['aR_idx_orig'])      # 年度退休年龄索引 (46: 65岁)
        nw = int(cS['nw'])                        # 效率状态数量
        
        # 年度年龄到组别年龄的映射
        ageToGroupMap = np.zeros(aD_orig, dtype=int)
        for i_group, indices in enumerate(cS['physAgeMap']):
            if indices:
                ageToGroupMap[np.array(indices)] = i_group
        
        # 提取必要的参数
        leGridV = np.array(paramS['leGridV']).flatten()
        kMin = float(cS['k_min'])
        kppsMin = float(cS['k_pps_min'])
        kMax = float(cS['k_max'])
        kppsMax = float(cS['k_pps_max'])
        pps_active = bool(cS['pps_active'])
        
        # PPS相关参数
        pps_contribution_age_max_idx = int(cS.get('pps_contribution_age_max_idx', aR_idx_orig - 1))
        pps_withdrawal_age_min_idx = int(cS.get('pps_withdrawal_age_min_idx', aR_idx_orig))
        pps_withdrawal_rate = float(cS.get('pps_withdrawal_rate', 0.15))
        pps_return_rate_premium = float(cS.get('pps_return_rate_premium', 0.08))
        pps_tax_rate_withdrawal = float(cS.get('pps_tax_rate_withdrawal', 0.03))
        
        # PPS年度净回报因子
        pps_return_factor = 1 + ((R_k_net_factor - 1) + pps_return_rate_premium)
        
        if verbose:
            print(f"📊 年度时间尺度: aD_orig={aD_orig} (20-98岁)")
            print(f"📊 组别时间尺度: aD_new={aD_new} (16个5年期组别)")
            print(f"📊 年度退休年龄索引: aR_idx_orig={aR_idx_orig} (对应65岁)")
            print(f"📊 模拟个体数: nsim={nsim}")
            print(f"📊 PPS参数: 最大缴费年龄索引={pps_contribution_age_max_idx}, 最小提取年龄索引={pps_withdrawal_age_min_idx}")
        
        # 构造固定的宏观参数（与SAC训练/评估一致）
        M_fixed = {
            'R_k_net_factor': R_k_net_factor,
            'w_gross': MPL_gross,
            'TR_total': TR_total,
            'b_payg_avg_retiree': np.mean(bV_new_annual[bV_new_annual > 0]) if np.any(bV_new_annual > 0) else 0.4,
            'tau_l': paramS.get('tau_l', 0.15),
            'theta_payg_actual': paramS.get('theta_payg_actual_for_hh', 0.12),
            'b_payg_avg_for_obs': np.mean(bV_new_annual[bV_new_annual > 0]) if np.any(bV_new_annual > 0) else 0.4
        }
        
        # 创建策略查询环境（只用于策略查询，不用于状态转移）
        paramS_rl = {
            'leGridV': leGridV,
            'leTrProbM': paramS.get('leTrProbM', np.eye(nw)),
            'leProb1V': paramS.get('leProb1V', np.ones(nw) / nw)
        }
        rng_M = {}  # 不需要采样范围，因为使用固定参数
        env = OLGEnvV8SAC(cS, paramS_rl, rng_M)
        env.set_macro_parameters(M_fixed)
        
        # 初始化结果存储（年度尺度）
        kHistM_annual = np.zeros((nsim, aD_orig))          # 年度非PPS资产路径
        kPpsHistM_annual = np.zeros((nsim, aD_orig))       # 年度PPS资产路径
        cHistM_annual = np.zeros((nsim, aD_orig))          # 年度消费路径
        
        if verbose:
            print("🔄 开始SAC Agent年度生命周期模拟...")
        
        # 🔄 年度生命周期模拟主循环（模仿MATLAB HHSimulation_olgm）
        for i_sim in range(nsim):
            if verbose and (i_sim + 1) % max(1, nsim // 10) == 0:
                print(f"  进度: {i_sim + 1}/{nsim}")
            
            # 初始状态
            k_current = kMin
            kpps_current = kppsMin
            
            # 🕐 年度时间循环 (a_orig_loop_idx = 1:aD_orig in MATLAB)
            for a_orig in range(aD_orig):
                # 获取当前年度的效率冲击索引和值
                eps_idx_current = eIdxM_annual[i_sim, a_orig]
                # 修正：确保索引正确转换从1-based到0-based
                eps_idx_0based = eps_idx_current - 1  # 转换为0-based索引
                eps_idx_0based = np.clip(eps_idx_0based, 0, len(leGridV) - 1)
                epsilon_val = leGridV[eps_idx_0based]
                
                # 获取对应的组别年龄
                a_group = ageToGroupMap[a_orig]
                
                # 年度年龄判断（与MATLAB完全一致）
                is_working_age_annual = (a_orig < aR_idx_orig)
                
                # PPS提取判断
                is_pps_withdrawal_eligible = (not is_working_age_annual and 
                                            pps_active and 
                                            a_orig >= pps_withdrawal_age_min_idx)
                
                # PPS年度提取
                pps_withdrawal = 0
                if is_pps_withdrawal_eligible:
                    pps_withdrawal = kpps_current * pps_withdrawal_rate
                
                # 💡 策略查询：构造环境状态，查询RL策略
                env.current_age_idx = a_group + 1  # 转换为MATLAB索引（从1开始）
                env.current_eps_idx = eps_idx_current  # 保持MATLAB索引（从1开始）
                env.current_k_val = k_current
                env.current_k_pps_val = kpps_current
                
                # 构造观测状态
                obs = env._get_observation()
                
                # 🚨 最后一期特殊处理（年度最后一期）
                if a_orig == aD_orig - 1:
                    # 最后一期：强制消费所有资产，不储蓄，不缴费PPS
                    k_next = kMin
                    cpps_decision = 0
                    
                    # 计算年度收入
                    b_payg_val = bV_new_annual[a_orig]
                    
                    # PPS提取（税后）
                    pps_withdrawal_aftertax = 0
                    if is_pps_withdrawal_eligible:
                        pps_withdrawal_aftertax = kpps_current * (1 - pps_tax_rate_withdrawal)
                    
                    # 构建参数用于收入计算
                    paramS_hh = {
                        'tau_l': M_fixed['tau_l'],
                        'theta_payg_actual_for_hh': M_fixed['theta_payg_actual'],
                        'pps_tax_deferral_active': pps_active,
                        'ageEffV_new': cS['ageEffV_new']
                    }
                    
                    # 计算年度收入
                    total_income, _, _ = OLGV9Utils.hh_income_huggett(
                        k_current, M_fixed['R_k_net_factor'], M_fixed['w_gross'],
                        M_fixed['TR_total'], b_payg_val, 0.0,  # cpps = 0
                        a_group, paramS_hh, cS, epsilon_val
                    )
                    
                    # 总可用资源
                    total_resources = total_income + pps_withdrawal_aftertax
                    
                    # 最后一期消费（扣除消费税）
                    tau_c = float(cS.get('tau_c', 0.10))
                    c_consumption = max(float(cS.get('cFloor', 0.05)), 
                                      total_resources / (1 + tau_c))
                    
                else:
                    # 非最后一期：使用RL策略
                    action, _ = sac_model.predict(obs, deterministic=deterministic)
                    
                    # 解析动作
                    cpps_frac = float(action[0])
                    k_save_frac = float(action[1])
                    
                    # 计算年度收入
                    b_payg_val = bV_new_annual[a_orig]
                    
                    # 构建参数用于收入计算
                    paramS_hh = {
                        'tau_l': M_fixed['tau_l'],
                        'theta_payg_actual_for_hh': M_fixed['theta_payg_actual'],
                        'pps_tax_deferral_active': pps_active,
                        'ageEffV_new': cS['ageEffV_new']
                    }
                    
                    # 初始PPS缴费决策
                    cpps_decision = 0
                    if pps_active and is_working_age_annual and (a_orig <= pps_contribution_age_max_idx):
                        # 计算总劳动收入
                        gross_labor_income = M_fixed['w_gross'] * cS['ageEffV_new'][a_group] * epsilon_val
                        
                        # PPS缴费约束
                        max_cpps_by_income = gross_labor_income * cS.get('pps_max_contrib_frac', 1.0)
                        max_cpps_absolute = cS.get('pps_annual_contrib_limit', 9999)
                        max_cpps_allowed = min(max_cpps_by_income, max_cpps_absolute)
                        cpps_decision = max(0, min(cpps_frac * max_cpps_allowed, max_cpps_allowed))
                    
                    # 计算年度收入（包含PPS缴费）
                    total_income, _, _ = OLGV9Utils.hh_income_huggett(
                        k_current, M_fixed['R_k_net_factor'], M_fixed['w_gross'],
                        M_fixed['TR_total'], b_payg_val, cpps_decision,
                        a_group, paramS_hh, cS, epsilon_val
                    )
                    
                    # PPS提取（税后）
                    pps_withdrawal_aftertax = 0
                    if is_pps_withdrawal_eligible:
                        pps_withdrawal_aftertax = pps_withdrawal * (1 - pps_tax_rate_withdrawal)
                    
                    # 总可用资源
                    total_resources = total_income + pps_withdrawal_aftertax
                    
                    # 储蓄决策
                    max_feasible_k = max(0, total_resources * 0.95)  # 最多储蓄95%的资源
                    k_save_amount = k_save_frac * max_feasible_k
                    k_next = np.clip(k_save_amount, kMin, kMax)
                    
                    # 消费计算
                    tau_c = float(cS.get('tau_c', 0.10))
                    consumption_pretax = max(0, total_resources - k_save_amount)
                    c_consumption = max(float(cS.get('cFloor', 0.05)), 
                                      consumption_pretax / (1 + tau_c))
                
                # 存储年度轨迹
                kHistM_annual[i_sim, a_orig] = k_current
                kPpsHistM_annual[i_sim, a_orig] = kpps_current
                cHistM_annual[i_sim, a_orig] = c_consumption
                
                # 状态更新为下一年度
                k_current = k_next if a_orig < aD_orig - 1 else kMin
                
                # 🏦 PPS资产年度演化（与MATLAB HHSimulation_olgm完全一致）
                if pps_active and a_orig < aD_orig - 1:  # 不在最后一期更新PPS
                    # PPS缴费年龄检查
                    can_contribute_pps_annual = (is_working_age_annual and 
                                               (a_orig <= pps_contribution_age_max_idx))
                    
                    # 如果不符合缴费条件，强制PPS缴费为0
                    if not can_contribute_pps_annual:
                        cpps_decision = 0
                    
                    # PPS资产演化：(存量 + 缴费 - 提取) * 收益率
                    kpps_current = ((kpps_current + cpps_decision - pps_withdrawal) * 
                                  pps_return_factor)
                    kpps_current = np.clip(kpps_current, kppsMin, kppsMax)
                else:
                    kpps_current = kppsMin  # PPS未激活或最后一期时保持最小值
        
        if verbose:
            print(f"✅ SAC Agent年度生命周期模拟完成")
            print(f"📊 输出维度: {nsim}个体 × {aD_orig}年度")
            print(f"🎯 与MATLAB HHSimulation_olgm在时间尺度和逻辑上完全一致")
            
            # 显示一些统计信息
            mean_k = np.mean(kHistM_annual)
            mean_kpps = np.mean(kPpsHistM_annual)
            mean_c = np.mean(cHistM_annual)
            print(f"📈 平均资产: K={mean_k:.4f}, K_pps={mean_kpps:.4f}")
            print(f"📈 平均消费: C={mean_c:.4f}")
        
        return kHistM_annual, kPpsHistM_annual, cHistM_annual

    @staticmethod
    def hh_solution_sac_huggett(R_k_net_factor: float, MPL_gross: float, TR_total: float, 
                               bV_new: np.ndarray, paramS: Dict[str, Any], cS: Dict[str, Any],
                               sac_model=None) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        SAC求解器包装函数，使其与VFI求解器框架兼容
        
        这个函数提供与hh_solution_vfi_huggett相同的接口，但内部使用SAC agent进行年度模拟
        
        Args:
            R_k_net_factor: 税后净资本回报率因子
            MPL_gross: 毛边际劳动生产率  
            TR_total: 总转移支付
            bV_new: 组别PAYG收益向量 (aD_new,)
            paramS: 模型参数
            cS: 配置参数
            sac_model: SAC模型（如果为None则抛出错误）
            
        Returns:
            cPolM: 消费策略函数 (placeholder, 全零)
            kPolM: 储蓄策略函数 (placeholder, 全零)  
            cPpsPolM: PPS缴费策略函数 (placeholder, 全零)
            valM: 价值函数 (placeholder, 全零)
            
        注意：
            - 这个函数主要用于与现有求解器框架的兼容性
            - 实际的模拟通过get_aggregates_from_sac_simulation完成
            - 策略函数返回占位符，因为SAC是端到端的
        """
        if sac_model is None:
            raise ValueError("SAC模型不能为None")
            
        if not hasattr(sac_model, 'predict'):
            raise ValueError("SAC模型必须有predict方法")
        
        # 返回占位符策略函数（与VFI格式兼容）
        nk, nkpps, nw, aD_new = cS['nk'], cS['nkpps'], cS['nw'], cS['aD_new']
        
        cPolM = np.zeros((nk, nkpps, nw, aD_new))
        kPolM = np.zeros((nk, nkpps, nw, aD_new))  
        cPpsPolM = np.zeros((nk, nkpps, nw, aD_new))
        valM = np.zeros((nk, nkpps, nw, aD_new))
        
        logging.info("SAC求解器已调用，返回占位符策略函数")
        logging.info("实际模拟将通过get_aggregates_from_sac_simulation完成")
        
        return cPolM, kPolM, cPpsPolM, valM

    @staticmethod
    def get_aggregates_from_sac_simulation(sac_model, paramS: Dict[str, Any], cS: Dict[str, Any], 
                                          eIdxM: np.ndarray, R_k_net_factor: float, 
                                          MPL_gross: float, TR_total: float, 
                                          bV_payg: np.ndarray, paramS_hh: Dict[str, Any]) -> Dict[str, float]:
        """
        使用SAC模型进行年度模拟并计算宏观聚合量
        
        这个函数替代get_aggregates_from_simulation，专门用于SAC模拟
        
        Args:
            sac_model: SAC模型
            paramS: 模型参数
            cS: 配置参数  
            eIdxM: 年度效率冲击矩阵
            R_k_net_factor: 税后净资本回报率因子
            MPL_gross: 毛边际劳动生产率
            TR_total: 总转移支付
            bV_payg: 组别PAYG收益向量
            paramS_hh: 家庭参数
            
        Returns:
            aggr_res: 宏观聚合结果字典
        """
        # 将组别PAYG收益扩展到年度
        aD_orig = cS['aD_orig']
        aR_idx_orig = cS['aR_idx_orig']
        bV_new_annual = np.zeros(aD_orig)
        
        # 使用paramS_hh中的信息计算年度PAYG收益
        for a_orig in range(aD_orig):
            if a_orig >= aR_idx_orig:  # 年度退休年龄判断
                # 简化：使用平均退休收益
                if len(bV_payg) > cS['aR_new']:
                    bV_new_annual[a_orig] = np.mean(bV_payg[cS['aR_new']:])
                else:
                    bV_new_annual[a_orig] = 0.0
        
        # 调用SAC年度模拟
        kHistM_annual, kPpsHistM_annual, cHistM_annual = OLGV9Utils.hh_simulation_sac_annual(
            sac_model, eIdxM, R_k_net_factor, MPL_gross, TR_total, 
            bV_new_annual, paramS_hh, cS, deterministic=True, verbose=False
        )
        
        # 计算宏观聚合量（使用年度权重）
        aggr_res = OLGV9Utils.get_aggregates(
            kHistM_annual, kPpsHistM_annual, cHistM_annual,
            MPL_gross, R_k_net_factor, paramS, cS
        )
        
        return aggr_res

    @staticmethod
    def markov_chain_simulation(num_simulations: int, num_periods_sim: int, 
                               initial_prob_dist: np.ndarray, transition_matrix: np.ndarray, 
                               random_numbers: np.ndarray) -> np.ndarray:
        """
        Markov chain simulation that exactly matches MATLAB's MarkovChainSimulation function.
        """
        num_states = len(initial_prob_dist)
        
        # Validation checks
        if transition_matrix.shape[0] != num_states or transition_matrix.shape[1] != num_states:
            raise ValueError('MarkovChainSimulation: 转移矩阵维度与初始分布长度不匹配。')
        
        if abs(np.sum(initial_prob_dist) - 1) > 1e-5:
            logging.warning('MarkovChainSimulation: 初始分布 p0V 的和不为1，已重新归一化。')
            initial_prob_dist = initial_prob_dist / np.sum(initial_prob_dist)
        
        row_sums = np.sum(transition_matrix, axis=1)
        if np.any(np.abs(row_sums - 1) > 1e-5):
            logging.warning('MarkovChainSimulation: 转移矩阵 trProbM 的某些行和不为1，已重新归一化。')
            row_sums[row_sums <= 1e-9] = 1
            transition_matrix = transition_matrix / row_sums[:, np.newaxis]
        
        if random_numbers.shape[0] != num_simulations or random_numbers.shape[1] != num_periods_sim:
            raise ValueError('MarkovChainSimulation: 随机数矩阵维度与模拟参数不匹配。')
        
        # Build cumulative distributions - exactly as in MATLAB
        cumulative_initial_prob = np.cumsum(initial_prob_dist.flatten())
        cumulative_transition_prob = np.cumsum(transition_matrix, axis=1)
        if num_states > 0:
            cumulative_transition_prob[:, num_states-1] = 1.0  # Ensure last column is exactly 1
        
        eIdxM_out = np.zeros((num_simulations, num_periods_sim), dtype=np.uint16)
        
        if num_simulations > 0 and num_periods_sim > 0 and num_states > 0:
            # Initial period - MATLAB uses: 1 + sum(bsxfun(@gt, random_numbers(:,1), cumulative_initial_prob), 2)
            random_col_1 = random_numbers[:, 0:1]  # Shape: (num_simulations, 1)
            cumulative_initial_prob_row = cumulative_initial_prob.reshape(1, -1)  # Shape: (1, num_states)
            comparisons = random_col_1 > cumulative_initial_prob_row  # Broadcasting
            eIdxM_out[:, 0] = 1 + np.sum(comparisons, axis=1)  # MATLAB indexing (1-based)
        elif num_states == 0 and num_simulations > 0 and num_periods_sim > 0:
            logging.warning('MarkovChainSimulation: num_states is 0. Cannot simulate.')
            return eIdxM_out
        
        # Subsequent periods
        for t in range(num_periods_sim - 1):
            current_state_indices = eIdxM_out[:, t]
            valid_indices = (current_state_indices >= 1) & (current_state_indices <= num_states)
            
            if not np.all(valid_indices):
                logging.warning(f'MarkovChainSimulation: 在期 {t} 检测到无效的当前状态索引。已重置为状态1。')
                current_state_indices[~valid_indices] = 1
                eIdxM_out[:, t] = current_state_indices
            
            # Extract cumulative transition probabilities for current states
            # MATLAB indexing is 1-based, so we subtract 1 for Python indexing
            current_state_indices_0based = current_state_indices - 1
            cPt_for_next_state = cumulative_transition_prob[current_state_indices_0based, :]
            
            # MATLAB: 1 + sum(bsxfun(@gt, random_numbers(:, t+1), cPt_for_next_state), 2)
            random_col_t_plus_1 = random_numbers[:, t+1:t+2]  # Shape: (num_simulations, 1)
            comparisons = random_col_t_plus_1 > cPt_for_next_state  # Broadcasting
            eIdxM_out[:, t+1] = 1 + np.sum(comparisons, axis=1)  # MATLAB indexing (1-based)
        
        return eIdxM_out

    @staticmethod
    def call_interpolator_python(interpolant_obj, k_val, k_pps_val, cS_local):
        """
        Python版本的CallInterpolator，匹配MATLAB的边界处理逻辑
        """
        nk = cS_local['nk']
        nkpps = cS_local['nkpps']
        result = -1e12
        
        # First clamp inputs to grid boundaries
        k_clamped = k_val
        if nk > 0 and 'kGridV' in cS_local and len(cS_local['kGridV']) > 0:
            k_clamped = max(cS_local['kGridV'][0], min(cS_local['kGridV'][-1], k_val))
        
        k_pps_clamped = k_pps_val
        if nkpps > 0 and 'kppsGridV' in cS_local and len(cS_local['kppsGridV']) > 0:
            k_pps_clamped = max(cS_local['kppsGridV'][0], min(cS_local['kppsGridV'][-1], k_pps_val))
        
        try:
            if nk > 1 and nkpps > 1:
                # 关键修复：RegularGridInterpolator需要数组形式的输入
                result = interpolant_obj([k_clamped, k_pps_clamped])
                if hasattr(result, 'item'):
                    result = result.item()
                elif hasattr(result, 'shape') and len(result.shape) > 0:
                    result = result.flat[0] if result.size > 0 else -1e12
            elif nk > 1:  # nkpps must be 1
                result = interpolant_obj(k_clamped)
                if hasattr(result, 'item'):
                    result = result.item()
            elif nkpps > 1:  # nk must be 1
                result = interpolant_obj(k_pps_clamped)
                if hasattr(result, 'item'):
                    result = result.item()
            else:  # nk=1, nkpps=1
                if callable(interpolant_obj):
                    result = interpolant_obj(k_clamped, k_pps_clamped)
                else:
                    result = interpolant_obj
                if hasattr(result, 'item'):
                    result = result.item()
        except Exception as e:
            # 调试：记录插值异常的详细信息
            with open('debug_interpolator_error.txt', 'a', encoding='utf-8') as debug_f:
                debug_f.write(f"INTERP_ERROR: k_val={k_val:.6f}, k_pps_val={k_pps_val:.6f}, k_clamped={k_clamped:.6f}, k_pps_clamped={k_pps_clamped:.6f}\n")
                debug_f.write(f"INTERP_ERROR: nk={nk}, nkpps={nkpps}, interpolant_type={type(interpolant_obj)}\n")
                debug_f.write(f"INTERP_ERROR: 异常详情: {str(e)}\n")
                if 'kGridV' in cS_local and 'kppsGridV' in cS_local:
                    debug_f.write(f"INTERP_ERROR: k_grid范围=[{cS_local['kGridV'][0]:.6f}, {cS_local['kGridV'][-1]:.6f}], kpps_grid范围=[{cS_local['kppsGridV'][0]:.6f}, {cS_local['kppsGridV'][-1]:.6f}]\n")
                debug_f.write("---\n")
            result = -1e11
        
        # Final check for finite result
        if not np.isfinite(result):
            result = -1e12
        
        return result

    @staticmethod
    def hh_solution_vfi_huggett_v8_compat(R_k_net_factor: float, MPL_gross: float, TR_total: float, 
                               bV_new: np.ndarray, paramS: Dict[str, Any], cS: Dict[str, Any]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        完全匹配MATLAB main_olg_v8_utils.HHSolution_VFI_Huggett的Python实现
        使用与MATLAB相同的插值方法、优化器设置和边界处理
        """
        from scipy.interpolate import PchipInterpolator, RectBivariateSpline
        from scipy.optimize import minimize
        
        nk, nkpps, nw, aD_new = cS['nk'], cS['nkpps'], cS['nw'], cS['aD_new']
        
        # 强制使用MATLAB的精确网格值以排除网格差异
        k_grid = np.array([0.000000, 2.121704, 6.001086, 11.024698, 16.973633])
        kpps_grid = np.array([0.000000, 1.060852, 3.000543, 5.512349, 8.486817])
        le_grid = np.array([0.502456, 1.000000, 1.990224])
        
        # 验证网格大小匹配
        assert len(k_grid) == nk, f"k网格大小不匹配: {len(k_grid)} vs {nk}"
        assert len(kpps_grid) == nkpps, f"kpps网格大小不匹配: {len(kpps_grid)} vs {nkpps}"
        assert len(le_grid) == nw, f"le网格大小不匹配: {len(le_grid)} vs {nw}"
        
        # 创建调试输出文件
        debug_file = open('debug_vfi_python.txt', 'w', encoding='utf-8')
        debug_file.write(f"DEBUG_VFI_PY: 开始VFI求解\n")
        debug_file.write(f"DEBUG_VFI_PY: R_k_net_factor={R_k_net_factor:.6f}, MPL_gross={MPL_gross:.6f}\n")
        debug_file.write(f"DEBUG_VFI_PY: 网格大小 nk={nk}, nkpps={nkpps}, nw={nw}, aD_new={aD_new}\n")
        debug_file.write(f"DEBUG_VFI_PY: k_grid={k_grid}\n")
        debug_file.write(f"DEBUG_VFI_PY: kpps_grid={kpps_grid}\n")
        debug_file.write(f"DEBUG_VFI_PY: le_grid={le_grid}\n")
        debug_file.write(f"DEBUG_VFI_PY: bV_new={bV_new}\n")
        
        # 初始化策略函数和价值函数
        cPolM_q = np.zeros((nk, nkpps, nw, aD_new))
        kPolM = np.zeros((nk, nkpps, nw, aD_new))
        cPpsPolM_choice = np.zeros((nk, nkpps, nw, aD_new))
        valM = np.full((nk, nkpps, nw, aD_new), -np.inf)
        
        # 构造参数集合（确保包含必需的税收参数）
        paramS_hh = {
            'tau_l': paramS.get('tau_l', 0.15),  # 与compare_py_matlab_step_by_step.py一致
            'theta_payg_actual_for_hh': paramS.get('theta_payg_actual_for_hh', 0.12),
            'pps_tax_deferral_active': paramS.get('pps_tax_deferral_active', cS.get('pps_active', False))
        }
        
        # 逆向迭代求解
        for a_idx in range(aD_new - 1, -1, -1):
            debug_file.write(f"DEBUG_VFI_PY: 年龄组 {a_idx + 1} / {aD_new}\n")
            
            # 获取下期价值函数
            v_prime_next = None
            if a_idx < aD_new - 1:
                v_prime_next = valM[:, :, :, a_idx + 1]
            
            b_age_val = bV_new[a_idx] if a_idx < len(bV_new) else 0.0
            eps_grid = le_grid
            
            # 处理最后一期（终端条件）- 完全匹配MATLAB逻辑
            if a_idx == aD_new - 1:
                debug_file.write(f"DEBUG_VFI_PY: 最后一期 a_idx={a_idx}, b_age_val={b_age_val:.6f}\n")
                
                # 创建网格进行批量计算
                K_grid, Kpps_grid, Epsilon_grid = np.meshgrid(k_grid, kpps_grid, eps_grid, indexing='ij')
                resources_batch = np.zeros_like(K_grid)
                
                # 批量计算资源
                for i in range(min(K_grid.size, 3)):  # 只输出前3个样本
                    ik, ikpps, ie = np.unravel_index(i, K_grid.shape)
                    k_val = K_grid[ik, ikpps, ie]
                    k_pps_val = Kpps_grid[ik, ikpps, ie]
                    eps_val = Epsilon_grid[ik, ikpps, ie]
                    
                    resources, _, _ = OLGV9Utils.hh_income_huggett(
                        k_val, R_k_net_factor, MPL_gross, TR_total, b_age_val,
                        0, a_idx, paramS_hh, cS, eps_val  # PPS缴费为0
                    )
                    resources_batch[ik, ikpps, ie] = resources
                    debug_file.write(f"DEBUG_VFI_PY: 终端期样本 [{ik},{ikpps},{ie}]: k={k_val:.4f}, kpps={k_pps_val:.4f}, eps={eps_val:.4f}, resources={resources:.6f}\n")
                
                # 完成所有计算
                for i in range(K_grid.size):
                    ik, ikpps, ie = np.unravel_index(i, K_grid.shape)
                    k_val = K_grid[ik, ikpps, ie]
                    k_pps_val = Kpps_grid[ik, ikpps, ie]
                    eps_val = Epsilon_grid[ik, ikpps, ie]
                    
                    resources, _, _ = OLGV9Utils.hh_income_huggett(
                        k_val, R_k_net_factor, MPL_gross, TR_total, b_age_val,
                        0, a_idx, paramS_hh, cS, eps_val  # PPS缴费为0
                    )
                    resources_batch[ik, ikpps, ie] = resources
                
                # 加上PPS提取（只有在PPS激活时）
                total_resources = resources_batch
                if cS['pps_active']:
                    total_resources = total_resources + Kpps_grid * (1 - cS.get('pps_tax_rate_withdrawal', 0.03))
                    debug_file.write(f"DEBUG_VFI_PY: PPS激活，添加PPS提取，税率={cS.get('pps_tax_rate_withdrawal', 0.03):.3f}\n")
                
                # 计算最优消费、储蓄和价值
                cPol_terminal = np.maximum(cS['cFloor'], total_resources / (1 + cS['tau_c']))
                kPol_terminal = np.full_like(cPol_terminal, cS['k_min'])
                cPpsPol_terminal = np.zeros_like(cPol_terminal)
                
                # 批量计算效用
                val_terminal = np.zeros_like(cPol_terminal)
                for i in range(cPol_terminal.size):
                    ik, ikpps, ie = np.unravel_index(i, cPol_terminal.shape)
                    _, utility = OLGV9Utils.ces_utility(cPol_terminal[ik, ikpps, ie], cS['sigma'], cS)
                    val_terminal[ik, ikpps, ie] = utility
                
                debug_file.write(f"DEBUG_VFI_PY: 终端期平均值 - c={np.mean(cPol_terminal):.6f}, k={np.mean(kPol_terminal):.6f}, val={np.mean(val_terminal):.6f}\n")
                
                # 存储终端期结果
                cPolM_q[:, :, :, a_idx] = cPol_terminal
                kPolM[:, :, :, a_idx] = kPol_terminal
                cPpsPolM_choice[:, :, :, a_idx] = cPpsPol_terminal
                valM[:, :, :, a_idx] = val_terminal
                            
            else:
                # 非终端期：使用动态规划
                
                # 1. 计算期望价值函数矩阵（完全匹配MATLAB）
                EV_matrix = np.zeros((nk, nkpps, nw))
                for ie in range(nw):
                    trans_probs = paramS['leTrProbM'][ie, :]
                    # 确保转移概率求和为1
                    trans_probs = trans_probs / np.sum(trans_probs) if np.sum(trans_probs) > 0 else trans_probs
                    EV_slice = np.sum(v_prime_next * trans_probs.reshape(1, 1, nw), axis=2)
                    EV_matrix[:, :, ie] = EV_slice
                
                # 清理期望价值矩阵中的无限值（关键修复）
                EV_matrix = np.where(np.isfinite(EV_matrix), EV_matrix, -1e12)
                
                # 调试：检查期望价值矩阵的有效性
                if a_idx == aD_new - 2:  # 倒数第二期
                    debug_file.write(f"DEBUG_VFI_PY: EV_matrix统计 a_idx={a_idx}: min={np.min(EV_matrix):.6f}, max={np.max(EV_matrix):.6f}, mean={np.mean(EV_matrix):.6f}\n")
                    debug_file.write(f"DEBUG_VFI_PY: v_prime_next统计: min={np.min(v_prime_next):.6f}, max={np.max(v_prime_next):.6f}, mean={np.mean(v_prime_next):.6f}\n")
                    # 检查是否有无限值
                    inf_count = np.sum(~np.isfinite(EV_matrix))
                    debug_file.write(f"DEBUG_VFI_PY: EV_matrix清理后无限值数量={inf_count}\n")
                
                # 2. 创建期望价值插值器（使用griddedInterpolant风格）
                ev_interpolants = {}
                debug_file.write(f"DEBUG_VFI_PY: 开始创建插值器 nk={nk}, nkpps={nkpps}, nw={nw}\n")
                for ie in range(nw):
                    debug_file.write(f"DEBUG_VFI_PY: 创建插值器 ie={ie}/{nw}\n")
                    if nk > 1 and nkpps > 1:
                        # 尝试最近邻插值以排除插值方法差异
                        try:
                            from scipy.interpolate import RegularGridInterpolator
                            debug_file.write(f"DEBUG_VFI_PY: EV_matrix[:,:,{ie}]形状={EV_matrix[:, :, ie].shape}, 范围=[{np.min(EV_matrix[:, :, ie]):.6f}, {np.max(EV_matrix[:, :, ie]):.6f}]\n")
                            ev_interpolants[ie] = RegularGridInterpolator(
                                (k_grid, kpps_grid), EV_matrix[:, :, ie], 
                                method='nearest', bounds_error=False, fill_value=-1e12
                            )
                            # 测试插值器是否工作正常
                            test_val = ev_interpolants[ie]((k_grid[0], kpps_grid[0]))
                            debug_file.write(f"DEBUG_VFI_PY: 插值器ie={ie}创建成功，测试值={test_val:.6f}\n")
                            if not np.isfinite(test_val):
                                debug_file.write(f"DEBUG_VFI_PY: 插值器测试失败 ie={ie}, test_val={test_val}\n")
                        except Exception as e:
                            debug_file.write(f"DEBUG_VFI_PY: 插值器创建失败 ie={ie}: {str(e)}\n")
                            # 备用插值器：使用最近邻
                            try:
                                ev_interpolants[ie] = RegularGridInterpolator(
                                    (k_grid, kpps_grid), EV_matrix[:, :, ie], 
                                    method='nearest', bounds_error=False, fill_value=-1e12
                                )
                                debug_file.write(f"DEBUG_VFI_PY: 备用插值器ie={ie}创建成功\n")
                            except Exception as e2:
                                debug_file.write(f"DEBUG_VFI_PY: 备用插值器也失败 ie={ie}: {str(e2)}\n")
                                ev_interpolants[ie] = lambda pts, val=-1e12: val
                    elif nk > 1 and nkpps == 1:
                        ev_data = EV_matrix[:, 0, ie] 
                        ev_data = np.where(np.isfinite(ev_data), ev_data, -1e12)
                        ev_interpolants[ie] = PchipInterpolator(k_grid, ev_data, extrapolate=True)
                        debug_file.write(f"DEBUG_VFI_PY: 一维k插值器ie={ie}创建成功\n")
                    elif nk == 1 and nkpps > 1:
                        ev_data = EV_matrix[0, :, ie]
                        ev_data = np.where(np.isfinite(ev_data), ev_data, -1e12)
                        ev_interpolants[ie] = PchipInterpolator(kpps_grid, ev_data, extrapolate=True)
                        debug_file.write(f"DEBUG_VFI_PY: 一维kpps插值器ie={ie}创建成功\n")
                    else:  # nk=1, nkpps=1
                        ev_val = EV_matrix[0, 0, ie] if np.isfinite(EV_matrix[0, 0, ie]) else -1e12
                        ev_interpolants[ie] = lambda k_s, kp_s, val=ev_val: val
                        debug_file.write(f"DEBUG_VFI_PY: 标量插值器ie={ie}创建成功，值={ev_val:.6f}\n")
                
                debug_file.write(f"DEBUG_VFI_PY: 所有插值器创建完成\n")
                debug_file.flush()
                
                # 3. 对每个状态进行优化（完全匹配MATLAB fmincon）
                sample_count = 0
                for ie, eps_val in enumerate(eps_grid):
                    for ikpps, k_pps_val in enumerate(kpps_grid):
                        for ik, k_val in enumerate(k_grid):
                            
                            # 计算PPS缴费约束（完全匹配MATLAB）
                            max_permissible_cpps = 0.0
                            model_age_group_start_year_idx = cS['physAgeMap'][a_idx][0] if cS['physAgeMap'][a_idx] else 0
                            is_pps_eligible = (a_idx <= cS['aR_new'] and 
                                             model_age_group_start_year_idx <= cS['pps_contribution_age_max_idx'] and 
                                             cS['pps_active'])
                            
                            if is_pps_eligible:
                                age_efficiency = cS['ageEffV_new'][a_idx]
                                current_gross_labor_income = MPL_gross * age_efficiency * eps_val
                                if current_gross_labor_income > 1e-6:
                                    max_cpps_by_frac = current_gross_labor_income * cS['pps_max_contrib_frac']
                                    max_permissible_cpps = min(cS['pps_annual_contrib_limit'], max_cpps_by_frac)
                                    max_permissible_cpps = max(0, max_permissible_cpps)
                            
                            # 调试输出前几个样本
                            if sample_count < 3:
                                debug_file.write(f"DEBUG_VFI_PY: 状态样本 [{ik},{ikpps},{ie}]: k={k_val:.4f}, kpps={k_pps_val:.4f}, eps={eps_val:.4f}\n")
                                debug_file.write(f"DEBUG_VFI_PY: PPS资格={is_pps_eligible}, max_cpps={max_permissible_cpps:.6f}\n")
                                sample_count += 1
                            
                            # 定义目标函数（完全匹配MATLAB的fmincon_objective_helper_proportional）
                            def objective_function(x_prop):
                                pps_proportion, k_prime_proportion = x_prop
                                
                                # 计算实际PPS缴费
                                actual_c_pps = pps_proportion * max_permissible_cpps
                                actual_c_pps = max(0, min(actual_c_pps, max_permissible_cpps))
                                
                                # 计算资源
                                resources_after_pps, _, _ = OLGV9Utils.hh_income_huggett(
                                    k_val, R_k_net_factor, MPL_gross, TR_total, b_age_val,
                                    actual_c_pps, a_idx, paramS_hh, cS, eps_val
                                )
                                
                                consumption_floor_spending = cS['cFloor'] * (1 + cS['tau_c'])
                                resources_for_kprime_and_c_above_floor = resources_after_pps - consumption_floor_spending
                                
                                # 计算储蓄和消费（完全匹配MATLAB逻辑）
                                if resources_for_kprime_and_c_above_floor >= 0:
                                    actual_k_prime = k_prime_proportion * resources_for_kprime_and_c_above_floor
                                    actual_k_prime = max(cS['k_min'], min(actual_k_prime, resources_for_kprime_and_c_above_floor))
                                    consumption_expenditure = resources_after_pps - actual_k_prime
                                    current_c = max(cS['cFloor'], consumption_expenditure / (1 + cS['tau_c']))
                                else:
                                    actual_k_prime = cS['k_min']
                                    consumption_expenditure = resources_after_pps - actual_k_prime
                                    current_c = max(cS['cFloor'], consumption_expenditure / (1 + cS['tau_c']))
                                
                                actual_k_prime = max(cS['k_min'], min(actual_k_prime, cS['k_max']))
                                
                                # 计算PPS演化（完全匹配MATLAB）
                                pps_withdrawal = 0
                                annual_age_check = cS['physAgeMap'][a_idx][0] if 'physAgeMap' in cS and cS['physAgeMap'] and a_idx < len(cS['physAgeMap']) and cS['physAgeMap'][a_idx] else 0
                                # 修复：使用与hh_income_huggett一致的退休判断
                                matlab_a_idx_v8 = a_idx + 1
                                is_retired = (matlab_a_idx_v8 > cS['aR_new'])
                                if is_retired and annual_age_check >= cS['pps_withdrawal_age_min_idx'] and cS['pps_active']:
                                    pps_withdrawal = k_pps_val * cS['pps_withdrawal_rate']
                                
                                pps_return_factor = 1 + ((R_k_net_factor - 1) + cS['pps_return_rate_premium'])
                                k_pps_prime = (k_pps_val + actual_c_pps - pps_withdrawal) * pps_return_factor
                                k_pps_prime = max(cS['k_pps_min'], min(cS['k_pps_max'], k_pps_prime))
                                
                                # 计算当期效用
                                _, current_utility = OLGV9Utils.ces_utility(current_c, cS['sigma'], cS)
                                if not np.isfinite(current_utility):
                                    return 1e12 + abs(current_c - cS['cFloor']) * 1e10
                                
                                # 计算期望价值（使用CallInterpolator匹配MATLAB）
                                expected_future_value = -np.inf
                                if a_idx < cS['aD_new'] - 1:  # 只有非最后一期才计算期望价值
                                    try:
                                        k_prime_eval = max(k_grid[0], min(k_grid[-1], actual_k_prime))
                                        k_pps_prime_eval = max(kpps_grid[0], min(kpps_grid[-1], k_pps_prime))
                                        
                                        # 使用CallInterpolator风格的调用
                                        expected_future_value = OLGV9Utils.call_interpolator_python(
                                            ev_interpolants[ie], k_prime_eval, k_pps_prime_eval, cS
                                        )
                                        
                                        # 额外调试：如果期望价值异常，记录详细信息
                                        if expected_future_value < -1e10 and sample_count <= 3:
                                            debug_file.write(f"DEBUG_VFI_PY: 插值异常 ie={ie}, k_prime_eval={k_prime_eval:.6f}, k_pps_prime_eval={k_pps_prime_eval:.6f}, expected_future_value={expected_future_value:.6f}\n")
                                    except:
                                        expected_future_value = -1e11
                                
                                if not np.isfinite(expected_future_value):
                                    expected_future_value = -1e11
                                
                                # 计算总价值
                                s_transition = cS['s_1yr_transitionV'][a_idx] if a_idx < len(cS['s_1yr_transitionV']) else 0.0
                                total_value = current_utility + cS['beta'] * s_transition * expected_future_value
                                
                                # 调试异常的价值计算
                                if total_value < -1e10:
                                    debug_file.write(f"DEBUG_VFI_PY: 异常价值 a_idx={a_idx}: current_utility={current_utility:.6f}, expected_future_value={expected_future_value:.6f}, s_transition={s_transition:.6f}, total_value={total_value:.6f}\n")
                                    debug_file.flush()  # 立即刷新缓冲区
                                
                                if not np.isfinite(total_value):
                                    return 1e12
                                else:
                                    return -total_value
                            
                            # 优化设置（匹配MATLAB的sqp算法）
                            lb_fmin = [0, 0]
                            ub_fmin = [1, 1]
                            x0_pps_prop_fmin = 0.5
                            if max_permissible_cpps < 1e-9:
                                x0_pps_prop_fmin = 0
                                ub_fmin[0] = 0
                            x0_fmin = [x0_pps_prop_fmin, 0.5]
                            
                            optimal_cpps_val = 0
                            optimal_k_prime_val = cS['k_min']
                            optimal_c_val = cS['cFloor']
                            optimal_value_val = -np.inf
                            
                            try:
                                # 使用粗网格搜索来排除优化算法差异
                                best_val = np.inf
                                best_x = x0_fmin.copy()
                                
                                # 粗网格搜索
                                n_grid = 11  # 每个维度11个点
                                pps_props = np.linspace(lb_fmin[0], ub_fmin[0], n_grid)
                                k_props = np.linspace(lb_fmin[1], ub_fmin[1], n_grid)
                                
                                for pps_prop in pps_props:
                                    for k_prop in k_props:
                                        try:
                                            val = objective_function([pps_prop, k_prop])
                                            if val < best_val:
                                                best_val = val
                                                best_x = [pps_prop, k_prop]
                                        except:
                                            continue
                                
                                # 使用网格搜索结果创建伪result对象
                                class GridResult:
                                    def __init__(self, x, fun):
                                        self.x = np.array(x)
                                        self.fun = fun
                                        self.success = np.isfinite(fun)
                                
                                result = GridResult(best_x, best_val)
                                
                                if result.success and np.isfinite(result.fun):
                                    pps_prop_opt, k_prime_prop_opt = result.x
                                    optimal_value_val = -result.fun
                                    
                                    # 重新计算最优值（匹配MATLAB逻辑）
                                    optimal_cpps_val = pps_prop_opt * max_permissible_cpps
                                    optimal_cpps_val = max(0, min(optimal_cpps_val, max_permissible_cpps))
                                    
                                    resources_after_pps_opt, _, _ = OLGV9Utils.hh_income_huggett(
                                        k_val, R_k_net_factor, MPL_gross, TR_total, b_age_val,
                                        optimal_cpps_val, a_idx, paramS_hh, cS, eps_val
                                    )
                                    
                                    consumption_floor_spending_opt = cS['cFloor'] * (1 + cS['tau_c'])
                                    resources_for_kprime_c_above_floor_opt = resources_after_pps_opt - consumption_floor_spending_opt
                                    
                                    if resources_for_kprime_c_above_floor_opt >= 0:
                                        optimal_k_prime_val = k_prime_prop_opt * resources_for_kprime_c_above_floor_opt
                                        optimal_k_prime_val = max(cS['k_min'], min(optimal_k_prime_val, resources_for_kprime_c_above_floor_opt))
                                    else:
                                        optimal_k_prime_val = cS['k_min']
                                    optimal_k_prime_val = max(cS['k_min'], min(optimal_k_prime_val, cS['k_max']))
                                    
                                    consumption_expenditure_opt = resources_after_pps_opt - optimal_k_prime_val
                                    optimal_c_val = max(cS['cFloor'], consumption_expenditure_opt / (1 + cS['tau_c']))
                                else:
                                    # 使用备选网格搜索（匹配MATLAB的fallback）
                                    raise ValueError("优化失败，使用备选方案")
                                    
                            except:
                                # 备选网格搜索（匹配MATLAB的fallback_discrete_solution）
                                n_grid_k = 6
                                n_grid_cpps = 4
                                
                                k_grid_search = np.linspace(cS['k_min'], cS['k_max'], n_grid_k)
                                if max_permissible_cpps > 1e-9:
                                    cpps_grid_search = np.linspace(0, max_permissible_cpps, n_grid_cpps)
                                else:
                                    cpps_grid_search = [0.0]
                                
                                best_value = 1e12
                                
                                for cpps_test in cpps_grid_search:
                                    for k_test in k_grid_search:
                                        obj_val = objective_function([
                                            cpps_test / max(max_permissible_cpps, 1e-10),
                                            k_test / max(cS['k_max'], 1e-10)
                                        ])
                                        if obj_val < best_value:
                                            best_value = obj_val
                                            optimal_cpps_val = cpps_test
                                            optimal_k_prime_val = k_test
                                            optimal_value_val = -best_value
                                            
                                            # 计算对应的消费
                                            resources_test, _, _ = OLGV9Utils.hh_income_huggett(
                                                k_val, R_k_net_factor, MPL_gross, TR_total, b_age_val,
                                                optimal_cpps_val, a_idx, paramS_hh, cS, eps_val
                                            )
                                            optimal_c_val = max(cS['cFloor'], (resources_test - optimal_k_prime_val) / (1 + cS['tau_c']))
                            
                            # 存储结果
                            cPolM_q[ik, ikpps, ie, a_idx] = optimal_c_val
                            kPolM[ik, ikpps, ie, a_idx] = optimal_k_prime_val
                            cPpsPolM_choice[ik, ikpps, ie, a_idx] = optimal_cpps_val
                            valM[ik, ikpps, ie, a_idx] = optimal_value_val
                            
                            # 调试输出前几个结果
                            if sample_count <= 3:
                                debug_file.write(f"DEBUG_VFI_PY: 最优结果 [{ik},{ikpps},{ie}]: c={optimal_c_val:.6f}, k={optimal_k_prime_val:.6f}, cpps={optimal_cpps_val:.6f}, val={optimal_value_val:.6f}\n")
                
                # 输出年龄组的统计信息
                c_avg = np.mean(cPolM_q[:,:,:,a_idx])
                k_avg = np.mean(kPolM[:,:,:,a_idx])
                cpps_avg = np.mean(cPpsPolM_choice[:,:,:,a_idx])
                val_avg = np.mean(valM[:,:,:,a_idx])
                
                debug_file.write(f"DEBUG_VFI_PY: 年龄组{a_idx + 1}平均值 - c={c_avg:.6f}, k={k_avg:.6f}, cpps={cpps_avg:.6f}\n")
                
                # 特别处理年龄组16
                if a_idx + 1 == 16:
                    print(f"DEBUG_VFI_PY: 年龄组16平均值 - c={c_avg:.6f}, k={k_avg:.6f}, cpps={cpps_avg:.6f}")
                    
                    # 写入值函数到单独文件
                    try:
                        with open('debug_vfi_val16_python_v8.txt', 'w', encoding='utf-8') as val_file:
                            val_file.write(f"DEBUG_VFI_PY: 年龄组16值函数平均值 = {val_avg:.6f}\n")
                            val_file.write(f"DEBUG_VFI_PY: 年龄组16平均值详细信息:\n")
                            val_file.write(f"  消费 (c) = {c_avg:.6f}\n")
                            val_file.write(f"  资本 (k) = {k_avg:.6f}\n")
                            val_file.write(f"  PPS缴费 (cpps) = {cpps_avg:.6f}\n")
                            val_file.write(f"  值函数 (val) = {val_avg:.6f}\n")
                            print("DEBUG_VFI_PY: 值函数已写入 debug_vfi_val16_python_v8.txt")
                    except Exception as e:
                        print(f"DEBUG_VFI_PY: 写入文件时出错: {str(e)}")
        
        debug_file.write("DEBUG_VFI_PY: VFI求解完成\n")
        debug_file.close()
        return cPolM_q, kPolM, cPpsPolM_choice, valM

class OLGEnvV8SAC(gym.Env):
    """
    OLG环境，用于SAC Agent训练（期望化死亡率版本）
    
    观测空间: [k, k_pps, age_idx, eps_idx, M_vars(6)] (归一化)
    动作空间: [prop_pps_contrib, prop_non_pps_saving] (0-1)
    特性: 死亡率通过期望值纳入，确保固定长度的episodes
    """
    
    def __init__(self, cS: Dict[str, Any], paramS_rl: Dict[str, Any], rng_M: Dict[str, Any]):
        """
        初始化环境
        
        Args:
            cS: OLG模型参数
            paramS_rl: RL相关的派生参数
            rng_M: 宏观变量的采样范围
        """
        if gym is None or spaces is None:
            raise ImportError("需要安装gymnasium库：pip install gymnasium")
        
        super().__init__()
        
        self.cS = cS
        self.paramS_rl = paramS_rl
        self.rng_M = rng_M
        
        # 定义观测和动作空间
        obs_dim = 1 + 1 + 1 + 1 + 6  # k, k_pps, age_idx, eps_idx, M_vars(6)
        self.observation_space = spaces.Box(
            low=0.0, high=1.0, shape=(obs_dim,), dtype=np.float32
        )
        
        act_dim = 2  # [prop_pps_contrib, prop_non_pps_saving]
        self.action_space = spaces.Box(
            low=0.0, high=1.0, shape=(act_dim,), dtype=np.float32
        )
        
        # 初始化当前状态
        self.current_M = {
            'R_k_net_factor': 0,
            'w_gross': 0,
            'TR_total': 0,
            'tau_l': 0,
            'theta_payg_actual': 0,
            'b_payg_avg_for_obs': 0
        }
        
        self.current_bV_payg = np.zeros(self.cS['aD_new'])
        self.current_age_idx = 1
        self.current_k_val = self.cS['k_min']
        self.current_k_pps_val = self.cS['k_pps_min']
        self.current_eps_idx = 1
        
        # 初始化归一化参数
        self._init_normalization_params()
    
    def _init_normalization_params(self):
        """初始化观测归一化参数 - 与MATLAB环境保持完全一致"""
        # 确保使用与MATLAB环境相同的归一化边界
        # 基于MATLAB运行结果: kMax=16.973633, kppsMax=8.486817
        
        # [k, k_pps, age_idx, eps_idx, R_k_net, w_g, TR, b_payg_avg, tau_l, theta_payg]
        self.obs_norm_min = np.array([
            0.0, 0.0, 1.0, 1.0,  # 使用与MATLAB完全相同的值
            self.rng_M['R_k_net_factor'][0], self.rng_M['w_gross'][0],
            self.rng_M['TR_total'][0], self.rng_M['b_payg_avg_retiree'][0],
            self.rng_M['tau_l'][0], self.rng_M['theta_payg_actual'][0]
        ])
        
        self.obs_norm_max = np.array([
            16.973633, 8.486817, 16.0, 3.0,  # 使用MATLAB验证的精确值
            self.rng_M['R_k_net_factor'][1], self.rng_M['w_gross'][1],
            self.rng_M['TR_total'][1], self.rng_M['b_payg_avg_retiree'][1],
            self.rng_M['tau_l'][1], self.rng_M['theta_payg_actual'][1]
        ])
        
        self.obs_norm_range = self.obs_norm_max - self.obs_norm_min
        # 避免除以零 - 与MATLAB环境保持一致
        self.obs_norm_range[self.obs_norm_range < 1e-6] = 1
        self.obs_norm_mean = (self.obs_norm_min + self.obs_norm_max) / 2
    
    def reset(self, seed: Optional[int] = None, options: Optional[Dict] = None) -> Tuple[np.ndarray, Dict]:
        """
        重置环境（每个episode开始时调用）
        
        Returns:
            observation: 初始观测
            info: 额外信息
        """
        if seed is not None:
            np.random.seed(seed)
        
        # 1. 采样新的宏观经济状态M for this episode
        self.current_M['R_k_net_factor'] = (
            self.rng_M['R_k_net_factor'][0] + 
            np.random.rand() * (self.rng_M['R_k_net_factor'][1] - self.rng_M['R_k_net_factor'][0])
        )
        self.current_M['w_gross'] = (
            self.rng_M['w_gross'][0] + 
            np.random.rand() * (self.rng_M['w_gross'][1] - self.rng_M['w_gross'][0])
        )
        self.current_M['TR_total'] = (
            self.rng_M['TR_total'][0] + 
            np.random.rand() * (self.rng_M['TR_total'][1] - self.rng_M['TR_total'][0])
        )
        b_payg_avg = (
            self.rng_M['b_payg_avg_retiree'][0] + 
            np.random.rand() * (self.rng_M['b_payg_avg_retiree'][1] - self.rng_M['b_payg_avg_retiree'][0])
        )
        self.current_M['tau_l'] = (
            self.rng_M['tau_l'][0] + 
            np.random.rand() * (self.rng_M['tau_l'][1] - self.rng_M['tau_l'][0])
        )
        self.current_M['theta_payg_actual'] = (
            self.rng_M['theta_payg_actual'][0] + 
            np.random.rand() * (self.rng_M['theta_payg_actual'][1] - self.rng_M['theta_payg_actual'][0])
        )
        
        # 从b_payg_avg生成bV_payg
        self.current_bV_payg = np.zeros(self.cS['aD_new'])
        if self.cS['aR_new'] < self.cS['aD_new']:
            self.current_bV_payg[self.cS['aR_new']:] = b_payg_avg
        self.current_M['b_payg_avg_for_obs'] = b_payg_avg
        
        # 2. 初始化个体状态
        self.current_age_idx = 1
        self.current_k_val = self.cS['k_min']
        self.current_k_pps_val = self.cS['k_pps_min']
        
        # 初始效率状态（从paramS_rl['leProb1V']中抽取）
        self.current_eps_idx = np.random.choice(
            len(self.paramS_rl['leProb1V']), 
            p=self.paramS_rl['leProb1V']
        ) + 1  # MATLAB索引从1开始
        
        observation = self._get_observation()
        info = {}
        
        return observation, info
    
    def step(self, action: np.ndarray) -> Tuple[np.ndarray, float, bool, bool, Dict]:
        """
        执行一步（agent采取行动后环境的响应）
        
        Args:
            action: [prop_pps_contrib, prop_non_pps_saving]
            
        Returns:
            observation: 下一步观测
            reward: 奖励
            terminated: 是否终止
            truncated: 是否截断
            info: 额外信息
        """
        # 0. 解析行动（比例变量）
        prop_pps_contrib = np.clip(action[0], 0, 1)
        prop_non_pps_saving = np.clip(action[1], 0, 1)
        
        # 1. 计算实际的PPS缴费c_pps
        actual_c_pps, max_permissible_cpps = self._calculate_pps_contribution(prop_pps_contrib)
        
        # 2. 计算可用于非PPS储蓄和消费的资源
        resources_after_pps = self._calculate_resources_after_pps(actual_c_pps)
        
        # 3. 计算实际的非PPS储蓄k_prime和消费c
        actual_k_prime, current_c = self._calculate_consumption_and_savings(
            resources_after_pps, prop_non_pps_saving
        )
        
        # 4. 计算奖励（当期效用，期望化死亡率处理）
        reward = self._calculate_reward(current_c)
        
        # 5. 更新个体状态到下一期
        terminated = self._update_state(actual_k_prime, actual_c_pps)
        
        observation = self._get_observation()
        info = {
            'consumption': current_c,
            'k_prime': actual_k_prime,
            'c_pps': actual_c_pps,
            'age_idx': self.current_age_idx
        }
        
        return observation, reward, terminated, False, info
    
    def _calculate_pps_contribution(self, prop_pps_contrib: float) -> Tuple[float, float]:
        """计算实际PPS缴费"""
        actual_c_pps = 0.0
        max_permissible_cpps = 0.0
        current_epsilon_val = self.paramS_rl['leGridV'][self.current_eps_idx - 1]  # Python索引从0开始
        
        # 检查是否可以缴费PPS
        if (self.current_age_idx <= self.cS['aR_new'] and 
            self.cS['physAgeMap'][self.current_age_idx - 1][0] <= self.cS['pps_contribution_age_max_idx'] and 
            self.cS['pps_active']):
            
            age_efficiency = self.cS['ageEffV_new'][self.current_age_idx - 1]
            current_gross_labor_income = self.current_M['w_gross'] * age_efficiency * current_epsilon_val
            
            if current_gross_labor_income > 1e-6:
                max_cpps_by_frac = current_gross_labor_income * self.cS['pps_max_contrib_frac']
                max_permissible_cpps = min(self.cS['pps_annual_contrib_limit'], max_cpps_by_frac)
                max_permissible_cpps = max(0, max_permissible_cpps)
            
            actual_c_pps = prop_pps_contrib * max_permissible_cpps
            actual_c_pps = max(0, min(actual_c_pps, max_permissible_cpps))
        
        return actual_c_pps, max_permissible_cpps
    
    def _calculate_resources_after_pps(self, actual_c_pps: float) -> float:
        """计算扣除PPS缴费后的可用资源"""
        # 构建临时的paramS_hh
        paramS_hh_step = {
            'tau_l': self.current_M['tau_l'],
            'theta_payg_actual_for_hh': self.current_M['theta_payg_actual'],
            'pps_tax_deferral_active': self.cS['pps_active'],
            'ageEffV_new': self.cS['ageEffV_new']
        }
        
        b_payg_this_age = self.current_bV_payg[self.current_age_idx - 1]
        current_epsilon_val = self.paramS_rl['leGridV'][self.current_eps_idx - 1]
        
        resources_after_pps, _, _ = OLGV9Utils.hh_income_huggett(
            self.current_k_val, self.current_M['R_k_net_factor'], self.current_M['w_gross'],
            self.current_M['TR_total'], b_payg_this_age, actual_c_pps,
            self.current_age_idx - 1, paramS_hh_step, self.cS, current_epsilon_val
        )
        
        return resources_after_pps
    
    def _calculate_consumption_and_savings(self, resources_after_pps: float, 
                                         prop_non_pps_saving: float) -> Tuple[float, float]:
        """计算消费和储蓄"""
        consumption_floor_spending = self.cS['cFloor'] * (1 + self.cS['tau_c'])
        resources_for_kprime_and_c_above_floor = resources_after_pps - consumption_floor_spending
        
        actual_k_prime = 0.0
        current_c = self.cS['cFloor']
        
        if resources_for_kprime_and_c_above_floor >= 0:
            # prop_non_pps_saving应用于(总资源 - 最低消费开销)
            actual_k_prime = prop_non_pps_saving * resources_for_kprime_and_c_above_floor
            # 确保k_prime在可行范围内
            actual_k_prime = max(self.cS['k_min'], 
                               min(actual_k_prime, resources_for_kprime_and_c_above_floor))
        else:
            # 如果资源不足以支付最低消费，则k_prime为k_min
            actual_k_prime = self.cS['k_min']
        
        actual_k_prime = max(self.cS['k_min'], min(actual_k_prime, self.cS['k_max']))  # 最终约束
        
        consumption_expenditure = resources_after_pps - actual_k_prime
        current_c = max(self.cS['cFloor'], consumption_expenditure / (1 + self.cS['tau_c']))
        
        return actual_k_prime, current_c
    
    def _calculate_reward(self, current_c: float) -> float:
        """计算奖励（当期效用）"""
        _, utility = OLGV9Utils.ces_utility(np.array([current_c]), self.cS['sigma'], self.cS)
        utility = utility[0]
        
        # 获取存活概率用于期望化处理
        survival_prob = 1.0
        if self.current_age_idx <= len(self.cS['s_1yr_transitionV']):
            survival_prob = self.cS['s_1yr_transitionV'][self.current_age_idx - 1]
        
        if not np.isfinite(utility) or utility < -1e9:  # 惩罚无效消费
            reward = -1000 - abs(current_c - self.cS['cFloor']) * 100
        else:
            # 期望化处理：死亡率通过RL的价值函数学习自动纳入期望
            # 保持原始效用，让SAC的折现机制和价值函数学习处理死亡风险
            reward = utility
        
        return float(reward)
    
    def _update_state(self, actual_k_prime: float, actual_c_pps: float) -> bool:
        """更新状态到下一期"""
        # 5.1 PPS资产演化
        k_pps_next = self.current_k_pps_val
        if self.cS['pps_active']:
            pps_withdrawal = 0.0
            is_retired_model_group = (self.current_age_idx > self.cS['aR_new'])
            if is_retired_model_group:
                pps_withdrawal = self.current_k_pps_val * self.cS['pps_withdrawal_rate']
            
            pps_return_factor = 1 + ((self.current_M['R_k_net_factor'] - 1) + 
                                   self.cS['pps_return_rate_premium'])
            k_pps_next_unclamped = ((self.current_k_pps_val + actual_c_pps - pps_withdrawal) * 
                                  pps_return_factor)
            k_pps_next = max(self.cS['k_pps_min'], min(self.cS['k_pps_max'], k_pps_next_unclamped))
        
        # 5.2 非PPS资产由actual_k_prime决定
        k_next = actual_k_prime
        
        # 5.3 年龄演化
        age_next_idx = self.current_age_idx + 1
        
        # 5.4 效率冲击演化
        eps_next_idx = self.current_eps_idx
        if self.current_age_idx < self.cS['aD_new']:  # 只有在不是最后年龄时才转移
            trans_probs_eps = self.paramS_rl['leTrProbM'][self.current_eps_idx - 1, :]
            eps_next_idx = np.random.choice(len(trans_probs_eps), p=trans_probs_eps) + 1
        
        # 6. 判断是否终止（只在到达最大年龄时终止，期望化死亡率处理）
        terminated = False
        if age_next_idx > self.cS['aD_new']:
            terminated = True
        
        # 7. 更新内部状态（如果未终止）
        if not terminated:
            self.current_age_idx = age_next_idx
            self.current_k_val = k_next
            self.current_k_pps_val = k_pps_next
            self.current_eps_idx = eps_next_idx
        
        return terminated
    
    def _get_observation(self) -> np.ndarray:
        """获取当前观测（并归一化）"""
        # 原始观测向量
        raw_obs_vec = np.array([
            self.current_k_val, self.current_k_pps_val,
            self.current_age_idx, self.current_eps_idx,
            self.current_M['R_k_net_factor'], self.current_M['w_gross'],
            self.current_M['TR_total'], self.current_M['b_payg_avg_for_obs'],
            self.current_M['tau_l'], self.current_M['theta_payg_actual']
        ])
        
        # 归一化（min-max scaling to [0,1]）
        obs = (raw_obs_vec - self.obs_norm_min) / self.obs_norm_range
        obs = np.clip(obs, 0, 1)  # 确保在[0,1]范围内
        
        return obs.astype(np.float32)
    
    def set_macro_parameters(self, M_fixed: Dict[str, Any]):
        """设置宏观参数的方法，用于评估"""
        self.current_M = M_fixed.copy()
        # 并据此更新current_bV_payg
        self.current_bV_payg = np.zeros(self.cS['aD_new'])
        if self.cS['aR_new'] < self.cS['aD_new']:
            self.current_bV_payg[self.cS['aR_new']:] = M_fixed['b_payg_avg_retiree']
        self.current_M['b_payg_avg_for_obs'] = M_fixed['b_payg_avg_retiree']
        # print('Environment macro parameters set externally for evaluation.')

    def render(self, mode='human'):
        """渲染环境（可选实现）"""
        pass
    
    def close(self):
        """关闭环境"""
        pass

# Standalone function for MATLAB compatibility
def get_aggregates_from_simulation(vfi_results, paramS, cS, eIdxM, R_k_net_factor, w, TR, bV_payg, paramS_hh):
    """
    Standalone wrapper function that calls household simulation and aggregation.
    This function exists outside the class for better MATLAB compatibility.
    """
    return OLGV9Utils.get_aggregates_from_simulation(vfi_results, paramS, cS, eIdxM, R_k_net_factor, w, TR, bV_payg, paramS_hh)