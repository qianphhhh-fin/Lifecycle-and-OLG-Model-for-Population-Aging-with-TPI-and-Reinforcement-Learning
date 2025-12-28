"""
OLG Model V8 Utilities - Python Implementation
转换自MATLAB版本的main_olg_v8_utils.m

实现内生PPS缴费决策的数值方法
包括参数设定、效用函数、收入计算等核心功能
"""

import numpy as np
import scipy.optimize as optimize
from typing import Dict, Tuple, List, Any
import warnings

class OLGUtils:
    """OLG模型V8工具函数类"""
    
    @staticmethod
    def parameter_values_huggett_style() -> Dict[str, Any]:
        """
        设置OLG模型V8的所有参数
        对应MATLAB的ParameterValues_HuggettStyle函数
        
        Returns:
            cS: 包含所有模型参数的字典
        """
        cS = {}
        
        # --- 人口结构基础参数 ---
        cS['age1_orig'] = 20              # 模型起始年龄（岁）
        cS['ageLast_orig'] = 98           # 模型终止年龄（岁）
        cS['ageRetire_orig'] = 65         # 退休年龄（岁）
        cS['popGrowth_orig'] = 0.012      # 原始人口增长率
        cS['aD_orig'] = cS['ageLast_orig'] - cS['age1_orig'] + 1        # 年度年龄组数
        cS['aR_idx_orig'] = cS['ageRetire_orig'] - cS['age1_orig'] + 1  # 退休年度年龄索引
        cS['aW_orig'] = cS['aR_idx_orig'] - 1                           # 工作年数
        cS['physAgeV_orig'] = np.arange(cS['age1_orig'], cS['ageLast_orig'] + 1)  # 年度年龄向量
        
        # --- 年度死亡率数据（基于中国生命表） ---
        cS['d_orig'] = np.array([
            0.00159,0.00169,0.00174,0.00172,0.00165,0.00156,0.00149,0.00145,0.00145,0.00149,
            0.00156,0.00163,0.00171,0.00181,0.00193,0.00207,0.00225,0.00246,0.00270,0.00299,
            0.00332,0.00368,0.00409,0.00454,0.00504,0.00558,0.00617,0.00686,0.00766,0.00865,
            0.00955,0.01058,0.01162,0.01264,0.01368,0.01475,0.01593,0.01730,0.01891,0.02074,
            0.02271,0.02476,0.02690,0.02912,0.03143,0.03389,0.03652,0.03930,0.04225,0.04538,
            0.04871,0.05230,0.05623,0.06060,0.06542,0.07066,0.07636,0.08271,0.08986,0.09788,
            0.10732,0.11799,0.12895,0.13920,0.14861,0.16039,0.17303,0.18665,0.20194,0.21877,
            0.23601,0.25289,0.26973,0.28612,0.30128,0.31416,0.32915,0.34450,0.36018
        ])
        
        if len(cS['d_orig']) != cS['aD_orig']:
            raise ValueError('年度死亡率数据d_orig长度与年龄跨度不匹配')
        
        cS['s_orig'] = 1 - cS['d_orig']  # 年度存活率
        
        # --- 年龄组聚合参数 ---
        cS['yearStep'] = 5                                    # 每个年龄组跨度（年）
        cS['aD_new'] = int(np.ceil(cS['aD_orig'] / cS['yearStep']))  # 年龄组数量
        cS['aR_new'] = int(np.ceil(cS['aW_orig'] / cS['yearStep']))  # 工作年龄组数量
        
        # 建立年度年龄到年龄组的映射关系
        cS['physAgeMap'] = {}
        for a in range(1, cS['aD_new'] + 1):
            start_idx = (a - 1) * cS['yearStep']
            end_idx = min(a * cS['yearStep'] - 1, cS['aD_orig'] - 1)
            cS['physAgeMap'][a] = list(range(start_idx, end_idx + 1))
        
        # 计算各年龄组代表性年龄
        cS['physAgeV_new'] = np.zeros(cS['aD_new'])
        for a in range(1, cS['aD_new'] + 1):
            cS['physAgeV_new'][a-1] = cS['physAgeV_orig'][cS['physAgeMap'][a][0]]
        
        # 计算年龄组间转移存活率
        cS['s_1yr_transitionV'] = np.zeros(cS['aD_new'])
        for a in range(1, cS['aD_new']):
            last_year_idx_in_group = cS['physAgeMap'][a][-1]
            if last_year_idx_in_group < cS['aD_orig'] - 1:
                cS['s_1yr_transitionV'][a-1] = cS['s_orig'][last_year_idx_in_group]
            else:
                cS['s_1yr_transitionV'][a-1] = 0
        cS['s_1yr_transitionV'][-1] = 0
        
        # --- 初始人口分布 ---
        cS['initial_pop'] = np.array([76.2, 86.4, 113.8, 98.6, 86.6, 102.7, 112.0, 99.0, 
                                     64.0, 66.9, 44.1, 25.4, 14.9, 6.8, 1.7, 0.2])
        if len(cS['initial_pop']) != cS['aD_new']:
            cS['initial_pop'] = np.ones(cS['aD_new']) * (100 / cS['aD_new'])
            warnings.warn('initial_pop长度与年龄组数不匹配，已重设为均匀分布。')
        
        # --- 年龄组间存活率 ---
        beta_surv_pop = np.array([0.998,0.996,0.994,0.992,0.988,0.984,0.980,0.976,
                                 0.970,0.960,0.945,0.920,0.880,0.800,0.680])
        if len(beta_surv_pop) != cS['aD_new'] - 1:
            raise ValueError(f'年龄组间存活率beta_surv_pop的长度对于{cS["aD_new"]}个年龄组不正确。应为{cS["aD_new"]-1}。')
        cS['survivalProbV_popdyn'] = np.append(beta_surv_pop, 0)
        
        # --- 人口动态收敛参数 ---
        cS['bgp_tolerance'] = 0.001
        cS['bgp_window'] = 5
        cS['max_periods'] = 50
        
        # --- 家庭偏好参数 ---
        cS['sigma'] = 1.5           # 相对风险厌恶系数
        cS['beta'] = 0.97           # 主观贴现因子（与MATLAB版本一致）
        cS['cFloor'] = 0.05         # 最低消费约束
        cS['nSim'] = 1000           # 蒙特卡洛模拟个体数
        
        # --- 生产技术参数 ---
        cS['A'] = 0.895944          # 全要素生产率
        cS['alpha'] = 0.36          # 资本产出弹性
        cS['ddk'] = 0.06            # 资本折旧率
        
        # --- 政府财政参数 ---
        cS['tau_k'] = 0.20          # 资本所得税率
        cS['tau_c'] = 0.10          # 消费税率
        cS['gov_exp_frac_Y'] = 0.15 # 政府支出占GDP比例
        cS['gov_debt_frac_Y'] = 0.60 # 政府债务占GDP比例
        
        # --- 劳动效率冲击过程参数 ---
        cS['leSigma1'] = np.sqrt(0.38)       # 初期效率分布标准差
        cS['leShockStd'] = np.sqrt(0.045)    # 效率冲击标准差
        cS['lePersistence'] = 0.96           # AR(1)持续性参数
        cS['leWidth'] = 4                    # Tauchen方法的标准差倍数
        cS['nw'] = 3                         # 效率状态网格点数
        
        # --- 资产网格参数 ---
        cS['tgKY'] = 3                      # 目标资本产出比
        cS['tgWage'] = (1-cS['alpha'])*cS['A']*((cS['tgKY']/cS['A'])**(cS['alpha']/(1-cS['alpha'])))
        cS['nk'] = 30                       # 非PPS资产网格点数
        cS['kMin'] = 0                      # 非PPS资产下界
        cS['kMax'] = 15 * cS['tgWage']      # 非PPS资产上界
        
        power = 1.5
        k_grid_v = cS['kMin'] + (cS['kMax'] - cS['kMin']) * (np.linspace(0, 1, cS['nk']) ** power)
        if cS['nk'] > 0:
            k_grid_v[0] = cS['kMin']
        cS['kGridV'] = k_grid_v
        
        # --- 年龄效率剖面 ---
        age_eff_v_orig_temp = np.zeros(100)
        
        # 对应MATLAB: ageEffV_orig_temp(20:72) = [各个年龄段]
        # Python索引需要减1，且要注意切片的右边界
        
        # 年龄20-36: linspace(0.3,1.5,36-20+1) = 17个点
        age_range_20_36 = np.linspace(0.3, 1.5, 36-20+1)  # 17个点
        # 年龄37-47: 1.5.*ones(1,47-37+1) = 11个点  
        age_range_37_47 = np.ones(47-37+1) * 1.5          # 11个点
        # 年龄48-65: linspace(1.5,0.2,65-48+1) = 18个点
        age_range_48_65 = np.linspace(1.5, 0.2, 65-48+1)  # 18个点
        # 年龄66-72: linspace(0.18,0,72-66+1) = 7个点
        age_range_66_72 = np.linspace(0.18, 0, 72-66+1)   # 7个点
        
        # Python索引：年龄20对应索引19，年龄36对应索引35
        age_eff_v_orig_temp[19:36] = age_range_20_36      # 索引19-35，17个点
        # 年龄37对应索引36，年龄47对应索引46  
        age_eff_v_orig_temp[36:47] = age_range_37_47      # 索引36-46，11个点
        # 年龄48对应索引47，年龄65对应索引64
        age_eff_v_orig_temp[47:65] = age_range_48_65      # 索引47-64，18个点
        # 年龄66对应索引65，年龄72对应索引71
        age_eff_v_orig_temp[65:72] = age_range_66_72      # 索引65-71，7个点
        
        cS['ageEffV_orig'] = age_eff_v_orig_temp[cS['age1_orig']-1:cS['ageLast_orig']]
        if len(cS['ageEffV_orig']) != cS['aD_orig']:
            raise ValueError('ageEffV_orig年度年龄效率剖面长度不匹配')
        
        # 计算年龄组平均效率
        cS['ageEffV_new'] = np.zeros(cS['aD_new'])
        for a in range(1, cS['aD_new'] + 1):
            age_indices = cS['physAgeMap'][a]
            cS['ageEffV_new'][a-1] = np.mean(cS['ageEffV_orig'][age_indices])
        
        # --- PPS制度参数 ---
        cS['use_continuous_optimization'] = True
        cS['pps_active'] = True
        cS['pps_tax_rate_withdrawal'] = 0.03
        cS['pps_return_rate_premium'] = 0.08
        cS['pps_withdrawal_rate'] = 0.15
        cS['pps_in_K'] = True
        cS['pps_bequeathable'] = True
        
        # --- PPS缴费约束参数 ---
        cS['pps_annual_contrib_limit'] = 9999
        cS['pps_max_contrib_frac'] = 1
        cS['pps_contribution_age_max_idx'] = cS['aR_idx_orig'] - 1
        cS['pps_withdrawal_age_min_idx'] = cS['aR_idx_orig']
        cS['n_pps_choice_grid_points'] = 12
        cS['power_pps_choice_grid'] = 1.3
        
        # --- PPS资产网格 ---
        cS['nkpps'] = 20
        cS['kppsMin'] = 0
        cS['kppsMax'] = cS['kMax']/2
        
        power_pps = 1.3
        kpps_grid_v = cS['kppsMin'] + (cS['kppsMax'] - cS['kppsMin']) * (np.linspace(0, 1, cS['nkpps']) ** power_pps)
        if cS['nkpps'] > 0:
            kpps_grid_v[0] = cS['kppsMin']
        cS['kppsGridV'] = kpps_grid_v
        
        # --- 数值优化参数 ---
        cS['fminbnd_TolX'] = 1e-4
        cS['fminbnd_Display'] = 'off'
        cS['maxfunevals'] = 300
        cS['maxiters'] = 120
        
        return cS
    
    @staticmethod
    def ces_utility(c_quantity, sigma_crra: float, cS: Dict[str, Any]) -> Tuple[np.ndarray, np.ndarray]:
        """
        计算CRRA效用函数和边际效用
        对应MATLAB的CES_utility函数
        
        Args:
            c_quantity: 消费量（可为数组或标量）
            sigma_crra: 相对风险厌恶系数γ
            cS: 参数结构体
            
        Returns:
            mu_m: 边际效用 ∂U/∂c
            util_m: 效用水平 U(c)
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
        
        util_m = np.full_like(c_quantity, -np.inf, dtype=float)
        mu_m = np.full_like(c_quantity, np.inf, dtype=float)
        
        if abs(sigma_crra - 1) < 1e-6:  # 对数效用 (sigma = 1)
            util_m[is_valid_consumption] = np.log(c_adjusted_quantity[is_valid_consumption])
            mu_m[is_valid_consumption] = 1.0 / c_adjusted_quantity[is_valid_consumption]
        else:  # CRRA效用 (sigma != 1)
            util_m[is_valid_consumption] = (c_adjusted_quantity[is_valid_consumption] ** (1-sigma_crra)) / (1-sigma_crra)
            mu_m[is_valid_consumption] = c_adjusted_quantity[is_valid_consumption] ** (-sigma_crra)
        
        # 惩罚低于最低消费的情况
        util_m[~is_valid_consumption] = -1e10 - (min_c_quantity - c_quantity[~is_valid_consumption]) * 1e10
        
        # 如果输入是标量，返回标量
        if is_scalar_input:
            return float(mu_m[0]), float(util_m[0])
        else:
            return mu_m, util_m
    
    @staticmethod
    def earning_process_olgm(cS: Dict[str, Any]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        使用Tauchen方法离散化AR(1)劳动效率过程
        对应MATLAB的EarningProcess_olgm函数
        
        Args:
            cS: 参数结构体
            
        Returns:
            le_log_grid_v: 对数效率网格
            le_tr_prob_m: 转移概率矩阵
            le_prob1_v: 初始分布
        """
        nw = cS['nw']
        persistence = cS['lePersistence']
        shock_std = cS['leShockStd']
        width = cS['leWidth']
        sigma1 = cS['leSigma1']
        
        # Tauchen方法实现
        std_long_run = shock_std / np.sqrt(1 - persistence**2)
        
        # 构建网格点
        y_max = width * std_long_run
        y_min = -y_max
        le_log_grid_v = np.linspace(y_min, y_max, nw)
        
        # 网格步长
        h = le_log_grid_v[1] - le_log_grid_v[0]
        
        # 构建转移概率矩阵
        le_tr_prob_m = np.zeros((nw, nw))
        
        from scipy.stats import norm
        
        for i in range(nw):
            for j in range(nw):
                if j == 0:
                    le_tr_prob_m[i, j] = norm.cdf((le_log_grid_v[0] - persistence * le_log_grid_v[i] + h/2) / shock_std)
                elif j == nw - 1:
                    le_tr_prob_m[i, j] = 1 - norm.cdf((le_log_grid_v[-1] - persistence * le_log_grid_v[i] - h/2) / shock_std)
                else:
                    le_tr_prob_m[i, j] = (norm.cdf((le_log_grid_v[j] - persistence * le_log_grid_v[i] + h/2) / shock_std) - 
                                         norm.cdf((le_log_grid_v[j] - persistence * le_log_grid_v[i] - h/2) / shock_std))
        
        # 规范化转移概率
        for i in range(nw):
            row_sum = np.sum(le_tr_prob_m[i, :])
            if row_sum > 0:
                le_tr_prob_m[i, :] /= row_sum
        
        # 计算初始分布（不变分布）
        le_prob1_v = np.ones(nw) / nw  # 简化版本：均匀分布
        
        # 更精确的不变分布计算
        try:
            eigenvals, eigenvecs = np.linalg.eig(le_tr_prob_m.T)
            stationary_idx = np.argmax(np.real(eigenvals))
            le_prob1_v = np.real(eigenvecs[:, stationary_idx])
            le_prob1_v = le_prob1_v / np.sum(le_prob1_v)
            le_prob1_v = np.abs(le_prob1_v)  # 确保非负
        except:
            # 如果特征值计算失败，使用均匀分布
            le_prob1_v = np.ones(nw) / nw
        
        return le_log_grid_v, le_tr_prob_m, le_prob1_v
    
    @staticmethod
    def hh_income_huggett(k_now_val: float, R_k_net_factor: float, w_gross: float,
                         TR_total: float, b_payg_val: float, c_pps_input_val: float,
                         a_idx: int, paramS_hh: Dict[str, Any], cS: Dict[str, Any], 
                         epsilon_val: float) -> Tuple[float, float, float]:
        """
        计算家庭可支配资源
        对应MATLAB的HHIncome_Huggett函数
        
        Args:
            k_now_val: 当前非PPS资产
            R_k_net_factor: 税后资本回报因子
            w_gross: 市场毛工资率
            TR_total: 总转移支付
            b_payg_val: PAYG养老金
            c_pps_input_val: PPS缴费额
            a_idx: 年龄组索引
            paramS_hh: 家庭参数
            cS: 模型参数
            epsilon_val: 当前劳动效率
            
        Returns:
            resources_for_c_and_k_prime: 可用于消费和储蓄的资源
            labor_income_gross_state: 税前劳动收入
            pps_deduction_actual_state: PPS税前扣除额
        """
        labor_income_gross_state = 0.0
        pps_deduction_actual_state = 0.0
        non_capital_income = 0.0
        
        # 确保PPS缴费非负
        actual_pps_contribution_expenditure = max(0, c_pps_input_val)
        
        if a_idx <= cS['aR_new']:  # 工作年龄组
            age_efficiency = cS['ageEffV_new'][a_idx - 1]  # Python索引从0开始
            labor_income_gross_state = w_gross * age_efficiency * epsilon_val
            
            # PPS税收递延处理
            if paramS_hh.get('pps_tax_deferral_active', False):
                pps_deduction_actual_state = actual_pps_contribution_expenditure
            else:
                pps_deduction_actual_state = 0
            
            # 应税劳动收入
            labor_income_taxable_for_tau_l = labor_income_gross_state - pps_deduction_actual_state
            labor_income_taxable_for_tau_l = max(0, labor_income_taxable_for_tau_l)
            
            # 税收计算
            income_tax_tau_l = labor_income_taxable_for_tau_l * paramS_hh['tau_l']
            payg_tax_theta = labor_income_gross_state * paramS_hh['theta_payg_actual_for_hh']
            
            # 税后劳动收入
            labor_income_net_of_all_taxes = labor_income_gross_state - income_tax_tau_l - payg_tax_theta
            non_capital_income = labor_income_net_of_all_taxes + TR_total + b_payg_val
        else:  # 退休年龄组
            actual_pps_contribution_expenditure = 0
            pps_deduction_actual_state = 0
            non_capital_income = TR_total + b_payg_val
        
        # 资本收入
        capital_income_net_of_tax = (R_k_net_factor - 1) * k_now_val
        
        # 总可支配资源
        resources_for_c_and_k_prime = (k_now_val + capital_income_net_of_tax + 
                                     non_capital_income - actual_pps_contribution_expenditure)
        
        if not np.isfinite(resources_for_c_and_k_prime):
            resources_for_c_and_k_prime = -1e10
        
        return resources_for_c_and_k_prime, labor_income_gross_state, pps_deduction_actual_state 