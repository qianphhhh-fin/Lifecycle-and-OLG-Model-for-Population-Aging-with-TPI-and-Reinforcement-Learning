"""
OLG Environment V8 SAC - Python Implementation
转换自MATLAB版本的OLGEnv_v8_SAC.m

用于SAC Agent训练的OLG环境（期望化死亡率版本）
家庭在给定的宏观经济状态M下进行生命周期决策
"""

import numpy as np
import gymnasium as gym
from gymnasium import spaces
from typing import Dict, Tuple, Any, Optional
from olg_utils import OLGUtils
import warnings

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
        self.current_k_val = self.cS['kMin']
        self.current_k_pps_val = self.cS['kppsMin']
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
        self.current_k_val = self.cS['kMin']
        self.current_k_pps_val = self.cS['kppsMin']
        
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
            self.cS['physAgeMap'][self.current_age_idx][0] <= self.cS['pps_contribution_age_max_idx'] and 
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
            'pps_tax_deferral_active': self.cS['pps_active']
        }
        
        b_payg_this_age = self.current_bV_payg[self.current_age_idx - 1]
        current_epsilon_val = self.paramS_rl['leGridV'][self.current_eps_idx - 1]
        
        resources_after_pps, _, _ = OLGUtils.hh_income_huggett(
            self.current_k_val, self.current_M['R_k_net_factor'], self.current_M['w_gross'],
            self.current_M['TR_total'], b_payg_this_age, actual_c_pps,
            self.current_age_idx, paramS_hh_step, self.cS, current_epsilon_val
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
            actual_k_prime = max(self.cS['kMin'], 
                               min(actual_k_prime, resources_for_kprime_and_c_above_floor))
        else:
            # 如果资源不足以支付最低消费，则k_prime为kMin
            actual_k_prime = self.cS['kMin']
        
        actual_k_prime = max(self.cS['kMin'], min(actual_k_prime, self.cS['kMax']))  # 最终约束
        
        consumption_expenditure = resources_after_pps - actual_k_prime
        current_c = max(self.cS['cFloor'], consumption_expenditure / (1 + self.cS['tau_c']))
        
        return actual_k_prime, current_c
    
    def _calculate_reward(self, current_c: float) -> float:
        """计算奖励（当期效用）"""
        _, utility = OLGUtils.ces_utility(np.array([current_c]), self.cS['sigma'], self.cS)
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
            k_pps_next = max(self.cS['kppsMin'], min(self.cS['kppsMax'], k_pps_next_unclamped))
        
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