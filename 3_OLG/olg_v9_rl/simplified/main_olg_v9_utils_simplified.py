"""
OLG Model V9 Utils - Simplified Version with Fixed Macro Parameters
简化版OLG模型工具类 - 使用固定宏观参数

🎯 核心简化：
- 宏观参数在环境初始化时固定，不作为状态变量
- RL状态变量与VFI基本相同：(k, k_pps, age, ε)
- 保持累积存活概率方法的理论等价性
- 更公平地比较RL和VFI性能

主要变化：
1. 环境初始化时传入固定宏观参数
2. 观测空间降维：从10维降到4维
3. 动作空间改为：[PPS缴费比例, 消费比例]
4. 决策顺序：先PPS缴费，再消费，最后储蓄
5. 保持所有其他逻辑不变
"""

import numpy as np
import gymnasium as gym
from gymnasium import spaces
import pickle
import time
from typing import Dict, Any, Tuple, Optional, List
from scipy.optimize import minimize
from scipy.interpolate import RegularGridInterpolator

# 从原始文件导入必要的类和函数
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from main_olg_v9_utils import OLG_V9_Utils, TempParamSHH

class OLGEnvV9SACSimplified(gym.Env):
    """
    简化版OLG环境 - 固定宏观参数版本
    
    🎯 核心特性：
    - 宏观参数在初始化时固定
    - 状态空间：(k, k_pps, age, ε) - 与VFI基本相同
    - 动作空间：[PPS缴费比例, 消费比例] - 先PPS后消费再储蓄
    - 保持累积存活概率方法
    - 训练/评估模式支持
    """
    metadata = {'render.modes': ['human']}  # 定义环境支持的渲染模式

    def __init__(self, cS: Any, paramS_rl: Dict[str, Any], M_fixed: Dict[str, float], training_mode: bool = True):
        """
        初始化简化版OLG环境
        
        Args:
            cS: 模型参数对象
            paramS_rl: RL专用参数字典
            M_fixed: 固定的宏观参数字典
            training_mode: 训练模式（True）或评估模式（False）
        """
        super().__init__()  # 调用父类gym.Env的初始化方法
        
        self.cS = cS  # 保存模型参数对象，包含所有经济参数
        self.paramS_rl = paramS_rl  # 保存RL专用参数字典，包含效率冲击转移概率等
        self.M_fixed = M_fixed.copy()  # 深拷贝固定宏观参数字典，避免意外修改
        self.training_mode = False  # 设置训练模式标志，影响奖励函数计算方式
        
        # 打印环境初始化信息，显示当前模式
        print(f"🏋️ OLG环境初始化完成 - {'训练' if training_mode else '评估'}模式 (简化版)")
        print(f"   状态空间: (k, k_pps, age, ε) - 4维")  # 状态空间降维到4维
        print(f"   固定宏观参数: R={M_fixed['R_k_net_factor']:.3f}, w={M_fixed['w_gross']:.3f}")  # 显示关键宏观参数
        if training_mode:
            # 训练模式使用累积存活概率加权的奖励函数
            print(f"   奖励函数: r'(t) = u(c) * ∏_{{i=1}}^{{t}} s(i) (累积存活概率加权)")
        else:
            # 评估模式使用纯效用奖励函数
            print(f"   奖励函数: r(t) = u(c) (纯效用奖励)")
        
        # 设置动作空间和观测空间
        self._setup_spaces()
        
        # 初始化归一化参数
        self._init_normalization_params()
        
        # 初始化状态变量：个体资产
        self.current_k_val = self.cS.kMin  # 当前非PPS资产，从最小值开始
        self.current_k_pps_val = self.cS.kppsMin  # 当前PPS资产，从最小值开始
        self.current_age_idx = 1  # 当前年龄组索引（1-based），从第1个年龄组开始
        self.current_eps_idx = 1  # 当前效率冲击索引（1-based），从第1个冲击开始
        
        # 累积存活概率（用于训练模式的奖励函数）
        self.cumulative_survival_prob = 1.0  # 初始化为1.0，随时间衰减
        
        # 构建福利向量
        self._build_payg_benefits()
        
        # 随机数生成器
        self.np_random = np.random.RandomState()  # 创建独立的随机数生成器

    def _setup_spaces(self):
        """设置动作空间和观测空间"""
        # 动作空间：[PPS缴费比例, 消费比例]
        self.action_space = spaces.Box(
            low=np.array([0.0, 0.0]),  # 动作下界：两个比例都不能为负
            high=np.array([self.cS.pps_max_contrib_frac, 1.0]),  # 动作上界：PPS缴费受限制，消费比例最大为1
            dtype=np.float32  # 使用32位浮点数以提高计算效率
        )
        
        # 观测空间：[k, k_pps, age, ε] - 4维，与VFI基本相同
        self.observation_space = spaces.Box(
            low=np.array([0.0, 0.0, 0.0, 0.0]),  # 观测下界：归一化后都是0
            high=np.array([1.0, 1.0, 1.0, 1.0]),  # 观测上界：归一化后都是1
            dtype=np.float32  # 使用32位浮点数以提高计算效率
        )

    def _init_normalization_params(self):
        """初始化观测归一化参数"""
        # 原始观测范围 [k, k_pps, age, ε]
        obs_min = np.array([
            self.cS.kMin,  # 非PPS资产最小值
            self.cS.kppsMin,  # PPS资产最小值
            1,  # 年龄组索引最小值（1-based）
            1   # 效率冲击索引最小值（1-based）
        ])
        
        obs_max = np.array([
            self.cS.kMax,  # 非PPS资产最大值
            self.cS.kppsMax,  # PPS资产最大值
            self.cS.aD_new,  # 年龄组数量（最大年龄组索引）
            self.cS.nw  # 效率冲击状态数量（最大冲击索引）
        ])
        
        self.obs_norm_min = obs_min  # 保存归一化最小值
        self.obs_norm_range = obs_max - obs_min  # 计算归一化范围
        
        # 打印归一化范围信息
        print(f"📊 观测归一化范围:")
        print(f"   k: [{self.cS.kMin:.2f}, {self.cS.kMax:.2f}]")  # 非PPS资产范围
        print(f"   k_pps: [{self.cS.kppsMin:.2f}, {self.cS.kppsMax:.2f}]")  # PPS资产范围
        print(f"   age: [1, {self.cS.aD_new}]")  # 年龄组范围
        print(f"   ε: [1, {self.cS.nw}]")  # 效率冲击范围

    def _build_payg_benefits(self):
        """构建PAYG福利向量"""
        self.current_bV_payg = np.zeros(self.cS.aD_new)  # 初始化所有年龄组的PAYG福利为0
        if self.cS.aR_new < self.cS.aD_new:  # 如果退休年龄小于最大年龄
            # 为退休年龄及以后的所有年龄组设置固定的PAYG福利
            self.current_bV_payg[self.cS.aR_new:] = self.M_fixed['b_payg_avg_retiree']

    def set_training_mode(self, training_mode: bool):
        """设置训练/评估模式"""
        old_mode = self.training_mode  # 记录旧模式
        self.training_mode = training_mode  # 设置新模式
        if old_mode != training_mode:  # 如果模式确实发生了变化
            mode_name = "训练" if training_mode else "评估"  # 确定模式名称
            print(f"🔄 环境模式切换为: {mode_name}模式")  # 打印模式切换信息

    def reset(self, seed: Optional[int] = None, options: Optional[Dict] = None) -> Tuple[np.ndarray, Dict]:
        """重置环境"""
        if seed is not None:  # 如果提供了随机种子
            self.np_random = np.random.RandomState(seed)  # 用种子初始化随机数生成器
        
        # 重置状态变量到初始值
        self.current_k_val = self.cS.kMin  # 重置非PPS资产到最小值
        self.current_k_pps_val = self.cS.kppsMin  # 重置PPS资产到最小值
        self.current_age_idx = 1  # 重置年龄组索引到第1组
        
        # 重置累积存活概率
        self.cumulative_survival_prob = 1.0  # 重置累积存活概率为1.0
        
        # 根据初始分布随机选择初始效率冲击
        leProb1V = np.array(self.paramS_rl['leProb1V']).flatten()  # 获取初始效率冲击分布
        self.current_eps_idx = self.np_random.choice(len(leProb1V), p=leProb1V) + 1  # 随机选择并转换为1-based索引
        
        observation = self._get_observation()  # 获取初始观测
        # 构建信息字典
        info = {
            'k': self.current_k_val,  # 当前非PPS资产
            'k_pps': self.current_k_pps_val,  # 当前PPS资产
            'age_idx': self.current_age_idx,  # 当前年龄组索引
            'eps_idx': self.current_eps_idx,  # 当前效率冲击索引
            'cumulative_survival_prob': self.cumulative_survival_prob,  # 累积存活概率
            'training_mode': self.training_mode  # 当前模式
        }
        
        return observation, info  # 返回观测和信息

    def step(self, action: np.ndarray) -> Tuple[np.ndarray, float, bool, bool, Dict]:
        """执行一步动作"""
        # 解析动作向量
        prop_pps_contrib = float(action[0])  # 提取PPS缴费比例
        prop_consumption = float(action[1])  # 提取消费比例
        
        # 1. 计算PPS缴费
        actual_c_pps, max_permissible_cpps = self._calculate_pps_contribution(prop_pps_contrib)
        
        # 2. 计算扣除PPS缴费后的可用资源
        resources_after_pps = self._calculate_resources_after_pps(actual_c_pps)
        
        # 3. 计算消费和储蓄（改为基于消费比例的决策）
        actual_k_prime, current_c = self._calculate_consumption_and_savings(
            resources_after_pps, prop_consumption
        )
        
        # 4. 计算纯效用
        pure_utility = self._calculate_pure_utility(current_c)
        
        # 5. 获取存活概率
        survival_prob = 1.0  # 默认存活概率为1
        vfi_age_idx = self.current_age_idx - 1  # 转换为VFI的0-based索引
        
        if vfi_age_idx >= len(self.cS.s_1yr_transitionV):  # 如果超出存活概率向量长度
            survival_prob = 0.0  # 存活概率为0（死亡）
        elif vfi_age_idx >= 0:  # 如果索引有效
            survival_prob = self.cS.s_1yr_transitionV[vfi_age_idx]  # 获取对应的存活概率
        else:  # 如果索引无效（小于0）
            survival_prob = 1.0  # 设为1（不应该发生）

        # 6. 根据模式计算奖励
        if self.training_mode:
            # 训练模式：使用累积存活概率加权的奖励
            reward = pure_utility * self.cumulative_survival_prob
        else:
            # 评估模式：使用纯效用奖励
            reward = pure_utility

        # 7. 更新状态
        terminated = self._update_state(actual_k_prime, actual_c_pps)
        
        # 8. 更新累积存活概率（为下一步准备）
        if not terminated:  # 如果游戏没有结束
            self.cumulative_survival_prob *= survival_prob  # 累积存活概率乘以当期存活概率

        observation = self._get_observation()  # 获取新的观测

        # 9. 构建info字典
        info = {
            "survival_prob": survival_prob,  # 当期存活概率
            "beta": self.cS.beta,  # 贴现因子
            'vfi_equivalent_info': {  # VFI等价信息
                'survival_prob': survival_prob,  # 存活概率
                'beta': self.cS.beta,  # 贴现因子
                'vfi_age_idx': vfi_age_idx,  # VFI年龄索引
                'discount_factor': self.cS.beta * survival_prob,  # 有效贴现因子
                'cumulative_survival_prob': self.cumulative_survival_prob,  # 累积存活概率
                'pure_utility': pure_utility,  # 纯效用
                'training_mode': self.training_mode  # 训练模式
            },
            'consumption': current_c,  # 当期消费
            'k_prime': actual_k_prime,  # 下期非PPS资产
            'c_pps': actual_c_pps,  # 当期PPS缴费
            'age_idx': self.current_age_idx,  # 当前年龄组索引
            'vfi_age_idx': vfi_age_idx,  # VFI年龄索引
            'cumulative_survival_prob': self.cumulative_survival_prob,  # 累积存活概率
            'reward_type': 'cumulative_weighted' if self.training_mode else 'pure_utility',  # 奖励类型
            # 固定宏观参数信息
            'M_fixed': self.M_fixed  # 固定宏观参数字典
        }

        truncated = False  # 游戏没有被截断
        return observation, reward, terminated, truncated, info  # 返回标准的step返回值

    def _calculate_pure_utility(self, current_c: float) -> float:
        """计算纯效用 u(c)"""
        # 使用CES效用函数计算消费的效用值
        _, utility = OLG_V9_Utils.CES_utility(np.array([current_c]), self.cS.sigma, self.cS)

        if not np.isfinite(utility):  # 如果效用值不是有限数
            # 返回负的惩罚值，惩罚程度与消费偏离最低消费的程度相关
            return -1000.0

        return float(utility)  # 返回有效的效用值

    def _calculate_pps_contribution(self, prop_pps_contrib: float) -> Tuple[float, float]:
        """计算实际PPS缴费"""
        actual_c_pps = 0.0  # 初始化实际PPS缴费为0
        max_permissible_cpps = 0.0  # 初始化最大允许PPS缴费为0
        current_epsilon_val = self.paramS_rl['leGridV'][self.current_eps_idx - 1]  # 获取当前效率冲击值
        
        # 获取当前的年度年龄索引
        current_annual_age_idx = self.cS.physAgeMap[self.current_age_idx - 1][0]
        # 判断是否为工作年龄
        is_working_age_annual = (current_annual_age_idx + 1 < self.cS.aR_idx_orig)
        # 判断是否符合PPS缴费条件
        is_pps_contribution_eligible = (is_working_age_annual and
                                       current_annual_age_idx <= self.cS.pps_contribution_age_max_idx and
                                       self.cS.pps_active)
        
        if is_pps_contribution_eligible:  # 如果符合PPS缴费条件
            # 获取当前年龄的效率系数
            age_efficiency = self.cS.ageEffV_new[self.current_age_idx - 1]
            # 计算当前的总劳动收入
            current_gross_labor_income = self.M_fixed['w_gross'] * age_efficiency * current_epsilon_val
            
            if current_gross_labor_income > 1e-6:  # 如果劳动收入大于很小的阈值
                # 按收入比例计算的最大PPS缴费
                max_cpps_by_frac = current_gross_labor_income * self.cS.pps_max_contrib_frac
                # 取收入比例限制和绝对限额的较小值
                max_permissible_cpps = min(self.cS.pps_annual_contrib_limit, max_cpps_by_frac)
                # 确保不为负数
                max_permissible_cpps = max(0, max_permissible_cpps)
            
            # 根据动作比例计算实际PPS缴费
            actual_c_pps = prop_pps_contrib * max_permissible_cpps
            # 确保实际缴费在有效范围内
            actual_c_pps = max(0, min(actual_c_pps, max_permissible_cpps))
        
        return actual_c_pps, max_permissible_cpps  # 返回实际缴费和最大允许缴费

    def _calculate_resources_after_pps(self, actual_c_pps: float) -> float:
        """计算扣除PPS缴费后的可用资源"""
        # 创建临时参数对象
        paramS_hh_step = TempParamSHH(
            self.M_fixed['tau_l'],  # 劳动所得税率
            self.M_fixed['theta_payg_actual'],  # 实际PAYG缴费率
            self.cS.pps_active,  # PPS是否激活
            self.cS.ageEffV_new  # 年龄效率向量
        )
        
        # 获取当前年度年龄索引
        current_annual_age_idx = self.cS.physAgeMap[self.current_age_idx - 1][0]
        # 判断是否为工作年龄
        is_working_age = (current_annual_age_idx < self.cS.aR_idx_orig)
        # 判断是否符合PPS提取条件
        is_pps_withdrawal_eligible = (not is_working_age and
                                    self.cS.pps_active and
                                    current_annual_age_idx >= self.cS.pps_withdrawal_age_min_idx)
        
        if is_pps_withdrawal_eligible:  # 如果符合PPS提取条件
            # 计算PPS提取金额
            pps_withdrawal_amount = self.current_k_pps_val * self.cS.pps_withdrawal_rate
            paramS_hh_step.current_pps_withdrawal = pps_withdrawal_amount  # 设置当期PPS提取
        else:
            paramS_hh_step.current_pps_withdrawal = 0  # 否则PPS提取为0
        
        # 获取当前年龄的PAYG福利
        b_payg_this_age = self.current_bV_payg[self.current_age_idx - 1]
        # 获取当前效率冲击值
        current_epsilon_val = self.paramS_rl['leGridV'][self.current_eps_idx - 1]
        
        # 调用收入计算函数
        resources_after_pps, _, _ = OLG_V9_Utils.HHIncome_Huggett(
            self.current_k_val,  # 当前非PPS资产
            self.M_fixed['R_k_net_factor'],  # 净资本回报率
            self.M_fixed['w_gross'],  # 总工资率
            self.M_fixed['TR_total'],  # 总转移支付
            b_payg_this_age,  # 当前年龄的PAYG福利
            actual_c_pps,  # 实际PPS缴费
            self.current_age_idx - 1,  # 年龄组索引（0-based）
            paramS_hh_step,  # 临时参数对象
            self.cS,  # 模型参数
            current_epsilon_val  # 当前效率冲击值
        )
        
        return resources_after_pps  # 返回扣除PPS缴费后的可用资源

    def _calculate_consumption_and_savings(self, resources_after_pps: float,
                                            prop_consumption: float) -> Tuple[float, float]:
        """计算消费和储蓄（先消费后储蓄）"""
        # 计算消费底线支出（包含消费税）
        consumption_floor_spending = self.cS.cFloor * (1 + self.cS.tau_c)
        # 计算可用于消费和储蓄的总资源
        total_available_resources = max(0, resources_after_pps)
        
        # 1. 根据消费比例确定消费支出
        if total_available_resources >= consumption_floor_spending:
            # 资源充足时，可在底线消费之上进行选择
            resources_above_floor = total_available_resources - consumption_floor_spending
            # 根据消费比例分配：底线消费 + 比例消费
            consumption_expenditure = consumption_floor_spending + prop_consumption * resources_above_floor
        else:
            # 资源不足时，只能维持底线消费
            consumption_expenditure = total_available_resources
        
        # 确保消费支出不超过总资源
        consumption_expenditure = min(consumption_expenditure, total_available_resources)
        consumption_expenditure = max(consumption_expenditure, consumption_floor_spending)
        
        # 2. 计算实际消费（扣除消费税）
        current_c = max(self.cS.cFloor, consumption_expenditure / (1 + self.cS.tau_c))
        
        # 3. 剩余资源全部用于储蓄
        actual_k_prime = total_available_resources - consumption_expenditure
        
        # 确保储蓄在允许范围内
        actual_k_prime = max(self.cS.kMin, min(actual_k_prime, self.cS.kMax))
        
        return actual_k_prime, current_c  # 返回下期资产和当期消费

    def _update_state(self, actual_k_prime: float, actual_c_pps: float) -> bool:
        """更新个体状态到下一期"""
        self.current_k_val = actual_k_prime  # 更新非PPS资产
        
        if self.cS.pps_active:  # 如果PPS系统激活
            # 获取当前年度年龄索引
            current_annual_age_idx = self.cS.physAgeMap[self.current_age_idx - 1][0]
            # 判断是否为工作年龄
            is_working_age_annual = (current_annual_age_idx + 1 < self.cS.aR_idx_orig)
            # 判断是否符合PPS提取条件
            is_pps_withdrawal_eligible = (not is_working_age_annual and
                                        self.cS.pps_active and
                                        current_annual_age_idx >= self.cS.pps_withdrawal_age_min_idx)
            
            pps_withdrawal = 0  # 初始化PPS提取为0
            if is_pps_withdrawal_eligible:  # 如果符合提取条件
                pps_withdrawal = self.current_k_pps_val * self.cS.pps_withdrawal_rate  # 计算提取金额
            
            # 计算PPS收益率
            pps_return_factor = 1 + ((self.M_fixed['R_k_net_factor'] - 1) + self.cS.pps_return_rate_premium)
            # 计算下期PPS资产（未约束）
            k_pps_next_unclamped = (self.current_k_pps_val + actual_c_pps - pps_withdrawal) * pps_return_factor
            
            # 约束PPS资产在允许范围内
            self.current_k_pps_val = max(self.cS.kppsMin, min(self.cS.kppsMax, k_pps_next_unclamped))
        else:  # 如果PPS系统未激活
            self.current_k_pps_val = self.cS.kppsMin  # PPS资产设为最小值
        
        terminated = False  # 初始化游戏结束标志为False
        if self.current_age_idx < self.cS.aD_new:  # 如果还未到最后一个年龄组
            # 获取效率冲击转移概率
            trans_probs = self.paramS_rl['leTrProbM'][self.current_eps_idx - 1, :]
            # 根据转移概率随机选择下期效率冲击
            self.current_eps_idx = self.np_random.choice(len(trans_probs), p=trans_probs) + 1
            self.current_age_idx += 1  # 年龄组索引加1
        else:  # 如果已到最后一个年龄组
            terminated = True  # 设置游戏结束标志为True
        
        return terminated  # 返回游戏结束标志

    def _get_observation(self) -> np.ndarray:
        """获取当前观测（并归一化）"""
        # 简化的观测：[k, k_pps, age, ε]
        raw_obs_vec = np.array([
            self.current_k_val,  # 当前非PPS资产
            self.current_k_pps_val,  # 当前PPS资产
            self.current_age_idx,  # 当前年龄组索引
            self.current_eps_idx  # 当前效率冲击索引
        ])
        return self._normalize_observation(raw_obs_vec)  # 归一化并返回观测

    def _normalize_observation(self, raw_obs_vec: np.ndarray) -> np.ndarray:
        """归一化观测值"""
        obs = (raw_obs_vec - self.obs_norm_min) / self.obs_norm_range  # 线性归一化到[0,1]
        obs = np.clip(obs, 0, 1)  # 确保值在[0,1]范围内
        return obs.astype(np.float32)  # 转换为32位浮点数并返回

    def render(self, mode='human'):
        """渲染环境（可选实现）"""
        pass  # 暂不实现渲染功能

    def close(self):
        """关闭环境"""
        pass  # 暂无需清理资源


class OLGUtilsSimplified:
    """简化版OLG工具类 - 包装原始工具类的必要功能"""
    
    @staticmethod
    def HHSimulation_olgm_rl_simplified(rl_model, cS, paramS_rl, M_fixed, eIdxM_input):
        """
        简化版RL生命周期模拟（年龄组版本）
        
        重要更新：改为按年龄组模拟，而不是年度年龄模拟
        - 模拟维度：从aD_orig（年度）改为aD_new（年龄组）
        - 效率冲击：直接使用年龄组效率冲击序列
        - 年龄判断：直接使用年龄组索引，无需复杂映射
        - 动作空间：[PPS缴费比例, 消费比例]
        
        Args:
            rl_model: RL模型
            cS: 参数对象
            paramS_rl: RL参数字典
            M_fixed: 固定宏观参数
            eIdxM_input: 年龄组效率冲击序列（nSim × aD_new）
            
        Returns:
            k_path, kpps_path, c_path, cpps_path: 生命周期路径（年龄组尺度）
        """
        # 自动检测输入格式
        
            # 已经是年龄组格式，直接使用
        eIdxM_group = eIdxM_input.astype(int)
        print(f"🔍 检测到年龄组效率冲击序列: {eIdxM_input.shape} (nSim × aD_new)")

        
        n_sim = eIdxM_group.shape[0]
        aD_new = int(cS.aD_new)
        
        # 创建临时环境用于归一化（评估模式）
        temp_env = OLGEnvV9SACSimplified(cS, paramS_rl, M_fixed, training_mode=False)
        
        # 初始化结果矩阵（年龄组尺度）
        kHistM_rl = np.zeros((n_sim, aD_new))
        kPpsHistM_rl = np.zeros((n_sim, aD_new))
        cHistM_rl = np.zeros((n_sim, aD_new))
        cppsHistM_rl = np.zeros((n_sim, aD_new))
        
        # 生成年龄组PAYG福利向量
        bV_payg_group = np.zeros(aD_new)
        if cS.aR_new < aD_new:
            bV_payg_group[cS.aR_new:] = M_fixed['b_payg_avg_retiree']
        
        # PPS相关参数
        pps_return_factor = 1 + ((M_fixed['R_k_net_factor'] - 1) + cS.pps_return_rate_premium)
        
        # 创建参数对象
        paramS_hh = TempParamSHH(
            M_fixed['tau_l'],
            M_fixed['theta_payg_actual'],
            cS.pps_active,
            cS.ageEffV_new
        )
        
        # 生命周期模拟主循环（年龄组版本）
        for i_sim in range(n_sim):
            # 初始状态
            k_current = cS.kMin
            kpps_current = cS.kppsMin
            
            # 年龄组时间循环
            for a_group in range(aD_new):
                # 获取效率冲击
                eps_idx_current = eIdxM_group[i_sim, a_group]
                epsilon_val = paramS_rl['leGridV'][eps_idx_current]
                
                # 年龄组判断（简化）
                is_working_age_group = (a_group < cS.aR_new)
                is_pps_withdrawal_eligible = (not is_working_age_group and
                                            cS.pps_active)
                
                # PPS提取
                pps_withdrawal_pretax = 0
                if is_pps_withdrawal_eligible:
                    pps_withdrawal_pretax = kpps_current * cS.pps_withdrawal_rate
                    paramS_hh.current_pps_withdrawal = pps_withdrawal_pretax
                else:
                    paramS_hh.current_pps_withdrawal = 0
                
                if a_group == aD_new - 1:
                    # 最后一期：消费所有资产
                    k_next = cS.kMin
                    cpps_decision = 0
                    
                    if cS.pps_active and kpps_current > cS.kppsMin:
                        pps_withdrawal_pretax = kpps_current
                        paramS_hh.current_pps_withdrawal = pps_withdrawal_pretax
                    
                    b_payg_val = bV_payg_group[a_group]
                    total_income, _, _ = OLG_V9_Utils.HHIncome_Huggett(
                        k_current, M_fixed['R_k_net_factor'], M_fixed['w_gross'],
                        M_fixed['TR_total'], b_payg_val, 0.0,
                        a_group, paramS_hh, cS, epsilon_val
                    )
                    
                    consumption_resources = total_income - k_next
                    c_consumption = max(cS.cFloor, consumption_resources / (1 + cS.tau_c))
                    
                else:
                    # 使用RL策略
                    # 构造观测状态（简化版：4维）
                    raw_obs = np.array([
                        k_current, kpps_current,
                        a_group + 1, eps_idx_current + 1
                    ])
                    
                    # 归一化观测
                    normalized_obs = temp_env._normalize_observation(raw_obs)
                    
                    # RL策略预测
                    action, _ = rl_model.predict(normalized_obs, deterministic=True)
                    prop_pps_contrib = float(action[0])
                    prop_consumption = float(action[1])  # 改为消费比例
                    
                    # PPS缴费计算
                    cpps_decision = 0
                    max_permissible_cpps = 0
                    if (is_working_age_group and cS.pps_active):
                        age_efficiency = cS.ageEffV_new[a_group]
                        gross_labor_income = M_fixed['w_gross'] * age_efficiency * epsilon_val
                        
                        if gross_labor_income > 1e-6:
                            max_cpps_by_frac = gross_labor_income * cS.pps_max_contrib_frac
                            max_permissible_cpps = min(cS.pps_annual_contrib_limit, max_cpps_by_frac)
                            max_permissible_cpps = max(0, max_permissible_cpps)
                            
                            cpps_decision = prop_pps_contrib * max_permissible_cpps
                            cpps_decision = max(0, min(cpps_decision, max_permissible_cpps))
                    
                    # 收入计算
                    b_payg_val = bV_payg_group[a_group]
                    total_income, _, _ = OLG_V9_Utils.HHIncome_Huggett(
                        k_current, M_fixed['R_k_net_factor'], M_fixed['w_gross'],
                        M_fixed['TR_total'], b_payg_val, cpps_decision,
                        a_group, paramS_hh, cS, epsilon_val
                    )
                    
                    # 消费和储蓄决策（先消费后储蓄）
                    consumption_floor_spending = cS.cFloor * (1 + cS.tau_c)
                    total_available_resources = max(0, total_income)
                    
                    # 根据消费比例确定消费支出
                    if total_available_resources >= consumption_floor_spending:
                        resources_above_floor = total_available_resources - consumption_floor_spending
                        consumption_expenditure = consumption_floor_spending + prop_consumption * resources_above_floor
                    else:
                        consumption_expenditure = total_available_resources
                    
                    # 确保消费支出合理
                    consumption_expenditure = min(consumption_expenditure, total_available_resources)
                    consumption_expenditure = max(consumption_expenditure, consumption_floor_spending)
                    
                    # 计算实际消费和储蓄
                    c_consumption = max(cS.cFloor, consumption_expenditure / (1 + cS.tau_c))
                    k_next = total_available_resources - consumption_expenditure
                    k_next = max(cS.kMin, min(k_next, cS.kMax))
                
                # 记录结果
                kHistM_rl[i_sim, a_group] = k_current
                kPpsHistM_rl[i_sim, a_group] = kpps_current
                cHistM_rl[i_sim, a_group] = c_consumption
                cppsHistM_rl[i_sim, a_group] = cpps_decision
                
                # 更新PPS资产
                if cS.pps_active:
                    pps_withdrawal = 0
                    if is_pps_withdrawal_eligible:
                        pps_withdrawal = kpps_current * cS.pps_withdrawal_rate
                    
                    kpps_next = (kpps_current + cpps_decision - pps_withdrawal) * pps_return_factor
                    kpps_next = max(cS.kppsMin, min(kpps_next, cS.kppsMax))
                else:
                    kpps_next = cS.kppsMin
                
                # 更新状态
                k_current = k_next
                kpps_current = kpps_next
        
        return kHistM_rl, kPpsHistM_rl, cHistM_rl, cppsHistM_rl


# 🔧 新增：定义临时参数类（与主版本保持一致）
class TempParamSHH:
    """
    临时参数结构体类 - 用于模拟MATLAB中的结构体
    与HHIncome_Huggett函数兼容
    """
    def __init__(self, tau_l, theta_payg_actual_for_hh, pps_tax_deferral_active, ageEffV_new):
        self.tau_l = tau_l
        self.theta_payg_actual_for_hh = theta_payg_actual_for_hh
        self.pps_tax_deferral_active = pps_tax_deferral_active
        self.ageEffV_new = ageEffV_new
        
        # PPS提取税收相关参数
        self.tau_k = 0.2
        self.pps_tax_rate_withdrawal = 0.15
        self.current_pps_withdrawal = 0