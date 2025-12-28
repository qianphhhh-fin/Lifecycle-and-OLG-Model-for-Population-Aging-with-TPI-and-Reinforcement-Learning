# 文件名: v3_cocco_env.py

import gymnasium as gym
from gymnasium import spaces
import numpy as np
from typing import Optional

class CoccoEnvV3(gym.Env):
    """
    Cocco et al. (2005) Life-Cycle Model Environment - V3 (Robust Version).

    本版本集成了多层防御体系以确保数值稳定性，同时保持理论纯粹性：
    1.  **对数动作空间 (Log Action Space)**: 智能体选择一个潜在动作，
        该动作被解码为对数消费 c = log(C)，从根源上避免 C<=0。
    2.  **消费截断 (Consumption Clipping)**: 在计算效用前，强制消费 C
        不低于一个经过校准的绝对下限 C_min，防止梯度爆炸。
    3.  **仿射奖励变换 (Affine Reward Transformation)**: 将原始CRRA效用
        通过 U' = a*U + b 的形式进行线性变换，在不改变风险规避曲率的
        前提下，将奖励信号映射到对神经网络友好的数值区间。
    """
    metadata = {'render.modes': []}

    def __init__(self, disable_random_death: bool = False):
        super(CoccoEnvV3, self).__init__()
        
        self.disable_random_death = disable_random_death

        # --- A. 经济参数 ---
        self.tb, self.tr, self.td = 22, 61, 100
        self.tn = self.td - self.tb + 1
        self.beta, self.gamma = 0.95, 3.84
        self.rf, self.mu, self.sigr = 1.02, 0.04, 0.27
        self.ret_fac = 0.6827
        self.aa = -2.170042 + 2.700381
        self.b1, self.b2, self.b3 = 0.16818, -0.0323371/10, 0.0019704/100
        self.smay, self.smav = np.sqrt(0.169993), np.sqrt(0.112572)

        # --- B. 自动参数校准 ---
        self._calibrate_parameters()

        # --- C. 派生参数与预计算 ---
        self._precompute_helpers()

        # --- D. 状态变量 ---
        self.W, self.P, self.P_retirement, self.age, self.is_done = 0.0, 1.0, 1.0, self.tb, False

        # --- E. Gym 空间定义 ---
        # 动作1: 潜在的对数消费动作 (将被解码)
        # 动作2: 风险资产比例
        self.action_space = spaces.Box(
            low=np.array([-1.0, 0.0]), high=np.array([1.0, 1.0]),
            shape=(2,), dtype=np.float32
        )
        # 观测: [归一化财富, 归一化年龄]
        self.observation_space = spaces.Box(
            low=np.array([0.0, -1.0]), high=np.array([1000.0, 1.0]),
            shape=(2,), dtype=np.float32
        )

    def _calibrate_parameters(self):
        """自动校准截断和奖励变换所需的参数。"""
        # 使用工作期间的确定性收入分布作为消费分布的代理
        det_incomes = [np.exp(self.aa + self.b1*age + self.b2*age**2 + self.b3*age**3) for age in range(self.tb, self.tr)]
        
        # 1. 校准消费截断值 C_min
        MIN_C_PERCENTILE = 5.0
        self.MIN_CONSUMPTION = np.percentile(det_incomes, MIN_C_PERCENTILE)

        # 2. 校准仿射变换的锚点
        REF_C_PERCENTILE = 50.0
        self.REFERENCE_CONSUMPTION = np.percentile(det_incomes, REF_C_PERCENTILE)
        
        # 3. 校准对数动作空间的解码范围
        # 允许智能体选择的消费范围是 [C_min, C_ref * 5]
        self.C_DECODE_LOWER = self.MIN_CONSUMPTION
        self.C_DECODE_UPPER = self.REFERENCE_CONSUMPTION * 5.0
        self.c_decode_lower = np.log(self.C_DECODE_LOWER)
        self.c_decode_upper = np.log(self.C_DECODE_UPPER)

        # 4. 校准仿射变换参数 a 和 b
        REWARD_FLOOR, REFERENCE_REWARD = -2.0, 1.0
        u_min = self._raw_utility(self.MIN_CONSUMPTION)
        u_ref = self._raw_utility(self.REFERENCE_CONSUMPTION)
        
        self.reward_scale_a = (REFERENCE_REWARD - REWARD_FLOOR) / (u_ref - u_min)
        self.reward_shift_b = REWARD_FLOOR - self.reward_scale_a * u_min
        
        print("\n--- CoccoEnvV3 Initialized with Calibrated Parameters ---")
        print(f"Consumption Clipping Lower Bound (C_min): {self.MIN_CONSUMPTION:.4f}")
        print(f"Log-Action Decode Range (C): [{self.C_DECODE_LOWER:.2f}, {self.C_DECODE_UPPER:.2f}]")
        print(f"Affine Reward Params: a(scale)={self.reward_scale_a:.4f}, b(shift)={self.reward_shift_b:.4f}")
        print(f" -> Verification: R(C_min)={self.reward_scale_a * u_min + self.reward_shift_b:.2f} (target {REWARD_FLOOR})")
        print(f" -> Verification: R(C_ref)={self.reward_scale_a * u_ref + self.reward_shift_b:.2f} (target {REFERENCE_REWARD})")
        print("-------------------------------------------------------\n")

    def _precompute_helpers(self):
        """预计算生存概率和确定性收入剖面。"""
        survprob_data = [0.99975, 0.99974, 0.99973, 0.99972, 0.99971, 0.99969, 0.99968, 0.99966, 0.99963, 0.99961, 0.99958, 0.99954, 0.99951, 0.99947, 0.99943, 0.99938, 0.99934, 0.99928, 0.99922, 0.99916, 0.99908, 0.999, 0.99891, 0.99881, 0.9987, 0.99857, 0.99843, 0.99828, 0.99812, 0.99794, 0.99776, 0.99756, 0.99735, 0.99713, 0.9969, 0.99666, 0.9964, 0.99613, 0.99583, 0.99551, 0.99515, 0.99476, 0.99432, 0.99383, 0.9933, 0.9927, 0.99205, 0.99133, 0.99053, 0.98961, 0.98852, 0.98718, 0.98553, 0.98346, 0.98089, 0.97772, 0.97391, 0.96943, 0.96429, 0.95854, 0.95221, 0.94537, 0.93805, 0.93027, 0.92202, 0.91327, 0.90393, 0.89389, 0.88304, 0.87126, 0.85846, 0.84452, 0.82935, 0.81282, 0.79485, 0.77543, 0.75458, 0.7324, 0.70893, 0.68424, 0.68424, 0.68424]
        self.survprob = np.ones(self.tn - 1)
        copy_len = min(len(self.survprob), len(survprob_data))
        self.survprob[:copy_len] = survprob_data[:copy_len]
        self.f_y = np.zeros(self.tn)
        for i in range(self.tn):
            age = self.tb + i
            if age < self.tr: self.f_y[i] = np.exp(self.aa + self.b1*age + self.b2*age**2 + self.b3*age**3)

    def _raw_utility(self, C):
        if self.gamma == 1: return np.log(C)
        else: return (C**(1 - self.gamma)) / (1 - self.gamma)

    def _get_obs(self):
        normalized_age = 2 * (self.age - self.tb) / (self.td - self.tb) - 1
        normalized_wealth = self.W / self.P
        return np.array([normalized_wealth, normalized_age], dtype=np.float32)

    def _get_info(self, C, Y, alpha):
        return {"age": self.age, "absolute_wealth": self.W, "absolute_consumption": C,
                "absolute_income": Y, "risky_share": alpha, "permanent_income": self.P}

    def reset(self, *, seed: Optional[int] = None, options: Optional[dict] = None):
        super().reset(seed=seed)
        self.age = self.tb
        self.is_done = False
        self.W = 0.0
        self.P = np.exp(np.random.normal(0, self.smav))
        return self._get_obs(), self._get_info(0, 0, 0)

    def step(self, action):
        if self.is_done:
            obs, info = self.reset()
            return obs, 0.0, True, False, info

        # 1. 解码动作
        latent_c_action, alpha = action
        alpha = np.clip(alpha, 0.0, 1.0)
        
        # --- 层防御 1: 对数动作空间解码 ---
        # 将[-1, 1]的潜在动作映射到[c_lower, c_upper]的对数消费空间
        c = self.c_decode_lower + (latent_c_action + 1.0) * 0.5 * (self.c_decode_upper - self.c_decode_lower)
        intended_C = np.exp(c)

        # 2. 计算当期资源
        current_age_idx = self.age - self.tb
        if self.age < self.tr:
            u_shock = np.random.normal(0, self.smay)
            Y = self.f_y[current_age_idx] * self.P * np.exp(u_shock)
        else:
            Y = self.ret_fac * self.P_retirement
        X = self.W + Y

        # 3. 计算最终消费
        # 智能体意图消费intended_C，但不能超过总资源
        C = min(intended_C, X)
        
        # --- 第二层防御: 消费截断 ---
        C_clipped = max(C, self.MIN_CONSUMPTION)

        # 4. 计算奖励
        raw_utility = self._raw_utility(C_clipped)
        # --- 第三层防御: 仿射变换 ---
        reward = self.reward_scale_a * raw_utility + self.reward_shift_b

        # 5. 世界演化
        S = X - C_clipped
        r_shock = np.random.normal(0, self.sigr)
        portfolio_return = alpha * max(0, self.rf + self.mu + r_shock) + (1 - alpha) * self.rf
        self.W = S * portfolio_return

        # 6. 永久收入演化
        if self.age < self.tr:
            z_shock = np.random.normal(0, self.smav)
            g_growth = self.f_y[current_age_idx + 1] / self.f_y[current_age_idx] if self.age < self.tr - 1 else 1.0
            self.P *= g_growth * np.exp(z_shock)
            if self.age == self.tr - 1: self.P_retirement = self.P
        
        self.age += 1

        # 7. 终止检查
        died = False
        if self.age > self.td: died = True
        elif not self.disable_random_death and self.np_random.uniform() > self.survprob[current_age_idx]:
            died = True
        self.is_done = died
        
        return self._get_obs(), float(reward), self.is_done, False, self._get_info(C_clipped, Y, alpha)