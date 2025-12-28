import gymnasium as gym
from gymnasium import spaces
import numpy as np
from typing import Optional
from scipy.stats import norm # <-- 新增导入

import gymnasium as gym
from gymnasium import spaces
import numpy as np
from typing import Optional
from scipy.stats import norm

class CoccoEnvV2(gym.Env):
    """
    Cocco et al. (2005) Life-Cycle Model Environment - V2.
    """
    metadata = {'render.modes': []}

    def __init__(self, disable_random_death: bool = False):
        super(CoccoEnvV2, self).__init__()
        
        self.disable_random_death = disable_random_death

        # ... [其他参数保持不变] ...
        self.tb = 22
        self.tr = 25
        self.td = 30
        self.tn = self.td - self.tb + 1
        self.beta = 0.95
        self.gamma = 3.84
        self.aa = -2.170042 + 2.700381
        self.b1 = 0.16818
        self.b2 = -0.0323371 / 10
        self.b3 = 0.0019704 / 100
        self.smay = np.sqrt(0.169993)
        self.smav = np.sqrt(0.112572)
        self.rf = 1.02
        self.mu = 0.04
        self.sigr = 0.27
        self.ret_fac = 0.6827
        self.tau_y = 0.06
        survprob_data = [
            0.99975, 0.99974, 0.99973, 0.99972, 0.99971, 0.99969, 0.99968, 0.99966, 0.99963, 0.99961,
            0.99958, 0.99954, 0.99951, 0.99947, 0.99943, 0.99938, 0.99934, 0.99928, 0.99922, 0.99916,
            0.99908, 0.999, 0.99891, 0.99881, 0.9987, 0.99857, 0.99843, 0.99828, 0.99812, 0.99794,
            0.99776, 0.99756, 0.99735, 0.99713, 0.9969, 0.99666, 0.9964, 0.99613, 0.99583, 0.99551,
            0.99515, 0.99476, 0.99432, 0.99383, 0.9933, 0.9927, 0.99205, 0.99133, 0.99053, 0.98961,
            0.98852, 0.98718, 0.98553, 0.98346, 0.98089, 0.97772, 0.97391, 0.96943, 0.96429, 0.95854,
            0.95221, 0.94537, 0.93805, 0.93027, 0.92202, 0.91327, 0.90393, 0.89389, 0.88304, 0.87126,
            0.85846, 0.84452, 0.82935, 0.81282, 0.79485, 0.77543, 0.75458, 0.7324, 0.70893, 0.68424,
            0.68424, 0.68424
        ]
        self.survprob = np.ones(self.tn - 1)
        copy_len = min(len(self.survprob), len(survprob_data))
        self.survprob[:copy_len] = survprob_data[:copy_len]
        self.f_y = np.zeros(self.tn)
        for i in range(self.tn):
            age = self.tb + i
            if age <= self.tr:
                 self.f_y[i] = np.exp(self.aa + self.b1*age + self.b2*age**2 + self.b3*age**3)/10
        n_shocks = 5
        tauchen_q = 2.0 
        self.shock_grid = np.linspace(-tauchen_q, tauchen_q, n_shocks)
        omega = self.shock_grid[1] - self.shock_grid[0]
        self.shock_weights = np.zeros(n_shocks)
        for i in range(n_shocks):
            if i == 0:
                self.shock_weights[i] = norm.cdf(self.shock_grid[i] + omega/2)
            elif i == n_shocks - 1:
                self.shock_weights[i] = 1 - norm.cdf(self.shock_grid[i] - omega/2)
            else:
                self.shock_weights[i] = norm.cdf(self.shock_grid[i] + omega/2) - norm.cdf(self.shock_grid[i] - omega/2)
        self.shock_weights /= np.sum(self.shock_weights)

        self.W = 0.0
        self.X = 0.0
        self.P = 1.0
        self.P_retirement = 1.0
        self.age = self.tb
        self.is_done = False

        self.action_space = spaces.Box(
            low=np.array([0.05, 0.0]),
            high=np.array([1.0, 1.0]),
            shape=(2,),
            dtype=np.float32
        )

        # --- 核心修改: 观测空间增加 W/P ratio ---
        # [cash_on_hand, permanent_income, normalized_age, wealth_to_income_ratio]
        self.observation_space = spaces.Box(
            low=np.array([0.0, 0.0, -1.0, 0.0]),
            high=np.array([10000.0, 1000.0, 1.0, 1000.0]), 
            shape=(4,),
            dtype=np.float32
        )

    def _get_obs(self):
        normalized_age = 2 * (self.age - self.tb) / (self.td - self.tb) - 1
        # P 可能极小，避免除以零
        wealth_to_income_ratio = self.W / (self.P + 1e-6)
        return np.array([self.X, self.P, normalized_age, wealth_to_income_ratio], dtype=np.float32)

    def _get_info(self, C, Y, alpha, W_pre_shock):
        return {
            "age": self.age,
            "absolute_wealth": W_pre_shock,
            "cash_on_hand": W_pre_shock + Y,
            "absolute_consumption": C,
            "absolute_income": Y,
            "risky_share": alpha,
            "permanent_income": self.P
        }

    def reset(self, *, seed: Optional[int] = None, options: Optional[dict] = None):
        super().reset(seed=seed)
        self.age = self.tb
        self.is_done = False
        self.P = 1.0
        self.P_retirement = 1.0
        self.W = 0.0

        current_age_idx = self.age - self.tb
        u_shock_std = self.np_random.choice(self.shock_grid, p=self.shock_weights)
        u_shock = u_shock_std * self.smay
        det_Y = self.f_y[current_age_idx]
        Y_gross = det_Y * self.P * np.exp(u_shock)
        Y = Y_gross * (1 - self.tau_y)
        self.X = self.W + Y
        
        info = self._get_info(C=0, Y=Y, alpha=0, W_pre_shock=self.W) 
        return self._get_obs(), info

    def step(self, action):
        if self.is_done:
            obs, info = self.reset()
            return obs, 0.0, True, False, info

        # --- 记录决策前的财富 W_t ---
        W_before_step = self.W

        c_prop, alpha = action
        c_prop = np.clip(c_prop, 0.0, 1.0)
        alpha = np.clip(alpha, 0.0, 1.0)

        C = c_prop * self.X
        C = min(C, self.X * 0.99999)
        if C < 1e-6: C = 1e-6

        if self.gamma == 1:
            utility = np.log(C)
        else:
            utility = (C**(1 - self.gamma)) / (1 - self.gamma)
        reward = np.maximum(utility,-150)
        
        S = self.X - C
        r_shock_std = self.np_random.choice(self.shock_grid, p=self.shock_weights)
        r_shock = r_shock_std * self.sigr
        risky_return = self.rf + self.mu + r_shock
        portfolio_return = alpha * risky_return + (1 - alpha) * self.rf
        
        W_next = S * portfolio_return

        current_age = self.age
        next_age = current_age + 1
        
        if next_age < self.tr:
            z_shock_std = self.np_random.choice(self.shock_grid, p=self.shock_weights)
            z_shock = z_shock_std * self.smav
            self.P = self.P * np.exp(z_shock)
            if next_age == self.tr:
                self.P_retirement = self.P
        
        self.age = next_age

        died = False
        if self.age >= self.td:
            died = True
        else:
            if not self.disable_random_death:
                surv_prob_idx = current_age - self.tb
                if surv_prob_idx < len(self.survprob):
                    survival_prob = self.survprob[surv_prob_idx]
                    if self.np_random.uniform() > survival_prob:
                        died = True
                else:
                    died = True
        self.is_done = died

        Y_next = 0.0
        if not self.is_done:
            next_age_idx = self.age - self.tb
            if self.age < self.tr:
                u_shock_std = self.np_random.choice(self.shock_grid, p=self.shock_weights)
                u_shock = u_shock_std * self.smay
                det_Y = self.f_y[next_age_idx]
                Y_gross = det_Y * self.P * np.exp(u_shock)
                Y_next = Y_gross * (1 - self.tau_y)
            else:
                Y_next = self.ret_fac * self.f_y[self.tr-self.tb] * self.P_retirement
        
        self.W = W_next
        self.X = self.W + Y_next
        
        # --- info 现在使用决策前的财富 W_before_step ---
        info = self._get_info(C=C, Y=self.X - self.W, alpha=alpha, W_pre_shock=W_before_step)
        
        return self._get_obs(), float(reward), self.is_done, False, info
        
# class CoccoEnvV2(gym.Env):
#     """
#     Cocco et al. (2005) Life-Cycle Model Environment - V2.

#     版本 V2 特点:
#     - **使用绝对单位 (primitive units)** 进行奖励计算和状态观测。
#     - 观测 (to agent): [absolute_wealth, permanent_income, normalized_age]
#     - 动作 (from agent): [consumption_proportion, risky_share]
#     - 奖励 (reward): 基于绝对消费 U(C) 计算，并进行数值缩放。
    
#     此版本直接面对CRRA效用函数的数值挑战。
#     """
#     metadata = {'render.modes': []}

# # v2_cocco_env.py

#     def __init__(self, disable_random_death: bool = False):
#         super(CoccoEnvV2, self).__init__()
        
#         self.disable_random_death = disable_random_death

#         # --- A. 生命周期与人口结构 ---
#         self.tb = 22
#         self.tr = 25
#         self.td = 30
#         self.tn = self.td - self.tb + 1

#         # --- B. 经济主体偏好 (CRRA) ---
#         self.beta = 0.95
#         self.gamma = 3.84

#         # --- C. 收入过程 ---
#         self.aa = -2.170042 + 2.700381
#         self.b1 = 0.16818
#         self.b2 = -0.0323371 / 10
#         self.b3 = 0.0019704 / 100
#         self.smay = np.sqrt(0.169993)
#         self.smav = np.sqrt(0.112572)

#         # --- D. 资产与回报率 ---
#         self.rf = 1.02
#         self.mu = 0.04
#         self.sigr = 0.27

#         # --- E. 养老金与税收 ---
#         self.ret_fac = 0.6827
#         self.tau_y = 0.06

#         # --- G. 派生参数与预计算 ---
#         survprob_data = [
#             0.99975, 0.99974, 0.99973, 0.99972, 0.99971, 0.99969, 0.99968, 0.99966, 0.99963, 0.99961,
#             0.99958, 0.99954, 0.99951, 0.99947, 0.99943, 0.99938, 0.99934, 0.99928, 0.99922, 0.99916,
#             0.99908, 0.999, 0.99891, 0.99881, 0.9987, 0.99857, 0.99843, 0.99828, 0.99812, 0.99794,
#             0.99776, 0.99756, 0.99735, 0.99713, 0.9969, 0.99666, 0.9964, 0.99613, 0.99583, 0.99551,
#             0.99515, 0.99476, 0.99432, 0.99383, 0.9933, 0.9927, 0.99205, 0.99133, 0.99053, 0.98961,
#             0.98852, 0.98718, 0.98553, 0.98346, 0.98089, 0.97772, 0.97391, 0.96943, 0.96429, 0.95854,
#             0.95221, 0.94537, 0.93805, 0.93027, 0.92202, 0.91327, 0.90393, 0.89389, 0.88304, 0.87126,
#             0.85846, 0.84452, 0.82935, 0.81282, 0.79485, 0.77543, 0.75458, 0.7324, 0.70893, 0.68424,
#             0.68424, 0.68424
#         ]
#         self.survprob = np.ones(self.tn - 1)
#         copy_len = min(len(self.survprob), len(survprob_data))
#         self.survprob[:copy_len] = survprob_data[:copy_len]

#         self.f_y = np.zeros(self.tn)
#         for i in range(self.tn):
#             age = self.tb + i
#             if age <= self.tr:
#                  self.f_y[i] = np.exp(self.aa + self.b1*age + self.b2*age**2 + self.b3*age**3)/10
        
#         n_shocks = 5
#         tauchen_q = 2.0 
#         self.shock_grid = np.linspace(-tauchen_q, tauchen_q, n_shocks)
        
#         omega = self.shock_grid[1] - self.shock_grid[0]
#         self.shock_weights = np.zeros(n_shocks)
#         for i in range(n_shocks):
#             if i == 0:
#                 self.shock_weights[i] = norm.cdf(self.shock_grid[i] + omega/2)
#             elif i == n_shocks - 1:
#                 self.shock_weights[i] = 1 - norm.cdf(self.shock_grid[i] - omega/2)
#             else:
#                 self.shock_weights[i] = norm.cdf(self.shock_grid[i] + omega/2) - norm.cdf(self.shock_grid[i] - omega/2)
#         self.shock_weights /= np.sum(self.shock_weights)

#         # --- 核心修改：状态变量现在包含现金持有量X ---
#         self.W = 0.0  # 期初财富
#         self.X = 0.0  # 现金持有量 (W+Y)
#         self.P = 1.0
#         self.P_retirement = 1.0
#         self.age = self.tb
#         self.is_done = False

#         self.action_space = spaces.Box(
#             low=np.array([0.05, 0.0]),
#             high=np.array([1.0, 1.0]),
#             shape=(2,),
#             dtype=np.float32
#         )

#         # --- 核心修改: 观测空间现在是 [cash_on_hand, permanent_income, normalized_age] ---
#         self.observation_space = spaces.Box(
#             low=np.array([0.0, 0.0, -1.0]),
#             high=np.array([10000.0, 1000.0, 1.0]), 
#             shape=(3,),
#             dtype=np.float32
#         )


# # v2_cocco_env.py

#     def _get_obs(self):
#         # 观测值现在是现金持有量 X
#         normalized_age = 2 * (self.age - self.tb) / (self.td - self.tb) - 1
#         return np.array([self.X, self.P, normalized_age], dtype=np.float32)

#     def _get_info(self, C, Y, alpha, W_pre_shock):
#         # info中的 'absolute_wealth' 现在明确指期初财富
#         return {
#             "age": self.age,
#             "absolute_wealth": W_pre_shock,
#             "cash_on_hand": W_pre_shock + Y,
#             "absolute_consumption": C,
#             "absolute_income": Y,
#             "risky_share": alpha,
#             "permanent_income": self.P
#         }

# # v2_cocco_env.py

#     def reset(self, *, seed: Optional[int] = None, options: Optional[dict] = None):
#         super().reset(seed=seed)

#         # 1. 初始化基本状态
#         self.age = self.tb
#         self.is_done = False
#         self.P = 1.0
#         self.P_retirement = 1.0
#         self.W = 0.0  # 期初财富为0

#         # 2. **核心修改**: 为第一期决策准备信息
#         # 计算第一期的收入 Y_t
#         current_age_idx = self.age - self.tb
#         u_shock_std = self.np_random.choice(self.shock_grid, p=self.shock_weights)
#         u_shock = u_shock_std * self.smay
#         det_Y = self.f_y[current_age_idx]
#         Y_gross = det_Y * self.P * np.exp(u_shock)
#         Y = Y_gross * (1 - self.tau_y)

#         # 计算第一期的现金持有量 X_t，并存储
#         self.X = self.W + Y
        
#         # 3. 返回正确的初始观测值
#         # 初始info中的消费等值为0，因为决策尚未做出
#         info = self._get_info(C=0, Y=Y, alpha=0, W_pre_shock=self.W) 
#         return self._get_obs(), info

# # v2_cocco_env.py

# # v2_cocco_env.py

#     def step(self, action):
#         if self.is_done:
#             obs, info = self.reset()
#             return obs, 0.0, True, False, info

#         # --- 1. 应用决策 ---
#         # 此时 self.X 是智能体做决策时看到的现金持有量
#         c_prop, alpha = action
#         c_prop = np.clip(c_prop, 0.0, 1.0)
#         alpha = np.clip(alpha, 0.0, 1.0)

#         C = c_prop * self.X
#         C = min(C, self.X * 0.99999)
#         if C < 1e-6: C = 1e-6

#         # --- 2. 计算奖励 ---
#         if self.gamma == 1:
#             utility = np.log(C)
#         else:
#             utility = (C**(1 - self.gamma)) / (1 - self.gamma)
#         reward = np.maximum(utility,-150)
#         # reward = utility
        
#         # --- 3. 状态转移：计算下一期的期初财富 W_{t+1} ---
#         S = self.X - C
#         r_shock_std = self.np_random.choice(self.shock_grid, p=self.shock_weights)
#         r_shock = r_shock_std * self.sigr
#         risky_return = self.rf + self.mu + r_shock
#         portfolio_return = alpha * risky_return + (1 - alpha) * self.rf
        
#         W_next = S * portfolio_return # 这是下一期的期初财富 W_{t+1}

#         # --- 4. 状态转移：计算下一期的永久收入 P_{t+1} 和年龄 ---
#         current_age = self.age
#         next_age = current_age + 1
        
#         if next_age < self.tr:
#             z_shock_std = self.np_random.choice(self.shock_grid, p=self.shock_weights)
#             z_shock = z_shock_std * self.smav
#             self.P = self.P * np.exp(z_shock)
#             if next_age == self.tr:
#                 self.P_retirement = self.P
        
#         self.age = next_age

#         # --- 5. 检查是否终结 ---
#         died = False
#         if self.age >= self.td:
#             died = True
#         else:
#             if not self.disable_random_death:
#                 surv_prob_idx = current_age - self.tb
#                 if surv_prob_idx < len(self.survprob):
#                     survival_prob = self.survprob[surv_prob_idx]
#                     if self.np_random.uniform() > survival_prob:
#                         died = True
#                 else:
#                     died = True
#         self.is_done = died

#         # --- 6. 为下一期决策准备信息 ---
#         Y_next = 0.0
#         if not self.is_done:
#             # 计算下一期的收入 Y_{t+1}
#             next_age_idx = self.age - self.tb
#             if self.age < self.tr:
#                 u_shock_std = self.np_random.choice(self.shock_grid, p=self.shock_weights)
#                 u_shock = u_shock_std * self.smay
#                 det_Y = self.f_y[next_age_idx]
#                 Y_gross = det_Y * self.P * np.exp(u_shock)
#                 Y_next = Y_gross * (1 - self.tau_y)
#             else:
#                 Y_next = self.ret_fac * self.f_y[self.tr-self.tb] * self.P_retirement
        
#         # 更新下一期的现金持有量 X_{t+1}
#         self.W = W_next # W 更新为下一期的期初财富
#         self.X = self.W + Y_next # X 更新为下一期的现金持有量
        
#         # --- 7. 整理并返回结果 ---
#         # 注意: info 记录的是本期t的决策和结果
#         # Y_current 是 self.X - W_before_step, 但我们没有保存 W_before_step
#         # 为了info的准确性，我们需要调整get_info的参数
#         info = self._get_info(C=C, Y=self.X - self.W, alpha=alpha, W_pre_shock=W_next - Y_next)
        
#         return self._get_obs(), float(reward), self.is_done, False, info