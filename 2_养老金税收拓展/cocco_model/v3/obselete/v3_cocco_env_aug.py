import gymnasium as gym
from gymnasium import spaces
import numpy as np
from typing import Optional
from scipy.stats import norm

class CoccoEnvV3_Augmented(gym.Env):
    """
    Cocco et al. (2005) Life-Cycle Model Environment - V3 (Augmented State).

    此版本根据 Bisi et al. (2022) 的思想，通过“状态增强”修改了环境，
    以便标准的风险中性 RL 算法可以优化风险敏感的目标。

    核心改动:
    1.  **状态空间增强**: 观测空间增加了两个维度:
        - `v`: 累计折现效用 (cumulative discounted utility)。
        - `w`: 当前的折现因子 (discount factor tracker, beta^t)。
        因此，观测空间从 (3,) 变为 (5,)。

    2.  **稀疏奖励结构**: 每一步(step)的奖励都为 0。只有在 episode 结束时，
        环境才会返回一个等于整个生命周期累计折现效用的奖励。
        这使得智能体的目标直接就是最大化整个生命周期的总效用。

    3.  **内部状态追踪**: 环境内部需额外追踪 `self.v` 和 `self.w`。
    """
    metadata = {'render.modes': []}

    def __init__(self, disable_random_death: bool = False):
        super(CoccoEnvV3_Augmented, self).__init__()

        self.disable_random_death = disable_random_death

        # --- A. 生命周期与人口结构 (与 v2 一致) ---
        self.tb = 22
        self.tr = 26
        self.td = 30
        self.tn = self.td - self.tb + 1

        # --- B. 经济主体偏好 (与 v2 一致) ---
        self.beta = 0.95
        self.gamma = 3.84

        # --- C. 收入过程参数 (与 v2 一致) ---
        self.aa = -2.170042 + 2.700381
        self.b1 = 0.16818
        self.b2 = -0.0323371 / 10
        self.b3 = 0.0019704 / 100
        self.smay = 0.169993**0.5
        self.smav = 0.112572**0.5

        # --- D. 资产与回报率 (与 v2 一致) ---
        self.rf = 1.02
        self.mu = 0.04
        self.sigr = 0.27

        # --- E. 养老金与税收 (与 v2 一致) ---
        self.ret_fac = 0.6827
        self.tau_y = 0.0

        # --- F. 生存概率 (与 v2 一致) ---
        survprob_data = [
            0.99975, 0.99974, 0.99973, 0.99972, 0.99971, 0.99969, 0.99968, 0.99966, 0.99963, 0.99961, 0.99958,
            0.99954, 0.99951, 0.99947, 0.99943, 0.99938, 0.99934, 0.99928, 0.99922, 0.99916, 0.99908, 0.999,
            0.99891, 0.99881, 0.9987, 0.99857, 0.99843, 0.99828, 0.99812, 0.99794, 0.99776, 0.99756, 0.99735,
            0.99713, 0.9969, 0.99666, 0.9964, 0.99613, 0.99583, 0.99551, 0.99515, 0.99476, 0.99432, 0.99383,
            0.9933, 0.9927, 0.99205, 0.99133, 0.99053, 0.98961, 0.8852, 0.98718, 0.98553, 0.98346, 0.98089,
            0.97772, 0.97391, 0.96943, 0.96429, 0.95854, 0.95221, 0.94537, 0.93805, 0.93027, 0.92202, 0.91327,
            0.90393, 0.89389, 0.88304, 0.87126, 0.85846, 0.84452, 0.82935, 0.81282, 0.79485, 0.77543, 0.75458,
            0.7324, 0.70893, 0.68424, 0.68424, 0.68424
        ]
        self.survprob = np.ones(self.tn - 1)
        copy_len = min(len(self.survprob), len(survprob_data))
        self.survprob[:copy_len] = survprob_data[:copy_len]

        # --- G. 派生参数与预计算 (与 v2 一致) ---
        self.f_y = np.zeros(self.tn)
        for i in range(self.tn):
            age = self.tb + i
            if age <= self.tr:
                self.f_y[i] = (self.aa + self.b1 * age + self.b2 * age**2 + self.b3 * age**3) - 3

        n_shocks = 5
        self.shock_grid, self.shock_weights = self._discretize_shocks(n_shocks)

        # --- H. 状态变量 ---
        self.W = 0.0
        self.Y = 0.0
        self.X = 0.0
        self.P = 0.0
        self.P_retirement = 0.0
        self.age = self.tb
        self.is_done = False
        
        # **[V3 新增]** 增强状态的内部追踪器
        self.v = 0.0  # 累计折现效用
        self.w = 1.0  # 折现因子追踪器

        # --- I. Gym 空间定义 ---
        self.action_space = spaces.Box(
            low=np.array([0.2, 0.05]),
            high=np.array([1.1, 1.1]),
            shape=(2,),
            dtype=np.float32
        )
        
        # **[V3 修改]** 观测空间维度从 3 扩展到 5
        # 原始: [X_norm, P, normalized_age]
        # 增强: [X_norm, P, normalized_age, v, w]
        self.observation_space = spaces.Box(
            low=np.array([0.0, -np.inf, -1.0, -np.inf, 0.0]),
            high=np.array([np.inf, np.inf, 1.0, np.inf, 1.0]),
            shape=(5,),
            dtype=np.float32
        )

    def _discretize_shocks(self, n_shocks):
        if n_shocks == 1:
            return np.array([0.0]), np.array([1.0])
        tauchen_q = 2.0
        grid = np.linspace(-tauchen_q, tauchen_q, n_shocks)
        omega = grid[1] - grid[0]
        weights = np.zeros(n_shocks)
        for i in range(n_shocks):
            if i == 0: weights[i] = norm.cdf(grid[i] + omega / 2)
            elif i == n_shocks - 1: weights[i] = 1 - norm.cdf(grid[i] - omega / 2)
            else: weights[i] = norm.cdf(grid[i] + omega / 2) - norm.cdf(grid[i] - omega / 2)
        return grid, weights / np.sum(weights)

    def _get_obs(self):
        """
        **[V3 修改]** 返回包含 v 和 w 的增强状态。
        """
        normalized_age = 2 * (self.age - self.tb) / (self.td - self.tb) - 1
        return np.array([self.X_norm, self.P, normalized_age, self.v, self.w], dtype=np.float32)

    def reset(self, *, seed: Optional[int] = None, options: Optional[dict] = None):
        super().reset(seed=seed)
        
        self.age = self.tb
        self.is_done = False
        self.P_retirement = 0.0
        
        # **[V3 修改]** 重置增强状态追踪器
        self.v = 0.0
        self.w = 1.0
        
        z_shock = self.np_random.choice(self.shock_grid, p=self.shock_weights) * self.smav
        self.P = z_shock

        self.W = 0.0
        u_shock = self.np_random.choice(self.shock_grid, p=self.shock_weights) * self.smay
        self.Y = np.exp(self.f_y[0] + self.P + u_shock)
        
        self.X = self.W + self.Y
        self.gyp = np.exp(self.f_y[0] + self.P)
        self.X_norm = self.X/self.gyp

        info = {
            "age": self.age,
            "absolute_wealth": self.W,
            "raw_permant_shock": self.P,
            "absolute_consumption": 0.0,
            "absolute_income": self.Y,
            "c_prop": 0.0,
            "risky_share": 0.0,
            "normalized_coh": self.X_norm,
            "period_utility": 0.0 # 新增用于调试
        }
        return self._get_obs(), info

    def step(self, action):
        if self.is_done:
            obs, info = self.reset()
            return obs, 0.0, True, False, info

        # 保存 t 期状态用于构建 info 字典
        W_t = self.W
        Y_t = self.Y
        age_t = self.age
        current_age_idx = age_t - self.tb

        # 0. 剪裁动作
        c_prop = np.clip(action[0], 0.0, 1.0)
        alpha = np.clip(action[1], 0.0, 1.0)

        # 1. 决策
        C = c_prop * self.X
        C = np.minimum(C, self.X * 0.99999)
        # 如果 X 接近于0，c_prop 可能为 nan，进行修正
        c_prop = C / self.X if self.X > 1e-9 else 0.0 
        if C < 1e-6: C = 1e-6
        S = self.X - C

        # 2. **[V3 修改]** 计算当期效用，但不作为奖励返回
        if self.gamma == 1: 
            period_utility = np.log(C)
        else: 
            period_utility = ((C)**(1 - self.gamma)) / (1 - self.gamma) + 1
        period_utility = np.maximum(period_utility, -30)    
        
        # **[V3 修改]** 更新增强状态
        self.v += self.w * period_utility
        self.w *= self.beta

        # 3. 计算 W_{t+1}
        r_shock = self.np_random.choice(self.shock_grid, p=self.shock_weights) * self.sigr
        portfolio_return = alpha * (self.rf + self.mu + r_shock) + (1 - alpha) * self.rf
        W_next = S * portfolio_return

        # 4. 演进 age 和 P
        next_age = age_t + 1
        z_shock = self.np_random.choice(self.shock_grid, p=self.shock_weights) * self.smav
        if next_age < self.tr:
            self.P = self.P + z_shock
        elif next_age == self.tr:
            self.P = self.P
            self.P_retirement = self.P
        
        self.age = next_age

        # 5. 实现 Y_{t+1}
        Y_next = 0.0
        next_age_idx = self.age - self.tb
        if self.age < self.tr:
            u_shock = self.np_random.choice(self.shock_grid, p=self.shock_weights) * self.smay
            Y_next = np.exp(self.f_y[next_age_idx] + self.P + u_shock)
            self.gyp = np.exp(self.f_y[next_age_idx] + self.P)
        else:
            if self.P_retirement == 0.0:
                 self.P_retirement = self.P
            Y_next = self.ret_fac * np.exp(self.f_y[self.tr - 1 - self.tb] + self.P_retirement)
            self.gyp = np.exp(self.f_y[self.tr - 1 - self.tb] + self.P_retirement)
        
        # 6. 更新 t+1 的状态
        self.W = W_next
        self.Y = Y_next
        self.X = self.W + self.Y
        if self.gyp > 1e-9:
            self.X_norm = self.X / self.gyp
        else:
            self.X_norm = self.X

        # 7. 检查终止条件
        died = False
        if self.age > self.td:
            died = True
        elif not self.disable_random_death:
            if current_age_idx < len(self.survprob) and self.np_random.uniform() > self.survprob[current_age_idx]:
                died = True
        self.is_done = died
        
        # 8. **[V3 修改]** 构造稀疏奖励
        final_reward = 0.0
        if self.is_done:
            final_reward = self.v # 在终点返回累计折现效用

        # 9. 构建 info 字典 (保持与 v2 相似，用于外部评估和调试)
        info = {
            "age": age_t,
            "absolute_wealth": self.W,
            "raw_permant_shock": self.P,
            "absolute_income": self.Y,
            "absolute_consumption": C,
            "c_prop": c_prop,
            "risky_share": alpha,
            "normalized_coh": self.X_norm,
            "period_utility": period_utility # 记录当期真实效用，便于分析
        }
        
        return self._get_obs(), float(final_reward), self.is_done, False, info
    # def step(self, action):
    #     if self.is_done:
    #         obs, info = self.reset()
    #         return obs, 0.0, True, False, info

    #     # 保存 t 期状态用于构建 info 字典
    #     age_t = self.age
    #     current_age_idx = age_t - self.tb
        
    #     # 获取 t 期的持久收入平减因子 P_t
    #     # 在我们的代码中，self.gyp 就是 t 期的 P_t 的确定性部分
    #     # 为了与理论一致，我们使用 self.gyp 作为 P_t
    #     P_t = self.gyp 

    #     # 0. 剪裁动作
    #     c_prop = np.clip(action[0], 0.0, 1.0)
    #     alpha = np.clip(action[1], 0.0, 1.0)

    #     # 1. 决策 (计算绝对消费和归一化消费)
    #     C_t = c_prop * self.X
    #     C_t = np.minimum(C_t, self.X * 0.99999)
    #     if C_t < 1e-6: C_t = 1e-6
        
    #     # 计算归一化消费 c_t
    #     c_t = C_t / P_t if P_t > 1e-9 else C_t

    #     c_prop = C_t / self.X if self.X > 1e-9 else 0.0
    #     S = self.X - C_t

    #     # 2. **[V4 核心修改]** 基于归一化消费计算效用，并应用修正因子
    #     if self.gamma == 1:
    #         # 对数效用: u(C) = log(C) = log(c * P) = log(c) + log(P)
    #         # 这里 u(c) = log(c)
    #         utility_from_normalized_c = np.log(c_t) if c_t > 1e-9 else -np.inf
    #         correction_factor = np.log(P_t) if P_t > 1e-9 else 0.0
    #         period_utility = utility_from_normalized_c + correction_factor
    #     else:
    #         # CRRA效用: u(C) = C^(1-g)/(1-g) = (c*P)^(1-g)/(1-g) = u(c) * P^(1-g)
    #         # 这里 u(c) = c^(1-g)/(1-g)
    #         utility_from_normalized_c = ((c_t)**(1 - self.gamma)) / (1 - self.gamma)
    #         correction_factor = P_t**(1 - self.gamma)
    #         period_utility = utility_from_normalized_c * correction_factor + 1

    #     # 更新增强状态 v (累计折现效用)
    #     self.v += self.w * period_utility
    #     self.w *= self.beta

    #     # 3. 计算 W_{t+1} (与之前一致)
    #     r_shock = self.np_random.choice(self.shock_grid, p=self.shock_weights) * self.sigr
    #     portfolio_return = alpha * (self.rf + self.mu + r_shock) + (1 - alpha) * self.rf
    #     W_next = S * portfolio_return

    #     # 4. 演进 age 和 P (与之前一致)
    #     next_age = age_t + 1
    #     z_shock = self.np_random.choice(self.shock_grid, p=self.shock_weights) * self.smav
    #     if next_age < self.tr:
    #         self.P = self.P + z_shock
    #     elif next_age == self.tr:
    #         self.P = self.P
    #         self.P_retirement = self.P
        
    #     self.age = next_age

    #     # 5. 实现 Y_{t+1} (与之前一致)
    #     Y_next = 0.0
    #     next_age_idx = self.age - self.tb
    #     if self.age < self.tr:
    #         u_shock = self.np_random.choice(self.shock_grid, p=self.shock_weights) * self.smay
    #         Y_next = np.exp(self.f_y[next_age_idx] + self.P + u_shock)
    #         self.gyp = np.exp(self.f_y[next_age_idx] + self.P)
    #     else:
    #         if self.P_retirement == 0.0:
    #             self.P_retirement = self.P
    #         Y_next = self.ret_fac * np.exp(self.f_y[self.tr - 1 - self.tb] + self.P_retirement)
    #         self.gyp = np.exp(self.f_y[self.tr - 1 - self.tb] + self.P_retirement)
        
    #     # 6. 更新 t+1 的状态 (与之前一致)
    #     self.W = W_next
    #     self.Y = Y_next
    #     self.X = self.W + self.Y
    #     if self.gyp > 1e-9:
    #         self.X_norm = self.X / self.gyp
    #     else:
    #         self.X_norm = self.X

    #     # 7. 检查终止条件 (与之前一致)
    #     died = False
    #     if self.age > self.td:
    #         died = True
    #     elif not self.disable_random_death:
    #         if current_age_idx < len(self.survprob) and self.np_random.uniform() > self.survprob[current_age_idx]:
    #             died = True
    #     self.is_done = died
        
    #     # 8. 构造稀疏奖励 (与之前一致)
    #     final_reward = 0.0
    #     if self.is_done:
    #         final_reward = self.v 

    #     # 9. 构建 info 字典 (与之前一致)
    #     info = {
    #         "age": age_t,
    #         "absolute_wealth": self.W,
    #         "raw_permant_shock": self.P,
    #         "absolute_income": self.Y,
    #         "absolute_consumption": C_t,
    #         "c_prop": c_prop,
    #         "risky_share": alpha,
    #         "normalized_coh": self.X_norm,
    #         "period_utility": period_utility 
    #     }
        
    #     return self._get_obs(), float(final_reward), self.is_done, False, info