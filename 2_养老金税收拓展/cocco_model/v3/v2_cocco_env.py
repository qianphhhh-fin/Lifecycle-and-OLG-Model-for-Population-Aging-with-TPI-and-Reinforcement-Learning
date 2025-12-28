import gymnasium as gym
from gymnasium import spaces
import numpy as np
from typing import Optional
from scipy.stats import norm

class CoccoEnvV2(gym.Env):
    """
    Cocco et al. (2005) Life-Cycle Model Environment - V2.
    **[修正版 v2]** 此版本修复了 step 函数返回的 info 字典，以支持评估回调。
    - 在环境内部显式追踪 W_t 和 Y_t，以确保 info 字典可以报告当期的准确值。
    """
    metadata = {'render.modes': []}

    def __init__(self, disable_random_death: bool = False):
        super(CoccoEnvV2, self).__init__()

        self.disable_random_death = disable_random_death

        # --- A. 生命周期与人口结构 (与 CoccoEnv.m 对齐) ---
        self.tb = 22
        self.tr = 26
        self.td = 30
        self.tn = self.td - self.tb + 1

        # --- B. 经济主体偏好 (与 CoccoEnv.m 对齐) ---
        self.beta = 0.95
        self.gamma = 3.84

        # --- C. 收入过程参数 (与 CoccoEnv.m 对齐) ---
        self.aa = -2.170042 + 2.700381
        self.b1 = 0.16818
        self.b2 = -0.0323371 / 10
        self.b3 = 0.0019704 / 100
        self.smay = 0.169993**0.5
        self.smav = 0.112572**0.5

        # --- D. 资产与回报率 (与 CoccoEnv.m 对齐) ---
        self.rf = 1.02
        self.mu = 0.04
        self.sigr = 0.27

        # --- E. 养老金与税收 (与 CoccoEnv.m 对齐) ---
        self.ret_fac = 0.6827
        self.tau_y = 0.0

        # --- F. 生存概率 (从 CoccoEnv.m 复制) ---
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

        # --- G. 派生参数与预计算 ---
        self.f_y = np.zeros(self.tn)
        for i in range(self.tn):
            age = self.tb + i
            if age <= self.tr:
                self.f_y[i] = (self.aa + self.b1 * age + self.b2 * age**2 + self.b3 * age**3) - 3

        n_shocks = 5
        self.shock_grid, self.shock_weights = self._discretize_shocks(n_shocks)

        # --- H. 状态变量 ---
        # **核心修改**: 显式追踪 W 和 Y
        self.W = 0.0  # 当期期初财富
        self.Y = 0.0  # 当期收入
        self.X = 0.0  # 当期现金持有量 (W+Y)
        self.P = 0.0
        self.P_retirement = 0.0
        self.age = self.tb
        self.is_done = False

        # --- I. Gym 空间定义 ---
        self.action_space = spaces.Box(
            low=np.array([0.2, 0.05]),
            high=np.array([1.1, 1.1]),
            shape=(2,),
            dtype=np.float32
        )
        self.observation_space = spaces.Box(
            low=np.array([0.0, -np.inf,-1.0]),
            high=np.array([np.inf, np.inf,1.0]),
            shape=(3,),
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
        normalized_age = 2 * (self.age - self.tb) / (self.td - self.tb) - 1
        return np.array([self.X_norm, self.P, normalized_age], dtype=np.float32)

    def reset(self, *, seed: Optional[int] = None, options: Optional[dict] = None):
        super().reset(seed=seed)
        
        self.age = self.tb
        self.is_done = False
        self.P_retirement = 0.0
        
        z_shock = self.np_random.choice(self.shock_grid, p=self.shock_weights) * self.smav
        self.P = z_shock

        # **核心修改**: 初始化 W_1 和 Y_1
        self.W = 0.0
        u_shock = self.np_random.choice(self.shock_grid, p=self.shock_weights) * self.smay
        self.Y = np.exp(self.f_y[0] + self.P + u_shock)
        
        self.X = self.W + self.Y
        self.gyp = np.exp(self.f_y[0] + self.P)
        self.X_norm = self.X/self.gyp

        # reset 返回的 info 结构与 step 保持一致，尽管在 SB3 中通常被忽略
        info = {
            "age": self.age,
            "absolute_wealth": self.W,
            "raw_permant_shock": self.P,
            "absolute_consumption": 0.0,
            "absolute_income": self.Y,
            "risky_share": 0.0,
            "normalized_coh": self.X_norm,
        }
        return self._get_obs(), info

    def step(self, action):
        if self.is_done:
            obs, info = self.reset()
            return obs, 0.0, True, False, info

        # **核心修改**: 保存 t 期状态用于构建 info 字典
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
        c_prop = C/self.X
        if C < 1e-6: C = 1e-6
        S = self.X - C

        # 2. 计算奖励
        if self.gamma == 1: utility = np.log(C)
        else: utility = ((C)**(1 - self.gamma)) / (1 - self.gamma) + 1
        reward = np.maximum(utility, -20)

        # 3. 计算 W_{t+1}
        r_shock = self.np_random.choice(self.shock_grid, p=self.shock_weights) * self.sigr
        portfolio_return = alpha * (self.rf + self.mu + r_shock) + (1 - alpha) * self.rf
        W_next = S * portfolio_return

        # 4. 演进 age 和 P
        next_age = age_t + 1
        z_shock = self.np_random.choice(self.shock_grid, p=self.shock_weights) * self.smav
        # [MODIFIED] 退休年龄逻辑与 MATLAB 对齐
        if next_age < self.tr:
            self.P = self.P + z_shock
        elif next_age == self.tr: # 刚刚进入退休年龄
            self.P = self.P
            self.P_retirement = self.P
        
        self.age = next_age

        # 5. 实现 Y_{t+1}
        Y_next = 0.0
        next_age_idx = self.age - self.tb
        # [MODIFIED] 退休年龄逻辑与 MATLAB 对齐
        if self.age < self.tr:
            u_shock = self.np_random.choice(self.shock_grid, p=self.shock_weights) * self.smay
            Y_next = np.exp(self.f_y[next_age_idx] + self.P + u_shock)
            self.gyp = np.exp(self.f_y[next_age_idx] + self.P)
        else:
            if self.P_retirement == 0.0: # 确保P_retirement在第一次退休时被设置
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
            self.X_norm = self.X # 避免除以零

        # 7. 检查终止条件
        died = False
        if self.age > self.td:
            died = True
        elif not self.disable_random_death:
            if current_age_idx < len(self.survprob) and self.np_random.uniform() > self.survprob[current_age_idx]:
                died = True
        self.is_done = died
        
        # 8. **核心修改**: 构建包含 t 期完整信息的 info 字典
        info = {
            "age": age_t,
            "absolute_wealth": self.W,
            "raw_permant_shock": self.P,
            "absolute_income": self.Y,
            "absolute_consumption": C,
            "c_prop": c_prop,
            "risky_share": alpha,
            "normalized_coh": self.X_norm
        }
        
        return self._get_obs(), float(reward), self.is_done, False, info


    # def step(self, action):
    #     if self.is_done:
    #         obs, info = self.reset()
    #         return obs, 0.0, True, False, info

    #     # **核心修改**: 保存 t 期状态用于构建 info 字典
    #     W_t = self.W
    #     Y_t = self.Y
    #     age_t = self.age
    #     current_age_idx = age_t - self.tb

    #     # 0. 剪裁动作
    #     c_prop = np.clip(action[0], 0.0, 1.0)
    #     alpha = np.clip(action[1], 0.0, 1.0)

    #     # 1. 决策
    #     C = c_prop * self.X
    #     C = np.minimum(C, self.X * 0.99999)
    #     c_prop = C/self.X
    #     if C < 1e-6: C = 1e-6
    #     S = self.X - C

    #     # 2. 计算奖励
    #     if self.gamma == 1: utility = np.log(C)
    #     else: utility = ((C)**(1 - self.gamma)) / (1 - self.gamma) + 1
    #     reward = np.maximum(utility, -20.0)

    #     # 3. 计算 W_{t+1}
    #     r_shock = self.np_random.choice(self.shock_grid, p=self.shock_weights) * self.sigr
    #     portfolio_return = alpha * (self.rf + self.mu + r_shock) + (1 - alpha) * self.rf
    #     W_next = S * portfolio_return

    #     # 4. 演进 age 和 P
    #     next_age = age_t + 1
    #     z_shock = self.np_random.choice(self.shock_grid, p=self.shock_weights) * self.smav
    #     if next_age <= self.tr:
    #         self.P = self.P + z_shock
    #         if next_age == self.tr:
    #             self.P_retirement = self.P
        
    #     self.age = next_age

    #     # 5. 实现 Y_{t+1}
    #     Y_next = 0.0
    #     next_age_idx = self.age - self.tb
    #     if self.age <= self.tr:
    #         u_shock = self.np_random.choice(self.shock_grid, p=self.shock_weights) * self.smay
    #         Y_next = np.exp(self.f_y[next_age_idx] + self.P + u_shock)
    #         self.gyp = np.exp(self.f_y[next_age_idx] + self.P)
    #     else:
    #         Y_next = self.ret_fac * np.exp(self.f_y[self.tr - self.tb] + self.P_retirement)
    #         self.gyp = np.exp(self.f_y[self.tr - self.tb] + self.P_retirement)
        
    #     # 6. 更新 t+1 的状态
    #     self.W = W_next
    #     self.Y = Y_next
    #     self.X = self.W + self.Y
    #     self.X_norm = self.X/self.gyp
        
    #     # 7. 检查终止条件
    #     died = False
    #     if self.age > self.td:
    #         died = True
    #     elif not self.disable_random_death:
    #         if current_age_idx < len(self.survprob) and self.np_random.uniform() > self.survprob[current_age_idx]:
    #             died = True
    #     self.is_done = died
        
    #     # 8. **核心修改**: 构建包含 t 期完整信息的 info 字典
    #     info = {
    #         "age": age_t,
    #         "absolute_wealth": self.W,
    #         "raw_permant_shock": self.P,
    #         "absolute_income": self.Y,
    #         "absolute_consumption": C,
    #         "c_prop": c_prop,
    #         "risky_share": alpha,
    #         "normalized_coh": self.X_norm
    #     }
        
    #     return self._get_obs(), float(reward), self.is_done, False, info