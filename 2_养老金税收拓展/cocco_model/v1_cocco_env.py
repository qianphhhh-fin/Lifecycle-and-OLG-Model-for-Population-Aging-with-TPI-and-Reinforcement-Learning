import gymnasium as gym
from gymnasium import spaces
import numpy as np
from typing import Optional

class CoccoEnv(gym.Env):
    """
    Cocco et al. (2005) Life-Cycle Model Environment for Reinforcement Learning.

    This environment implements economic normalization internally.
    - Observations (to agent): [normalized_wealth, normalized_age]
    - Actions (from agent): [normalized_consumption, risky_share]
    - Internal State: Tracks absolute wealth and absolute permanent income.
    """
    metadata = {'render.modes': []}

    def __init__(self, disable_random_death: bool = False): # 增加一个参数
        super(CoccoEnv, self).__init__()
        
        self.disable_random_death = disable_random_death # 存储这个参数

        # --- A. 生命周期与人口结构 ---
        self.tb = 22       # 初始工作年龄
        self.tr = 61       # 退休年龄
        self.td = 100      # 最高寿命
        self.tn = self.td - self.tb + 1 # 总期数

        # --- B. 经济主体偏好 (CRRA) ---
        self.beta = 0.95     # 折现因子
        self.gamma = 3.84    # 相对风险规避系数

        # --- C. 收入过程 ---
        # 确定性部分 g(t) = exp(aa + b1*t + b2*t^2 + b3*t^3)
        self.aa = -2.170042 + 2.700381
        self.b1 = 0.16818
        self.b2 = -0.0323371 / 10
        self.b3 = 0.0019704 / 100
        # 随机部分
        self.smay = np.sqrt(0.169993) # 暂时冲击标准差 (u_t)
        self.smav = np.sqrt(0.112572) # 持久冲击标准差 (z_t)

        # --- D. 资产与回报率 ---
        self.rf = 1.02     # 无风险总回报率
        self.mu = 0.04     # 风险资产超额回报率
        self.sigr = 0.27   # 风险资产回报率标准差

        # --- E. 养老金 ---
        self.ret_fac = 0.6827      # 退休后基础养老金 (替代率)

        # --- G. 派生参数与预计算 ---
        # 生存概率 (p_t 是从 age t 存活到 age t+1 的概率)
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

        # 确定性收入剖面
        self.f_y = np.zeros(self.tn)
        for i in range(self.tn):
            age = self.tb + i
            if age <= self.tr:
                 self.f_y[i] = np.exp(self.aa + self.b1*age + self.b2*age**2 + self.b3*age**3)

        # 状态变量和内部变量
        self.W = 0.0
        self.P = 1.0
        self.P_retirement = 1.0
        self.age = self.tb
        self.is_done = False

        # --- 定义动作和状态空间 ---
        # ** 修改: 消费动作现在是占总资源的比例 **
        # Action: [consumption_proportion, risky_share]
        self.action_space = spaces.Box(
            low=np.array([0.0001, 0.0]),
            high=np.array([1.0, 1.0]),
            shape=(2,),
            dtype=np.float32
        )

        # Observation: [normalized_wealth, normalized_age]
        self.observation_space = spaces.Box(
            low=np.array([0.0, -1.0]),
            high=np.array([1000.0, 1.0]), # 保持一个足够大的上界
            shape=(2,),
            dtype=np.float32
        )

    def _get_obs(self):
        # 年龄归一化到 [-1, 1]
        normalized_age = 2 * (self.age - self.tb) / (self.td - self.tb) - 1
        normalized_wealth = self.W / self.P
        return np.array([normalized_wealth, normalized_age], dtype=np.float32)

    def _get_info(self, C, Y, alpha):
        return {
            "age": self.age,
            "absolute_wealth": self.W,
            "absolute_consumption": C,
            "absolute_income": Y,
            "risky_share": alpha,
            "permanent_income": self.P
        }

    def reset(self, *, seed: Optional[int] = None, options: Optional[dict] = None):
        super().reset(seed=seed)

        self.age = self.tb
        self.is_done = False
        self.W = 0.0

        # 初始化永久收入
        z_shock = np.random.normal(0, self.smav)
        self.P = np.exp(z_shock)

        info = self._get_info(0, 0, 0) # Placeholder info for reset
        return self._get_obs(), info

# v1_cocco_env.py

    def step(self, action):
        if self.is_done:
            obs, info = self.reset()
            # 当一个VecEnv中的某个环境结束后，SB3会期望返回一个有效的 (obs, reward, terminated, truncated, info) 元组
            # 在这种情况下，我们返回重置后的观测，奖励为0，并将done标志设为True。
            return obs, 0.0, True, False, info

        # 1. 解包并约束动作 (比例)
        c_prop, alpha = action
        # 依然保留此关键修正，从根源上防止 C=0
        c_prop = np.clip(c_prop, 0.05, 1.0)
        alpha = np.clip(alpha, 0.0, 1.0)

        # 2. 计算当期绝对资源
        Y = 0.0
        current_age = self.age
        current_age_idx = current_age - self.tb
        
        if current_age < self.tr: # 工作期
            u_shock = np.random.normal(0, self.smay)
            det_Y = self.f_y[current_age_idx]
            Y_gross = det_Y * self.P * np.exp(u_shock)
            
            # 加入工资税，与VFI模型匹配
            tau_y = 0.06 
            Y = Y_gross * (1 - tau_y)

        else: # 退休期
            Y = self.ret_fac * self.P_retirement

        X = self.W + Y # 绝对现金持有量

        # 3. 将消费比例动作转换为绝对消费 C
        C = c_prop * X
        # 确保消费不会超过可用资源
        C = min(C, X)
        
        # --- 核心修改: 使用归一化消费计算奖励 ---
        # 4. 计算归一化消费 c = C / P
        # 在退休期, P 保持为 P_retirement，此逻辑依然有效
        # 增加一个极小值防止 self.P 为0（虽然理论上不会）
        c_normalized = C / (self.P + 1e-8)

        # 5. 基于归一化消费 c 计算效用和奖励
        if self.gamma == 1:
            # 确保 c_normalized > 0, c_prop的裁剪已经保证了这一点
            utility = np.log(c_normalized)
        else:
            utility = (c_normalized**(1 - self.gamma)) / (1 - self.gamma)

        # 奖励整形策略
        reward_scale =  0.00411
        reward_shift =  0.00002
        scaled_reward = utility * reward_scale + reward_shift
        reward_lower_bound = -99999
        reward = max(scaled_reward, reward_lower_bound)

        # 6. 模拟世界演化 (此部分必须使用绝对值)
        S = X - C # 绝对储蓄

        # 投资组合回报
        r_shock = np.random.normal(0, self.sigr)
        risky_return = self.rf + self.mu + r_shock
        risky_return = max(0, risky_return) # 风险资产投资收益率不可能小于0
        portfolio_return = alpha * risky_return + (1 - alpha) * self.rf

        # 更新下一期绝对财富
        self.W = S * portfolio_return

        # 7. 更新下一期永久收入 (P)
        if current_age < self.tr:
            z_shock = np.random.normal(0, self.smav)
            
            if current_age < self.tr - 1:
                g_growth = self.f_y[current_age_idx + 1] / self.f_y[current_age_idx]
                self.P = self.P * g_growth * np.exp(z_shock)
            elif current_age == self.tr - 1: # 工作期的最后一年
                self.P = self.P * np.exp(z_shock)
                self.P_retirement = self.P
        # 退休后 P 保持不变

        self.age += 1

        # 8. 随机死亡检查
        died = False
        if self.age > self.td:
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
        info = self._get_info(C, Y, alpha)

        return self._get_obs(), float(reward), self.is_done, False, info