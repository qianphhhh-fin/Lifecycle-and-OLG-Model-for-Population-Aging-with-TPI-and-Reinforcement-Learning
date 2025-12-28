import gymnasium as gym
from gymnasium import spaces
from gymnasium.envs.registration import register
import numpy as np
from typing import Optional
import warnings
warnings.filterwarnings('ignore')

register(
    id='cocco-v0',
    entry_point='cocco_env:CoccoEnv',
)

class CoccoEnv(gym.Env):
    def __init__(self, params=None):
        super(CoccoEnv, self).__init__()  # 调用父类(gym.Env)的初始化方法,确保环境正确初始化
        
        # 设置默认参数
        self.tb = 20  # 初始年龄
        self.tr = 65  # 退休年龄
        self.td = 100  # 最大年龄
        self.tn = self.td - self.tb + 1
        self.r = 1.015  # 无风险收益
        self.mu = 0.04  # 超额收益
        self.sigr = 0.5  # 风险资产收益率标准差
        self.riskaversion = 3  # CRRA风险规避系数
        self.delta = 0.97  # 贴现因子
        self.smay = 0.1  # 收入冲击标准差
        self.aa = -2.170042 + 2.700381
        self.b1, self.b2, self.b3 = 0.16818, -0.0323371/10, 0.0019704/100
        self.initial_cash = 2  # 初始现金
        
        # 如果传入参数则覆盖默认值
        if params:
            for key, value in params.items():
                setattr(self, key, value)
                
        # 动作空间: [消费比例, 风险资产比例]
        self.action_space = spaces.Box(
            low=np.array([0, 0]), 
            high=np.array([1, 1]),
            shape=(2,),
            dtype=float
        )
        
        # 状态空间: [现金, 年龄]
        self.observation_space = spaces.Box(
            low=np.array([0.25, -self.tb]),
            high=np.array([50.0, self.td]),
            shape=(2,),
            dtype=float
        )
        
        # 生存概率表
        self.survprob = [
        0.99975, 0.99974, 0.99973, 0.99972, 0.99971, 0.99969, 0.99968, 0.99966,
        0.99963, 0.99961, 0.99958, 0.99954, 0.99951, 0.99947, 0.99943, 0.99938,
        0.99934, 0.99928, 0.99922, 0.99916, 0.99908, 0.999, 0.99891, 0.99881,
        0.9987, 0.99857, 0.99843, 0.99828, 0.99812, 0.99794, 0.99776, 0.99756,
        0.99735, 0.99713, 0.9969, 0.99666, 0.9964, 0.99613, 0.99583, 0.99551,
        0.99515, 0.99476, 0.99432, 0.99383, 0.9933, 0.9927, 0.99205, 0.99133,
        0.99053, 0.98961, 0.98852, 0.98718, 0.98553, 0.98346, 0.98089, 0.97772,
        0.97391, 0.96943, 0.96429, 0.95854, 0.95221, 0.94537, 0.93805, 0.93027,
        0.92202, 0.91327, 0.90393, 0.89389, 0.88304, 0.87126, 0.85846, 0.84452,
        0.82935, 0.81282, 0.79485, 0.77543, 0.75458, 0.7324, 0.70893, 0.68424]
        
        # 初始化状态
        self.state = None
        self.done = False
        
    def step(self, action):
        # 解包动作
        consum_pct, risky_pct = action
        
        # 获取当前状态
        cash, age = self.state
        
        # 计算风险资产收益
        epsilon = np.random.normal(0, 1)
        rand_return = self.mu + self.r + epsilon * self.sigr
        
        # 计算下一期现金
        # 计算工资
        z_t = np.random.normal(0, self.smay)

        
        # 基础工资函数
        if age < self.tr:
            wage = np.clip(np.exp(self.aa + self.b1*age + self.b2*(age**2) + self.b3*(age**3))/100 + z_t, 0, np.inf)
        else:
            wage = self.final_wage*0.6 # 退休后的固定收入,设为退休前最后一年收入的60%

        if age == self.tr-1:
            # 退休前最后一年收入
            self.final_wage = wage
        
        # 计算下一期现金
        new_cash = cash * (1-consum_pct) * (risky_pct * (1 + rand_return) + 
                    (1-risky_pct) * (1 + self.r)) + wage
        # print(wage)
        # if new_cash<0:
        #     raise
        
        # 计算效用(奖励)
        if cash * consum_pct > 0:
            reward = ((cash * consum_pct)**(1-self.riskaversion))/(1-self.riskaversion)
        else:
            reward = -99
            
        # 更新年龄和状态
        new_age = age + 1
        self.state = np.array([new_cash, new_age])
        
        # 判断是否结束
        self.done = new_age >= self.td or \
            np.random.random() > self.survprob[int(age-self.tb)]
            
        info = {
            'age': age,
            'wage': wage,
            'consumption': cash * consum_pct
        }
        
        return self.state, reward, self.done, False, info
        
    def reset(self, *, seed: Optional[int] = None, options: Optional[dict] = None):
        super().reset(seed=seed)
        
        # 初始状态
        init_age = self.tb  # 初始年龄
        
        self.state = np.array([self.initial_cash, init_age])
        self.done = False
        
        info = {
            'age': init_age,
            'initial_cash': self.initial_cash,
        }
        
        return self.state, info


if __name__ == '__main__':
    # 测试环境
     # 定义参数字典
    params = {
        'tb': 20, # 初始年龄
        'td': 100, # 最大年龄
        'tr': 65, # 退休年龄
        'r': 1.015, # 无风险收益
        'mu': 0.04, # 超额收益
        'sigr': 0.5, # 风险资产收益率标准差
        'riskaversion': 3, # CRRA风险规避系数
        'delta': 0.97, # 贴现因子
        'smay': 0.1, # 收入冲击标准差
        'initial_cash': 2, # 初始现金
    }
    
    env = gym.make('cocco-v0',params=params)
    # env = gym.wrappers.FlattenObservation(env)  # deal with dm_control's Dict observation space
    env = gym.wrappers.RecordEpisodeStatistics(env)
    
    # 重置环境
    obs, info = env.reset()
    print("初始状态:", obs)
    print("初始信息:", info)
    
    # 模拟100个agent的路径
    n_agents = 1
    
    # 存储所有agent的轨迹
    all_trajectories = []
    
    for agent in range(n_agents):
        # 重置环境获取初始状态
        obs, info = env.reset()
        done = False
        agent_trajectory = {
            'states': [obs],
            'actions': [],
            'rewards': [],
            'infos': [info],
            'consumptions': []
        }
        
        while not done:  
            # 随机动作 (消费比例, 风险资产比例)
            action = env.action_space.sample()
            
            # 执行动作
            obs, reward, done, truncated, info = env.step(action)
            
            # 记录轨迹
            agent_trajectory['states'].append(obs)
            # print(obs)
            agent_trajectory['actions'].append(action)
            agent_trajectory['rewards'].append(reward)
            agent_trajectory['consumptions'].append(info['consumption'])
            agent_trajectory['infos'].append(info)
            
            if done:
                break
                
        all_trajectories.append(agent_trajectory)
        
        # print(f"完成Agent {agent+1:1d}: "
        #       f"年龄为{int(info['age']):3d}, "
        #       f"现金为{obs[0]:8.2f}, "
        #       f"平均工资为{info['wage']:8.2f}, "
        #       f"消费为{info['consumption']:8.2f}")
    # 计算所有agent的reward统计信息
    all_rewards = [sum(traj['rewards']) for traj in all_trajectories]
    max_reward = max(all_rewards)
    min_reward = min(all_rewards)
    avg_reward = sum(all_rewards) / len(all_rewards)
    
    print(f'最大累积奖励: {max_reward:.2f}')
    print(f'最小累积奖励: {min_reward:.2f}') 
    print(f'平均累积奖励: {avg_reward:.2f}')
    
    # 计算所有agent的平均消费
    all_consumptions = [sum(traj['consumptions']) for traj in all_trajectories]
    avg_consumption = sum(all_consumptions) / len(all_consumptions)
    print(f'平均消费: {avg_consumption:.2f}')
    print(f'最大消费: {max(all_consumptions):.2f}')
    print(f'最小消费: {min(all_consumptions):.2f}')
    
    print("\n所有Agent模拟完成!")
            
    env.close()

