# -*- coding: utf-8 -*-
"""
基于gymnasium标准格式的env，随机死亡，累进税率，个人养老金购买上限
连续选择：消费比例，风险资产比例，购买个人养老金比例，个人养老金中购买风险资产比例
状态变量：财富，上期工资持续冲击，个人养老金账户余额，年龄
注：cash和pen_balance为万元/100，perm_shock为万元/10，norm_age为-1到1
action space为-1到2(在env中进行clip)
"""

import gymnasium as gym
from gymnasium import spaces
from gymnasium.envs.registration import register
# from gymnasium.utils.env_checker import check_env
import numpy as np
from typing import Optional
import warnings
warnings.filterwarnings('ignore')
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

# Register this module as a gym environment. Once registered, the id is usable in gym.make().
register(
    id='cocco-v0',                                # call it whatever you want
    entry_point='v0_cocco:CoccoEnv', # module_name:class_name
)


'''
v0: cocco基准模型
'''

#         #按年龄收入调整系数 2023人均年工资性收入2.21万元 
#         # (此参数平均收入:21岁1.766万元，49岁3.592万元，64岁3.329万元)
#         # desmos: \frac{e^{0.5+0.172x-0.00323371x^{2}+0.000019x^{3}}}{100}


class CoccoEnv(gym.Env):
    """Custom Environment that follows gym interface"""
    # metadata = {'render.modes': ['human']}
    
    def __init__(self,params=None): # 
        
        super(CoccoEnv, self).__init__()
                
         # 设置默认参数
        self.scale_inc = 2
        self.subtract_income = 1 # 相当于将以万为单位的收入除以np.exp(self.subtract_income), 若为0，则收入单位为万元
        self.scale_reward = 1 # 小于1有助于增强探索
        self.reward_lb = -np.inf # 单期reward的下限

        
        self.tb = 20 # 初始年龄
        self.tr = 61 # 退休年龄
        self.td = 100  # 最大年龄
        self.tn = self.td - self.tb + 1
        self.rf = 0.02  # 无风险收益
        self.mu = 0.04  # 超额收益
        self.sigr = 0.27  # 风险资产收益率标准差
        self.gamma = 3.84 # 3.84  # 风险厌恶系数
        self.beta = 0.95  # 贴现因子
        self.smay = np.sqrt(0.01)  # 暂时性shock的标准差
        self.smav = np.sqrt(0.01)  # 永久性shock的标准差
        self.aa = (-2.170042 + 2.700381)
        self.b1, self.b2, self.b3 = 0.16818, -0.0323371/10, 0.0019704/100
        # self.aa,self.b1, self.b2, self.b3 = 0.5,0.172,-0.00323371,0.000019
        self.ret_fac = 0.6827  # 退休后固定支付的基本养老金,为退休前工资的68.27%
        self.init_cash = 0  # 初始现金
        

        # 如果传入参数则覆盖默认值
        if params:
            for key, value in params.items():
                setattr(self, key, value)
        
        
        
        # 20-100岁的生存概率
        survprob_values = [
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
        self.survprob_dict = {age: value for age, value in zip(range(20, 101), survprob_values)}
        # self.survprob_dict = {age: 1 for age, value in zip(range(20, 101), survprob_values)} # 无随机死亡
    
        # ===================== 行动空间：=====================
        # 连续选择：消费比例（不允许小于零），风险资产比例
        self.action_space = spaces.Box(low=np.array([0.1, 0]),\
                                       high=np.array([1, 1]), shape=(2,), dtype=float)
        
        # ===================== 状态空间：=====================
        # 现金余额，永久收入冲击, 年龄
        self.observation_space = spaces.Box(
            low=np.array([0, -1e3, -0.5]),
            high=np.array([1e3,  1e3, 0.5]),
            shape=(3,),
            dtype=float
        )
        
        self.state = None
        



    def step(self, action):
        
        # if self.done:
        #     self.state,_ = self.reset()
        # print(self.pension_limit)
        # raise
        # Unpack the action
        consum_pct = np.clip(action[0],0,1)  # 消费比例
        risky_pct = np.clip(action[1],0,1) # 风险资产比例 


        # if any(np.isnan(action)):
        #     raise

        # Get current state values
        # cash = self.state['cash']  
        # previous_wage = self.wage # 上期工资
        # perm_shock = self.state['perm_shock']
        # pen_balance = self.state['pension_balance']
        # norm_age = self.state['norm_age']      
        cash,perm_shock,_ = self.state
        
        # age = self.state[3]


        #  === 风险资产随机性 ===
        # 随机项
        epsilon = np.random.normal(0,1)       
        # epsilon = np.random.randn(1000)  
        rand_return = self.mu + self.rf + epsilon * self.sigr # 风险资产收益率

        new_cash = 0 
# === 状态转移方程（退休前） ===       
        if self.age < self.tr:

            # 工资
            u_t = np.random.normal(0, self.smay)
            z_t = np.random.normal(0, self.smav)
            new_perm_shock = perm_shock + z_t # 工资持续性冲击
            # 如果下一期进入退休期,保存最后一期工作期的工资冲击   
            if self.age+1 ==self.tr:
                self.last_perm_shock = new_perm_shock
                self.last_f = self.f(self.age+1)
            # log_rand_wage = self.f(norm_age) + new_perm_shock + u_t   
            log_rand_wage = self.f(self.age+1) + new_perm_shock + u_t   
            new_wage = np.exp(log_rand_wage) # 下一期工资

            # 下一期到手现金 这里的处理与Gomes et al. (2009)相同, 
            new_cash = cash * (1-consum_pct) * (risky_pct * (1 + rand_return) \
                                                  + (1-risky_pct) * (1 + self.rf)) + new_wage
                  
            self.info = {'age':np.array([self.age]),'basic_income':self.wage*np.exp(self.subtract_income),\
                    'real_actions':np.array([consum_pct,risky_pct]), # 记录实际的action
                    'normalized_cash':(new_cash)/np.exp(self.f(self.age+1) + new_perm_shock),\
                    'real_cash':new_cash*np.exp(self.subtract_income),\
                    'normalized_income_working':(new_wage)/np.exp(self.f(self.age+1) + new_perm_shock),\
                    'state':self.state,'status':'working'}
                
# === 状态转移方程（退休后） ===   
        else:  
            new_perm_shock = self.last_perm_shock # 工资持续性冲击为最后一期工作期的工资冲击保持不变
            log_rand_wage = self.last_perm_shock + self.last_f + np.log(self.ret_fac)   # ret_fac是退休后固定支付的工资,为最后一期工作期工资的比例
            new_wage = np.exp(log_rand_wage) # 下一期工资
            # 当期到手现金
            new_cash = cash * (1-consum_pct) * (risky_pct * (1 + rand_return) \
                                                  + (1-risky_pct) * (1 + self.rf)) \
                + new_wage
                    
            # print([new_cash,self.age,'retired'])
            self.info = {'age':np.array([self.age]),'basic_income':self.wage*np.exp(self.subtract_income),\
                    'real_actions':np.array([consum_pct,risky_pct]), # 记录实际的action
                    'normalized_cash':(new_cash)/np.exp(self.last_perm_shock + self.last_f),\
                    'real_cash':new_cash*np.exp(self.subtract_income),\
                    'normalized_income_retired':(new_wage)/np.exp(self.last_perm_shock + self.last_f),\
                    'state':self.state,'status':'retired'}

        reward = -(consum_pct*cash) ** (1 - self.gamma)
        # 奖励一定要设置合理的下bound, 否则训练不稳定     
        reward = np.maximum(reward, self.reward_lb) 

        reward *= self.scale_reward # 缩放reward
        # Bounded below 
        # if reward<-150:
        #     reward = -150   
        # info.update({'消费':norm_consum})
        self.info.update({'实际消费':cash * consum_pct *np.exp(self.subtract_income)})    
        # print(norm_consum)

        self.age += 1
        self.wage = new_wage    
        # Determine if the episode is done
        # 如果完成，则返回同样的state，在计算时game结束的下一期的值函数无论如何都为0（不影响
        # 无法存活至下一期or超过最大年龄       
        if (self.age > self.td):
            self.done = True
            self.info.update({'PeriodEnd':'最大年龄死亡'})
            # print(self.state['pension_balance'],'万元')
            return self.state, reward, self.done, False,  self.info
        
        # 随机存活
        if  np.random.choice([True,False], 1,  p=[1-self.survprob_dict[self.age-1],self.survprob_dict[self.age-1]]): 
            self.done = True       
            self.info.update({'PeriodEnd':'提前死亡'})
            # print(self.state['pension_balance'],'万元')
            return self.state, reward, self.done, False, self.info
            # print([self.age])
            
        # else:
        self.done = False
        
        # 下一期状态
        self.state[0] = new_cash  
        self.state[1] = new_perm_shock
        self.state[2] = -0.5 + (self.age - self.tb)/self.tn

        return self.state, reward, self.done, False, self.info

    def reset(self, *, seed: Optional[int] = None, options: Optional[dict] = None):
        super().reset(seed=seed)
        np.random.seed(seed)
        
        self.age = self.tb # 起始年龄
        self.index = 0

        # 第一份随机工资收入
        u_t = np.random.normal(0, self.smay)
        z_t = np.random.normal(0, self.smav)
        perm_shock = z_t # 工资持续性冲击
        
        log_rand_wage = self.f(self.age) + perm_shock  + u_t  
        rand_wage = np.exp(log_rand_wage)
        self.wage = rand_wage
        new_cash = self.init_cash*np.exp(self.f(self.age) + perm_shock) # raw 
        
        # Take a sample from the observation_space to modify the values of
        self.state = self.observation_space.sample()

        self.state[0] = np.array([new_cash + self.wage],dtype=float)
        self.state[1] = np.array([perm_shock],dtype=float)
        self.state[2] = np.array([-0.5],dtype=float)
        
        # self.state = np.array([rand_wage,perm_shock,-1],dtype=np.float32)
        self.info = {'age':self.age,'basic_income':self.wage,\
                          'state':self.state,
                          'normalized_cash':(new_cash + self.wage)/np.exp(self.f(self.age) + perm_shock),
                          'normalized_income':self.wage/np.exp(self.f(self.age) + perm_shock),
                          'status':'wokring'}
        # if np.isnan(self.state).any(): # 检测是否有nan，返回的state确保没有nan，以免弄坏vec_env的rms
        #     self.state = np.nan_to_num(self.state, 0)
            
        return self.state, self.info

    def render(self, mode='human', close=False):
         # Render the environment to the screen
         print(f"State: {self.state}")
    
    # Define the function f(age)
    def f(self, age):
        wage = self.aa +self. b1 * age + self.b2 * age**2 + self.b3 * age**3
        return wage/self.scale_inc - self.subtract_income # 将收入缩放至效用函数斜率较大的范围内
        
    

            
if __name__=="__main__":

    
    env = gym.make('cocco-v0')
    env = gym.wrappers.RecordEpisodeStatistics(env)
    
    # 模拟100个agent的路径
    n_agents = 10000
    
    # 存储所有agent的轨迹
    all_trajectories = []
    
    # 创建进度条
    from tqdm import tqdm
    pbar = tqdm(total=n_agents, desc='模拟进度')
    for agent in range(n_agents):
        # 重置环境获取初始状态
        obs, info = env.reset()
        done = False
        agent_trajectory = {
            'normalized_income_working':[],
            'normalized_income_retired':[],
            'states': [obs],
            'actions': [],
            'rewards': [],
            'infos': [info],
            'consumption': [],
            'norm_age':[obs[-1]]
        }
        while not done:
            # 固定动作策略
            action = np.clip(env.action_space.sample(),-1,2)
            
            # 执行动作
            obs, reward, done, _, info = env.step(action)
            
            # 记录轨迹
            try:
                agent_trajectory['normalized_income_working'].append(info['normalized_income_working'])
            except:
                agent_trajectory['normalized_income_retired'].append(info['normalized_income_retired'])
            agent_trajectory['states'].append(obs)
            agent_trajectory['actions'].append(action)
            agent_trajectory['rewards'].append(reward)
            agent_trajectory['consumption'].append(info['消费'])
            agent_trajectory['infos'].append(info)
            agent_trajectory['norm_age'].append(obs[-1])
            
        all_trajectories.append(agent_trajectory)
        pbar.update(1)
    pbar.close()
        
    # 计算所有agent的reward统计信息
    all_income_working = np.concatenate([traj['normalized_income_working'] for traj in all_trajectories])
    all_income_retired = np.concatenate([traj['normalized_income_retired'] for traj in all_trajectories])
    max_income_working = max(all_income_working)
    min_income_working = min(all_income_working)
    avg_income_working = sum(all_income_working) / len(all_income_working)
    std_income_working = np.std(all_income_working)
    max_income_retired = max(all_income_retired)
    min_income_retired = min(all_income_retired)
    avg_income_retired = sum(all_income_retired) / len(all_income_retired)
    std_income_retired = np.std(all_income_retired)

    print(f'最大工作单期收入: {max_income_working:.2f} 万元')
    print(f'最小工作单期收入: {min_income_working:.2f} 万元')
    print(f'平均工作单期收入: {avg_income_working:.2f} 万元')
    print(f'工作单期收入std: {std_income_working:.2f} 万元')
    
    print(f'最大退休单期收入: {max_income_retired:.2f} 万元')
    print(f'最小退休单期收入: {min_income_retired:.2f} 万元')
    print(f'平均退休单期收入: {avg_income_retired:.2f} 万元')
    print(f'退休单期收入std: {std_income_retired:.2f} 万元')
    



    # all_age = np.concatenate([traj['norm_age'] for traj in all_trajectories])
    # print(np.unique(all_age))
    # print(np.shape(np.unique(all_age)))

    
    env.close()
    # import matplotlib.pyplot as plt
    # plt.plot(score_history)
    
    
    
            
