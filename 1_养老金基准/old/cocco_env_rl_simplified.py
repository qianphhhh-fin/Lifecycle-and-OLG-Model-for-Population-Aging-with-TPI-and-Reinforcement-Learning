# -*- coding: utf-8 -*-
"""
基于gymnasium标准格式的env，随机死亡，累进税率，个人养老金购买上限
连续选择：消费比例，风险资产比例，购买个人养老金比例，个人养老金中购买风险资产比例
状态变量：可用于购买个人养老金的现金，剩余现金，上期工资持续冲击，个人养老金账户余额，年龄
注：没有normalized状态

"""

import gymnasium as gym
from gymnasium import spaces
from gymnasium.envs.registration import register
# from gymnasium.utils.env_checker import check_env
from stable_baselines3.common.env_checker import check_env
import matplotlib.pyplot as plt
import numpy as np
from typing import Optional
import warnings
warnings.filterwarnings('ignore')
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

# Register this module as a gym environment. Once registered, the id is usable in gym.make().
register(
    id='cocco-rl',
    entry_point='cocco_env_rl_simplified:CoccoEnv',
)
  


class CoccoEnv(gym.Env):
    """Custom Environment that follows gym interface"""
    # metadata = {'render.modes': ['human']}
    
    def __init__(self, params): # 
        
        super(CoccoEnv, self).__init__()
        
         # 设置默认参数
        self.tb = 96 # 96#20  # 初始年龄
        self.tr = 98 # 98#61  # 退休年龄
        self.td = 100  # 最大年龄
        self.tn = self.td - self.tb + 1
        self.rf = 0.02  # 无风险收益
        self.mu = 0.04  # 超额收益
        self.sigr = 0.27  # 风险资产收益率标准差
        self.gamma = 3.84  # 风险厌恶系数
        self.beta = 0.95  # 贴现因子
        self.smay = np.sqrt(0.01)  # 暂时性shock的标准差
        self.smav = np.sqrt(0.01)  # 永久性shock的标准差
        self.aa = (-2.170042 + 2.700381)
        self.b1, self.b2, self.b3 = 0.16818, -0.0323371/10, 0.0019704/100
        self.ret_fac = 0.6827  # 退休后固定支付的工资
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
        
    
        # ===================== 行动空间：=====================
        # 连续选择：消费比例（不允许小于零），风险资产比例，购买个人养老金比例，个人养老金中购买风险资产比例
        self.action_space = spaces.Box(low=np.array([0.01, 0]),\
                                       high=np.array([1.5, 1.5]), shape=(2,), dtype=float)
        
        # 状态空间: [现金, 年龄, 永久收入]
        self.observation_space = spaces.Box(
            low=np.array([0.25, self.tb]),
            high=np.array([50.0, self.td]),
            shape=(2,),
            dtype=float
        )             
        self.state = None
        

    def step(self, action):
        
        # if self.done:
        #     self.state,_ = self.reset()
        
        # Unpack the action
        consum = np.clip(action[0],1e-3,1) # 消费比例
        risky_pct = np.clip(action[1],0,1) # 风险资产比例 

 
        cash, age= self.state 


        #  === 风险资产随机性 ===
        # 随机项
        epsilon = np.random.normal(0,1)        
        # epsilon = np.random.randn(1000)  
        rand_return = self.rf + self.mu + epsilon * self.sigr # 风险资产收益率

        new_cash = 0 
# === 状态转移方程（退休前） ===       
        if self.age < self.tr:

            # 工资
            u_t = np.random.normal(0, self.smay)
            z_t = np.random.normal(0, self.smav)
            new_perm_shock = self.perm_shock + z_t # 工资持续性冲击
             # 如果下一期进入退休期,保存最后一期工作期的工资冲击   
            if self.age+1 ==self.tr:
                self.last_perm_shock = new_perm_shock
                self.last_f = self.f(self.age+1)
            # log_rand_wage = self.f(norm_age) + new_perm_shock + u_t   
            log_rand_wage = self.f(self.age+1) + new_perm_shock + u_t   
            self.wage = np.exp(log_rand_wage)  # 下一期工资
            
            # next state cash normalized by exp(f(age+1)+M_t+1)
            self.cash = (cash* (1 - consum)*np.exp(self.f(self.age) + self.perm_shock)) * (risky_pct * (1 + rand_return) \
                                                  + (1-risky_pct) * (1 + self.rf)) + self.wage #   raw cash
            new_cash = self.cash/np.exp(self.f(self.age+1) + new_perm_shock) # normalized cash
                 
                  
            info = {'age':np.array([self.age]),'normalized_income':self.wage/np.exp(self.f(self.age+1) + new_perm_shock),\
                    'real_actions':np.array([consum,risky_pct]),\
                    'raw_consumption':consum*np.exp(self.f(self.age) + self.perm_shock),\
                    'raw_income':self.wage,'state':self.state,'status':'working'}
                
# === 状态转移方程（退休后） ===   
        else:  
            
            new_perm_shock = self.last_perm_shock # 工资持续性冲击为最后一期工作期的工资冲击保持不变
            log_rand_wage = self.last_perm_shock + self.last_f + np.log(self.ret_fac)   # ret_fac是退休后固定支付的工资,为最后一期工作期工资的比例
            self.wage = np.exp(log_rand_wage)# 下一期工资
            # 当期到手现金
            self.cash = (cash* (1 - consum)*np.exp(self.last_f + self.last_perm_shock)) * (risky_pct * (1 + rand_return) \
                                                  + (1-risky_pct) * (1 + self.rf)) + self.wage #  raw cash
            new_cash = self.cash/np.exp(self.last_perm_shock + self.last_f) # normalized cash

                                               
            # print([new_cash,self.age,'retired'])
            info = {'age':np.array([self.age]),'normalized_income':self.ret_fac,'real_actions':np.array([consum,risky_pct]),\
                    'raw_consumption':consum*np.exp(self.last_f + self.last_perm_shock),\
                    'raw_income':self.wage,'state':self.state,'status':'retired'}    
            
        # 计算效用（奖励函数）        
        # norm_consum = 2 + (cash * consum_pct - 0.5)/(7 - 0.5)
        # reward = ((consum))**(1-self.gamma)/((1-self.gamma)) # normalize reward to [0,1]
        if self.age < self.tr:
            reward = -((consum*cash*np.exp(self.f(self.age) + self.perm_shock))**(1-self.gamma))*1000 + 10 
            # reward = -((consum*np.exp(self.f(self.age) + perm_shock))**(1-self.gamma))             
        else:
            reward = -((consum*cash*np.exp(self.last_perm_shock + self.last_f))**(1-self.gamma))*1000 + 10 
            # reward = -((consum*np.exp(self.last_perm_shock + self.last_f))**(1-self.gamma))

        
        info.update({'消费':consum})    
        # print(norm_consum)
        self.perm_shock = new_perm_shock
        self.age += 1
            
        # Determine if the episode is done
        # 如果完成，则返回同样的state，在计算时game结束的下一期的值函数无论如何都为0（不影响
        # 无法存活至下一期or超过最大年龄       
        if (self.age > self.td):
            self.done = True
            info.update({'PeriodEnd':'最大年龄死亡'})
            return self.state, reward, self.done, False, info
        
        # 随机存活
        # if  np.random.choice([True,False], 1,  p=[1-self.survprob_dict[self.age-1],self.survprob_dict[self.age-1]]): 
        #     self.done = True       
        #     info.update({'PeriodEnd':'提前死亡'})
        #     return self.state, reward, self.done, False, info
            
        # else:
        self.done = False
        
        # 下一期状态
        self.state[0] = new_cash
        self.state[1] = np.array([self.age],dtype=np.float32)   
           
        
        self.index  += 1
        
        return self.state, reward, self.done, False, info

    def reset(self, *, seed: Optional[int] = None, options: Optional[dict] = None):
        super().reset(seed=seed)
        np.random.seed(seed)
        
        self.age = self.tb # 起始年龄
        
        self.index = 0
        self.cash = 0

        # 第一份随机工资收入
        u_t = np.random.normal(0, self.smay)
        z_t = np.random.normal(0, self.smav)
        # perm_shock = z_t # 工资持续性冲击
        self.perm_shock = z_t # 工资持续性冲击
        
        # log_rand_wage = self.f(norm_age) + new_perm_shock + u_t   
        log_rand_wage = self.f(self.age) + self.perm_shock + u_t  
        rand_wage = np.exp(log_rand_wage)
        self.wage = rand_wage
        self.cash = self.init_cash*np.exp(self.f(self.age) + self.perm_shock) + self.wage # raw cash
        # Take a sample from the observation_space to modify the values of
        self.state = self.observation_space.sample()
        self.state[0] = np.array([self.cash/np.exp(self.f(self.age) + self.perm_shock)],dtype=float)  # normalized cash
        self.state[1] = np.array([self.age],dtype=float)
        
        # self.state = np.array([rand_wage,perm_shock,-1],dtype=np.float32)
        info = {'age':self.age,'normalized_income':self.state[0],\
                'raw_income':self.cash,'state':self.state,'status':'working'}
        return self.state, info

    def render(self, mode='human', close=False):
         # Render the environment to the screen
         print(f"State: {self.state}")
    
    # Define the function f(age)
    def f(self, age):
        wage = self.aa +self. b1 * age + self.b2 * age**2 + self.b3 * age**3
        return wage/10
    
    
            
