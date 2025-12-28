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
    id='cocco-dp',
    entry_point='v3_cocco:CoccoEnv',
)
  


class CoccoEnv(gym.Env):
    """Custom Environment that follows gym interface"""
    # metadata = {'render.modes': ['human']}
    
    def __init__(self, params): # 
        
        super(CoccoEnv, self).__init__()
        
         # 设置默认参数
        self.tb = 20 # 96#20  # 初始年龄
        self.tr = 61 # 98#61  # 退休年龄
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
        # 连续选择：消费比例（不允许小于零），风险资产比例
        self.action_space = spaces.Box(low=np.array([0.05, -0.2,0,0]),\
                                       high=np.array([1.5, 1.5,1,1]), shape=(4,), dtype=float)
        
        # 状态空间: [现金, 年龄, 永久收入]
        self.observation_space = spaces.Box(
            low=np.array([0, -1e3,  -0.5]),
            high=np.array([1e3,  1e3, 0.5]),
            shape=(3,),
            dtype=float
        )          
        self.state = None
        

    def step(self, action):
        
        # if self.done:
        #     self.state,_ = self.reset()
        
        # Unpack the action
        consum_pct = np.clip(action[0],0,1) # 消费数量
        risky_pct = np.clip(action[1],0,1) # 风险资产比例 

 
        cash, perm_shock, _ = self.state 

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
            new_perm_shock = perm_shock + z_t # 工资持续性冲击
             # 如果下一期进入退休期,保存最后一期工作期的工资冲击   
            if self.age+1 ==self.tr:
                self.last_perm_shock = new_perm_shock
                self.last_f = self.f(self.age+1)
            # log_rand_wage = self.f(norm_age) + new_perm_shock + u_t   
            log_rand_wage = self.f(self.age+1) + new_perm_shock + u_t   
            new_wage = np.exp(log_rand_wage)# 下一期工资
            
            # next state cash normalized by exp(f(age+1)+M_t+1) 
            new_cash = cash * (1-consum_pct) * (risky_pct * (1 + rand_return) \
                                        + (1-risky_pct) * (1 + self.rf)) + new_wage
                  
            info = {'age':np.array([self.age]),
                    'basic_income':self.wage,\
                    'normalized_cash':new_cash/np.exp(self.f(self.age+1) + new_perm_shock),
                    'raw_income':self.wage,'state':self.state,'status':'working'}
                
# === 状态转移方程（退休后） ===   
        else:  
            
            new_perm_shock = self.last_perm_shock # 工资持续性冲击为最后一期工作期的工资冲击保持不变
            log_rand_wage = self.last_perm_shock + self.last_f + np.log(self.ret_fac)   # ret_fac是退休后固定支付的工资,为最后一期工作期工资的比例
            new_wage = np.exp(log_rand_wage)# 下一期工资
            # 当期到手现金
            new_cash = cash * (1-consum_pct) * (risky_pct * (1 + rand_return) \
                                                  + (1-risky_pct) * (1 + self.rf)) \
                + new_wage  #  raw cash

                                               
            # print([new_cash,self.age,'retired'])
            info = {'age':np.array([self.age]),'basic_income':self.wage,\
                    'normalized_income':self.ret_fac,\
                    'normalized_cash':(new_cash)/np.exp(self.last_perm_shock + self.last_f),\
                    'raw_income':self.wage,'state':self.state,'status':'retired'}    
              
        # 计算效用（奖励函数）        
        # norm_consum = 2 + (cash * consum_pct - 0.5)/(7 - 0.5)
        # reward = ((consum))**(1-self.gamma)/((1-self.gamma)) # normalize reward to [0,1]
        reward = ((consum_pct*cash))**(1-self.gamma)/(1-self.gamma) # normalize reward to [0,1]      
        # 奖励一定要设置合理的下bound, 否则训练不稳定     
        # print(reward)
        if (reward <-5) or (consum_pct==0):
            reward = -5
        
        info.update({'消费':cash * consum_pct})    
        # print(norm_consum)

        self.age += 1
        self.wage = new_wage     
        # Determine if the episode is done
        # 如果完成，则返回同样的state，在计算时game结束的下一期的值函数无论如何都为0（不影响
        # 无法存活至下一期or超过最大年龄       
        if (self.age > self.td):
            self.done = True
            info.update({'PeriodEnd':'最大年龄死亡'})
            return self.state, reward, self.done, False, info
        
        # 随机存活
        if  np.random.choice([True,False], 1,  p=[1-self.survprob_dict[self.age-1],self.survprob_dict[self.age-1]]): 
            self.done = True       
            info.update({'PeriodEnd':'提前死亡'})
            return self.state, reward, self.done, False, info
            
        # else:
        self.done = False
        
        # 下一期状态
        self.state[0] = new_cash
        self.state[1] = new_perm_shock
        self.state[2] = -0.5 + (self.age - self.tb)/self.tn
           
        
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
        perm_shock = z_t # 工资持续性冲击
        
        # log_rand_wage = self.f(norm_age) + new_perm_shock + u_t   
        log_rand_wage = self.f(self.age) + perm_shock + u_t  
        rand_wage = np.exp(log_rand_wage)
        self.wage = rand_wage
        new_cash = self.init_cash*np.exp(self.f(self.age) + perm_shock) # raw cash
        # Take a sample from the observation_space to modify the values of
        self.state = self.observation_space.sample()
        self.state[0] = np.array([new_cash + self.wage],dtype=float)
        self.state[1] = np.array([perm_shock],dtype=float)
        self.state[2] = np.array([-0.5],dtype=float)
        
        # self.state = np.array([rand_wage,perm_shock,-1],dtype=np.float32)
        info = {'age':self.age,'normalized_income':self.state[0],\
                'normalized_cash':(new_cash + self.wage)/np.exp(self.f(self.age) + perm_shock),\
                'raw_income':self.cash,'state':self.state,'status':'working'}
        return self.state, info

    def render(self, mode='human', close=False):
         # Render the environment to the screen
         print(f"State: {self.state}")
    
    # Define the function f(age)
    def f(self, age):
        wage = self.aa +self. b1 * age + self.b2 * age**2 + self.b3 * age**3
        return wage/10
    
            
if __name__=="__main__":
    # 获取当前文件所在路径
    import os
    current_path = os.path.dirname(os.path.abspath(__file__))
    # 读取当前文件所在路径下的A,C,V矩阵, 没有列名和索引
    import pandas as pd
    import numpy as np
    df_A = pd.read_excel(os.path.join(current_path,'result_cocco_matlab','A.xlsx'),header=None)
    df_C = pd.read_excel(os.path.join(current_path,'result_cocco_matlab','C.xlsx'),header=None) 
    df_V = pd.read_excel(os.path.join(current_path,'result_cocco_matlab','V.xlsx'),header=None)  
    df_gcash = pd.read_excel(os.path.join(current_path,'result_cocco_matlab','gcash.xlsx'),header=None)
    # 转换为numpy数组
    A = df_A.to_numpy()
    C = df_C.to_numpy()
    V = df_V.to_numpy()
    gcash = df_gcash.to_numpy().squeeze()

    # 测试环境   
    env = gym.make('cocco-dp',params=None)
    # env = gym.wrappers.RecordEpisodeStatistics(env)
    
    # 模拟100个agent的路径
    n_agents = 1000
    
    # 存储所有agent的轨迹
    all_trajectories = []
    score_history = []

    # 创建进度条
    from tqdm import tqdm
    pbar = tqdm(total=n_agents, desc='模拟进度')
    for agent in range(n_agents):
        # 重置环境获取初始状态
        obs, info = env.reset()
        done = False
        agent_trajectory = {
            'states': [obs.copy()],
            'actions': [],
            'rewards': [],
            'infos': [],
            'raw_consumption': [],
            'normalized_income': [],
            'raw_income': [],
            'age': [obs[1].copy()],
        }
        score = 0
        while not done:

            # 根据state[0]获取现金
            cash = np.squeeze(obs[0])        
            age_loc =  int(obs[1]) - 98# 获取年龄

            # 根据现金和年龄插值从A中获取action
            action = np.array([np.interp(cash, gcash, C[:, age_loc]),np.interp(cash, gcash, A[:, age_loc])])
            if age_loc == np.shape(C)[1]-1: # 如果到达最后一期,则固定消费比例为1（全部消费干净）
                action[0] = cash
            else:
                action[0] =  np.clip(action[0],0,cash)
            action[1] = np.clip(action[1],0,1)

            # if int(obs[1]) == 100:  
            #     raise
            # 固定动作策略
            # action = np.clip(env.action_space.sample(),0,1)
            # print(obs[1],obs[0],obs[2],action[0],action[1])
            # 执行动作
            obs, reward, done, _, info = env.step(action)
            action[0] = action[0]/cash # 消费比例
            score += reward * (env.beta ** len(agent_trajectory['actions']))
            agent_trajectory['actions'].append(action)            
            # 记录轨迹
            agent_trajectory['rewards'].append(reward)
            agent_trajectory['raw_consumption'].append(info['raw_consumption'])
            agent_trajectory['normalized_income'].append(info['normalized_income'])
            agent_trajectory['raw_income'].append(info['raw_income'])
            agent_trajectory['infos'].append(info)
            
            
            if not done:
                agent_trajectory['states'].append(obs.copy()) #                  
                agent_trajectory['age'].append(obs[1].copy())
        score_history.append(score)  
        all_trajectories.append(agent_trajectory)
        pbar.update(1)
    pbar.close()
        # 输出平均reward
    all_rewards = np.concatenate([traj['rewards'] for traj in all_trajectories])
    avg_reward = np.nanmean(all_rewards)
    print(f'平均单期reward: {avg_reward:.2f}')
    print(f'平均终身效用: {np.nanmean(score_history):.2f} +- {np.nanstd(score_history):.2f}')






    # 计算每个年龄的平均收入
    # mean_income_by_age = np.nanmean(normalized_income_by_age, axis=1)
    # print("每个年龄的平均收入:")
    # for age, mean_income in zip(unique_age, mean_income_by_age):
    #     print(f"年龄 {age}: {mean_income:.6f}")




    # all_rewards = np.concatenate([traj['rewards'] for traj in all_trajectories])
    # max_reward = max(all_rewards)
    # min_reward = min(all_rewards)
    # avg_reward = sum(all_rewards) / len(all_rewards)
    
    # print(f'最大单期效用: {max_reward:.2f}')
    # print(f'最小单期效用: {min_reward:.2f}')
    # print(f'平均单期效用: {avg_reward:.2f}')
    
    # # 计算所有agent的平均消费
    # all_consumption = np.concatenate([traj['consumption'] for traj in all_trajectories])
    # avg_consumption = sum(all_consumption) / len(all_consumption)
    # print(f'平均消费: {avg_consumption:.2f}')
    # print(f'最大消费: {max(all_consumption):.2f}')
    # print(f'最小消费: {min(all_consumption):.2f}')
    # print(f'标准差: {np.std(all_consumption):.2f}')
    # print(f'25分位数: {np.percentile(all_consumption, 25):.2f}')
    # print(f'中位数: {np.median(all_consumption):.2f}')
    # print(f'75分位数: {np.percentile(all_consumption, 75):.2f}')

    
    env.close()
    # import matplotlib.pyplot as plt
    # plt.plot(score_history)
    
    
    
            
