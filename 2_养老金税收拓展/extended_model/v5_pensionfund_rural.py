# -*- coding: utf-8 -*-
"""
v5-extend: 城乡居民养老保险+个人养老金
基于gymnasium标准格式的env，随机死亡，累进税率，个人养老金购买上限
连续选择：消费比例，风险资产比例，购买城乡居民养老保险比例，购买个人养老金比例(若基本养老金个人账户余额不为0，则允许购买个人养老金)，个人养老金中购买风险资产比例
状态变量：财富，上期工资持续冲击，基本养老金个人账户余额, 城乡保缴费年限, 个人养老金账户余额，年龄，
注：cash*np.exp(self.subtract_income)和pen_balance*np.exp(self.subtract_income)的单位万元/100，norm_age为-0.5到0.5
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
    id='pensionfund-v5-rural',                                # call it whatever you want
    entry_point='v5_pensionfund_rural:PensionFundEnv', # module_name:class_name
)


'''
v5: 养老金的约束为比例
'''

#         #按年龄收入调整系数 2023人均年工资性收入2.21万元 
#         # (此参数平均收入:21岁1.766万元，49岁3.592万元，64岁3.329万元)
#         # desmos: \frac{e^{0.5+0.172x-0.00323371x^{2}+0.000019x^{3}}}{100}


class PensionFundEnv(gym.Env):
    """Custom Environment that follows gym interface"""
    # metadata = {'render.modes': ['human']}
    
    def __init__(self,params=None): # 
        
        super(PensionFundEnv, self).__init__()
                
         # 设置默认参数
        self.scale_inc = 1 # 2
        self.subtract_income = 0 # 相当于将以万为单位的收入除以np.exp(self.subtract_income), 若为0，则收入单位为万元
        self.scale_reward = 1 # 小于1有助于增强探索
        self.reward_lb = -0.1 # 单期reward的下限

        self.delta = 0.2
        
        self.tb = 18 # 初始年龄
        self.tr = 61 # 退休年龄
        self.td = 100  # 最大年龄
        self.tn = self.td - self.tb + 1
        self.rf = 0.02  # 无风险收益
        self.mu = 0.04  # 超额收益
        self.sigr = 0.27  # 风险资产收益率标准差
        self.gamma = 3.84 # 3.84  # 风险厌恶系数
        self.beta = 0.95  # 贴现因子
        self.aa = -1.09533849e+01
        self.b1 = 1.19983877e+00
        self.b2 = -4.18823630e-02
        self.b3 = 6.23094506e-04
        self.b4 = -3.40219380e-06
        self.ar = 0.838679
        self.smay = np.sqrt(0.269354) # 暂时性shock的标准差
        self.smav = np.sqrt(0.155491) # 永久性shock的标准差
        # self.aa,self.b1, self.b2, self.b3 = 0.5,0.172,-0.00323371,0.000019
        self.ret_fac = 0.6827  # 退休后固定支付的基本养老金,为退休前工资的68.27%
        self.init_cash = 0  # 初始现金
        
        self.basic_pension = 200/10000*12 # 基础养老金，单位：万元 各地平均100-200元/月, 上海最多，2024年1490元/月

        self.tda_tax = 0.03 # 个人养老金税率
        self.fund_excess_return = 0.04  # 养老金风险投资超额收益
        self.fund_riskfree_return = 0.04 # 0.02 无差异  # 养老金无风险收益
        self.sigma_eta = 0.27  # 养老金收益随机性的标准差
        self.corr_fund_risky = 0.5  # 养老金收益与股市收益的相关性
        self.pension_limit = 1 # 缴费上限 100%的工资

        # 如果传入参数则覆盖默认值
        if params:
            for key, value in params.items():
                setattr(self, key, value)
        
        
        
        # 20-100岁的生存概率
        survprob_values = [0.9999,0.9998,
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
        self.survprob_dict = {age: value for age, value in zip(range(self.tb, self.td+1), survprob_values)}
        # self.survprob_dict = {age: 1 for age, value in zip(range(20, 101), survprob_values)} # 无随机死亡
    
        # ===================== 行动空间：=====================
        # 连续选择：消费比例（不允许小于零），风险资产比例，城乡居民养老保险购买比例（最多为工资的50%），个人养老金购买比例(最多为工资的30%)，个人养老金中购买风险资产比例
        self.action_space = spaces.Box(low=np.array([0.1, 0, -0.2, -0.2, 0]),\
                                       high=np.array([1.2, 1.2, 0.5, 0.3, 1.2]), shape=(5,), dtype=float)
        
        # ===================== 状态空间：=====================
        # 现金余额，永久收入冲击, 基本养老金个人账户余额，城乡保缴费年限，个人养老金账户余额，年龄
        self.observation_space = spaces.Box(
            low =np.array([0,   -1e20,   0, -0.5,   0,  -0.5]),
            high=np.array([1e20,  1e20, 1e20,  0.5, 1e20,   0.5]),
            shape=(6,),
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
        basic_pen_pct = np.clip(action[2],0,1) # 购买城乡居民养老保险比例
        tda_pct = np.clip(action[3],0,1) # 购买个人养老金比例
        tda_risky_pct = np.clip(action[4],0,1) # 个人养老金中购买风险资产比例
        

  
        cash,perm_shock,basic_pen_balance,basic_pen_years,pen_balance,_ = self.state
        
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
            new_perm_shock = self.ar*perm_shock + z_t # 工资持续性冲击
            # 如果下一期进入退休期,保存最后一期工作期的工资冲击   
            if self.age+1 ==self.tr:
                self.last_perm_shock = new_perm_shock
                self.last_f = self.f(self.age+1)
            # log_rand_wage = self.f(norm_age) + new_perm_shock + u_t   
            log_rand_wage = self.f(self.age+1) + new_perm_shock + u_t   
            new_wage = max(np.exp(log_rand_wage),0.6) # 下一期工资
            
            # 养老金风险收益
            # 生成独立的标准正态随机变量   
            zeta = np.random.normal(0, 1)  # zeta = np.random.randn(1000)            
            eta = self.corr_fund_risky * epsilon + (1 - self.corr_fund_risky**2)**0.5 * zeta
            fund_return = self.fund_riskfree_return  + self.fund_excess_return + \
                 self.sigma_eta * eta
            # 计算相关性
            # np.corrcoef(np.vstack([[fund_return],[rand_return]]))
            
            '''
            购买城乡居民养老保险
            '''
            basic_pen_amount = self.wage * basic_pen_pct # 城乡居民养老保险购买额(实际)
            basic_pen_years += 1/(self.tr - self.tb) if basic_pen_amount>0 else 0 # 城乡保缴费年限

            '''
            购买个人养老金
            '''
            actual_tda_pct = tda_pct.copy()
            if basic_pen_amount>0 or basic_pen_balance>0:
                
                if hasattr(self, 'pension_limit_yuan'): # 如果存在这个参数
                    tda_pct = min(tda_pct, self.pension_limit_yuan/(self.wage*np.exp(self.subtract_income))) if cash>0 else 0 # 这里penlimit的单位是万元，要 除以np.exp(self.subtract_income)       
                else:
                    tda_pct = min(tda_pct, self.pension_limit) if cash>0 else 0     
                pen_amount = (self.wage - basic_pen_amount) * tda_pct # 个人养老金购买额(实际)
                opt_pen_amount = (self.wage - basic_pen_amount) * actual_tda_pct # 个人养老金购买额(最优)
                self.info = {'个人养老金购买资格':True}
            else:
                pen_amount = 0
                opt_pen_amount = self.wage * actual_tda_pct
                self.info = {'个人养老金购买资格':False}

            
            # 下一期到手现金 这里的处理与Gomes et al. (2009)相同, 
            '''
            纳税
            '''    
             # 在缴纳阶段购买基本养老金和个人养老金均免税      
            tax = self.taxation(self.wage) # 当期已纳税额
            tax_refund = tax - self.taxation(self.wage - pen_amount - basic_pen_amount) # 购买养老金的节税额
            # if tax>0:
            #     raise
            cash = cash + tax_refund - pen_amount - basic_pen_amount #  tax_refund - pen_amount 对应 Gomes et al. (2009)的 -k_tY_t(1-\tau)
            # if tax_refund>tax:
            #     raise
            '''
            现金
            '''
            new_cash = cash * (1-consum_pct) * (risky_pct * (1 + rand_return) \
                                                  + (1-risky_pct) * (1 + self.rf)) + new_wage - self.taxation(new_wage)
            
            
            '''
            城乡保个人账户余额和个人养老金账户余额
            '''            
            new_basic_pen_balance = (basic_pen_balance + basic_pen_amount) * (1 + self.rf)                   
            new_pen_balance  = (pen_balance + pen_amount) *  (tda_risky_pct * (1 + fund_return) +
                                      (1 - tda_risky_pct) * (1 + self.rf)) 
            self.real_pension = new_pen_balance
                  
            self.info.update({'age':np.array([self.age]),
                    'basic_income':self.wage*np.exp(self.subtract_income),\
                    # '每年个人养老金购买额(实际)':np.array([pen_amount*np.exp(self.subtract_income)]),\
                    # '每年个人养老金购买额(最优)':np.array([opt_pen_amount*np.exp(self.subtract_income)]),\
                    # '每年个人养老金购买比例(实际)':np.array([tda_pct]),\
                    # '每年个人养老金购买比例(最优)':np.array([actual_tda_pct]),\
                    # '每年个人养老金退税额':np.array([tax_refund*np.exp(self.subtract_income)]),\
                    # '每年收入应缴税额(万元)':np.array([tax*np.exp(self.subtract_income)]),\
                    'real_actions':np.array([consum_pct,risky_pct,basic_pen_pct,tda_pct,tda_risky_pct]), # 记录实际的action
                    'normalized_cash':(new_cash)/np.exp(self.f(self.age+1) + new_perm_shock),\
                    'state':self.state,'status':'working'})
                
# === 状态转移方程（退休后仍然有随机性的收入，但仅有退休前的ret_fac比例） ===   
        else:  
            # 工资
            u_t = np.random.normal(0, self.smay)
            z_t = np.random.normal(0, self.smav)
            new_perm_shock = self.ar*perm_shock + z_t # 工资持续性冲击
            # log_rand_wage = self.last_f + new_perm_shock + u_t #+ np.log(self.ret_fac)   
            log_rand_wage = self.f(self.age+1) + new_perm_shock + u_t #+ np.log(self.ret_fac)   
            new_wage = max(np.exp(log_rand_wage),0.6) # 最低年收入 6000元 # 下一期工资
            
            '''
            领取城乡保养老金
            个人账户年度领取额: 余额/101*12, 101为65岁的计发月数
            基本账户: self.basic_pension = 200*12/10000万元
            必须缴满15年才能领取
            '''
            basic_pension = (basic_pen_balance/101*12 + self.basic_pension*np.exp(self.subtract_income)) \
                if basic_pen_years>=-0.5+15/(self.tr - self.tb) else 0 
            '''
            领取个人养老金
            '''
            pension = pen_balance  * (1 - self.tda_tax) / (self.td - self.tr + 1) # 领取个人养老金
            # pension = np.exp(np.log(np.exp(np.log(pen_balance)*self.income_rescale)* (1 - self.tda_tax))/self.income_rescale)/ (self.td - self.tr + 1)
            new_pen_balance  = pen_balance  # 养老金余额不再变化 
            new_basic_pen_balance = basic_pen_balance

            # 当期到手现金
            new_cash = cash * (1-consum_pct) * (risky_pct * (1 + rand_return) \
                                                  + (1-risky_pct) * (1 + self.rf)) \
                + pension + basic_pension + new_wage - self.taxation(new_wage)
            self.real_pension -= pension # 个人养老金账户实际余额  

            # print([new_cash,self.age,'retired'])
            self.info = {'age':np.array([self.age]),
                    'basic_income':(self.wage)*np.exp(self.subtract_income),\
                    'pension_income':(pension + basic_pension)*np.exp(self.subtract_income),\
                    # '每年个人养老金购买比例(实际)':np.array([0]),\
                    'real_actions':np.array([consum_pct,risky_pct,0,0,0]), # 记录实际的action
                    'normalized_cash':(new_cash)/np.exp(new_perm_shock + self.last_f),\
                    # '当期养老金缴税额':np.array([pension*self.tda_tax*np.exp(self.subtract_income)]),\
                    # '养老金实际余额':np.array([self.real_pension*np.exp(self.subtract_income)]),\
                    'state':self.state,'status':'retired'}
           


        
        # 计算效用（奖励函数）
        
        # reward = -((consum_pct*cash))**(1-self.gamma)/10 # normalize reward to [0,1]      
        # reward = np.log((consum_pct*cash))

        # if consum_pct*cash < self.delta:
        #     C = self.gamma * (1-self.gamma) *(1-self.delta)**(-self.gamma-1)/2
        #     B = (1-self.gamma)*(1-self.delta)**(-self.gamma)- 2 * C * self.delta
        #     A =-(1-self.delta)**(1-self.gamma) - B * self.delta - C * self.delta ** 2
        #     reward = A + B * consum_pct*cash + C * (consum_pct*cash) ** 2
        # else:
        reward = -(consum_pct*cash*10) ** (1 - self.gamma)
 


        # 奖励一定要设置合理的下bound, 否则训练不稳定     
        # print(reward)
        reward = np.maximum(reward, self.reward_lb) 

        reward *= self.scale_reward # 缩放reward
        # Bounded below 
        # if reward<-150:
        #     reward = -150   
        # info.update({'消费':norm_consum})
        self.info.update({'消费':cash * consum_pct})    
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
        
        # 下一期状态 # 现金余额，永久收入冲击, 基本养老金个人账户余额，城乡保缴费年限，个人养老金账户余额，年龄
        self.state[0] = new_cash  
        self.state[1] = new_perm_shock
        self.state[2] = new_basic_pen_balance
        self.state[3] = basic_pen_years
        self.state[4] = new_pen_balance  
        self.state[5] = -0.5 + (self.age - self.tb)/self.tn


        if np.isnan(self.state).any():
            self.state = np.nan_to_num(self.state, 0)
        
        return self.state, reward, self.done, False, self.info

    def reset(self, *, seed: Optional[int] = None, options: Optional[dict] = None):
        super().reset(seed=seed)
        np.random.seed(seed)
        
        self.age = self.tb # 起始年龄
        self.index = 0
        self.real_pension = 0 # 养老金实际余额

        # 第一份随机工资收入
        u_t = np.random.normal(0, self.smay)
        z_t = np.random.normal(0, self.smav)
        perm_shock = z_t # 工资持续性冲击
        # perm_shock = max(z_t,0) # 工资持续性冲击
        
        log_rand_wage = self.f(self.age) + perm_shock  + u_t  
        rand_wage = np.exp(log_rand_wage)
        self.wage = rand_wage
        new_cash = self.init_cash*np.exp(self.f(self.age) + perm_shock) # raw cash
        # print(rand_wage)
        # self.wage = np.array([rand_wage],dtype=float) # 第一份工资
        
        # Take a sample from the observation_space to modify the values of
        self.state = self.observation_space.sample()

        self.state[0] = np.array([new_cash + self.wage-self.taxation(self.wage)],dtype=float) # 现金余额
        self.state[1] = np.array([perm_shock],dtype=float) # 永久收入冲击
        self.state[2] = np.array([0],dtype=float)  # 基本养老金个人账户余额
        self.state[3] = np.array([-0.5],dtype=float) # 城乡保缴费年限
        self.state[4] = np.array([0],dtype=float) # 个人养老金账户余额
        self.state[5] = np.array([-0.5],dtype=float) # 年龄

        # self.state = np.array([rand_wage,perm_shock,-1],dtype=np.float32)
        self.info = {'age':self.age,'basic_income':self.wage*np.exp(self.subtract_income),\
                          'state':self.state,
                          'normalized_cash':(new_cash + self.wage-self.taxation(self.wage))/np.exp(self.f(self.age) + perm_shock),
                          'status':'wokring'}
        if np.isnan(self.state).any(): # 检测是否有nan，返回的state确保没有nan，以免弄坏vec_env的rms
            self.state = np.nan_to_num(self.state, 0)
            
        return self.state, self.info

    def render(self, mode='human', close=False):
         # Render the environment to the screen
         print(f"State: {self.state}")
    
    # Define the function f(age)
    def f(self, age):
        wage = self.aa +self. b1 * age + self.b2 * (age)**2 + self.b3 * (age)**3 + self.b4 *(age)**4
        return wage/self.scale_inc - self.subtract_income # 将收入缩放至效用函数斜率较大的范围内
    
    # 累进制税率  (单位：万元)  
    def taxation(self, income_scaled):
        # 将收入从万元转换为元
        income_in_ten_thousand = income_scaled*np.exp(self.subtract_income)
        income = income_in_ten_thousand * 10000 - 60000 # 60000为免税额
        
        if income <= 0:
            return 0
        
        if income <= 0:
            return 0
        
        # 根据税率表确定税率和速算扣除数
        if income <= 36000:
            tax_rate = 0.03
            deduction = 0
        elif income <= 144000:
            tax_rate = 0.10
            deduction = 2520
        elif income <= 300000:
            tax_rate = 0.20
            deduction = 16920
        elif income <= 420000:
            tax_rate = 0.25
            deduction = 31920
        elif income <= 660000:
            tax_rate = 0.30
            deduction = 52920
        elif income <= 960000:
            tax_rate = 0.35
            deduction = 85920
        else:
            tax_rate = 0.45
            deduction = 181920

        # 计算税款
        # 计算税款
        tax = (income * tax_rate - deduction)/10000
        return tax/np.exp(self.subtract_income)
        
    

            
if __name__=="__main__":

    # 测试环境
    # 定义参数字典
    params = {
        'pension_limit': 1 # 缴费上限
    }
    
    env = gym.make('pensionfund-v5-rural', params=params)
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
            'raw_income':[info['basic_income']],
            'states': [obs],
            'actions': [],
            'rewards': [],
            'infos': [info],
            'consumption': [],
            'norm_age':[obs[-1]],
            'age':[info['age']]
        }
        while not done:
            # 固定动作策略
            # 连续选择：消费比例，风险资产比例，购买城乡居民养老保险比例，购买个人养老金比例(若基本养老金个人账户余额不为0，则允许购买个人养老金)，个人养老金中购买风险资产比例
            action = np.clip(env.action_space.sample(),-1,2)
            action[0] = np.clip(action[0],1,1)
            action[2] = 0.1
            action[3] = 0.1
            
            # 执行动作
            obs, reward, done, _, info = env.step(action)
            
            # 记录轨迹
            if info['status'] == 'working':
                agent_trajectory['raw_income'].append(info['basic_income'])
                agent_trajectory['age'].append(info['age'][0])

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
    all_income = np.concatenate([traj['raw_income'] for traj in all_trajectories])
    all_age = np.concatenate([traj['age'] for traj in all_trajectories])
    # import pandas as pd
    # age_income = pd.DataFrame({'age':all_age,'income':all_income})
    # mean_income_by_age = age_income.groupby('age')['income'].mean()
    max_income = max(all_income)
    min_income = min(all_income)
    avg_income = sum(all_income) / len(all_income)
    std_income = np.std(all_income)

    print(f'最大单期收入: {max_income:.2f} 万元')
    print(f'最小单期收入: {min_income:.2f} 万元')
    print(f'平均单期收入: {avg_income:.2f} 万元')
    print(f'单期收入std: {std_income:.2f} 万元')

    all_consumption = np.concatenate([traj['consumption'] for traj in all_trajectories])
    max_consumption = max(all_consumption)
    min_consumption = min(all_consumption)
    avg_consumption = sum(all_consumption) / len(all_consumption)
    std_consumption = np.std(all_consumption) 

    print(f'最大单期消费: {max_consumption:.2f}')
    print(f'最小单期消费: {min_consumption:.2f}')
    print(f'平均单期消费: {avg_consumption:.2f}')
    print(f'单期消费std: {std_consumption:.2f}')

    all_reward = np.concatenate([traj['rewards'] for traj in all_trajectories])
    max_reward = max(all_reward)
    min_reward = min(all_reward)
    avg_reward = sum(all_reward) / len(all_reward)
    std_reward = np.std(all_reward)

    print(f'最大单期奖励: {max_reward:.10f}')
    print(f'最小单期奖励: {min_reward:.10f}')
    print(f'平均单期奖励: {avg_reward:.10f}')
    print(f'单期奖励中位数: {np.median(all_reward):.10f}')
    print(f'单期奖励std: {std_reward:.10f}')
    



    # all_age = np.concatenate([traj['norm_age'] for traj in all_trajectories])
    # print(np.unique(all_age))
    # print(np.shape(np.unique(all_age)))

    
    env.close()
    # import matplotlib.pyplot as plt
    # plt.plot(score_history)
    
    
    
            
