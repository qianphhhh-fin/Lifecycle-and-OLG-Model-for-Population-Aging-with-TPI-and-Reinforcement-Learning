import csv
import os
from tqdm import tqdm
import pandas as pd
import numpy as np
import gymnasium as gym
from statsmodels.tools.tools import add_constant
from linearmodels.panel.model import PanelOLS
from arch import arch_model
import numpy as np
import matplotlib.pyplot as plt
# Register this module as a gym environment. Once registered, the id is usable in gym.make().
from gymnasium.envs.registration import register
register(
    id='pensionfund-v5-rural',                                # call it whatever you want
    entry_point='v5_pensionfund_rural:PensionFundEnv', # module_name:class_name
)


'''
模拟数据
'''
def simulate_data(params,n_agents = 10000,age_range=(18,60)):
    # 测试环境
    # 定义参数字典


    env = gym.make('pensionfund-v5-rural', params=params)
    env = gym.wrappers.RecordEpisodeStatistics(env)

    # 模拟100个agent的路径
    

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
            # if info['status'] == 'working':
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
    import pandas as pd
    age_income = pd.DataFrame({'age':all_age,'income':all_income})
    mean_income_by_age_simu = age_income.groupby('age')['income'].mean()
    mean_income_by_age_simu.name = 'simu'
    mean_income_by_age_simu = mean_income_by_age_simu.loc[age_range[0]:age_range[1]]
    mean_income_by_age_f = pd.Series([np.exp(env.f(i)) for i in range(age_range[0],age_range[1]+1)], index=mean_income_by_age_simu.index)
    mean_income_by_age_f.name = 'f'
    return mean_income_by_age_simu,mean_income_by_age_f



'''
变量选择
pid: 个体唯一ID
age:年龄
qa301:户口类型
emp_income: 过去12个月所有工作（主要工作+一般工作）的税后工资性收入
retire: 是否退休 （“已退休”/“不适用”）
'''
n_sims = 10000
var_list = ['pid','age','qa301','emp_income','retire']


# 读取数据
da_rural = pd.read_csv('./cfps_data/da_rural.csv', delimiter=",", encoding='gbk')
da_non_rural = pd.read_csv('./cfps_data/da_non_rural.csv', delimiter=",", encoding='gbk')

    # 设置字体路径
plt.rcParams['font.sans-serif'] = ['SimHei']  # 使用简体中文字体
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题
# 使用subplot，左图为城镇，右图为农村
fig, axs = plt.subplots(1, 2, figsize=(8, 3))
for da in [da_rural,da_non_rural]:
    '''
    真实数据
    '''
    data = da.set_index(['pid', 'year'])
    data = data.dropna() # 删除包含 NaN 的行
    # 删除'emp_income_before_tax'小于10的行
    data = data.loc[data['emp_income_before_tax'] >= 10,:]
    data['age'] = data['age'].apply(lambda x: int(x) if x.isdigit() else np.nan)
    '''
    城乡69岁仍未退休的样本有32个，70岁仍未退休的样本有17个
    因此城乡居民取到69岁,此后为拟合数据
    '''
    if da is da_rural:
        data = data.loc[(data['age']>=18) & (data['age']<=75),:] 
    else:
        data = data.loc[(data['age']>=18) & (data['age']<=60),:] 

    data = data.groupby('pid').filter(lambda x: x['age'].nunique() >= 2)

    mean_income_by_age_real = data.groupby('age')['emp_income_before_tax'].mean()/10000
    mean_income_by_age_real.name = 'real'
    '''
    模拟数据
    '''
    if da is da_non_rural:
        # 城镇参数
        params = {
            'aa': -1.32157880e+01,
            'b1': 1.34951371e+00,
            'b2': -4.33632834e-02,
            'b3': 5.85163503e-04,
            'b4':-2.81009579e-06,
            'ar': 0.917415,
            'smay': np.sqrt(0.167829), # 暂时性shock的标准差
            'smav': np.sqrt(0.084805), # 永久性shock的标准差
        }
    else:
        # 城乡参数
        params = {
            'aa': -1.09533849e+01,
            'b1': 1.19983877e+00,
            'b2': -4.18823630e-02,
            'b3': 6.23094506e-04,
            'b4':-3.40219380e-06,
            'ar': 0.838679,
            'smay': np.sqrt(0.269354), # 暂时性shock的标准差
            'smav': np.sqrt(0.155491), # 永久性shock的标准差
        }
    mean_income_by_age_simu,mean_income_by_age_f = simulate_data(params,n_agents = n_sims, age_range=(18,100) if da is da_rural else (18,60))

    merge_income = pd.merge(mean_income_by_age_real, mean_income_by_age_simu, left_index=True, right_index=True, how='outer')
    merge_income = pd.merge(merge_income, mean_income_by_age_f, left_index=True, right_index=True, how='outer')
    merge_income.columns = ['real','simu','f']
    print(merge_income)


    if da is da_non_rural:
        ax = axs[0]
        ax.set_title('城镇')
    else:
        ax = axs[1]
        ax.set_title('城乡')
    ax.plot(merge_income.index, merge_income['real'], label='CFPS', marker='o',markersize=3)
    ax.plot(merge_income.index, merge_income['simu'], label='模拟', marker='o',markersize=3)
    ax.plot(merge_income.index, merge_income['f'], label='$\hat f_t$', marker='o',markersize=3)
    if da is da_rural:
        ax.legend()
    ax.set_xlabel('年龄')
    if da is da_non_rural:
        ax.set_ylabel('收入(万元)')
    # 添加grid
    ax.grid(True, linestyle='--', alpha=1)
# 保存为png

plt.tight_layout()
plt.savefig('fig/p6-income_path_comparison.png',dpi=300)

