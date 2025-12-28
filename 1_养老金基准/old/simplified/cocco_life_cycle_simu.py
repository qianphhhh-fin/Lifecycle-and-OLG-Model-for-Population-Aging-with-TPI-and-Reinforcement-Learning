import numpy as np
import gymnasium as gym
from gymnasium import spaces
from gymnasium.envs.registration import register
import numpy as np
import scipy as sp
import scipy.stats
from scipy.interpolate import CubicSpline
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize
import pandas as pd
import os
import matplotlib.pyplot as plt

cdir = os.path.dirname(os.path.abspath(__file__)) # 获取当前py文件所在目录

from gymnasium.envs.registration import register
register(
    id='cocco-v0',
    entry_point='cocco_env:CoccoEnv',
)
# 加载最优策略数据
data = np.load(cdir + '\\datas\\cocco.npz')
A = data['A']  # 最优风险资产配置比例
V = data['V']  # 值函数
C = data['C']  # 最优消费比例
gcash = data['gcash']  # 状态变量的grid

params = {
        'tb': 20, # 初始年龄
        'td': 101, # 最大年龄
        'r': 1.015, # 无风险收益
        'mu': 0.04, # 超额收益
        'sigr': 0.2, # 风险资产收益率标准差
        'riskaversion': 3, # CRRA风险规避系数
        'delta': 0.97, # 贴现因子
        'smay': 0.1, # 收入冲击标准差
        'initial_cash': 3, # 初始现金
    }

env = gym.make('cocco-v0',params=params)
# env = gym.wrappers.FlattenObservation(env)  # deal with dm_control's Dict observation space
env = gym.wrappers.RecordEpisodeStatistics(env)

# 模拟100个agent的路径
n_agents = 1000

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
        # print(obs)
        # 随机动作 (消费比例, 风险资产比例)
        # action = env.action_space.sample()
        # 根据obs所属的grid，得到策略函数（C,A) 的插值
        action = np.zeros(2)
        action[0] = np.interp(obs[0], gcash.squeeze(), C[:,int(obs[1]-params['tb'])])
        action[1] = np.interp(obs[0], gcash.squeeze(), A[:,int(obs[1]-params['tb'])])
        # if obs[1]-params['tb']==80:
        #     print(obs[1])
        #     raise
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

# 首先获取所有trajectory中所包含最大的年龄
max_length = int(max([max([state[1] for state in traj['states']]) for traj in all_trajectories]) - params['tb'])
# 将all_trajectories中每个trajectory的action扩展至max_age
for traj in all_trajectories:
    traj['actions'] = np.pad(traj['actions'], ((0, int(max_length) - len(traj['actions'])), (0, 0)), mode='constant', constant_values=np.nan)
    
avg_actions = np.zeros((int(max_length), 2))

# 将所有轨迹的actions堆叠起来
actions = np.array([traj['actions'] for traj in all_trajectories])
# 计算每个时间点的平均动作
avg_actions = np.nanmean(actions, axis=0)

# 画图：x轴为tb:td, y轴为消费比例和风险资产比例
plt.figure(figsize=(10, 6))
plt.plot(np.arange(params['tb'], params['tb'] + max_length), avg_actions[:, 0], label='Consumption')
plt.plot(np.arange(params['tb'], params['tb'] + max_length), avg_actions[:, 1], label='Risk Asset')
plt.xlabel('Age')
plt.ylabel('Action')
plt.legend()
plt.show()

# # 计算所有agent的reward统计信息
# all_rewards = [sum(traj['rewards']) for traj in all_trajectories]
# max_reward = max(all_rewards)
# min_reward = min(all_rewards)
# avg_reward = sum(all_rewards) / len(all_rewards)

# print(f'最大累积奖励: {max_reward:.2f}')
# print(f'最小累积奖励: {min_reward:.2f}') 
# print(f'平均累积奖励: {avg_reward:.2f}')

# # 计算所有agent的平均消费
# all_consumptions = [sum(traj['consumptions']) for traj in all_trajectories]
# avg_consumption = sum(all_consumptions) / len(all_consumptions)
# print(f'平均消费: {avg_consumption:.2f}')
# print(f'最大消费: {max(all_consumptions):.2f}')
# print(f'最小消费: {min(all_consumptions):.2f}')

# print("\n所有Agent模拟完成!")
        
env.close()