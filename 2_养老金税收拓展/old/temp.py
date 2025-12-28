import gymnasium as gym
from gymnasium import spaces
from gymnasium.envs.registration import register
# from gymnasium.utils.env_checker import check_env
from stable_baselines3.common.env_checker import check_env
import numpy as np
import pandas as pd
from typing import Optional
import warnings
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
# Register this module as a gym environment. Once registered, the id is usable in gym.make().
register(
    id='pensionfund-v0',                                # call it whatever you want
    entry_point='v0_pensionfund:PensionFundEnv', # module_name:class_name
)
# 测试环境
# 定义参数字典
params = {
    'pension_limit': 999 # 缴费上限
}

env = gym.make('pensionfund-v0', params=params)

# 模拟100个agent的路径
n_agents = 1000

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
        'raw_income':[],
        'states': [obs],
        'actions': [],
        'rewards': [],
        'infos': [info],
        'consumption': []
    }
    while not done:
        # 固定动作策略
        action = np.clip(env.action_space.sample(),-1,2)
        
        # 执行动作
        obs, reward, done, _, info = env.step(action)
        
        # 记录轨迹
        agent_trajectory['raw_income'].append(info['basic_income'])
        agent_trajectory['states'].append(obs)
        agent_trajectory['actions'].append(action)
        agent_trajectory['rewards'].append(reward)
        agent_trajectory['consumption'].append(info['消费'])
        agent_trajectory['infos'].append(info)
        
    all_trajectories.append(agent_trajectory)
    pbar.update(1)
pbar.close()
    

# 计算所有agent的收入
max_length = max(len(traj['raw_income']) for traj in all_trajectories)

income_path_array = np.full((n_agents, max_length), np.nan)
for i, traj in enumerate(all_trajectories):
    income_path_array[i, :len(traj['raw_income'])] = traj['raw_income']
age_array = np.arange(env.get_wrapper_attr('tb'), env.get_wrapper_attr('tb') + max_length)
full_income_path = np.column_stack((age_array, income_path_array.T))
income_path_df = pd.DataFrame(full_income_path,
                                columns=['年龄'] + [f'个体{i+1}' for i in range(n_agents)])
income_path_df.to_excel("fig\\income_path.xlsx", index=False)

# 画出每个年龄的平均收入，横轴年龄，纵轴收入
age_income_mean = np.nanmean(income_path_array.T, axis=1)
plt.plot(age_array, age_income_mean)
plt.show()

print(f'平均收入: {np.nanmean(income_path_array):.2f}')
print(f'最大收入: {np.nanmax(income_path_array):.2f}')
print(f'最小收入: {np.nanmin(income_path_array):.2f}')
print(f'标准差: {np.nanstd(income_path_array):.2f}')
print(f'25分位数: {np.nanpercentile(income_path_array, 25):.2f}')
print(f'中位数: {np.nanmedian(income_path_array):.2f}')
print(f'75分位数: {np.nanpercentile(income_path_array, 75):.2f}')


env.close()
# import matplotlib.pyplot as plt
# plt.plot(score_history)