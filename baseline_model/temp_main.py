from typing import Callable
import json
import gymnasium as gym
import numpy as np
import os
# from stable_baselines3 import DDPG
# from stable_baselines3 import TD3
# from stable_baselines3 import PPO
from sbx import SAC
# from utils.sac import SAC
from utils.env_util import make_vec_env
from stable_baselines3.common.vec_env import VecNormalize
from gymnasium.envs.registration import register
import pandas as pd
import matplotlib.pyplot as plt
# import matplotlib
# print(matplotlib.matplotlib_fname())
# raise

import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

# get current directory of the current .ipynb 
cdir = os.path.abspath('.') # 获取当前目录
pardir = os.path.abspath(os.path.join(os.getcwd(), "../")) # 获取上级目录
dbdir = os.path.abspath(os.path.join(os.getcwd(), "../..")) # 获取上上级目录(为了读取美股数据库数据)


# import torch 
# torch.autograd.set_detect_anomaly(True)
import warnings
warnings.filterwarnings('ignore')

# 将变量处理为终身路径
def var_path(vars, eval_episodes, env):
    max_length = max(len(path) for path in vars)
    var_path_array = np.full((eval_episodes, max_length), np.nan)
    for i, path in enumerate(vars):
        var_path_array[i, :len(path)] = path
    age_array = np.arange(env.get_attr('tb')[0], env.get_attr('tb')[0] + max_length)
    full_var_path = np.column_stack((age_array, var_path_array.T))
    var_path_df = pd.DataFrame(full_var_path,
                                    columns=['年龄'] + [f'个体{i+1}' for i in range(eval_episodes)])
    return var_path_df

import sys
cdir = os.path.abspath('.') # 获取当前目录
# print(sys.path)
# 状态空间：现金余额，个人养老金账户余额，永久收入冲击，年龄
register(
    id='pensionfund-v5',                                # call it whatever you want
    entry_point='v5_pensionfund:PensionFundEnv', # module_name:class_name
)
log_dir = '..//autodl-tmp//models//sac'

eval_dir = 'eval//'
# model_name = ['v13_pensionfund_penlim0_run7',
#               'v13_pensionfund_penlim1p2_run15',
#               'v13_pensionfund_nopenlim_run7'
#               ]
rl_model_name = 'pensionfund-v5_run50//RT//RT_3//RT_1//best_model' # 
rl_vec_env_name = 'pensionfund-v5_run50//RT//RT_3//RT_1//best_model_vec_env.pkl'


# model_name = [d for d in os.listdir(log_dir) if os.path.isdir(os.path.join(log_dir, d))]
eval_episodes = 10000
eval_seed =   3687851522
distf = 0.95
import multiprocessing as mp
from functools import partial

# def evalute_model_erl(idx, eval_episodes=1000, eval_seed=3687851522):


def evaluate_model_rl(model_name,vec_name, eval_episodes=1000, eval_seed=3687851522):
    current_path = os.path.dirname(os.path.abspath(__file__))    
    # try:
    #     with open(os.path.join(current_path, log_dir, model_name.split('//')[0]) +"//params.json", 'r') as file:
    #         params = json.load(file)   
    # except:
    params = {'pension_limit':1}
    # if 'cocco' in idx:
    #     env = make_vec_env("cocco-rl-norm",seed=1, n_envs=1, env_kwargs={'params': params})
    # else:
    # env = make_vec_env("pensionfund-v3",seed=1, n_envs=1, env_kwargs={'params': params})````````
    env = make_vec_env("pensionfund-v5",seed=1, n_envs=1, env_kwargs={'params': params},\
                            monitor_kwargs={'discount_factor':distf})
    

    env = VecNormalize.load(os.path.join(current_path, log_dir,vec_name), env)
    # env.obs_rms.mean
    # env.obs_rms.var
    
    model = SAC.load(os.path.join(current_path, log_dir, model_name), env=env)

    

    # env = make_vec_env("pensionfund-v1",seed=eval_seed, n_envs=1, env_kwargs={'params': None})
    
    # env = VecNormalize.load(os.path.join(current_path, log_dir, idx) + "//eval_env.pkl",env)

    env = make_vec_env("pensionfund-v5",seed=eval_seed, n_envs=1, env_kwargs={'params':params},\
                            monitor_kwargs={'discount_factor': distf})         
    env = VecNormalize(env, norm_obs=False, norm_reward=False) # 评估不需要norm_obs


    idx = model_name.split('//')[0]
    # model = DDPG.load(os.path.join(log_dir, idx) + "//latest_model", env)

    action_set = []
    age = []
    pen_amount = []
    pen_amount_temp = []
    income = []
    tmp_income = []
    
    # from tqdm import tqdm
    pbar = tqdm(total=eval_episodes, desc=f'模拟进度 {idx}')
    
    for episode in range(eval_episodes):       
        pbar.update(1)
        env.seed(seed=eval_seed+episode)
        obs = env.reset()    
        done = False
        score = 0
        episode_actions = []
        while not done:
            action, _ = model.predict(model.env.normalize_obs(obs), deterministic=True) 
            # action[0,2] = 0

            obs, reward, done, info = env.step(action)               
            # reward = -(((cash*action[0][0]))**(1-env.get_attr('gamma')[0]))*10 + 10       
            score += reward * (distf ** len(episode_actions))
            episode_actions.append(info[0]['real_actions'])
            tmp_income.append(info[0]['basic_income'])
            try:
                pen_amount_temp.append(info[0]['每年个人养老金购买额(实际)'][0])
            except:
                pen_amount_temp.append(np.nan)
            # if info[0]['status']=='working':
        action_set.append(episode_actions)
        age.append(info[0]['age'][0])
        pen_amount.append(pen_amount_temp)
        pen_amount_temp = []
        income.append(tmp_income)
        tmp_income = []
    pen_amount_df = var_path(pen_amount, eval_episodes, env)
    income_history_df = var_path(income, eval_episodes, env)

    max_length = max(len(path) for path in action_set)
    unique_age = np.arange(env.get_attr('tb')[0], env.get_attr('td')[0] + 1)[:max_length]
    num_action = np.shape(action_set[0])[1]
    aligned_actions = np.full((max_length, num_action, eval_episodes), np.nan)
    for i, actions in enumerate(action_set):
        aligned_actions[:len(actions), :, i] = actions
        

    return unique_age, aligned_actions, income_history_df,pen_amount_df

if __name__ == '__main__':


    # 创建进程池
    # pool = mp.Pool()
    
    # # 循环执行evaluate_model
    from tqdm import tqdm
    unique_age, rl_actions, income_history_df,pen_amount_df = evaluate_model_rl(rl_model_name,rl_vec_env_name,eval_episodes=eval_episodes, eval_seed=eval_seed)
    
    ''' 
    （1）根据收入和年龄分组计算平均养老金购买比例
    '''
    # 定义收入分位数和年龄段的分组
    income_quantile_bins = np.array([
        [0, 0.1],
        [0.1, 0.2],    
        [0.2, 0.3],
        [0.3, 0.4],
        [0.4, 0.5],
        [0.5, 0.6],
        [0.6, 0.7],
        [0.7, 0.8],
        [0.8, 0.9],
        [0.9, 0.95],
        [0.95, 1]
    ])
    age_bins = np.array([
        [20, 25],
        [26, 30],
        [31, 35],
        [36, 40],
        [41, 45],
        [46, 50],
        [51, 55],
        [56, 60],
        [61, 65],
    ])
    # 计算每个收入分组和年龄段的平均行动
    mean_actions_by_income_and_age = []
    temp_income = income_history_df.iloc[:, 1:]
    pension_action = rl_actions[:, 2, :]
    pensin_amount = pen_amount_df.iloc[:, 1:].values
    mean_pension_by_income_age = np.full((len(income_quantile_bins),len(age_bins)), np.nan) 
    idx_age = 0
    for start_age, end_age in age_bins:
        age_mask = (unique_age >= start_age) & (unique_age <= end_age)
        mean_actions_by_income = []
        idx_income = 0
        for start_income, end_income in income_quantile_bins:
            
            agent_meanincome_by_age = temp_income.iloc[age_mask].mean(axis=0)
            agent_mask = (agent_meanincome_by_age.quantile(start_income) <= agent_meanincome_by_age) & (agent_meanincome_by_age < agent_meanincome_by_age.quantile(end_income))
            mean_pension_by_income_age[idx_income, idx_age] = np.nanmean(pensin_amount[age_mask][:, agent_mask])
            idx_income+=1
        idx_age+=1

    # 将mean_actions_by_age转换为DataFrame
    mean_pension_by_income_age = pd.DataFrame(mean_pension_by_income_age*100, columns=[f'{int(start)}-{int(end)}' for start, end in age_bins], index=[f'{int(start*100)}-{int(end*100)}' for start, end in income_quantile_bins])
    mean_pension_by_income_age.index.name = '收入分位数'

    ''' 
    （2）根据收入和年龄分组计算平均养老金风险资产购买比例
    '''
    mean_pension_risk_by_income_age = np.full((len(income_quantile_bins),len(age_bins)), np.nan) 
    pension_risk = rl_actions[:, 3, :]
    pension_risky_amount = pensin_amount * pension_risk
    idx_age = 0
    for start_age, end_age in age_bins:
        age_mask = (unique_age >= start_age) & (unique_age <= end_age)
        mean_actions_by_income = []
        idx_income = 0
        for start_income, end_income in income_quantile_bins:
            agent_meanincome_by_age = temp_income.iloc[age_mask].mean(axis=0)
            agent_mask = (agent_meanincome_by_age.quantile(start_income) <= agent_meanincome_by_age) & (agent_meanincome_by_age < agent_meanincome_by_age.quantile(end_income))
            mean_pension_risk_by_income_age[idx_income, idx_age] = np.nanmean(pension_risky_amount[age_mask][:, agent_mask])
            idx_income+=1
        idx_age+=1
    mean_pension_risk_by_income_age = pd.DataFrame(mean_pension_risk_by_income_age*100, columns=[f'{int(start)}-{int(end)}' for start, end in age_bins], index=[f'{int(start*100)}-{int(end*100)}' for start, end in income_quantile_bins])
    mean_pension_risk_by_income_age.index.name = '收入分位数'

    # 保存为excel

    '''
    将表格以图的形式展现，并加上根据值的热力图
    '''
    import seaborn as sns
    import matplotlib.pyplot as plt

    plt.rcParams['font.sans-serif'] = ['SimHei']  # 使用简体中文字体
    plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题


    # 创建热力图，第一个子图为养老金购买比例，第二个子图为养老金风险资产购买比例, 两个子图用不一样的颜色
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    sns.heatmap(mean_pension_by_income_age, annot=True, cmap='Greys', fmt='.2f', cbar=False, ax=axes[0])
    sns.heatmap(mean_pension_risk_by_income_age, annot=True, cmap='Greys', fmt='.2f', cbar=False, ax=axes[1])
    # 行为年龄段为横坐标，收入分位数为纵坐标
    axes[0].set_xlabel('工作期年龄段')
    axes[0].set_ylabel('收入分位数')
    axes[0].set_title(r'(a) 养老金购买额(万元)')
    axes[1].set_xlabel('工作期年龄段')
    axes[1].set_ylabel('')  # 第二个子图不显示ylabel和坐标
    axes[1].set_yticks([])
    axes[1].set_title(r'(b) 养老金风险资产购买额(万元)')
    # 将两个子图的水平距离拉近一些
    plt.subplots_adjust(wspace=0.01)
    plt.tight_layout()
    plt.savefig('figs/p5-mean_pension_amount_by_income_age.png', dpi=300)
    
    



