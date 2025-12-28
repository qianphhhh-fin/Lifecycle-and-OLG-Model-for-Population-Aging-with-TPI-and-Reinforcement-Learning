from typing import Callable
import json
import gymnasium as gym
import numpy as np
import os
from tqdm import tqdm
# from stable_baselines3 import DDPG
# from stable_baselines3 import TD3
# from stable_baselines3 import PPO
# from stable_baselines3 import SAC
from sbx import DDPG, DQN, PPO, SAC, TD3, TQC, CrossQ
from utils.env_util import make_vec_env
from stable_baselines3.common.vec_env import VecNormalize
from gymnasium.envs.registration import register
import pandas as pd
import matplotlib.pyplot as plt
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
# get current directory of the current .ipynb 
cdir = os.path.abspath('.') # 获取当前目录
pardir = os.path.abspath(os.path.join(os.getcwd(), "../")) # 获取上级目录
dbdir = os.path.abspath(os.path.join(os.getcwd(), "../..")) # 获取上上级目录(为了读取美股数据库数据)

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

# import torch 
# torch.autograd.set_detect_anomaly(True)
import warnings
warnings.filterwarnings('ignore')


import sys
cdir = os.path.abspath('.') # 获取当前目录
# print(sys.path)
# 状态空间：现金余额，个人养老金账户余额，永久收入冲击，年龄
register(
    id='cocco-rl-norm',
    entry_point='cocco_env_rl_norm:CoccoEnv',
)
register(
    id='cocco-dp',
    entry_point='cocco_env_dp:CoccoEnv',
)
log_dir = 'models//cocco_sac'

eval_dir = 'eval//'
# model_name = ['v13_pensionfund_penlim0_run7',
#               'v13_pensionfund_penlim1p2_run15',
#               'v13_pensionfund_nopenlim_run7'
#               ]
rl_model_name = [
            #   'cocco-rl_cocco_run225//best_model_67', # 无随机死亡最佳
            # 'cocco-rl_cocco_run231//best_model_67', # 有随机死亡最佳 
            # 'cocco-rl_cocco_run248//best_model_567'
            'cocco-rl-norm_cocco_run112//best_model',
            # 'cocco-rl-norm_cocco_run11//best_model'
              ]
# model_name = [d for d in os.listdir(log_dir) if os.path.isdir(os.path.join(log_dir, d))]
eval_episodes = 200
eval_seed =  3687851522
# eval_seed = 1000
params = {'distf': 0.95}



def evaluate_model_dp(eval_episodes=10000, eval_seed=3687851522):

    current_path = os.path.dirname(os.path.abspath(__file__))
    df_A = pd.read_excel(os.path.join(current_path,'result_cocco_matlab','A.xlsx'),header=None)
    df_C = pd.read_excel(os.path.join(current_path,'result_cocco_matlab','C.xlsx'),header=None) 
    df_V = pd.read_excel(os.path.join(current_path,'result_cocco_matlab','V.xlsx'),header=None)  
    df_gcash = pd.read_excel(os.path.join(current_path,'result_cocco_matlab','gcash.xlsx'),header=None)
    # 转换为numpy数组
    A = df_A.to_numpy()
    C = df_C.to_numpy()
    V = df_V.to_numpy()
    gcash = df_gcash.to_numpy().squeeze()

    env = make_vec_env("cocco-dp",seed=1, n_envs=1,  env_kwargs={'params': None})

    action_set = []
    score_history = []
    age = []
    income = []
    consum = []
    tmp_income = []
    temp_consum = []
    
    reward_history = []

    # from tqdm import tqdm
    pbar = tqdm(total=eval_episodes, desc=f'dp模拟进度')    
    for episode in range(eval_episodes):
        pbar.update(1)
        env.seed(seed=eval_seed+episode)
        obs = env.reset()    
        done = False
        score = 0
        episode_actions = np.empty((0,2))
        temp_reward = []
        while not done:
            # 根据state[0]获取现金
            cash = np.squeeze(obs[0][0])        
            age_loc =  int(obs[0][1]) - env.get_attr('tb')[0] # 获取年龄

            # 根据现金和年龄插值从A中获取action
            action = np.array([np.interp(cash, gcash, C[:, age_loc]),np.interp(cash, gcash, A[:, age_loc])])
            if age_loc == np.shape(C)[1]-1: # 如果到达最后一期,则固定消费比例为1（全部消费干净）
                action[0] = cash
            else:
                action[0] = np.clip(action[0],0,cash)
            # action[0] = action[0]/cash # 转换为消费比例
            action[1] = np.clip(action[1],0,1)


            obs, reward, done, info = env.step(np.array([action]))
            score += reward * (params['distf'] ** len(episode_actions))
            action[0] = action[0]/cash # 转换为消费比例
            episode_actions = np.vstack((episode_actions, action))
            temp_consum.append(info[0]['raw_consumption'])
            temp_reward.append(reward[0])
            # if info[0]['status']=='working':
            tmp_income.append(info[0]['raw_income'])

        age.append(info[0]['state'][1])
        score_history.append(score[0])
        action_set.append(episode_actions)
        reward_history.append(temp_reward)
        # norm_age = info[0]['state'][4]
        # age.append(np.round((norm_age + 1) * (params['max_age'] - params['born_age']) / 2 + params['born_age'], 1) - 1)
        

        consum.append(temp_consum)
        income.append(tmp_income)


        tmp_income = []
        temp_consum = []


    life_path_detailed = pd.DataFrame({
        '死亡年龄': age,
        '工作期平均收入(万)': income,
        '终身效用': score_history,
    })
    # 保存评估文件
    # 从模型名中提取标识符（第二个_之后的字符串）
    model_identifier = 'dp'

    os.makedirs(eval_dir + model_identifier + '//', exist_ok=True)

    # =====保存消费路径=====
    consum_path_df = var_path(consum, eval_episodes, env)
    consum_path_df.to_excel(eval_dir+ model_identifier + "//" + model_identifier + "_consum_history.xlsx", index=False)
   # =====保存奖励历史=====
    reward_history_df = var_path(reward_history, eval_episodes, env)
    reward_history_df.to_excel(eval_dir+ model_identifier + "//" + model_identifier + "_reward_history.xlsx", index=False)
    # 保存收入历史
    income_history_df = var_path(income, eval_episodes, env)
    income_history_df.to_excel(eval_dir+ model_identifier + "//" + model_identifier + "_income_history.xlsx", index=False)
    # 保存风险资产比例决策历史
    risky_pct = []
    for episode in range(eval_episodes):
        risky_pct.append(action_set[episode][:,1])
    risky_pct_history_df = var_path(risky_pct, eval_episodes, env)
    risky_pct_history_df.to_excel(eval_dir+ model_identifier + "//" + model_identifier + "_risky_pct_history.xlsx", index=False)

        # 保存终身效用
    utility_path = eval_dir + model_identifier + "//" + model_identifier + "_score_history.xlsx"
    utility_path_df = pd.DataFrame(score_history, columns=['终身实现效用'])
    utility_path_df.insert(0, 'agent', np.arange(1, eval_episodes+1)) # 增加第一列为agent
    utility_path_df.to_excel(eval_dir+ model_identifier + "//" + model_identifier + "_score_history.xlsx", index=False)

    life_summary = np.array([
        np.mean(age),
        np.mean(score_history),
        np.std(score_history),
    ])
    life_summary_df = pd.DataFrame([life_summary], 
                                columns=[ '死亡年龄', '终身效用','终身效用标准差'])
    print(model_identifier,': ',life_summary_df)


    

def evaluate_model_rl(idx, eval_episodes=1000, eval_seed=3687851522):
    # # with open(os.path.join(log_dir, idx) +"//params.json", 'r') as file:
    # #     params = json.load(file)
    # current_path = os.path.dirname(os.path.abspath(__file__))    
    # # env = make_vec_env("cocco-rl",seed=1, n_envs=16, env_kwargs={'params': None})
    # # env = VecNormalize.load(os.path.join(current_path, log_dir,idx) + "//vec_env.pkl", env)
    # model = SAC.load(os.path.join(current_path, log_dir, idx), env=None) # , print_system_info=True

    # env = make_vec_env("cocco-rl",seed=eval_seed, n_envs=1, env_kwargs={'params': None})
    # idx = idx.split('//')[0]
    # env = VecNormalize.load(os.path.join(current_path, log_dir, idx) + "//eval_env.pkl",env)
    # # env = make_vec_env("cocco-rl",seed=eval_seed, n_envs=1, env_kwargs={'params':None},\
    # #                         monitor_kwargs={'discount_factor': params['distf']})         
    # # env = VecNormalize(env, norm_obs=True, norm_reward=False)
    # # model = DDPG.load(os.path.join(log_dir, idx) + "//latest_model", env)

    current_path = os.path.dirname(os.path.abspath(__file__))     
   
    env = make_vec_env("cocco-rl-norm",seed=1, n_envs=1, env_kwargs={'params': None})
    env = VecNormalize.load(os.path.join(current_path, log_dir,idx.split('//')[0]) + "//best_model_vec_env.pkl", env)
    model = SAC.load(os.path.join(current_path, log_dir, idx), env=env)

    env = make_vec_env("cocco-rl-norm",seed=eval_seed, n_envs=1, env_kwargs={'params':None},\
                            monitor_kwargs={'discount_factor':params['distf']})         
    env = VecNormalize(env, norm_obs=False, norm_reward=False) # 评估不需要norm_obs


    idx = idx.split('//')[0]
    

    action_set = []
    score_history = []
    age = []
    income = []
    consum = []
    tmp_income = []
    temp_consum = []
    reward_history = []
    # from tqdm import tqdm
    pbar = tqdm(total=eval_episodes, desc=f'rl模拟进度')  
    for episode in range(eval_episodes):  
        pbar.update(1)
        env.seed(seed=eval_seed+episode)
        obs = env.reset()    
        done = False
        score = 0
        episode_actions = np.empty((0,2))
        temp_reward = []
        while not done:
            action, _states = model.predict(model.env.normalize_obs(obs), deterministic=True)
            # 如果为最后一期，则固定消费比例为1（全部消费干净）
            if len(episode_actions) == env.get_attr('tn')[0]-1:
                action[0] = 1
            # cash = env.unnormalize_obs(obs[0])[0] # 反归一化
            obs, reward, done, info = env.step(action)               
            # reward = -(((cash*action[0][0]))**(1-env.get_attr('gamma')[0]))*10 + 10       
            score += reward * (params['distf'] ** len(episode_actions))
            episode_actions = np.vstack((episode_actions, info[0]['real_actions']))
            temp_consum.append(info[0]['消费'])
            temp_reward.append(reward[0])
            tmp_income.append(info[0]['raw_income'])
        
        score_history.append(score)
        reward_history.append(temp_reward)
        action_set.append(episode_actions)
        # norm_age = info[0]['state'][4]
        # age.append(np.round((norm_age + 1) * (params['max_age'] - params['born_age']) / 2 + params['born_age'], 1) - 1)
        age.append(info[0]['state'][1])

        consum.append(temp_consum)
        income.append(tmp_income)


        tmp_income = []
        temp_consum = []


    life_path_detailed = pd.DataFrame({
        '死亡年龄': age,
        '工作期平均收入(万)': income,
        '终身效用': score_history,
    })
    # 保存评估文件
    # 从模型名中提取标识符（第二个_之后的字符串）
    model_identifier = 'rl_' +idx.split('_', 2)[2] 

    os.makedirs(eval_dir + model_identifier + '//', exist_ok=True)

# =====保存消费路径=====
    consum_path_df = var_path(consum, eval_episodes, env)
    consum_path_df.to_excel(eval_dir+ model_identifier + "//" + model_identifier + "_consum_history.xlsx", index=False)
    # 计算各个分位数
    # consum_path_df = consum_path_df.iloc[:,]
    # consum_stats = pd.DataFrame({
    #     '均值': consum_path_df.mean(),
    #     '0.1%': consum_path_df.quantile(0.001),
    #     '1%': consum_path_df.quantile(0.01),
    #     '5%': consum_path_df.quantile(0.05),
    #     '10%': consum_path_df.quantile(0.1),
    #     '25%': consum_path_df.quantile(0.25),
    #     '50%': consum_path_df.quantile(0.5),
    #     '75%': consum_path_df.quantile(0.75),
    #     '90%': consum_path_df.quantile(0.9),
    #     '95%': consum_path_df.quantile(0.95),
    #     '99%': consum_path_df.quantile(0.99),
    #     '99.9%': consum_path_df.quantile(0.999)
    # })
    # # 保存统计结果
    # print(consum_stats)
    
   # =====保存奖励历史=====
    reward_history_df = var_path(reward_history, eval_episodes, env)
    reward_history_df.to_excel(eval_dir+ model_identifier + "//" + model_identifier + "_reward_history.xlsx", index=False)
    # 保存收入历史
    income_history_df = var_path(income, eval_episodes, env)
    income_history_df.to_excel(eval_dir+ model_identifier + "//" + model_identifier + "_income_history.xlsx", index=False)
    # 保存风险资产比例决策历史
    risky_pct = []
    for episode in range(eval_episodes):
        risky_pct.append(action_set[episode][:,1])
    risky_pct_history_df = var_path(risky_pct, eval_episodes, env)
    risky_pct_history_df.to_excel(eval_dir+ model_identifier + "//" + model_identifier + "_risky_pct_history.xlsx", index=False)
     # 保存终身效用
    utility_path = eval_dir + model_identifier + "//" + model_identifier + "_score_history.xlsx"
    utility_path_df = pd.DataFrame(score_history, columns=['终身实现效用'])
    utility_path_df.insert(0, 'agent', np.arange(1, eval_episodes+1)) # 增加第一列为agent
    utility_path_df.to_excel(eval_dir+ model_identifier + "//" + model_identifier + "_score_history.xlsx", index=False)
    

    # # =====保存生命总结=====
    life_summary = np.array([
        np.mean(age),
        np.mean(score_history),
        np.std(score_history),
    ])
    

    life_summary_df = pd.DataFrame([life_summary], 
                                columns=[ '死亡年龄', '终身效用','终身效用标准差'])
    print(idx,': ',life_summary_df)
 

if __name__ == '__main__':
 
    # # 循环执行evaluate_model
    evaluate_model_dp(eval_episodes=eval_episodes, eval_seed=eval_seed)
    for model_name in rl_model_name:
        evaluate_model_rl(model_name, eval_episodes=eval_episodes, eval_seed=eval_seed)
    # evaluate_model_rl(rl_model_name[0], eval_episodes=eval_episodes, eval_seed=eval_seed)
