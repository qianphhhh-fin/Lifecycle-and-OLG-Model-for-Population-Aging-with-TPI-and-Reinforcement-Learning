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
from tqdm import tqdm
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
    id='pensionfund-v5-rural',                                # call it whatever you want
    entry_point='v5_pensionfund_rural:PensionFundEnv', # module_name:class_name
)
log_dir = '../models//sac//'

eval_dir = 'eval//'
# model_name = ['v13_pensionfund_penlim0_run7',
#               'v13_pensionfund_penlim1p2_run15',
#               'v13_pensionfund_nopenlim_run7'
#               ]
rl_model_name = 'pensionfund-v5_run50//best_model' # 
rl_vec_env_name = 'pensionfund-v5_run50//best_model_vec_env.pkl'


# model_name = [d for d in os.listdir(log_dir) if os.path.isdir(os.path.join(log_dir, d))]
eval_episodes = 200
eval_seed =   3687851522
distf = 0.95
import multiprocessing as mp
from functools import partial

def evaluate_model_naive(eval_episodes=10000, eval_seed=3687851522):
    # 读取cocco模型的A,C,V,gcash


    params = {'pension_limit':1}
    params = None

    env = make_vec_env("pensionfund-v5-rural",seed=eval_seed, n_envs=1,  env_kwargs={'params': params})

    # 读取matlab dp解
    current_path = os.path.dirname(os.path.abspath(__file__))
    age_group = str(env.get_attr('tb')[0])+'-'+str(env.get_attr('td')[0])
    df_A = pd.read_excel(os.path.join(current_path,'result_cocco_matlab',age_group,'A.xlsx'),header=None)
    df_C = pd.read_excel(os.path.join(current_path,'result_cocco_matlab',age_group,'C.xlsx'),header=None) 
    df_V = pd.read_excel(os.path.join(current_path,'result_cocco_matlab',age_group,'V.xlsx'),header=None)  
    df_gcash = pd.read_excel(os.path.join(current_path,'result_cocco_matlab',age_group,'gcash.xlsx'),header=None)
    # 转换为numpy数组
    A = df_A.to_numpy()
    C = df_C.to_numpy()
    V = df_V.to_numpy()
    gcash = df_gcash.to_numpy().squeeze()



    # env = make_vec_env("cocco-dp",seed=1, n_envs=1,  env_kwargs={'params': None})

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
        info = env.reset_infos
        done = False
        score = 0
        episode_actions = []
        temp_reward = []
        '''
        5个连续选择：消费比例，风险资产比例，购买城乡居民养老保险比例，购买个人养老金比例(若基本养老金个人账户余额不为0\则允许购买个人养老金)，个人养老金中购买风险资产比例
        '''
        while not done:
            # 根据state[0]获取现金
            normalized_cash = np.squeeze(info[0]['normalized_cash']) # 归一化的现金
            cash = obs[0][0]  # 现金
            age_loc =  env.get_attr('age')[0] - env.get_attr('tb')[0] # 获取年龄

            # 根据现金和年龄插值从A中获取action
            action = np.array([np.interp(normalized_cash, gcash, C[:, age_loc]),
                               np.interp(normalized_cash, gcash, A[:, age_loc]),0.2,0.1,0.5])
            if age_loc == np.shape(C)[1]-1: # 如果到达最后一期,则固定消费比例为1（全部消费干净）
                action[0] = normalized_cash
            else:
                action[0] = np.clip(action[0],0,normalized_cash)
            # action[0] = action[0]/cash # 转换为消费比例
            action[0] = action[0]/(normalized_cash) # 消费比例  
            # action[0] = 1 #          
            action[1] = np.clip(action[1],0,1) # 风险资产比例
            action[2] = np.clip(action[2],0,0) # 城乡居民养老保险购买比例
            action[3] = np.clip(action[3],0,0) # 个人养老金购买比例
            action[4] = np.clip(action[4],0,0) # 个人养老金风险资产比例 



            obs, reward, done, info = env.step(np.array([action]))
            score += reward * (distf ** len(episode_actions))

            episode_actions.append(action)
            temp_consum.append(cash*action[0])
            temp_reward.append(reward[0])
            # if info[0]['status']=='working':
            tmp_income.append(info[0]['basic_income'])

        age.append(info[0]['age'][0])
        score_history.append(score[0])
        action_set.append(episode_actions)
        reward_history.append(temp_reward)
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
    model_identifier = 'naive'

    os.makedirs(eval_dir + model_identifier + '//', exist_ok=True)


    reward_history_df = var_path(reward_history, eval_episodes, env)
    reward_history_df.to_excel(eval_dir+ model_identifier + "//" + model_identifier + "_reward_history.xlsx", index=False)

    life_summary = np.array([
        np.mean(age),
        np.mean(score_history),
        np.std(score_history),
        np.min(reward_history_df.iloc[:,1:]),
        np.max(reward_history_df.iloc[:,1:]),
        np.mean(reward_history_df.iloc[:,1:]),
    ])


    life_summary_df = pd.DataFrame([life_summary], 
                                columns=['死亡年龄', '终身效用','终身效用标准差','最小单期奖励','最大单期奖励','平均单期奖励'])
    # 在life_summary_df前面插入模型标识符列
    life_summary_df.insert(0, '模型', model_identifier)
    return life_summary_df

# def evalute_model_erl(idx, eval_episodes=1000, eval_seed=3687851522):


def evaluate_model_rl(model_name,vec_name, eval_episodes=1000, eval_seed=3687851522):
    current_path = os.path.dirname(os.path.abspath(__file__))    
    # try:
    with open(os.path.join(current_path, log_dir, model_name.split('//best_model')[0]) +"//params.json", 'r') as file:
        params = json.load(file)   
    # except:
    params['pension_limit'] = 1
    # params['sigma_eta'] = 0.27

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

    env = make_vec_env("pensionfund-v5-rural",seed=eval_seed, n_envs=1, env_kwargs={'params':params},\
                            monitor_kwargs={'discount_factor': distf})         
    env = VecNormalize(env, norm_obs=False, norm_reward=False) # 评估不需要norm_obs


    idx = model_name.split('//')[0]
    # model = DDPG.load(os.path.join(log_dir, idx) + "//latest_model", env)

    action_set = []
    score_history = []
    age = []
    income = []
    consum = []
    pension_pct = []
    tmp_income = []
    temp_consum = []
    temp_pension_pct = []
    reward_history = []
    # from tqdm import tqdm
    pbar = tqdm(total=eval_episodes, desc=f'模拟进度 {idx}')
    
    for episode in range(eval_episodes):       
        pbar.update(1)
        env.seed(seed=eval_seed+episode)
        obs = env.reset()    
        done = False
        score = 0
        episode_actions = []
        temp_reward = []
        while not done:
            action, _ = model.predict(model.env.normalize_obs(obs), deterministic=True) 
            # action[0,2] = 0
            cash = obs[0][0] #

            obs, reward, done, info = env.step(action)               
            # reward = -(((cash*action[0][0]))**(1-env.get_attr('gamma')[0]))*10 + 10       
            score += reward * (distf ** len(episode_actions))
            episode_actions.append(info[0]['real_actions'])
            temp_consum.append(cash*action[0][0])
            temp_pension_pct.append(info[0]['每年个人养老金购买比例(实际)'][0])
            temp_reward.append(reward[0])
            # if info[0]['status']=='working':
            tmp_income.append(info[0]['basic_income'])
        
        score_history.append(score[0])
        reward_history.append(temp_reward)
        action_set.append(episode_actions)
        # norm_age = info[0]['state'][4]
        # age.append(np.round((norm_age + 1) * (params['max_age'] - params['born_age']) / 2 + params['born_age'], 1) - 1)
        age.append(info[0]['age'][0])

        consum.append(temp_consum)
        pension_pct.append(temp_pension_pct)
        income.append(tmp_income)

        tmp_income = []
        temp_consum = []
        temp_pension_pct = []


    # 保存评估文件
    # 从模型名中提取标识符（第二个_之后的字符串）
    model_identifier = 'rl_' + idx.split('-', 1)[1] 

    os.makedirs(eval_dir + model_identifier + '//', exist_ok=True)


    life_summary = np.array([
        np.mean(age),
        np.mean(score_history),
        np.std(score_history),
    ])

    life_summary_df = pd.DataFrame([life_summary], 
                                columns=[ '死亡年龄', '终身效用','终身效用标准差'])

    # =====保存消费路径=====
    consum_path_df = var_path(consum, eval_episodes, env)
    consum_path_df.to_excel(eval_dir+ model_identifier + "//" + model_identifier + "_consum_history.xlsx", index=False)


    pension_pct_path_df = var_path(pension_pct, eval_episodes, env)
    pension_pct_path_df.to_excel(eval_dir+ model_identifier + "//" + model_identifier + "_pension_pct_history.xlsx", index=False)

    
        # 保存收入历史
    income_history_df = var_path(income, eval_episodes, env)
    income_history_df.to_excel(eval_dir+ model_identifier + "//" + model_identifier + "_income_history.xlsx", index=False)

    life_path_detailed = pd.DataFrame({
        '死亡年龄': age,
        '工作期平均收入(万)': np.mean(income_history_df.iloc[:,1:],axis=0),
        '终身效用': score_history,
    })
    
    life_path_detailed_path = eval_dir+ model_identifier + "//" + model_identifier + "_life_path_detailed.xlsx"
    life_path_detailed.to_excel(life_path_detailed_path, index=True)

    # income = income_history_df.iloc[:env.get_attr('tr')[0],1:]
    # income_summary_df = pd.DataFrame([np.array([
    #     np.mean(income),
    #     np.nanstd(income.values.flatten()),
    #     np.min(income),
    #     np.max(income),
    #     np.quantile(income,)
    # ])],columns=['工作期平均收入(w)', '标准差','最小值','最大值'])
    # print(income_summary_df)
       # =====保存奖励历史=====
    reward_history_df = var_path(reward_history, eval_episodes, env)
    reward_history_df.to_excel(eval_dir+ model_identifier + "//" + model_identifier + "_reward_history.xlsx", index=False)

    life_summary = np.array([
        np.mean(age),
        np.mean(score_history),
        np.std(score_history),
    ])

    life_path_detailed_path = eval_dir+ model_identifier + "//" + model_identifier + "_life_path_detailed.xlsx"
    life_path_detailed.to_excel(life_path_detailed_path, index=False)


    life_summary_df = pd.DataFrame([life_summary], 
                                columns=[ '死亡年龄', '终身效用','终身效用标准差'])
    life_summary_path = eval_dir+ model_identifier + "//" + model_identifier + "_life_summary.xlsx"
    life_summary_df.to_excel(life_summary_path, index=False)

    max_length = max(len(path) for path in consum)
    unique_age = np.arange(env.get_attr('tb')[0], env.get_attr('td')[0] + 1)[:max_length]
    max_length = max(len(arr) for arr in action_set)
    num_action = np.shape(action_set[0])[1]
    aligned_actions = np.full((max_length, num_action, eval_episodes), np.nan)
    for i, actions in enumerate(action_set):
        aligned_actions[:len(actions), :, i] = actions

    mean_action_by_age = np.nanmean(aligned_actions, axis=2)
    mean_actions = np.column_stack((unique_age, mean_action_by_age))
    mean_actions_df = pd.DataFrame(mean_actions, 
                                    columns=['年龄', '消费', '风险资产','个人养老金购买比例','个人养老金风险资产比例'])
    mean_actions_path = eval_dir+ model_identifier + "//" + model_identifier + "_mean_actions.xlsx"
    mean_actions_df.to_excel(mean_actions_path, index=False)
  
    # 在life_summary_df前面插入模型标识符列
    life_summary_df.insert(0, '模型', model_identifier)
    return life_summary_df
    # print(f"模型{idx}评估完毕, 评估次数为{eval_episodes}")


if __name__ == '__main__':


    # 创建进程池
    # pool = mp.Pool()
    
    # # 循环执行evaluate_model
    # from tqdm import tqdm
    life_summary_naive = evaluate_model_naive(eval_episodes=eval_episodes, eval_seed=eval_seed)
    print(life_summary_naive.to_string(index=False))
    # life_summary_rl = evaluate_model_rl(rl_model_name,rl_vec_env_name,eval_episodes=eval_episodes, eval_seed=eval_seed)
    # print(life_summary_rl.to_string(index=False))
    
    # life_summary = pd.concat([life_summary_naive, life_summary_rl], ignore_index=True)
    # pd.set_option('display.float_format', lambda x: '{:.8f}'.format(x))
    # pd.set_option('display.colheader_justify', 'center')
    # pd.set_option('display.unicode.ambiguous_as_wide', True)
    # pd.set_option('display.unicode.east_asian_width', True)
    # print(life_summary.to_string(index=False))


