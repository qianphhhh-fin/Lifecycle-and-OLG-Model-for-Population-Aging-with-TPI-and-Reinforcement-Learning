from typing import Callable
import json
import gymnasium as gym
import numpy as np
import os
from sbx import SAC
# from utils.sac import SAC
from utils.env_util import make_vec_env
from stable_baselines3.common.vec_env import VecNormalize
from gymnasium.envs.registration import register
import pandas as pd
import matplotlib.pyplot as plt
import pickle
from tqdm import tqdm
# import matplotlib
# print(matplotlib.matplotlib_fname())
# raise
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# import torch 
# torch.autograd.set_detect_anomaly(True)
import warnings
warnings.filterwarnings('ignore')
'''
pensionfund-v5-rural: 状态变量：财富，上期工资持续冲击，基本养老金个人账户余额, 城乡保缴费年限, 个人养老金账户余额，年龄，
'''
register(
id='pensionfund-v5-rural',                                # call it whatever you want
entry_point='v5_pensionfund_rural:PensionFundEnv', # module_name:class_name
)

''''
pensionfund-v5-nonrural: 状态变量：财富，上期工资持续冲击，个人养老金账户余额，年龄
'''
register(
    id='pensionfund-v5-nonrural',                                # call it whatever you want
    entry_point='v5_pensionfund_nonrural:PensionFundEnv', # module_name:class_name
)

distf = 0.95
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

def evaluate_model_rl(model_path,env_name,eval_episodes=1000, eval_seed=3687851522,params=None):


    # try:
    #     with open(os.path.join(current_path, log_dir, model_name.split('//')[0]) +"//params.json", 'r') as file:
    #         params = json.load(file)   
    # except:
    # if 'cocco' in idx:
    #     env = make_vec_env("cocco-rl-norm",seed=1, n_envs=1, env_kwargs={'params': params})
    # else:
    # env = make_vec_env("pensionfund-v3",seed=1, n_envs=1, env_kwargs={'params': params})````````
    env = make_vec_env(env_name,seed=1, n_envs=1, env_kwargs={'params':params},\
                            monitor_kwargs={'discount_factor':distf})
    

    env = VecNormalize.load(os.path.join(model_path,"best_model_vec_env.pkl"), env)
    # env.obs_rms.mean
    # env.obs_rms.var
    
    model = SAC.load(os.path.join(model_path, "best_model"), env=env)

    

    # env = make_vec_env("pensionfund-v1",seed=eval_seed, n_envs=1, env_kwargs={'params': None})
    
    # env = VecNormalize.load(os.path.join(current_path, log_dir, idx) + "//eval_env.pkl",env)

    env = make_vec_env(env_name,seed=eval_seed, n_envs=1, env_kwargs={'params':params},\
                            monitor_kwargs={'discount_factor': distf})         
    env = VecNormalize(env, norm_obs=False, norm_reward=False) # 评估不需要norm_obs


    idx = model_path.split('//')[0]
    # model = DDPG.load(os.path.join(log_dir, idx) + "//latest_model", env)



    # from tqdm import tqdm
    pbar = tqdm(total=eval_episodes, desc=f'模拟进度 {idx}')
    all_trajectories = []
    for episode in range(eval_episodes):       
        pbar.update(1)
        env.seed(seed=eval_seed+episode)
        obs = env.reset()    
        info = env.envs[0].env.env.env.info
        done = False
        agent_trajectory = {
            'raw_income':[],
            'states': [obs],
            'actions': [],
            'rewards': [],
            'infos': [info],
            'consumption': [],
            'norm_age':[obs[-1]],
            'tda_pen_amount':[],
            'basic_pen_amount':[],
            'tax':[],
            'tda_tax_savings':[],
            'basic_tax_savings':[],
            'tda_tax':[],
            'age':[info['age']],
            'cash':[obs[0][0]],
            'tda_pen_balance':[info['个人养老金实际余额'][0]],
            'basic_pen_balance':[info['基本养老金个人账户余额'][0] if '基本养老金个人账户余额' in info else np.nan]
        }
        while not done:
            action, _ = model.predict(model.env.normalize_obs(obs), deterministic=True) 
            # action[0,2] = 0

            obs, reward, done, info = env.step(action)               
            # 记录轨迹
            # if info[0]['status'] == 'working':
            agent_trajectory['raw_income'].append(info[0]['basic_income'])
            agent_trajectory['age'].append(info[0]['age'][0])

            agent_trajectory['states'].append(obs)
            agent_trajectory['actions'].append(info[0]['real_actions'])
            agent_trajectory['rewards'].append(reward)
            agent_trajectory['consumption'].append(info[0]['消费'])
            agent_trajectory['infos'].append(info[0])
            agent_trajectory['norm_age'].append(obs[0][-1])
            agent_trajectory['cash'].append(obs[0][0])
            agent_trajectory['tda_pen_amount'].append(info[0]['每年个人养老金购买额(实际)'][0] if '每年个人养老金购买额(实际)' in info[0] else np.nan)
            agent_trajectory['basic_pen_amount'].append(info[0]['每年基本养老金购买额(实际)'][0] if '每年基本养老金购买额(实际)' in info[0] else np.nan)
            agent_trajectory['tda_pen_balance'].append(info[0]['个人养老金实际余额'][0])
            agent_trajectory['basic_pen_balance'].append(info[0]['基本养老金个人账户余额'][0] if '基本养老金个人账户余额' in info[0] else np.nan)
            agent_trajectory['tax'].append(info[0]['每年收入应缴税额(万元)'][0] if '每年收入应缴税额(万元)' in info[0] else np.nan)
            agent_trajectory['tda_tax_savings'].append(info[0]['个人养老金退税额'][0] if '个人养老金退税额' in info[0] else np.nan)
            agent_trajectory['basic_tax_savings'].append(info[0]['基本养老金退税额'][0] if '基本养老金退税额' in info[0] else np.nan)
            agent_trajectory['tda_tax'].append(info[0]['当期个人养老金缴税额'][0] if '当期个人养老金缴税额' in info[0] else np.nan)

        all_trajectories.append(agent_trajectory)
    pbar.close()
        
    # 计算所有agent的reward统计信息
    all_age = np.concatenate([traj['age'] for traj in all_trajectories])
    # import pandas as pd
    # age_income = pd.DataFrame({'age':all_age,'income':all_income})
    # mean_income_by_age = age_income.groupby('age')['income'].mean()
    action_set = [traj['actions'] for traj in all_trajectories]
    max_length = max(len(path) for path in action_set)
    unique_age = np.arange(env.get_attr('tb')[0], env.get_attr('td')[0] + 1)[:max_length]
    num_action = np.shape(action_set[0])[1]
    aligned_actions = np.full((max_length, num_action, eval_episodes), np.nan)
    for i, actions in enumerate(action_set):
        aligned_actions[:len(actions), :, i] = actions


    data = {'income_history_df':var_path(list([traj['raw_income'] for traj in all_trajectories]), eval_episodes, env),
            'wealth_history_df':var_path(list([traj['cash'] for traj in all_trajectories]), eval_episodes, env),
            'tda_pen_balance_history_df':var_path(list([traj['tda_pen_balance'] for traj in all_trajectories]), eval_episodes, env),
            'basic_pen_balance_history_df':var_path(list([traj['basic_pen_balance'] for traj in all_trajectories]), eval_episodes, env),
          'aligned_actions':aligned_actions,
          'consum_history_df':var_path(list([traj['consumption'] for traj in all_trajectories]), eval_episodes, env),
          'tda_pen_amount_history_df':var_path(list([traj['tda_pen_amount'] for traj in all_trajectories]), eval_episodes, env),
          'basic_pen_amount_history_df':var_path(list([traj['basic_pen_amount'] for traj in all_trajectories]), eval_episodes, env),
          'tax_history_df':var_path(list([traj['tax'] for traj in all_trajectories]), eval_episodes, env),
          'tda_tax_savings_history_df':var_path(list([traj['tda_tax_savings'] for traj in all_trajectories]), eval_episodes, env),
          'basic_tax_savings_history_df':var_path(list([traj['basic_tax_savings'] for traj in all_trajectories]), eval_episodes, env),
          'tda_tax_history_df':var_path(list([traj['tda_tax'] for traj in all_trajectories]), eval_episodes, env),
          'unique_age':unique_age}
    # 保存pension_balance_history_df
    # model_identifier = 'rl_' + idx.split('-', 1)[1] 
    # os.makedirs(eval_path + model_identifier + '//', exist_ok=True)
    
    # pension_balance_history_df.to_excel(eval_path + model_identifier + '//' + model_identifier + "_pension_balance_history.xlsx",index=False)
    
    return data








if __name__ == '__main__':

    env_model_path_set = {'pensionfund-v5-rural':'models//rural//pensionfund-v5-rural_run13//',
                        'pensionfund-v5-nonrural':'models//nonrural//pensionfund-v5-nonrural_run83//'}

    # 创建进程池
    # pool = mp.Pool()
    
    # # 循环执行evaluate_model
    eval_episodes = 10000
    eval_seed =   3687851522


    rural_data = evaluate_model_rl(env_model_path_set['pensionfund-v5-rural'],'pensionfund-v5-rural',eval_episodes=eval_episodes, eval_seed=eval_seed)
    nonrural_data = evaluate_model_rl(env_model_path_set['pensionfund-v5-nonrural'],'pensionfund-v5-nonrural',eval_episodes=eval_episodes, eval_seed=eval_seed)

    
    save_data = {'rural_data':rural_data,
                 'nonrural_data':nonrural_data}
    with open('simu_data/1w_data.pkl', 'wb') as f:
        pickle.dump(save_data, f)