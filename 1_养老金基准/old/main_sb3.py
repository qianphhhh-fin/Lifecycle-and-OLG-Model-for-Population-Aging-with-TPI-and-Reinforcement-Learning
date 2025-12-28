from typing import Callable
import json
import gymnasium as gym
import argparse  # 用于解析命令行参数
import numpy as np
import os
# from stable_baselines3 import DDPG
from stable_baselines3 import TD3
from stable_baselines3 import PPO
from stable_baselines3 import DDPG
from stable_baselines3 import SAC
from sb3_contrib import ARS
from sb3_contrib.common.vec_env import AsyncEval
from stable_baselines3.common.noise import NormalActionNoise, OrnsteinUhlenbeckActionNoise
# from stable_baselines3.common.logger import TensorBoardOutputFormat
from utils.env_util import make_vec_env
from stable_baselines3.common.vec_env import VecNormalize
# from stable_baselines3.common.callbacks import EvalCallback,BaseCallback
from utils.cocco_AllocationEvalCallback import EvalCallback
from gymnasium.envs.registration import register
import pandas as pd
import matplotlib.pyplot as plt
import sys
# import torch 
# torch.autograd.set_detect_anomaly(True)
import warnings
warnings.filterwarnings('ignore')
# 加载PyTorch并设置线程数
import torch
# torch.set_num_threads(2)  # 不设置这个多进程会变得非常慢
# os.environ['OPENBLAS_NUM_THREADS'] = '2'  # 在某些架构上不设置这个可能会导致Numpy导入错误


cdir = os.path.abspath('.') # 获取当前目录
# print(sys.path)
# 状态空间：现金余额，个人养老金账户余额，永久收入冲击，年龄
register(
    id='cocco-rl',
    entry_point='cocco_env_rl:CoccoEnv',
)

def nonlin_schedule(initial_value: float) -> Callable[[float], float]:
    def func(progress_remaining: float) -> float:
        """
        Progress will decrease from 1 (beginning) to 0.
        使用指数函数实现学习率的非线性衰减，开始时衰减缓慢，后期衰减加快。

        :param progress_remaining: 剩余进度，从1减小到0
        :return: 当前的学习率
        """
        # 使用10次幂函数实现非线性衰减
        # 当progress_remaining=1时,学习率为initial_value
        # 当progress_remaining接近0时,学习率迅速衰减到0
        return -initial_value * (progress_remaining - 1)**10 + initial_value
    return func



def run(env_id):
    
    log_dir = 'models//cocco_sac//'
    # log_dir = 'models//sb3_td3_//'
    # log_dir = 'models//sb3_ddpg_//'
    # log_dir = 'models//sb3_ppo_//'
    log_env_id = env_id.replace('cocco-v', 'v') + '_cocco'
    model_name = log_env_id  + '_run'
   
    # 检查并创建日志目录
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    # 获取目录中所有文件夹
    existing_dirs = [d for d in os.listdir(log_dir) if os.path.isdir(os.path.join(log_dir, d)) and d.startswith(model_name)] 
    # 获取现有run编号
    run_numbers = []
    for d in existing_dirs:
        try:
            run_num = int(d.split('run')[-1])
            run_numbers.append(run_num) 
        except:
            continue    
    # 确定新的run编号
    next_run = max(run_numbers) + 1 if run_numbers else 0
    model_name = f"{model_name}{next_run}"  
    # 创建模型目录
    model_dir = os.path.join(log_dir, model_name)
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)  
    
    params = {'distf': 0.95}
    vec_env = make_vec_env("cocco-rl",seed=1, n_envs=1, env_kwargs={'params': None},\
                            monitor_kwargs={'discount_factor': params['distf']})
    vec_env = VecNormalize(vec_env, norm_obs=True, norm_reward=True)

    eval_env = make_vec_env("cocco-rl",seed=3687851522, n_envs=1, env_kwargs={'params': None},\
                            monitor_kwargs={'discount_factor': params['distf']})         
    eval_env = VecNormalize(eval_env, norm_obs=True, norm_reward=False)

    eval_callback = EvalCallback(eval_env, best_model_save_path=os.path.join(log_dir, model_name),\
                              eval_freq=5000, n_eval_episodes = 1000, params = params,\
                              deterministic=True, render=False, discount_factor=params['distf'])

    model = SAC("MlpPolicy", vec_env , verbose=1, gamma=params['distf'], device="cpu",\
                tensorboard_log=".//runs//",
                    # gradient_steps=4,
                    # learning_rate=nonlin_schedule(0.003),
                    # learning_rate=0.003,
                    # policy_kwargs={'activation_fn':torch.nn.Tanh},
                    action_noise = OrnsteinUhlenbeckActionNoise(mean=np.zeros(vec_env.action_space.shape[0]), sigma=0.5 * np.ones(vec_env.action_space.shape[0])),
                    stats_window_size=1000,
                    # ent_coef=0.3,
                    # batch_size=512,
                    # use_sde=True,
                    # batch_size=1024,
                    learning_starts=1000)
    # model = DDPG("MlpPolicy", vec_env , verbose=1, gamma=params['distf'], device="cpu",\
    #         tensorboard_log=".//runs//",
    #         action_noise = OrnsteinUhlenbeckActionNoise(mean=np.zeros(vec_env.action_space.shape[0]), sigma=0.5 * np.ones(vec_env.action_space.shape[0])),
    #         learning_starts=1000)
    model.learn(total_timesteps=1_000_000, log_interval=10,callback=eval_callback,progress_bar=True)

    # model = ARS("MlpPolicy", vec_env,tensorboard_log=".//runs//",delta_std=0.01)
    # model.learn(total_timesteps=1_000_000, log_interval=4,callback=eval_callback)

    # model = TD3("MlpPolicy", vec_env , verbose=1, gamma=params['distf'], device="cpu",\
    #         tensorboard_log=".//runs//",
    #             # gradient_steps=4,
    #             action_noise = OrnsteinUhlenbeckActionNoise(mean=np.zeros(vec_env.action_space.shape[0]), sigma=0.1 * np.ones(vec_env.action_space.shape[0])),
    #             stats_window_size=1000,
    #             learning_starts=1000)






    
    

    vec_env.save(os.path.join(log_dir, model_name)+"//vec_env.pkl") 
    # 保存params
    with open(os.path.join(log_dir, model_name)+"//params.json", 'w') as file:
        json.dump(params, file)

def parse_args(args=None):
    # 解析命令行参数
    parser = argparse.ArgumentParser("Experiments")
    # 通用参数:
    parser.add_argument('--envid', default="cocco-rl", type=str, metavar='ENVID', help="环境 (默认: %(default)s)")
    
    params = parser.parse_args(args)

    return params

if __name__ == "__main__":
    params = parse_args()
    # 将参数解包传递给run函数
    run( params.envid)