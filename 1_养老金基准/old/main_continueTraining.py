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
torch.set_num_threads(2)  # 不设置这个多进程会变得非常慢
os.environ['OPENBLAS_NUM_THREADS'] = '2'  # 在某些架构上不设置这个可能会导致Numpy导入错误

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

num_episodes_to_continue = 100000
rl_model_name = 'cocco-rl_cocco_run47'
log_dir = 'models//cocco_sac'
eval_episodes = 1000
eval_seed = 3687851522
params = {'distf': 0.95}

current_path = os.path.dirname(os.path.abspath(__file__))    
env = make_vec_env("cocco-rl",seed=1, n_envs=1, env_kwargs={'params': None},\
                        monitor_kwargs={'discount_factor': params['distf']})         
env = VecNormalize.load(os.path.join(current_path, log_dir,rl_model_name) + "//vec_env.pkl", env)
model = SAC.load(os.path.join(current_path, log_dir, rl_model_name) + "//latest_model", env)

eval_env = make_vec_env("cocco-rl",seed=eval_seed, n_envs=1, env_kwargs={'params': None},\
                        monitor_kwargs={'discount_factor': params['distf']})         
eval_env = VecNormalize(eval_env, norm_obs=True, norm_reward=False)
eval_callback = EvalCallback(eval_env, best_model_save_path=os.path.join(log_dir, rl_model_name),\
                            eval_freq=5000, n_eval_episodes = 1000, params = params,\
                            deterministic=True, render=False, discount_factor=params['distf'])

model.learn(total_timesteps=num_episodes_to_continue, callback=eval_callback)