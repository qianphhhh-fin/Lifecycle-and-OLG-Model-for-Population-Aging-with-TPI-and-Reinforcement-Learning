from typing import Callable
import json
import gymnasium as gym
import argparse  # 用于解析命令行参数
import numpy as np
import os
# from stable_baselines3 import DDPG
import gymnasium as gym
from sbx import DDPG, DQN, PPO, SAC, TD3, TQC, CrossQ
from stable_baselines3.common.type_aliases import Schedule
from stable_baselines3.common.noise import NormalActionNoise, OrnsteinUhlenbeckActionNoise
from stable_baselines3.common.vec_env import VecCheckNan
# from stable_baselines3.common.logger import TensorBoardOutputFormat
from utils.env_util import make_vec_env
from utils.vec_normalize import VecNormalize
# from stable_baselines3.common.vec_env import VecNormalize, VecCheckNan
# from stable_baselines3.common.callbacks import EvalCallback,BaseCallback
from utils.AllocationEvalCallback import EvalCallback
from stable_baselines3.common.callbacks import StopTrainingOnNoModelImprovement

from gymnasium.envs.registration import register
import pandas as pd
import matplotlib.pyplot as plt
import sys
import torch as th
th.autograd.set_detect_anomaly(True)
import warnings
warnings.filterwarnings('ignore')
# 加载PyTorch并设置线程数
import torch
torch.set_num_threads(2)  # 不设置这个多进程会变得非常慢
os.environ['OPENBLAS_NUM_THREADS'] = '2'  # 在某些架构上不设置这个可能会导致Numpy导入错误


cdir = os.path.abspath('.') # 获取当前目录
distf = 0.95
# print(sys.path)
# 状态空间：现金余额，个人养老金账户余额，永久收入冲击，年龄
register(
    id='pensionfund-v3-simplified',                                # call it whatever you want
    entry_point='v3_pensionfund_simplified:PensionFundEnv', # module_name:class_name
)

model_name = 'pensionfund-v3_run66'
log_dir = 'models//sac//'

# load params
# with open(os.path.join(log_dir, model_name)+"//params.json", 'r') as file:
#     params = json.load(file)   
vec_env = make_vec_env("pensionfund-v3-simplified",seed=1, n_envs=22, env_kwargs={'params': None},\
                   monitor_kwargs={'discount_factor': distf})
vec_env = VecNormalize.load(os.path.join(log_dir, model_name)+ "//best_model_vec_env.pkl", vec_env) # load env
model = SAC.load(os.path.join(log_dir, model_name)+'//best_model', env=vec_env) # load model
model.load_replay_buffer(os.path.join(log_dir, model_name)+"//best_model_replay_buffer") # load rb

# 重新训练的参数设置
model.learning_starts = 0 # 重新训练不需要
model.learning_rate = 1e-06 
model.batch_size = 256

# 检查并创建RT目录
rt_dir = os.path.join(log_dir, model_name, 'RT')
if not os.path.exists(rt_dir):
    os.makedirs(rt_dir)
    next_rt = 0
else:
    # 获取目录中所有RT_X文件夹
    existing_dirs = [d for d in os.listdir(rt_dir) if os.path.isdir(os.path.join(rt_dir, d)) and d.startswith('RT_')]
    # 获取现有RT编号
    rt_numbers = []
    for d in existing_dirs:
        try:
            rt_num = int(d.split('RT_')[-1])
            rt_numbers.append(rt_num)
        except:
            continue
    # 确定新的RT编号
    next_rt = max(rt_numbers) + 1 if rt_numbers else 0

# 创建新的RT_X目录
rt_name = f"RT_{next_rt}"
new_rt_dir = os.path.join(rt_dir, rt_name)
if not os.path.exists(new_rt_dir):
    os.makedirs(new_rt_dir)


# create eval env
eval_env = make_vec_env("pensionfund-v3-simplified",seed=3687851522, n_envs=1, env_kwargs={'params': None},\
                        monitor_kwargs={'discount_factor': distf})       
eval_env = VecNormalize(eval_env, norm_obs=False, norm_reward=False)
eval_callback = EvalCallback(eval_env, seed=3687851522, best_model_save_path=new_rt_dir,\
                            eval_freq=5000, n_eval_episodes = 200, params =None,\
                            deterministic=True, render=False, discount_factor=distf ,score_target=30.5)

model.learn(total_timesteps=20_000_000, log_interval=100,callback=eval_callback,progress_bar=True,\
            ) #  reset_num_timesteps=False 这个可以继续使用原来的log


# 保存replay buffer
model.save_replay_buffer(new_rt_dir+"//best_model_replay_buffer")

vec_env.save(new_rt_dir+"//best_model_vec_env.pkl") 