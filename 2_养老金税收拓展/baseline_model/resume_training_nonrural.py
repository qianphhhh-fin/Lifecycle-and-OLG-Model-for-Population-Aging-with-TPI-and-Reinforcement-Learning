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
from optax import cosine_decay_schedule, chain, adamw, GradientTransformation, constant_schedule
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

os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"  # add this
os.environ["JAX_TRACEBACK_FILTERING"]= "off"

cdir = os.path.abspath('.') # 获取当前目录
distf = 0.95
# print(sys.path)
# 状态空间：现金余额，个人养老金账户余额，永久收入冲击，年龄
register(
    id='pensionfund-v5-nonrural',                                # call it whatever you want
    entry_point='v5_pensionfund_nonrural:PensionFundEnv', # module_name:class_name
)


model_name = 'pensionfund-v5-nonrural_run52//'
log_dir = '../models//nonrural//'

# load params
# with open(os.path.join(log_dir, model_name)+"//params.json", 'r') as file:
#     params = json.load(file)   
# params['pension_limit'] = 0.3
vec_env = make_vec_env("pensionfund-v5-nonrural",seed=1, n_envs=88, env_kwargs={'params':None},\
                   monitor_kwargs={'discount_factor': distf})
# vec_env = VecNormalize(vec_env, norm_obs=True, norm_reward=True) #  新定义env , norm_obs_index_mask=[3]
# vec_env = VecNormalize.load('./vec_env/vec_env_stoch.pkl', vec_env) # load env
vec_env = VecNormalize.load(os.path.join(log_dir, model_name)+ "//best_model_vec_env.pkl", vec_env) # load env
# 虽然模型不更新，但是norm_obs和norm_reward但是会更新，导致learning_start阶段eval结果产生微小变化
# vec_env.norm_obs = False
# vec_env.norm_reward = False

# obs = vec_env.reset()
# for _ in range(10000):
#     action = np.array([vec_env.action_space.sample() for _ in range(88)])
#     obs, reward, done, info = vec_env.step(action)
# print('sample vec_env completed.')

# vec_env = VecNormalize.load(os.path.join(log_dir, model_name)+ "//best_model_vec_env.pkl", vec_env) # load env
model = SAC.load(os.path.join(log_dir, model_name)+'//best_model', env=vec_env) # load model rb
model.load_replay_buffer(os.path.join(log_dir, model_name)+"//best_model_replay_buffer") # load rb

# 重新训练的参数设置
model.learning_starts = 0 # 重新训练不需要/使用已有policy生成

initial_learning_rate = 3e-7
final_learning_rate = 1e-8
total_steps = 1_000_000  # 总训练步数

# 余弦退火调度 (从 initial_learning_rate 衰减到 final_learning_rate)
lr_schedule = cosine_decay_schedule(
    init_value=initial_learning_rate,
    decay_steps=total_steps,
    alpha=final_learning_rate / initial_learning_rate  # 最终学习率比例
)

model.learning_rate = constant_schedule(initial_learning_rate)
model.batch_size = 256 
model.tau = 0.005
# model.train_freq = 1
model.gradient_steps = 6
model.buffer_size = 2_000_000
model.ent_coef = 0.1
# model.action_noise = OrnsteinUhlenbeckActionNoise(mean=np.zeros(vec_env.action_space.shape[0]), sigma=0.1* np.ones(vec_env.action_space.shape[0]))


# 检查并创建RT目录
# 检测是否包含RT
if 'RT' in model_name:
    rt_dir = os.path.join(log_dir, model_name) # 如果目录中包含RT, 则无需创建RT目录
else:
    rt_dir = os.path.join(log_dir, model_name, 'RT') # 如果目录中不包含RT, 则创建RT目录

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
eval_env = make_vec_env("pensionfund-v5-nonrural",seed=3687851522, n_envs=1, env_kwargs={'params': None},\
                        monitor_kwargs={'discount_factor': distf})       
eval_env = VecNormalize(eval_env, norm_obs=False, norm_reward=False)

eval_callback = EvalCallback(eval_env, seed=3687851522, best_model_save_path=new_rt_dir,\
                            eval_freq=10, n_eval_episodes = 200, params = None,\
                            deterministic=True, render=False, discount_factor=distf)

model.learn(total_timesteps=total_steps, log_interval=100,callback=eval_callback,progress_bar=True,\
            ) #  reset_num_timesteps=False 这个可以继续使用原来的log


# 保存replay buffer
model.save_replay_buffer(new_rt_dir+"//best_model_replay_buffer")
model.save(new_rt_dir+"//latest_model.zip")
vec_env.save(new_rt_dir+"//vec_env.pkl") 
# 保存params
# with open(new_rt_dir+"//params.json", 'w') as file:
#     json.dump(params, file)
