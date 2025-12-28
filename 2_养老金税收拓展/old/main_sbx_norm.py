from typing import Callable
import json
import gymnasium as gym
import argparse  # 用于解析命令行参数
import numpy as np
import os
# from stable_baselines3 import DDPG
import gymnasium as gym
from sbx import DDPG, DQN, PPO, SAC, TD3, TQC, CrossQ
# from stable_baselines3 import SAC
from stable_baselines3.common.type_aliases import Schedule
from stable_baselines3.common.noise import NormalActionNoise, OrnsteinUhlenbeckActionNoise
from stable_baselines3.common.vec_env import VecCheckNan
# from stable_baselines3.common.logger import TensorBoardOutputFormat
from utils.env_util import make_vec_env
from utils.vec_normalize import VecNormalize
# from stable_baselines3.common.vec_env import VecNormalize, VecCheckNan
# from stable_baselines3.common.callbacks import EvalCallback,BaseCallback
# from utils.AllocationEvalCallback_PPO import EvalCallback
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
# import torch
# torch.set_num_threads(2)  # 不设置这个多进程会变得非常慢
# os.environ['OPENBLAS_NUM_THREADS'] = '2'  # 在某些架构上不设置这个可能会导致Numpy导入错误


cdir = os.path.abspath('.') # 获取当前目录
penlimit = 999 # 默认无养老金购买上限
distf = 0.95
# print(sys.path)
# 状态空间：现金余额，个人养老金账户余额，永久收入冲击，年龄

register(
    id='pensionfund-v3-norm',                                # call it whatever you want
    entry_point='v3_pensionfund_norm:PensionFundEnv', # module_name:class_name
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



def run(penlimit,env_id):
    
    log_dir = 'models//sac//'
    # log_dir = 'models//ppo//'
    # log_dir = 'models//td3//'
    log_env_id = env_id
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
    
    params = {'pension_limit': penlimit}
    vec_env = make_vec_env("pensionfund-v3-norm",seed=1, n_envs=16, # env_kwargs={'params':params},\
                            monitor_kwargs={'discount_factor': distf})
    vec_env = VecNormalize(vec_env, norm_obs=False, norm_reward=False, norm_obs_index_mask=[3])
    # vec_env =VecCheckNan(vec_env, raise_exception=True)
    # vec_env = VecCheckNan(vec_env, raise_exception=True)

    eval_env = make_vec_env("pensionfund-v3-norm",seed=3687851522, n_envs=1, # env_kwargs={'params': params},\
                            monitor_kwargs={'discount_factor': distf})         
    eval_env = VecNormalize(eval_env, norm_obs=False, norm_reward=False)

    # 让eval_env先运行一段时间以便于normalize
    # obs = eval_env.reset()
    # for _ in range(10000):
    #     action = np.array([[1,1,0,0]])
    #     obs, reward, done, info = eval_env.step(action)

    eval_callback = EvalCallback(eval_env, seed=3687851522, best_model_save_path=os.path.join(log_dir, model_name),\
                              eval_freq=10_000, n_eval_episodes = 200, params = params,\
                              deterministic=True, render=False, discount_factor=distf,score_target=26)

    # SAC
    model = SAC("MlpPolicy", vec_env , verbose=1, gamma=distf, device="cpu",\
                tensorboard_log=".//runs//",
                gradient_steps=4,
                learning_rate=3e-04, # 4e-04, # 0.00001
                policy_kwargs={'net_arch':[64,64]}, # 5priod: [32,32] ,'log_std_init':-3.67
                # buffer_size=100_000,
                # action_noise = OrnsteinUhlenbeckActionNoise(mean=np.zeros(vec_env.action_space.shape[0]), sigma=0.3* np.ones(vec_env.action_space.shape[0])),
                # tau=0.01,
                # ent_coef=0.3,
                batch_size=256,
                # 5priod: batch_size=64,
                # use_sde=True, 
                # buffer_size= 50_000,
                train_freq=4,
                # batch_size=2048,
                learning_starts=10_000)
    # model = TD3("MlpPolicy", vec_env , verbose=1, gamma=distf, device="cpu",\
    #         tensorboard_log=".//runs//",
    #         # gradient_steps=4,
    #         # learning_rate=3e-04, # 4e-04, # 0.00001
    #         policy_kwargs={'net_arch':[64,64]}, # 5priod: [32,32] ,'log_std_init':-3.67
    #         # buffer_size=100_000,
    #         # action_noise = OrnsteinUhlenbeckActionNoise(mean=np.zeros(vec_env.action_space.shape[0]), sigma=0.3* np.ones(vec_env.action_space.shape[0])),
    #         # tau=0.01,
    #         # ent_coef=0.3,
    #         # batch_size=256,
    #         # 5priod: batch_size=64,
    #         # use_sde=True, 
    #         # buffer_size= 50_000,
    #         # train_freq=4,
    #         # batch_size=2048,
    #         learning_starts=10_000)
    # # model = PPO("MlpPolicy", vec_env , verbose=1, gamma=distf, device="cpu",\
    # #         tensorboard_log=".//runs//",
    # #         learning_rate=3e-04, # 4e-04, # 0.00001
    # #         policy_kwargs={'net_arch':[64,64]},ent_coef=10,use_sde=True)

    model.learn(total_timesteps=50_000_000, log_interval=100,callback=eval_callback,progress_bar=True)

    model.save_replay_buffer(os.path.join(log_dir, model_name)+"//replay_buffer")
    vec_env.save(os.path.join(log_dir, model_name)+"//vec_env.pkl") 
    # 保存评估环境
    eval_env.save(os.path.join(log_dir, model_name)+"//eval_env.pkl")
    # 保存params
    with open(os.path.join(log_dir, model_name)+"//params.json", 'w') as file:
        json.dump(params, file)

def parse_args(args=None):
    # 解析命令行参数
    parser = argparse.ArgumentParser("Experiments")
    # 通用参数:
    parser.add_argument('--envid', default="pensionfund-v3-norm", type=str, metavar='ENVID', help="环境 (默认: %(default)s)")
    parser.add_argument('--penlimit', default=penlimit, type=float, metavar='PENLIMIT', help="养老金购买上限 (默认: %(default)s)")
    params = parser.parse_args(args)

    return params

if __name__ == "__main__":
    params = parse_args()
    # 将参数解包传递给run函数
    run(params.penlimit,params.envid)