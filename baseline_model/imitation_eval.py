import numpy as np
import gymnasium as gym
from imitation.util.util import make_vec_env
from imitation.policies.serialize import load_policy
import utils.rollout as rollout
from imitation.algorithms import bc
from imitation.data import serialize
from stable_baselines3 import SAC
from imitation.algorithms.adversarial.airl import AIRL
import imitation
from utils import rollout
from imitation.rewards.reward_nets import BasicShapedRewardNet
from utils.evaluate_policy import evaluate_policy,evaluate_policy_dp
from imitation.util.networks import RunningNorm
import pandas as pd
import os
from tqdm import tqdm
from utils.env_util import make_vec_env
from utils.vec_normalize import VecNormalize
import matplotlib.pyplot as plt
os.environ['KMP_DUPLICATE_LIB_OK']='True'
from gymnasium.envs.registration import register
register(
    id='pensionfund-v3',                                # call it whatever you want
    entry_point='v3_pensionfund:PensionFundEnv', # module_name:class_name
)

def scale_action(env, action: np.ndarray) -> np.ndarray:
    """
    Rescale the action from [low, high] to [-1, 1]
    (no need for symmetric action space)

    :param action: Action to scale
    :return: Scaled action
    """
    low, high = env.action_space.low, env.action_space.high
    return 2.0 * ((action - low) / (high - low)) - 1.0

class dp_policy():
    def __init__(self,vec_env,A,C,gcash): # 读取policy
        self.A = A
        self.C = C
        self.gcash= gcash
        self.f = vec_env.envs[0].env.env.f # 获取年龄收入函数
        # self.aa = self.env.
        self.n_envs = np.shape(vec_env.envs)[0]

    def predict(self,obs,vec_env):
        normalized_cash = np.array([i.env.env.env.info['normalized_cash'] for i in vec_env.envs])
        age_loc =  np.array(vec_env.get_attr('age')) - np.array(vec_env.get_attr('tb'))# 获取年龄
        action = np.zeros((self.n_envs,vec_env.action_space.shape[0]))
        # 根据现金和年龄插值从A中获取action
        # 消费
        action[:,0] = np.array([np.interp(normalized_cash[i], self.gcash, np.squeeze(self.C[:, age_loc[i]])) for i in range(self.n_envs)])
        action[:,0] = action[:,0]/(normalized_cash) # 消费比例 
        # risky
        action[:,1] = np.array([np.interp(normalized_cash[i], self.gcash, np.squeeze(self.A[:, age_loc[i]])) for i in range(self.n_envs)])

        action[age_loc == np.shape(self.C)[1]-1,0] = 1 # 如果到达最后一期,则固定消费比例为1（全部消费干净）
        action[age_loc == np.shape(self.C)[1]-1,1] = 0 # 如果到达最后一期,则riksy比例为0    
        return action
         


# 生成rb
distf = 0.95
n_envs = 88


vec_env = make_vec_env("pensionfund-v3",seed=1, n_envs=n_envs, # vec_env_kwargs={'params':params},\
                        monitor_kwargs={'discount_factor': distf})
# vec_env = VecNormalize(vec_env, norm_obs=True, norm_reward=False) # 
vec_env = VecNormalize.load('./im_models/bc_ppo_256_vec_env.pkl',vec_env)


# 让eval_env先运行一段时间以便于normalize
# obs = vec_env.reset()
# pbar = tqdm(total=5000, desc=f'模拟进度')
# for _ in range(5000):
#     pbar.update(1)
#     action = np.array([vec_env.action_space.sample() for _ in range(n_envs)])
#     obs, reward, done, info = vec_env.step(action)



eval_env = make_vec_env("pensionfund-v3",seed=3687851522, n_envs=1, # env_kwargs={'params': params},\
                        monitor_kwargs={'discount_factor': distf})         
eval_env = VecNormalize(eval_env, norm_obs=False, norm_reward=False)

# 读取matlab dp解
current_path = os.path.dirname(os.path.abspath(__file__))
age_group = str(eval_env.get_attr('tb')[0])+'-'+str(eval_env.get_attr('td')[0])
df_A = pd.read_excel(os.path.join(current_path,'result_cocco_matlab',age_group,'A.xlsx'),header=None)
df_C = pd.read_excel(os.path.join(current_path,'result_cocco_matlab',age_group,'C.xlsx'),header=None) 
df_gcash = pd.read_excel(os.path.join(current_path,'result_cocco_matlab',age_group,'gcash.xlsx'),header=None)
# 转换为numpy数组
df_A = df_A.to_numpy()
df_C = df_C.to_numpy()
df_gcash = df_gcash.to_numpy().squeeze()
dpp = dp_policy(eval_env,df_A,df_C,df_gcash)

reward_dp, _ ,_ ,mean_actions = evaluate_policy_dp(dpp, vec_env, eval_env,seed=3687851522, n_eval_episodes=200,deterministic=True, render=False, discount_factor=distf)
print(f"Reward DP: {reward_dp}")

plt.figure(figsize=(10, 6))
plt.plot(mean_actions[:,0], mean_actions[:,1], label='consum')
plt.plot(mean_actions[:,0], mean_actions[:,2], label='riksy')   
plt.plot(mean_actions[:,0], mean_actions[:,3], label='tda_pct')  
plt.plot(mean_actions[:,0], mean_actions[:,4], label='riksy_tda')  
plt.savefig('figs/dp.png', dpi=300, bbox_inches='tight')  # 保存为PNG文件

# policy = bc.reconstruct_policy('./im_models/bc.zip')
# policy = SAC.load('./im_models/airl_sac.zip', print_system_info=True)
from stable_baselines3.sac.policies import MlpPolicy
from stable_baselines3 import SAC,A2C,PPO
# model = PPO("MlpPolicy", vec_env , verbose=1, gamma=distf, device="cuda",
#             policy_kwargs={"net_arch": [256,256],})
policy = PPO.load('./im_models/bc_ppo_256.zip',vec_env).policy

# policy = load_policy(policy_type='sac',venv=vec_env,path='./im_models/airl_sac.zip')
reward_after_training, _ ,_, mean_actions = evaluate_policy(policy, vec_env, eval_env,seed=3687851522, n_eval_episodes=200,deterministic=True, render=False, discount_factor=distf)
print(f"Reward after training: {reward_after_training}")



plt.figure(figsize=(10, 6))
plt.plot(mean_actions[:,0], mean_actions[:,1], label='consum')
plt.plot(mean_actions[:,0], mean_actions[:,2], label='riksy')   
plt.plot(mean_actions[:,0], mean_actions[:,3], label='tda_pct')  
plt.plot(mean_actions[:,0], mean_actions[:,4], label='riksy_tda')  
plt.savefig('figs/bc.png', dpi=300, bbox_inches='tight')  # 保存为PNG文件


