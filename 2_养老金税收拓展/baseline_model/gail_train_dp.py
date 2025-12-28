import numpy as np
import gymnasium as gym
from imitation.util.util import make_vec_env
from imitation.util.util import save_policy
import utils.rollout as rollout
from imitation.algorithms import bc
from imitation.data import serialize
from stable_baselines3 import SAC
from imitation.algorithms.adversarial.gail import GAIL
from utils import rollout
from imitation.rewards.reward_nets import BasicShapedRewardNet
from utils.evaluate_policy import evaluate_policy
from imitation.util.networks import RunningNorm
import pandas as pd
import os
from tqdm import tqdm
from utils.env_util import make_vec_env
from utils.vec_normalize import VecNormalize
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

    def policy(self,obs,vec_env):
        normalized_cash = np.array([i.env.env.env.info['normalized_cash'] for i in vec_env.envs])
        age_loc =  np.array(vec_env.get_attr('age')) - np.array(vec_env.get_attr('tb'))# 获取年龄
        action = np.zeros((n_envs,vec_env.action_space.shape[0]))
        # 根据现金和年龄插值从A中获取action
        # 消费
        action[:,0] = np.array([np.interp(normalized_cash[i], self.gcash, np.squeeze(self.C[:, age_loc[i]])) for i in range(n_envs)])
        action[:,0] = action[:,0]/(normalized_cash) # 消费比例 
        # risky
        action[:,1] = np.array([np.interp(normalized_cash[i], self.gcash, np.squeeze(self.A[:, age_loc[i]])) for i in range(n_envs)])

        action[age_loc == np.shape(self.C)[1]-1,0] = 1 # 如果到达最后一期,则固定消费比例为1（全部消费干净）
        action[age_loc == np.shape(self.C)[1]-1,1] = 0 # 如果到达最后一期,则riksy比例为0    
        return action
         


# 生成rb
distf = 0.95
n_envs = 88


vec_env = make_vec_env("pensionfund-v3",seed=1, n_envs=n_envs, # vec_env_kwargs={'params':params},\
                        monitor_kwargs={'discount_factor': distf})
vec_env = VecNormalize(vec_env, norm_obs=True, norm_reward=False) # 

# 让eval_env先运行一段时间以便于normalize
obs = vec_env.reset()
pbar = tqdm(total=1000, desc=f'模拟进度')
for _ in range(1000):
    pbar.update(1)
    action = np.array([vec_env.action_space.sample() for _ in range(n_envs)])
    obs, reward, done, info = vec_env.step(action)

# print(vec_env.obs_rms.mean)


# 读取matlab dp解
# current_path = os.path.dirname(os.path.abspath(__file__))
# age_group = str(vec_env.get_attr('tb')[0])+'-'+str(vec_env.get_attr('td')[0])
# df_A = pd.read_excel(os.path.join(current_path,'result_cocco_matlab',age_group,'A.xlsx'),header=None)
# df_C = pd.read_excel(os.path.join(current_path,'result_cocco_matlab',age_group,'C.xlsx'),header=None) 
# df_gcash = pd.read_excel(os.path.join(current_path,'result_cocco_matlab',age_group,'gcash.xlsx'),header=None)
# # 转换为numpy数组
# df_A = df_A.to_numpy()
# df_C = df_C.to_numpy()
# df_gcash = df_gcash.to_numpy().squeeze()
# dpp = dp_policy(vec_env,df_A,df_C,df_gcash)



# rng = np.random.default_rng()
# rollouts = rollout.rollout(
#     dpp.policy,
#     vec_env,
#     rollout.make_sample_until(min_timesteps=100_000),
#     rng=rng,
#     unwrap=False
# )


# serialize.save('./traj/', rollouts)

print('loading trajs....')
rollouts = serialize.load('./traj/')
transitions = rollout.flatten_trajectories(rollouts)

print('loading trajs completed.')
# from stable_baselines3.common.evaluation import evaluate_policy

rng = np.random.default_rng()
learner = SAC("MlpPolicy", vec_env , verbose=1, gamma=distf, device="cuda",\
              tensorboard_log="../tf-logs/",
            policy_kwargs={"net_arch": [32,32],})

reward_net = BasicShapedRewardNet(
    observation_space=vec_env.observation_space,
    action_space=vec_env.action_space,
    normalize_input_layer=RunningNorm,
    discount_factor=distf,
    # use_next_state = True,
    # use_done = True,    
    # reward_hid_sizes = (64,),
    # potential_hid_sizes = (64, 64),
)
trainer = GAIL(
    demonstrations=rollouts,
    demo_batch_size=2048,
    gen_replay_buffer_capacity=512,
    n_disc_updates_per_round=16,
    venv=vec_env,
    gen_algo=learner,
    reward_net=reward_net,
)

# sqil_trainer = bc.BC(
#     observation_space=vec_env.observation_space,
#     action_space=vec_env.action_space,
#     demonstrations=transitions,
#     rng=rng,
# )

eval_env = make_vec_env("pensionfund-v3",seed=3687851522, n_envs=1, # env_kwargs={'params': params},\
                        monitor_kwargs={'discount_factor': distf})         
eval_env = VecNormalize(eval_env, norm_obs=False, norm_reward=False)


reward_before_training, _ ,_, _ = evaluate_policy(trainer.policy, vec_env, eval_env,seed=3687851522, n_eval_episodes=200,deterministic=True, render=False, discount_factor=distf)


trainer.train(10_000) 
# bc_trainer.train(expert_traj)
reward_after_training, _ ,_, _ = evaluate_policy(trainer.policy, vec_env, eval_env,seed=3687851522, n_eval_episodes=200,deterministic=True, render=False, discount_factor=distf)
print(f"Reward before training: {reward_before_training}")
print(f"Reward after training: {reward_after_training}")
# trainer.policy.save('./im_models/gail_sac_32_800k.pkl')
# vec_env.save('./im_models/gail_sac_32_800k_vec_env.pkl') 
# save_policy(trainer.policy,'./im_models/airl_sac.zip')
