from sbx import SAC
from stable_baselines3.common.buffers import ReplayBuffer
from stable_baselines3.common.policies import BasePolicy
from sbx.common.type_aliases import ReplayBufferSamplesNp
import pandas as pd
import os
from tqdm import tqdm
from utils.env_util import make_vec_env
from utils.vec_normalize import VecNormalize
from stable_baselines3.common.save_util import load_from_pkl, save_to_pkl
import numpy as np
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

distf = 0.95

eval_seed = 1

buffer_size = 5_000_000
n_eval_episodes = buffer_size//81 + 81
n_envs = 88


vec_env = make_vec_env("pensionfund-v3",seed=eval_seed, n_envs=n_envs, # vec_env_kwargs={'params':params},\
                        monitor_kwargs={'discount_factor': distf})
vec_env = VecNormalize(vec_env, norm_obs=True, norm_reward=True, norm_obs_index_mask=[3]) # 

# model = SAC("MlpPolicy", vec_env , verbose=1, gamma=distf, device="cuda",\
#             tensorboard_log="../tf-logs/",
#             gradient_steps=22,
#             learning_rate=lr_schedule,
#             # action_noise = OrnsteinUhlenbeckActionNoise(mean=np.zeros(vec_env.action_space.shape[0]), sigma=0.5* np.ones(vec_env.action_space.shape[0])),
#             # learning_rate=1e-03,
#             policy_kwargs={
#                 "net_arch": [512, 512],
#             },
#             # learning_rate=1.0,  # 此处设为1.0，因为调度器已控制实际学习率
#             batch_size=128,
#             tau=0.005,
#             # 5priod: batch_size=64,
#             # use_sde=True, 
#             buffer_size= 2_000_000,
#             train_freq=6,
#             # batch_size=2048,
#             learning_starts=100_000)

# 生成rb
rb_class = ReplayBuffer
rb = rb_class(  # type: ignore[misc]
    buffer_size,
    observation_space = vec_env.observation_space,
    action_space = vec_env.action_space,
    device="cuda",  # force cpu device to easy torch -> numpy conversion
    n_envs=n_envs,
)

# 读取matlab dp解
current_path = os.path.dirname(os.path.abspath(__file__))
age_group = str(vec_env.get_attr('tb')[0])+'-'+str(vec_env.get_attr('td')[0])
df_A = pd.read_excel(os.path.join(current_path,'result_cocco_matlab',age_group,'A.xlsx'),header=None)
df_C = pd.read_excel(os.path.join(current_path,'result_cocco_matlab',age_group,'C.xlsx'),header=None) 
df_V = pd.read_excel(os.path.join(current_path,'result_cocco_matlab',age_group,'V.xlsx'),header=None)  
df_gcash = pd.read_excel(os.path.join(current_path,'result_cocco_matlab',age_group,'gcash.xlsx'),header=None)
# 转换为numpy数组
A = df_A.to_numpy()
C = df_C.to_numpy()
V = df_V.to_numpy()
gcash = df_gcash.to_numpy().squeeze()




episode_rewards = []
episode_lengths = []

episode_counts = np.zeros(n_envs, dtype="int")
# Divides episodes among different sub environments in the vector as evenly as possible
episode_count_targets = np.array([(n_eval_episodes + i) // n_envs for i in range(n_envs)], dtype="int")

current_rewards = np.zeros(n_envs)
current_lengths = np.zeros(n_envs, dtype="int")

states = None
episode_starts = np.ones((vec_env.num_envs,), dtype=bool)

# 初始化seed
# seed_list = [seed+i for i in range(n_eval_episodes)] # 所有episode的seed
vec_env.seed(seed=eval_seed)
observations = vec_env.reset()   
# for idx in range(n_envs):
#     observations[idx,:] = np.array([env.envs[idx].reset(seed=seed_list[idx])[0]])

from tqdm import tqdm

# 计算总episode数
total_episodes = sum(episode_count_targets)
# 创建进度条
pbar = tqdm(total=total_episodes, desc="Generating Replay Buffer", unit="episode")

while (episode_counts < episode_count_targets).any():
    # 更新进度条
    pbar.n = sum(episode_counts)
    pbar.refresh()
    # env.seed(seed=seed+env.num_envs)
    normalized_cash = np.array([i.env.env.env.info['normalized_cash'] for i in vec_env.envs])
    age_loc =  np.array(vec_env.get_attr('age')) - np.array(vec_env.get_attr('tb'))# 获取年龄
    action = np.zeros((n_envs,vec_env.action_space.shape[0]))
    # 根据现金和年龄插值从A中获取action
    # 消费
    action[:,0] = np.array([np.interp(normalized_cash[i], gcash, np.squeeze(C[:, age_loc[i]])) for i in range(n_envs)])
    action[:,0] = action[:,0]/(normalized_cash) # 消费比例 
    # risky
    action[:,1] = np.array([np.interp(normalized_cash[i], gcash, np.squeeze(A[:, age_loc[i]])) for i in range(n_envs)])

    action[age_loc == np.shape(C)[1]-1,0] = 1 # 如果到达最后一期,则固定消费比例为1（全部消费干净）
    action[age_loc == np.shape(C)[1]-1,1] = 0 # 如果到达最后一期,则riksy比例为0

    new_observations, rewards, dones, infos = vec_env.step(action)



    for i in range(n_envs):
        if episode_counts[i] < episode_count_targets[i]:
            # unpack values so that the callback can access the local variables
            reward = rewards[i]
            done = dones[i]
            info = infos[i]
            episode_starts[i] = done

            if dones[i]:
                    episode_rewards.append(current_rewards[i])
                    episode_lengths.append(current_lengths[i])
                    episode_counts[i] += 1
            current_rewards[i] = 0 # 重新开始累积reward
            current_lengths[i] = 0


    rb.add(
        obs=vec_env.unnormalize_obs(observations),
        action=scale_action(vec_env,action),
        reward=vec_env.get_original_reward(),
        next_obs=vec_env.unnormalize_obs(new_observations),
        done=done,
        infos=infos,
    )
    observations = new_observations

# data = rb.sample(100) #,env=vec_env
# data = ReplayBufferSamplesNp(  # type: ignore[assignment]
#             data.observations.numpy(),
#             data.actions.numpy(),
#             data.next_observations.numpy(),
#             data.dones.numpy().flatten(),
#             data.rewards.numpy().flatten(),
#         )

save_to_pkl('./replay_buffer/v3_no_stoch_death_rb.pkl',rb)
vec_env.save('./replay_buffer/v3_no_stoch_death_env.pkl')