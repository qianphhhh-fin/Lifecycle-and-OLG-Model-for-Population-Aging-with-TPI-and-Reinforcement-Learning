import gymnasium as gym
from tqdm import tqdm
import numpy as np

print(f"{gym.__version__}")

import torch as th
import torch.nn as nn
import torch.optim as optim
from torch.optim.lr_scheduler import StepLR
from torch.utils.data.dataset import Dataset, random_split
import numpy as np
import gymnasium as gym
from stable_baselines3 import SAC,A2C,PPO,DDPG,TD3
from utils.evaluate_policy import evaluate_policy
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
        action = np.zeros((self.n_envs,vec_env.action_space.shape[0]))
        # 根据现金和年龄插值从A中获取action
        # 消费
        action[:,0] = np.array([np.interp(normalized_cash[i], self.gcash, np.squeeze(self.C[:, age_loc[i]])) for i in range(n_envs)])
        action[:,0] = action[:,0]/(normalized_cash) # 消费比例 
        # risky
        action[:,1] = np.array([np.interp(normalized_cash[i], self.gcash, np.squeeze(self.A[:, age_loc[i]])) for i in range(n_envs)])

        action[age_loc == np.shape(self.C)[1]-1,0] = 1 # 如果到达最后一期,则固定消费比例为1（全部消费干净）
        action[age_loc == np.shape(self.C)[1]-1,1] = 0 # 如果到达最后一期,则riksy比例为0    
        return action
    

class ExpertDataSet(Dataset):
    def __init__(self, expert_observations, expert_actions):
        self.observations = expert_observations
        self.actions = expert_actions

    def __getitem__(self, index):
        return (self.observations[index], self.actions[index])

    def __len__(self):
        return len(self.observations)

def pretrain_agent(
    student,
    batch_size=64,
    epochs=1000,
    scheduler_gamma=0.7,
    learning_rate=1.0,
    log_interval=100,
    cuda=True,
    seed=1,
    test_batch_size=64,
):
    use_cuda = cuda and th.cuda.is_available()
    th.manual_seed(seed)
    device = th.device("cuda" if use_cuda else "cpu")
    kwargs = {"num_workers": 10, "pin_memory": True} if use_cuda else {}

    if isinstance(env.action_space, gym.spaces.Box):
        criterion = nn.MSELoss()
    else:
        criterion = nn.CrossEntropyLoss()

    # Extract initial policy
    model = student.policy.to(device)

    def train(model, device, train_loader, optimizer):
        model.train()

        for batch_idx, (data, target) in enumerate(train_loader):
            data, target = data.to(device), target.to(device)
            optimizer.zero_grad()

            if isinstance(env.action_space, gym.spaces.Box):
                # A2C/PPO policy outputs actions, values, log_prob
                # SAC/TD3 policy outputs actions only
                if isinstance(student, (A2C, PPO)):
                    action, _, _ = model(data)
                else:
                    # SAC/TD3:
                    action = model(data)
                action_prediction = action.double()
            else:
                # Retrieve the logits for A2C/PPO when using discrete actions
                dist = model.get_distribution(data)
                action_prediction = dist.distribution.logits
                target = target.long()

            loss = criterion(action_prediction, target)
            loss.backward()
            optimizer.step()
            # if batch_idx % log_interval == 0:
            #     print(
            #         "Train Epoch: {} [{}/{} ({:.0f}%)]\tLoss: {:.6f}".format(
            #             epoch,
            #             batch_idx * len(data),
            #             len(train_loader.dataset),
            #             100.0 * batch_idx / len(train_loader),
            #             loss.item(),
            #         )
            #     )

    def test(model, device, test_loader):
        model.eval()
        test_loss = 0
        with th.no_grad():
            for data, target in test_loader:
                data, target = data.to(device), target.to(device)

                if isinstance(env.action_space, gym.spaces.Box):
                    # A2C/PPO policy outputs actions, values, log_prob
                    # SAC/TD3 policy outputs actions only
                    if isinstance(student, (A2C, PPO)):
                        action, _, _ = model(data)
                    else:
                        # SAC/TD3:
                        action = model(data)
                    action_prediction = action.double()
                else:
                    # Retrieve the logits for A2C/PPO when using discrete actions
                    dist = model.get_distribution(data)
                    action_prediction = dist.distribution.logits
                    target = target.long()

                test_loss = criterion(action_prediction, target)
        test_loss /= len(test_loader.dataset)
        print(f"Test set: Average loss: {test_loss:.8f}")

    # Here, we use PyTorch `DataLoader` to our load previously created `ExpertDataset` for training
    # and testing
    train_loader = th.utils.data.DataLoader(
        dataset=train_expert_dataset, batch_size=batch_size, shuffle=True, **kwargs
    )
    test_loader = th.utils.data.DataLoader(
        dataset=test_expert_dataset,
        batch_size=test_batch_size,
        shuffle=True,
        **kwargs,
    )

    # Define an Optimizer and a learning rate schedule.
    optimizer = optim.Adadelta(model.parameters(), lr=learning_rate)
    scheduler = StepLR(optimizer, step_size=1, gamma=scheduler_gamma)

    # Now we are finally ready to train the policy model.
    pbar = tqdm(total=epochs, desc=f'training..')
    for epoch in range(1, epochs + 1):
        pbar.update(1)
        train(model, device, train_loader, optimizer)
        # test(model, device, test_loader)
        # if epoch % 1 == 0:
        #     mean_reward, std_reward, _, _ = evaluate_policy(model, env, eval_env,seed=3687851522, n_eval_episodes=200,deterministic=True, render=False, discount_factor=distf)
        #     print(f"Mean reward = {mean_reward}")
        scheduler.step()

    # Implant the trained policy network back into the RL student agent
    student.policy = model

# Example for continuous actions
# env_id = "LunarLanderContinuous-v2"

# Example for discrete actions
# env_id = "CartPole-v1"

# env = gym.make(env_id)

## Train Expert Model

# We create an expert RL agent and let it learn to solve a task by interacting with the evironment.

# print('creating expert....')
# ppo_expert = PPO("MlpPolicy", env_id, verbose=1)
# ppo_expert.learn(total_timesteps=3e4,progress_bar=True)
# ppo_expert.save("ppo_expert")

# check the performance of the trained agent

# mean_reward, std_reward = evaluate_policy(ppo_expert, env, n_eval_episodes=10)

# print(f"Mean reward = {mean_reward} +/- {std_reward}")

## Create Student

# We also create a student RL agent, which will later be trained with the expert dataset


# a2c_student = A2C("MlpPolicy", env_id, verbose=1)

# only valid for continuous actions
# sac_student = SAC('MlpPolicy', env_id, verbose=1, policy_kwargs=dict(net_arch=[64, 64]))


# We now let our expert interact with the environment (except we already have expert data) and store resultant expert observations and actions to build an expert dataset.

# 生成rb
distf = 0.95
n_envs = 1
num_interactions = 1_000_000

env = make_vec_env("pensionfund-v3",seed=1, n_envs=n_envs, # vec_env_kwargs={'params':params},\
                        monitor_kwargs={'discount_factor': distf})
env = VecNormalize(env, norm_obs=True, norm_reward=False) # 

obs = env.reset()
pbar = tqdm(total=10000, desc=f'模拟进度')
for _ in range(10000):
    pbar.update(1)
    action = np.array([env.action_space.sample() for _ in range(n_envs)])
    obs, reward, done, info = env.step(action)

# 生成traj
current_path = os.path.dirname(os.path.abspath(__file__))
age_group = str(env.get_attr('tb')[0])+'-'+str(env.get_attr('td')[0])
df_A = pd.read_excel(os.path.join(current_path,'result_cocco_matlab',age_group,'A.xlsx'),header=None)
df_C = pd.read_excel(os.path.join(current_path,'result_cocco_matlab',age_group,'C.xlsx'),header=None) 
df_gcash = pd.read_excel(os.path.join(current_path,'result_cocco_matlab',age_group,'gcash.xlsx'),header=None)
# 转换为numpy数组
df_A = df_A.to_numpy()
df_C = df_C.to_numpy()
df_gcash = df_gcash.to_numpy().squeeze()
dpp = dp_policy(env,df_A,df_C,df_gcash)



if isinstance(env.action_space, gym.spaces.Box):
    expert_observations = np.empty((num_interactions,) + env.observation_space.shape)
    expert_actions = np.empty((num_interactions,) + (env.action_space.shape[0],))

else:
    expert_observations = np.empty((num_interactions,) + env.observation_space.shape)
    expert_actions = np.empty((num_interactions,) + env.action_space.shape)

obs = env.reset()

for i in tqdm(range(num_interactions)):
    action = dpp.policy(obs, env)
    expert_observations[i] = obs
    expert_actions[i] = action
    obs, reward, done, info = env.step(action)
    if done:
        obs = env.reset()

np.savez_compressed(
    "/traj/1M_trivial_expert_data",
    expert_actions=expert_actions,
    expert_observations=expert_observations,
)

data = np.load( "./traj/1M_trivial_expert_data.npz")
expert_actions = data['expert_actions']
expert_observations = data['expert_observations']

# - To seamlessly use PyTorch in the training process, we subclass an `ExpertDataset` from PyTorch's base `Dataset`.
# - Note that we initialize the dataset with the previously generated expert observations and actions.
# - We further implement Python's `__getitem__` and `__len__` magic functions to allow PyTorch's dataset-handling to access arbitrary rows in the dataset and inform it about the length of the dataset.
# - For more information about PyTorch's datasets, you can read: https://pytorch.org/docs/stable/data.html.

# We now instantiate the `ExpertDataSet` and split it into training and test datasets.

expert_dataset = ExpertDataSet(expert_observations, expert_actions)

train_size = int(0.95 * len(expert_dataset))

test_size = len(expert_dataset) - train_size

train_expert_dataset, test_expert_dataset = random_split(
    expert_dataset, [train_size, test_size]
)



# NOTE: The supervised learning section of this code is adapted from: https://github.com/pytorch/examples/blob/master/mnist/main.py
# 1. We extract the policy network of our RL student agent.
# 2. We load the (labeled) expert dataset containing expert observations as inputs and expert actions as targets.
# 3. We perform supervised learning, that is, we adjust the policy network's parameters such that given expert observations as inputs to the network, its outputs match the targets (expert actions).
# By training the policy network in this way the corresponding RL student agent is taught to behave like the expert agent that was used to created the expert dataset (Behavior Cloning).


# Evaluate the agent before pretraining, it should be random
student = SAC("MlpPolicy", env , verbose=1, gamma=distf, device="cuda",\
            policy_kwargs={"net_arch": [512,512],})
# 冻结除actor之外的所有参数
# for name, param in student.policy.named_parameters():
#     if "actor" not in name:
#         param.requires_grad = False  
# student = PPO("MlpPolicy", env , verbose=1, gamma=distf, device="cuda",
#             policy_kwargs={"net_arch": [256,256],})
# student = TD3("MlpPolicy", env , verbose=1, gamma=distf, device="cuda",
#             policy_kwargs={"net_arch": [32,32],})

eval_env = make_vec_env("pensionfund-v3",seed=3687851522, n_envs=1, # env_kwargs={'params': params},\
                        monitor_kwargs={'discount_factor': distf})         
eval_env = VecNormalize(eval_env, norm_obs=False, norm_reward=False)


before_mean_reward, before_std_reward, _, _ = evaluate_policy(student.policy, env, eval_env ,seed=3687851522, n_eval_episodes=200,deterministic=True, render=False, discount_factor=distf)

# Having defined the training procedure we can now run the pretraining!

print('start pretrain....')
pretrain_agent(
    student,
    epochs=300,
    scheduler_gamma=1,
    learning_rate=0.0001,
    log_interval=1000,
    cuda=True,
    seed=1,
    batch_size=256,
    test_batch_size=1000,
)
# a2c_student.save("a2c_student")

# Finally, let us test how well our RL agent student learned to mimic the behavior of the expert
mean_reward, std_reward, _, _ = evaluate_policy(student.policy, env, eval_env,seed=3687851522, n_eval_episodes=200,deterministic=True, render=False, discount_factor=distf)

print(f"before BC Mean reward = {before_mean_reward}")
print(f"after BC Mean reward = {mean_reward}")
# student.policy.save('./im_models/trivial_sac_256.pkl')
# env.save('./im_models/trivial_sac_256_vec_env.pkl') 
