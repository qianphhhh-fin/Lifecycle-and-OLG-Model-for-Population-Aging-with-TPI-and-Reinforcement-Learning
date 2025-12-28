import gymnasium as gym
from elegantrl.train.config import Config
from elegantrl.agents.AgentSAC import AgentModSAC
from elegantrl.train.config import build_env
from elegantrl.train.evaluator import get_rewards_actions_and_steps
from gymnasium.envs.registration import register
import torch as th
import numpy as np
import os
import matplotlib.pyplot as plt

register(
    id='pensionfund-v3-norm',                                # call it whatever you want
    entry_point='v3_pensionfund_norm:PensionFundEnv', # module_name:class_name
)

if __name__ == '__main__':
    env = gym.make


    env_args = {
        "env_name": "pensionfund-v3-norm",
        "num_envs": 1,
        "max_step": 200,
        "state_dim": 4,
        "action_dim": 4,
        "if_discrete": False,
        'if_build_vec_env': False,
    }
    args = Config(AgentModSAC, env_class=env, env_args=env_args)
    args.max_step = 1000
    args.reward_scale = 1  # RewardRange: -1800 < -200 < -50 < 0
    args.gamma = 0.95
    args.target_step = args.max_step
    args.eval_per_step = int(1000) 
    args.eval_times = 200
    args.gpu_id = 0
    args.eval_seed = 3687851522
    


    env = build_env(env_class=env, env_args=env_args)

    eval_times = 200
    net_dim = [128,128]
    actor_path = 'E:/OneDrive - pku.edu.cn/工作台/博士论文/程序/养老金基准模型/pensionfund-v3-norm_ModSAC_0/actor__000000088064_00011.325.pt'

    state_dim = env_args['state_dim']
    action_dim = env_args['action_dim']
    agent = AgentModSAC(net_dim, state_dim, action_dim, gpu_id=args.gpu_id)
    act = agent.act
    act = th.load(actor_path, map_location=lambda storage, loc: storage, weights_only=False)

    r_s_ary = []
    steps = []
    action_set = []
    age_set = []
    for i in range(eval_times):
        score, n_steps, action , age = get_rewards_actions_and_steps(env, act, args.gamma, args.eval_seed + i)
        r_s_ary.append(score)
        steps.append(n_steps)
        action_set.extend(action)  # 将每个episode的action序列加入action_set
        age_set.extend(age)  # 将每个episode的age序列加入age_set

    age_set =  np.vstack(age_set)
    age = age_set 
    action_set = np.vstack(action_set)
    unique_age = np.unique(age)
    mean_actions = np.zeros((len(unique_age),np.shape(action_set)[1]+1))    # mean_bound = np.zeros((len(unique_age),1))
    for idx in range(len(unique_age)):
        mean_actions[idx,0] = unique_age[idx]
        mean_actions[idx,1:] = np.mean(action_set[(age==unique_age[idx]).squeeze(),:],0)

    r_s_ary = np.array(r_s_ary, dtype=np.float32)
    r_avg = r_s_ary.mean()  # average of episode return and episode step
    r_std = r_s_ary.std()  # average of episode return and episode step
    steps = np.array(steps, dtype=np.float32)
    steps_avg = steps.mean()  # average of episode return and episode step


    # 画图
    plt.figure(figsize=(10, 6))
    plt.plot(mean_actions[:,0], mean_actions[:,1], label='consum')
    plt.plot(mean_actions[:,0], mean_actions[:,2], label='riksy')   
    plt.plot(mean_actions[:,0], mean_actions[:,3], label='tda_pct')  
    plt.plot(mean_actions[:,0], mean_actions[:,4], label='riksy_tda')  
    plt.title( f"episode_reward={r_avg:.4f} ± {r_std:.4f}")
    plt.xlabel('Age')
    plt.legend()
    plt.grid(True)
    #
    # 保存图片
    save_path = actor_path.replace('.pt', '.png')
    plt.savefig(save_path)
    plt.show()
    print(f"图片已保存至: {save_path}")
