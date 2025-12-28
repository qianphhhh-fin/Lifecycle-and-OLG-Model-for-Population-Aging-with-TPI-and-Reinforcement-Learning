import gymnasium as gym
import time
import numpy as np

from gymnasium.envs.registration import register

register(
    id='cocco-rl',
    entry_point='cocco_env_rl:CoccoEnv',
)

def linear(obs, vec):
    control = np.sum(obs*vec)
    action = 1 if control>0 else 0
    return(action)

env = gym.make('cocco-rl',params=None)

best_vec, best_score = np.zeros(2), 0
num_draws = 50
for k in range(num_draws):
    vec = np.random.uniform(low=0, high=1, size=2)
    avg_reward = 0
    num_eval_eps = 1
    for i in range(num_eval_eps):
        ep_reward = 0
        obs = env.reset()[0]
        done = False
        while not done:
            action = linear(obs, vec)
            obs, reward, done, info = env.step(action)
            ep_reward += reward
            if done:
                avg_reward += ep_reward/num_eval_eps
                ep_reward = 0
    if avg_reward > best_score:
        best_score, best_vec = avg_reward, vec

env.close()

print('Best score {}\nBest vec {}'.format(best_score, best_vec))