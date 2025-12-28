import os
import warnings
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Optional, Union
import gymnasium as gym
import numpy as np
from stable_baselines3.common.logger import Logger
from stable_baselines3.common.vec_env import DummyVecEnv, VecEnv, sync_envs_normalization
from stable_baselines3.common.callbacks import EventCallback,BaseCallback
from stable_baselines3.common.logger import Figure
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import gymnasium as gym
from stable_baselines3.common import type_aliases
from stable_baselines3.common.vec_env import VecEnv, VecMonitor, is_vecenv_wrapped

# 解决OMP: Error #15:
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"

def evaluate_policy(
    model: "type_aliases.PolicyPredictor",
    vec_env: Union[gym.Env, VecEnv],
    env: Union[gym.Env, VecEnv],
    seed: int = None,
    n_eval_episodes: int = 10,
    deterministic: bool = True,
    render: bool = False,
    callback: Optional[Callable[[Dict[str, Any], Dict[str, Any]], None]] = None,
    reward_threshold: Optional[float] = None,
    return_episode_rewards: bool = False,
    warn: bool = True,
    discount_factor: float = 1,
    params: dict = None
) -> Union[Tuple[float, float], Tuple[List[float], List[int]]]:
    """
    Runs policy for ``n_eval_episodes`` episodes and returns average reward.
    If a vector env is passed in, this divides the episodes to evaluate onto the
    different elements of the vector env. This static division of work is done to
    remove bias. See https://github.com/DLR-RM/stable-baselines3/issues/402 for more
    details and discussion.

    .. note::
        If environment has not been wrapped with ``Monitor`` wrapper, reward and
        episode lengths are counted as it appears with ``env.step`` calls. If
        the environment contains wrappers that modify rewards or episode lengths
        (e.g. reward scaling, early episode reset), these will affect the evaluation
        results as well. You can avoid this by wrapping environment with ``Monitor``
        wrapper before anything else.

    :param model: The RL agent you want to evaluate. This can be any object
        that implements a `predict` method, such as an RL algorithm (``BaseAlgorithm``)
        or policy (``BasePolicy``).
    :param env: The gym environment or ``VecEnv`` environment.
    :param n_eval_episodes: Number of episode to evaluate the agent
    :param deterministic: Whether to use deterministic or stochastic actions
    :param render: Whether to render the environment or not
    :param callback: callback function to do additional checks,
        called after each step. Gets locals() and globals() passed as parameters.
    :param reward_threshold: Minimum expected reward per episode,
        this will raise an error if the performance is not met
    :param return_episode_rewards: If True, a list of rewards and episode lengths
        per episode will be returned instead of the mean.
    :param warn: If True (default), warns user about lack of a Monitor wrapper in the
        evaluation environment.
    :return: Mean reward per episode, std of reward per episode.
        Returns ([float], [int]) when ``return_episode_rewards`` is True, first
        list containing per-episode rewards and second containing per-episode lengths
        (in number of steps).
    """
    is_monitor_wrapped = False
    # Avoid circular import
    from stable_baselines3.common.monitor import Monitor

    if not isinstance(env, VecEnv):
        env = DummyVecEnv([lambda: env])  # type: ignore[list-item, return-value]

    is_monitor_wrapped = is_vecenv_wrapped(env, VecMonitor) or env.env_is_wrapped(Monitor)[0]

    if not is_monitor_wrapped and warn:
        warnings.warn(
            "Evaluation environment is not wrapped with a ``Monitor`` wrapper. "
            "This may result in reporting modified episode lengths and rewards, if other wrappers happen to modify these. "
            "Consider wrapping environment first with ``Monitor`` wrapper.",
            UserWarning,
        )
    
    n_envs = env.num_envs
    episode_rewards = []
    episode_lengths = []

    # Divides episodes among different sub environments in the vector as evenly as possible
    episode_count_targets = np.array([(n_eval_episodes + i) // n_envs for i in range(n_envs)], dtype="int")

    current_rewards = np.zeros(n_envs)
    current_lengths = np.zeros(n_envs, dtype="int")

    scale_reward = env.get_attr('scale_reward')[0]
    age_set = []
    action_set = []
    # 存储所有agent的轨迹
    # bound_set = []

    for episode in range(episode_count_targets[0]):       
        env.seed(seed=seed+episode)
        observations = env.reset()    
        dones = False
        current_rewards = 0
        current_lengths = 0
        while not dones:
            actions, states = model.predict(
            vec_env.normalize_obs(observations),  # type: ignore[arg-type]
            # state=states,
            # episode_start=episode_starts,
            deterministic=deterministic,
        )
            # unnorm_obs = env.unnormalize_obs(observations) # 反归一化
            # age_set += [round(unnorm_obs[0][-1],1)] # 取整为整数
            age_set += [env.get_attr('age')[0]] 
            new_observations, rewards, dones, infos = env.step(actions)
            action_set += [infos[0]['real_actions']]
            current_rewards += rewards/scale_reward * discount_factor**current_lengths
            current_lengths += 1
            observations = new_observations
        episode_rewards.append(current_rewards)
        episode_lengths.append(current_lengths)

    age_set =  np.vstack(age_set)
    # age = np.round((age_set + 1 )*(params['max_age']-params['born_age'])/2+params['born_age'],1)
    age = age_set 
    action_set = np.vstack(action_set)
    
    # bound_set = np.vstack(bound_set)
    unique_age = np.unique(age)
    mean_actions = np.zeros((len(unique_age),np.shape(action_set)[1]+1))    # mean_bound = np.zeros((len(unique_age),1))
    for idx in range(len(unique_age)):
        mean_actions[idx,0] = unique_age[idx]
        mean_actions[idx,1:] = np.mean(action_set[(age==unique_age[idx]).squeeze(),:],0)
        
    # 计算折现后的平均奖励:
    # 1. discount_factor**np.arange(len(episode_rewards)) 生成折现因子序列,如[1, df, df^2, df^3, ...]
    # 2. episode_rewards 乘以折现因子序列,得到折现后的奖励序列
    # 3. np.mean() 计算折现后奖励序列的平均值
    mean_reward = np.mean(episode_rewards)
    std_reward = np.std(episode_rewards)
    if reward_threshold is not None:
        assert mean_reward > reward_threshold, "Mean reward below threshold: " f"{mean_reward:.2f} < {reward_threshold:.2f}"
    if return_episode_rewards:
        return episode_rewards, episode_lengths, std_reward 
    return mean_reward, episode_lengths, std_reward,  mean_actions


def evaluate_policy_dp(
    model: "type_aliases.PolicyPredictor",
    vec_env: Union[gym.Env, VecEnv],
    env: Union[gym.Env, VecEnv],
    seed: int = None,
    n_eval_episodes: int = 10,
    deterministic: bool = True,
    render: bool = False,
    callback: Optional[Callable[[Dict[str, Any], Dict[str, Any]], None]] = None,
    reward_threshold: Optional[float] = None,
    return_episode_rewards: bool = False,
    warn: bool = True,
    discount_factor: float = 1,
    params: dict = None
) -> Union[Tuple[float, float], Tuple[List[float], List[int]]]:
    """
    Runs policy for ``n_eval_episodes`` episodes and returns average reward.
    If a vector env is passed in, this divides the episodes to evaluate onto the
    different elements of the vector env. This static division of work is done to
    remove bias. See https://github.com/DLR-RM/stable-baselines3/issues/402 for more
    details and discussion.

    .. note::
        If environment has not been wrapped with ``Monitor`` wrapper, reward and
        episode lengths are counted as it appears with ``env.step`` calls. If
        the environment contains wrappers that modify rewards or episode lengths
        (e.g. reward scaling, early episode reset), these will affect the evaluation
        results as well. You can avoid this by wrapping environment with ``Monitor``
        wrapper before anything else.

    :param model: The RL agent you want to evaluate. This can be any object
        that implements a `predict` method, such as an RL algorithm (``BaseAlgorithm``)
        or policy (``BasePolicy``).
    :param env: The gym environment or ``VecEnv`` environment.
    :param n_eval_episodes: Number of episode to evaluate the agent
    :param deterministic: Whether to use deterministic or stochastic actions
    :param render: Whether to render the environment or not
    :param callback: callback function to do additional checks,
        called after each step. Gets locals() and globals() passed as parameters.
    :param reward_threshold: Minimum expected reward per episode,
        this will raise an error if the performance is not met
    :param return_episode_rewards: If True, a list of rewards and episode lengths
        per episode will be returned instead of the mean.
    :param warn: If True (default), warns user about lack of a Monitor wrapper in the
        evaluation environment.
    :return: Mean reward per episode, std of reward per episode.
        Returns ([float], [int]) when ``return_episode_rewards`` is True, first
        list containing per-episode rewards and second containing per-episode lengths
        (in number of steps).
    """
    is_monitor_wrapped = False
    # Avoid circular import
    from stable_baselines3.common.monitor import Monitor

    if not isinstance(env, VecEnv):
        env = DummyVecEnv([lambda: env])  # type: ignore[list-item, return-value]

    is_monitor_wrapped = is_vecenv_wrapped(env, VecMonitor) or env.env_is_wrapped(Monitor)[0]

    if not is_monitor_wrapped and warn:
        warnings.warn(
            "Evaluation environment is not wrapped with a ``Monitor`` wrapper. "
            "This may result in reporting modified episode lengths and rewards, if other wrappers happen to modify these. "
            "Consider wrapping environment first with ``Monitor`` wrapper.",
            UserWarning,
        )
    
    n_envs = env.num_envs
    episode_rewards = []
    episode_lengths = []

    # Divides episodes among different sub environments in the vector as evenly as possible
    episode_count_targets = np.array([(n_eval_episodes + i) // n_envs for i in range(n_envs)], dtype="int")

    current_rewards = np.zeros(n_envs)
    current_lengths = np.zeros(n_envs, dtype="int")

    scale_reward = env.get_attr('scale_reward')[0]
    age_set = []
    action_set = []
    # 存储所有agent的轨迹
    # bound_set = []

    for episode in range(episode_count_targets[0]):       
        env.seed(seed=seed+episode)
        observations = env.reset()    
        dones = False
        current_rewards = 0
        current_lengths = 0
        while not dones:
            actions = model.predict(observations,env)  # type: ignore[arg-type]
            # unnorm_obs = env.unnormalize_obs(observations) # 反归一化
            # age_set += [round(unnorm_obs[0][-1],1)] # 取整为整数
            age_set += [env.get_attr('age')[0]] 
            new_observations, rewards, dones, infos = env.step(actions)
            action_set += [infos[0]['real_actions']]
            current_rewards += rewards/scale_reward * discount_factor**current_lengths
            current_lengths += 1
            observations = new_observations
        episode_rewards.append(current_rewards)
        episode_lengths.append(current_lengths)

    age_set =  np.vstack(age_set)
    # age = np.round((age_set + 1 )*(params['max_age']-params['born_age'])/2+params['born_age'],1)
    age = age_set 
    action_set = np.vstack(action_set)
    
    # bound_set = np.vstack(bound_set)
    unique_age = np.unique(age)
    mean_actions = np.zeros((len(unique_age),np.shape(action_set)[1]+1))    # mean_bound = np.zeros((len(unique_age),1))
    for idx in range(len(unique_age)):
        mean_actions[idx,0] = unique_age[idx]
        mean_actions[idx,1:] = np.mean(action_set[(age==unique_age[idx]).squeeze(),:],0)


        
    # 计算折现后的平均奖励:
    # 1. discount_factor**np.arange(len(episode_rewards)) 生成折现因子序列,如[1, df, df^2, df^3, ...]
    # 2. episode_rewards 乘以折现因子序列,得到折现后的奖励序列
    # 3. np.mean() 计算折现后奖励序列的平均值
    mean_reward = np.mean(episode_rewards)
    std_reward = np.std(episode_rewards)
    if reward_threshold is not None:
        assert mean_reward > reward_threshold, "Mean reward below threshold: " f"{mean_reward:.2f} < {reward_threshold:.2f}"
    if return_episode_rewards:
        return episode_rewards, episode_lengths, std_reward 
    return mean_reward, episode_lengths, std_reward, mean_actions