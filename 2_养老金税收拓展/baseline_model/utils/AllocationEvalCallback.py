import os
import warnings
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Optional, Union

import gymnasium as gym
import numpy as np

from stable_baselines3.common.logger import Logger

from stable_baselines3.common.vec_env import DummyVecEnv, VecEnv, sync_envs_normalization
from stable_baselines3.common.callbacks import EventCallback,BaseCallback
from stable_baselines3.common.logger import Figure

import json

from typing import Any, Callable, Dict, List, Optional, Tuple, Union
import matplotlib.pyplot as plt

import gymnasium as gym

from stable_baselines3.common import type_aliases
from stable_baselines3.common.vec_env import VecEnv, VecMonitor, is_vecenv_wrapped

# 解决OMP: Error #15:
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"

def nonlin_schedule(initial_value,progress_remaining):
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

def evaluate_policy(
    model: "type_aliases.PolicyPredictor",
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

    episode_counts = np.zeros(n_envs, dtype="int")
    # Divides episodes among different sub environments in the vector as evenly as possible
    episode_count_targets = np.array([(n_eval_episodes + i) // n_envs for i in range(n_envs)], dtype="int")

    current_rewards = np.zeros(n_envs)
    current_lengths = np.zeros(n_envs, dtype="int")

    # scale_reward = env.get_attr('scale_reward')[0]
    states = None
    episode_starts = np.ones((env.num_envs,), dtype=bool)
    action_set = []
    age_set = []
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
            model.env.normalize_obs(observations),  # type: ignore[arg-type]
            # state=states,
            # episode_start=episode_starts,
            deterministic=deterministic,
        )
            # unnorm_obs = env.unnormalize_obs(observations) # 反归一化
            # age_set += [round(unnorm_obs[0][-1],1)] # 取整为整数
            age_set += [env.get_attr('age')[0]] 
            new_observations, rewards, dones, infos = env.step(actions)
            action_set += [infos[0]['real_actions']]
            current_rewards += rewards * discount_factor**current_lengths
            current_lengths += 1
            observations = new_observations
        episode_rewards.append(current_rewards)
        episode_lengths.append(current_lengths)
    # while (episode_counts < episode_count_targets).any():
    #     env.seed(seed=seed+episode_counts)
    #     observations = env.reset()
    #     actions, states = model.predict(
    #         observations,  # type: ignore[arg-type]
    #         state=states,
    #         episode_start=episode_starts,
    #         deterministic=deterministic,
    #     )
        
        # # age_set+= [observations[0][3]]
        
        # unnorm_obs = env.unnormalize_obs(observations) # 反归一化
        # age_set += [round(unnorm_obs[0][1],1)] # 取整为整数
        
        # new_observations, rewards, dones, infos = env.step(actions)
        # action_set += [infos[0]['real_actions']]
        # current_rewards += rewards * discount_factor**current_lengths
        # current_lengths += 1
        # for i in range(n_envs):
        #     if episode_counts[i] < episode_count_targets[i]:
        #         # unpack values so that the callback can access the local variables
        #         reward = rewards[i]
        #         done = dones[i]
        #         info = infos[i]
        #         episode_starts[i] = done

        #         if callback is not None:
        #             callback(locals(), globals())

        #         if dones[i]:
        #             if is_monitor_wrapped:
        #                 # Atari wrapper can send a "done" signal when
        #                 # the agent loses a life, but it does not correspond
        #                 # to the true end of episode
        #                 if "episode" in info.keys():
        #                     # Do not trust "done" with episode endings.
        #                     # Monitor wrapper includes "episode" key in info if environment
        #                     # has been wrapped with it. Use those rewards instead.
        #                     episode_rewards.append(info["episode"]["r"])
        #                     episode_lengths.append(info["episode"]["l"])
        #                     # Only increment at the real end of an episode
        #                     episode_counts[i] += 1
        #             else:
        #                 episode_rewards.append(current_rewards[i])
        #                 episode_lengths.append(current_lengths[i])
        #                 episode_counts[i] += 1
        #             current_rewards[i] = 0
        #             current_lengths[i] = 0

        # observations = new_observations

        # if render:
        #     env.render()
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
        # mean_bound[idx] = np.mean(bound_set[(age==unique_age[idx]).squeeze()])

        
        
    # 计算折现后的平均奖励:
    # 1. discount_factor**np.arange(len(episode_rewards)) 生成折现因子序列,如[1, df, df^2, df^3, ...]
    # 2. episode_rewards 乘以折现因子序列,得到折现后的奖励序列
    # 3. np.mean() 计算折现后奖励序列的平均值
    mean_reward = np.mean(episode_rewards)
    std_reward = np.std(episode_rewards)
    if reward_threshold is not None:
        assert mean_reward > reward_threshold, "Mean reward below threshold: " f"{mean_reward:.2f} < {reward_threshold:.2f}"
    if return_episode_rewards:
        return episode_rewards, episode_lengths, std_reward , mean_actions
    return mean_reward, episode_lengths, std_reward, mean_actions


class EvalCallback(EventCallback):
    """
    Callback for evaluating an agent.

    .. warning::

      When using multiple environments, each call to  ``env.step()``
      will effectively correspond to ``n_envs`` steps.
      To account for that, you can use ``eval_freq = max(eval_freq // n_envs, 1)``

    :param eval_env: The environment used for initialization
    :param callback_on_new_best: Callback to trigger
        when there is a new best model according to the ``mean_reward``
    :param callback_after_eval: Callback to trigger after every evaluation
    :param n_eval_episodes: The number of episodes to test the agent
    :param eval_freq: Evaluate the agent every ``eval_freq`` call of the callback.
    :param log_path: Path to a folder where the evaluations (``evaluations.npz``)
        will be saved. It will be updated at each evaluation.
    :param best_model_save_path: Path to a folder where the best model
        according to performance on the eval env will be saved.
    :param deterministic: Whether the evaluation should
        use a stochastic or deterministic actions.
    :param render: Whether to render or not the environment during evaluation
    :param verbose: Verbosity level: 0 for no output, 1 for indicating information about evaluation results
    :param warn: Passed to ``evaluate_policy`` (warns if ``eval_env`` has not been
        wrapped with a Monitor wrapper)
    """

    def __init__(
        self,
        eval_env: Union[gym.Env, VecEnv],
        seed: int = None,
        callback_on_new_best: Optional[BaseCallback] = None,
        callback_after_eval: Optional[BaseCallback] = None,
        n_eval_episodes: int = 5,
        eval_freq: int = 10000,
        log_path: Optional[str] = None,
        best_model_save_path: Optional[str] = None,
        deterministic: bool = True,
        render: bool = False,
        verbose: int = 1,
        warn: bool = True,
        discount_factor: float = 1,
        params: dict = None,
        score_target: float = 999,
    ):
        super().__init__(callback_after_eval, verbose=verbose)
        self.n_evals = 0
        self.callback_on_new_best = callback_on_new_best
        if self.callback_on_new_best is not None:
            # Give access to the parent
            self.callback_on_new_best.parent = self

        self.n_eval_episodes = n_eval_episodes
        self.eval_freq = eval_freq
        self.best_mean_reward = -np.inf
        self.last_mean_reward = -np.inf
        self.deterministic = deterministic
        self.render = render
        self.warn = warn
        self.discount_factor = discount_factor
        self.params = params
        self.seed = seed
        self.n_info = 0

        self.score_target = score_target
        self.last_best_model_n = 0 # 最近一次最佳模型出现的序号
        # Convert to VecEnv for consistency
        if not isinstance(eval_env, VecEnv):
            eval_env = DummyVecEnv([lambda: eval_env])  # type: ignore[list-item, return-value]

        self.eval_env = eval_env
        self.best_model_save_path = best_model_save_path
        # Logs will be written in ``evaluations.npz``
        if log_path is not None:
            log_path = os.path.join(log_path, "evaluations")
        self.log_path = log_path
        self.evaluations_results: List[List[float]] = []
        self.evaluations_timesteps: List[int] = []
        self.evaluations_length: List[List[int]] = []
        # For computing success rate
        self._is_success_buffer: List[bool] = []
        self.evaluations_successes: List[List[bool]] = []

    def _init_callback(self) -> None:
        # Does not work in some corner cases, where the wrapper is not the same
        if not isinstance(self.training_env, type(self.eval_env)):
            warnings.warn("Training and eval env are not of the same type" f"{self.training_env} != {self.eval_env}")

        # Create folders if needed
        if self.best_model_save_path is not None:
            os.makedirs(self.best_model_save_path, exist_ok=True)
        if self.log_path is not None:
            os.makedirs(os.path.dirname(self.log_path), exist_ok=True)

        # Init callback called on new best model
        if self.callback_on_new_best is not None:
            self.callback_on_new_best.init_callback(self.model)

    def _log_success_callback(self, locals_: Dict[str, Any], globals_: Dict[str, Any]) -> None:
        """
        Callback passed to the  ``evaluate_policy`` function
        in order to log the success rate (when applicable),
        for instance when using HER.

        :param locals_:
        :param globals_:
        """
        info = locals_["info"]

        if locals_["done"]:
            maybe_is_success = info.get("is_success")
            if maybe_is_success is not None:
                self._is_success_buffer.append(maybe_is_success)

    # def _on_training_start(self) -> None:
        
    #     metric_dict = {
    #             "eval/mean_reward": 0.0,
    #         }
    #     # 将环境名称和参数一起记录
    #     env_name = self.eval_env.unwrapped.spec.id if hasattr(self.eval_env.unwrapped, 'spec') else 'Unknown'
    #     self.logger.record("env_name", env_name)
        
    #     for key, value in self.params.items():
    #         self.logger.record(f"param/{key}", value)
            
        
    def _on_step(self) -> bool:
        continue_training = True
         
        # 学习率指数衰减
        # print(f"{self.model.learning_rate:.10f}")
        if self.n_calls == 1:
            self.lr_init = self.model.learning_rate

        # self.logger.record("train/lr", float(self.model.learning_rate))
        try:
            self.logger.record("train/lr", float(self.model.learning_rate(self.num_timesteps)))
        except:
            self.logger.record("train/lr", float(self.model.learning_rate))


        if np.isnan(self.model.predict(self.eval_env.reset())[0]).any(): # 如果model输出任何nan，提前结束
            print('nan detected, early stop!!!')
            return False
      
        if self.eval_freq > 0 and self.n_calls % self.eval_freq == 0 and self.num_timesteps>self.model.learning_starts:
            # Sync training and eval env if there is VecNormalize
            if self.model.get_vec_normalize_env() is not None:
                try:
                    sync_envs_normalization(self.training_env, self.eval_env)
                except AttributeError as e:
                    raise AssertionError(
                        "Training and eval env are not wrapped the same way, "
                        "see https://stable-baselines3.readthedocs.io/en/master/guide/callbacks.html#evalcallback "
                        "and warning above."
                    ) from e

            # Reset success rate buffer
            self._is_success_buffer = []

            episode_rewards, episode_lengths, _, mean_actions = evaluate_policy(
                self.model,
                self.eval_env,
                seed=self.seed,
                n_eval_episodes=self.n_eval_episodes,
                render=self.render,
                deterministic=self.deterministic,
                return_episode_rewards=True,
                warn=self.warn,
                callback=self._log_success_callback,
                discount_factor = self.discount_factor,
                params=self.params
            )
            

            
            if self.log_path is not None:
                assert isinstance(episode_rewards, list)
                assert isinstance(episode_lengths, list)
                self.evaluations_timesteps.append(self.num_timesteps)
                self.evaluations_results.append(episode_rewards)
                self.evaluations_length.append(episode_lengths)

                kwargs = {}
                # Save success log if present
                if len(self._is_success_buffer) > 0:
                    self.evaluations_successes.append(self._is_success_buffer)
                    kwargs = dict(successes=self.evaluations_successes)

                np.savez(
                    self.log_path,
                    timesteps=self.evaluations_timesteps,
                    results=self.evaluations_results,
                    ep_lengths=self.evaluations_length,
                    **kwargs,
                )
            mean_reward = np.mean(episode_rewards)
            std_reward = np.std(episode_rewards)
            mean_ep_length, std_ep_length = np.mean(episode_lengths), np.std(episode_lengths)
            self.last_mean_reward = float(mean_reward)

            if self.verbose >= 1:
                print(f"Eval num_timesteps={self.num_timesteps}, " f"episode_reward={mean_reward:.2f} +/- {std_reward:.2f}")
                print(f"Episode length: {mean_ep_length:.2f} +/- {std_ep_length:.2f}")
                
            # Add to current Logger
            self.logger.record("eval/mean_reward", float(mean_reward))
            self.logger.record("eval/mean_length", mean_ep_length)

            if len(self._is_success_buffer) > 0:
                success_rate = np.mean(self._is_success_buffer)
                if self.verbose >= 1:
                    print(f"Success rate: {100 * success_rate:.2f}%")
                self.logger.record("eval/success_rate", success_rate)

            # Dump log so the evaluation results are printed with the correct timestep
            self.logger.record("time/total_timesteps", self.num_timesteps, exclude="tensorboard")
            self.logger.dump(self.num_timesteps)
            
            if self.n_info == 0:
                self.logger.record("models/policy_kwargs", str(self.model.policy_kwargs))
                self.logger.record("models/lr", str(self.model.learning_rate))
                self.logger.record("models/gamma", str(self.model.gamma))
                self.logger.record("models/batch_size", str(self.model.batch_size))
                self.logger.record("models/buffer_size", str(self.model.buffer_size))
                self.logger.record("models/ent_coef", str(self.model.ent_coef_init))
                self.logger.record("models/tau", str(self.model.tau))
                self.logger.record("models/use_sde", str(self.model.use_sde))
                self.logger.record("models/gradient_steps", str(self.model.gradient_steps))
                self.logger.record("models/action_noise", str(self.model.action_noise)) 
                self.logger.record("models/runs",str(self.best_model_save_path.split('_')[-1]))
                self.n_info = 1 
            # 如果连续超过10次评估没有出现更好的模型且已经超过基准得分,那么开始降低学习率,每次降低90%
            # if self.model.num_timesteps>self.model.learning_starts \
            #     and self.n_evals - self.last_best_model_n>10 \
            #         and self.best_mean_reward >= self.score_target:    
            #     # self.model.learning_rate = nonlin_schedule(self.lr_init,self.model._current_progress_remaining)
            #     self.last_best_model_n = self.n_evals # 重新累积
            #     self.model.learning_rate = self.model.learning_rate * 0.1
                # self.model.batch_size = self.model.batch_size * 0.5
            
            if mean_reward >= self.best_mean_reward:
                self.last_best_model_n = self.n_evals
                if self.verbose >= 1:
                    print("New best mean reward!")
                if self.best_model_save_path is not None:
                    self.model.save(os.path.join(self.best_model_save_path, f"best_model"))
                    self.model.save_replay_buffer(os.path.join(self.best_model_save_path, f"best_model_replay_buffer"))
                    self.model.env.save(self.best_model_save_path+"//best_model_vec_env.pkl") 
                    # 保存params
                    if self.params is not None:
                        with open(self.best_model_save_path+"//params.json", 'w') as file:
                            json.dump(self.params, file)
                self.best_mean_reward = float(mean_reward)
                # Trigger callback on new best model, if needed
                if self.callback_on_new_best is not None:
                    continue_training = self.callback_on_new_best.on_step()  
                    
                # 增加action走势图
                runs = self.best_model_save_path.split('_')[-1]
                figure = plt.figure(figsize=(10, 6))
                figure.add_subplot()
                plt.plot(mean_actions[:,0], mean_actions[:,1], label='consum')
                plt.plot(mean_actions[:,0], mean_actions[:,2], label='riksy')   
                plt.plot(mean_actions[:,0], mean_actions[:,3], label='tda_pct')  
                plt.plot(mean_actions[:,0], mean_actions[:,4], label='riksy_tda')  
                plt.title(f"{runs}: Best Model {self.n_evals}: Timesteps={self.num_timesteps}, " f"episode_reward={mean_reward:.7f} ± {std_reward:.7f}")
                plt.xlabel('Age')
                plt.legend()
                plt.grid(True)
                self.logger.record("trajectory/figure", Figure(figure, close=True), exclude=("stdout", "log", "json", "csv"))
                plt.close()
                
            # else:
            #     # self.model.save(os.path.join(self.best_model_save_path, f"{self.n_evals}_model"))
            #     runs = self.best_model_save_path.split('_')[-1]
            #     # 增加action走勢圖
            #     figure = plt.figure(figsize=(10, 6))
            #     figure.add_subplot()
            #     plt.plot(mean_actions[:,0], mean_actions[:,1], label='consum')
            #     plt.plot(mean_actions[:,0], mean_actions[:,2], label='riksy') 
            #     plt.plot(mean_actions[:,0], mean_actions[:,3], label='tda_pct')  
            #     plt.plot(mean_actions[:,0], mean_actions[:,4], label='riksy_tda')    
            #     plt.title(f"{runs}: Model {self.n_evals}: Timesteps={self.num_timesteps}, " f"reward={mean_reward:.2f} ± {std_reward:.2f},  Best={self.best_mean_reward:.2f}")
            #     plt.xlabel('Age')
            #     plt.legend()
            #     plt.grid(True)
            #     self.logger.record("trajectory/figure", Figure(figure, close=True), exclude=("stdout", "log", "json", "csv"))
            #     plt.close()
            self.n_evals += 1
                
               
                

            # Trigger callback after every evaluation, if needed
            if self.callback is not None:
                continue_training = continue_training and self._on_event()

        return continue_training

    def update_child_locals(self, locals_: Dict[str, Any]) -> None:
        """
        Update the references to the local variables.

        :param locals_: the local variables during rollout collection
        """
        if self.callback:
            self.callback.update_locals(locals_)


