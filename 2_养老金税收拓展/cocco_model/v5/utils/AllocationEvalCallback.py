import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict, Union
import gymnasium as gym

from stable_baselines3.common.callbacks import EventCallback
from stable_baselines3.common.vec_env import VecEnv, sync_envs_normalization, DummyVecEnv
from stable_baselines3.common.logger import Figure

# 解决OMP: Error #15:
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

# def evaluate_policy(
#     model: "type_aliases.PolicyPredictor",
#     env: Union[gym.Env, VecEnv],
#     n_eval_episodes: int = 100,
#     deterministic: bool = True,
#     seed: int=42
# ) -> Dict[str, any]:
#     """
#     通过直接与Gym环境交互（reset/step）来评估RL智能体。
#     **修正版: 确保使用原始的、未归一化的奖励来计算真实终身效用。**
#     """
#     if not isinstance(env, VecEnv):
#         env = DummyVecEnv([lambda: env])
        
#     all_episode_data = []
#     episode_rewards = []
#     beta = env.get_attr("beta")[0]
#     np.random.seed(seed)

#     for ii in range(n_eval_episodes):
#         obs = env.reset()
#         done = False
        
#         trajectory_data = {
#             'age': [], 'wealth': [], 'consumption': [], 
#             'income': [], 'alpha': [],'c_prop':[]
#         }
        
#         current_discounted_reward = 0.0
#         time_step = 0
        
#         while not done:
#             action, _ = model.predict(obs, deterministic=deterministic) # 
#             # 这里的 obs 是 VecNormalize 自动归一化后的
#             obs, reward, done, info = env.step(action)
            
#             # **核心修正: 使用 env.get_original_reward() 获取真实的、未归一化的效用值**
#             # 这对于被 VecNormalize 封装的环境至关重要。
#             original_reward = env.get_original_reward()
            
#             info_dict = info[0]
            
#             trajectory_data['age'].append(info_dict['age'])
#             trajectory_data['wealth'].append(info_dict['absolute_wealth'])
#             trajectory_data['consumption'].append(info_dict['absolute_consumption'])
#             trajectory_data['income'].append(info_dict['absolute_income'])
#             trajectory_data['alpha'].append(info_dict['risky_share'])
#             trajectory_data['c_prop'].append(info_dict['c_prop'])
            
#             # 使用真实的效用值来累加终身回报
#             current_discounted_reward +=  (beta ** time_step) * original_reward[0] # 
#             time_step += 1

#         episode_rewards.append(current_discounted_reward)
#         episode_df = pd.DataFrame(trajectory_data)
#         all_episode_data.append(episode_df)
        
#     combined_df = pd.concat(all_episode_data)
#     mean_profiles_df = combined_df.groupby('age').mean().reset_index()

#     return {
#         "mean_reward": np.mean(episode_rewards),
#         "std_reward": np.std(episode_rewards),
#         "mean_profiles": mean_profiles_df,
#     }

def evaluate_policy(
    model: "type_aliases.PolicyPredictor",
    env: Union[gym.Env, VecEnv],
    n_eval_episodes: int = 100,
    deterministic: bool = True,
    seed: int=42
) -> Dict[str, any]:
    """
    通过直接与Gym环境交互（reset/step）来评估RL智能体。
    **修正版: 支持并行评估 (n_envs > 1) 并保证可复现性。**
    """
    if not isinstance(env, VecEnv):
        env = DummyVecEnv([lambda: env])

    n_envs = env.num_envs
    all_episode_data = []
    episode_rewards = []
    beta = env.get_attr("beta")[0]
    gamma = env.get_attr("gamma")[0]

    # --- [MODIFIED FOR PARALLEL EVALUATION] ---
    # 1. 为每个并行环境设置可复现的种子
    env.seed(seed)
    obs = env.reset()

    # 2. 为每个环境初始化独立的数据追踪器
    trajectories_data = [
        {'age': [], 'wealth': [], 'consumption': [], 'income': [], 'alpha': [], 'c_prop': []}
        for _ in range(n_envs)
    ]
    current_discounted_rewards = np.zeros(n_envs)
    time_steps = np.zeros(n_envs, dtype=int)

    # 3. 循环直到收集到足够数量的 episodes
    while len(episode_rewards) < n_eval_episodes:
        action, _ = model.predict(obs, deterministic=deterministic)
        obs, _, dones, infos = env.step(action)
        
        # 核心修正: 使用 env.get_original_reward() 获取真实的、未归一化的效用值
        original_rewards = env.get_original_reward()

        for i in range(n_envs):
            info_dict = infos[i]

            # 记录当前步的数据
            trajectories_data[i]['age'].append(info_dict['age'])
            trajectories_data[i]['wealth'].append(info_dict['absolute_wealth'])
            trajectories_data[i]['consumption'].append(info_dict['absolute_consumption'])
            trajectories_data[i]['income'].append(info_dict['absolute_income'])
            trajectories_data[i]['alpha'].append(info_dict['risky_share'])
            trajectories_data[i]['c_prop'].append(info_dict['c_prop'])

            # 使用真实的效用值来累加终身回报
            current_discounted_rewards[i] += (beta ** time_steps[i]) * (info_dict['absolute_consumption'])**(1-gamma)/(1-gamma)
            time_steps[i] += 1

            # 如果一个 episode 结束
            if dones[i]:
                # 仅在还需要收集更多episodes时才保存结果，防止超出 n_eval_episodes
                if len(episode_rewards) < n_eval_episodes:
                    episode_rewards.append(current_discounted_rewards[i])
                    episode_df = pd.DataFrame(trajectories_data[i])
                    all_episode_data.append(episode_df)

                # 重置该特定环境的追踪器，为其新开始的 episode 做准备
                trajectories_data[i] = {
                    'age': [], 'wealth': [], 'consumption': [], 
                    'income': [], 'alpha': [], 'c_prop': []
                }
                current_discounted_rewards[i] = 0.0
                time_steps[i] = 0
    # --- [END MODIFICATION] ---

    combined_df = pd.concat(all_episode_data)
    mean_profiles_df = combined_df.groupby('age').mean().reset_index()

    return {
        "mean_reward": np.mean(episode_rewards),
        "std_reward": np.std(episode_rewards),
        "mean_profiles": mean_profiles_df,
    }

# ==============================================================================
# 自定义评估回调类 (与之前版本逻辑一致)
# ==============================================================================
class EvalCallback(EventCallback):
    """
    一个专门用于评估生命周期模型的Callback。
    它会定期调用 `evaluate_policy` 函数，记录评估结果，
    并绘制和保存生命周期剖面图。
    """
    def __init__(
        self,
        eval_env: VecEnv,
        n_eval_episodes: int = 500,
        eval_freq: int = 10000,
        log_path: str = None,
        best_model_save_path: str = None,
        deterministic: bool = True,
        verbose: int = 1,
    ):
        super().__init__(None, verbose=verbose)
        self.eval_env = eval_env
        self.n_eval_episodes = n_eval_episodes
        self.eval_freq = eval_freq
        self.best_mean_reward = -np.inf
        self.deterministic = deterministic
        self.log_path = log_path
        self.best_model_save_path = best_model_save_path

    def _init_callback(self) -> None:
        if self.best_model_save_path is not None:
            os.makedirs(self.best_model_save_path, exist_ok=True)
        if self.log_path is not None:
            os.makedirs(self.log_path, exist_ok=True)

    def _on_step(self) -> bool:
        if self.n_calls > 0 and self.n_calls % self.eval_freq == 0:
            
            sync_envs_normalization(self.training_env, self.eval_env)
            
            eval_results = evaluate_policy(
                self.model, self.eval_env,
                n_eval_episodes=self.n_eval_episodes,
                deterministic=self.deterministic
            )
            
            mean_reward = eval_results["mean_reward"]
            std_reward = eval_results["std_reward"]
            rl_profiles_df = eval_results["mean_profiles"]

            if self.verbose > 0:
                print(f"\nEval num_timesteps={self.num_timesteps}, "
                      f"mean_lifetime_utility={mean_reward:.4f} +/- {std_reward:.4f}")

            self.logger.record("eval/mean_reward", mean_reward)
            self.logger.dump(step=self.num_timesteps)
            
            if mean_reward > self.best_mean_reward:
                if self.verbose > 0: print("New best mean lifetime utility!")
                if self.best_model_save_path is not None:
                    self.model.save(os.path.join(self.best_model_save_path, "best_model"))
                    stats_path = os.path.join(self.log_path, "vecnormalize_best.pkl")
                    # self.model.save_replay_buffer(os.path.join(self.best_model_save_path, "best_model_rb"))
                    self.eval_env.save(stats_path)
                self.best_mean_reward = mean_reward

                # --- [MODIFIED] 绘图逻辑与 MATLAB 的 AgentEvaluator.m 对齐 ---
                
                # 计算 c_prop
                # 添加一个小常数以避免除以零

                fig, axes = plt.subplots(2, 2, figsize=(14, 10))
                retirement_age = self.eval_env.get_attr("tr")[0] - 1 # 与MATLAB的retirement_age_line对齐

                # --- 1. 子图: 消费策略 (c_prop) ---
                ax1 = axes[0, 0]
                ax1.plot(rl_profiles_df['age'], rl_profiles_df['c_prop'], 'r--', linewidth=2, label='c-prop (RL)')
                ax1.axvline(x=retirement_age, color='grey', linestyle='--', label='Retirement')
                ax1.set_title('c-prop')
                ax1.set_xlabel('Age')
                ax1.set_ylabel('Proportion of Cash-on-Hand')
                ax1.legend(loc='best')
                ax1.grid(True)
                ax1.set_xlim(rl_profiles_df['age'].min(), rl_profiles_df['age'].max())
                ax1.set_ylim(-0.05, 1.05)
                
                # --- 2. 子图: 投资策略 (alpha) ---
                ax2 = axes[0, 1]
                ax2.plot(rl_profiles_df['age'], rl_profiles_df['alpha'], 'm--', linewidth=2, label='Risky Share (RL)')
                ax2.axvline(x=retirement_age, color='grey', linestyle='--', label='Retirement')
                ax2.set_title('Investment Policy (alpha)')
                ax2.set_xlabel('Age')
                ax2.set_ylabel('Share in Risky Asset')
                ax2.legend(loc='best')
                ax2.grid(True)
                ax2.set_xlim(rl_profiles_df['age'].min(), rl_profiles_df['age'].max())
                ax2.set_ylim(-0.05, 1.05)

                # --- 3. 子图: 宏观变量 (财富) ---
                ax3 = axes[1, 0]
                ax3.plot(rl_profiles_df['age'], rl_profiles_df['wealth'], 'b--', linewidth=2, label='Wealth (RL)')
                ax3.set_title('Mean Lifecycle Profiles (Absolute Units)')
                ax3.set_xlabel('Age')
                ax3.set_ylabel('Value')
                ax3.legend(loc='upper left') # MODIFIED: 'northwest' -> 'upper left'
                ax3.grid(True)
                ax3.set_xlim(rl_profiles_df['age'].min(), rl_profiles_df['age'].max())

                # --- 4. 隐藏未使用的子图 ---
                axes[1, 1].axis('off')
                
                # --- 总标题 ---
                base_title = "RL Agent Evaluation (Saved in " + self.best_model_save_path.split('/')[-1] +")"
                utility_title = f"RL Utility: {mean_reward:.4f}"
                fig.suptitle(f"{base_title}\n{utility_title}", fontsize=14)
                
                fig.tight_layout(rect=[0, 0.03, 1, 0.94])
                
                self.logger.record("eval/lifecycle_profiles", Figure(fig, close=True), exclude=("stdout", "log", "json", "csv"))
                plt.close(fig)

        return True