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

def evaluate_policy(
    model: "type_aliases.PolicyPredictor",
    env: Union[gym.Env, VecEnv],
    n_eval_episodes: int = 100,
    deterministic: bool = True,
    seed: int=42
) -> Dict[str, any]:
    """
    通过直接与Gym环境交互（reset/step）来评估RL智能体。
    **修正版: 确保使用原始的、未归一化的奖励来计算真实终身效用。**
    """
    if not isinstance(env, VecEnv):
        env = DummyVecEnv([lambda: env])
        
    all_episode_data = []
    episode_rewards = []
    beta = env.get_attr("beta")[0]
    np.random.seed(seed)

    for ii in range(n_eval_episodes):
        obs = env.reset()
        done = False
        
        trajectory_data = {
            'age': [], 'wealth': [], 'consumption': [], 
            'income': [], 'alpha': []
        }
        
        current_discounted_reward = 0.0
        time_step = 0
        
        while not done:
            action, _ = model.predict(obs, deterministic=deterministic) # 
            # 这里的 obs 是 VecNormalize 自动归一化后的
            obs, reward, done, info = env.step(action)
            
            # **核心修正: 使用 env.get_original_reward() 获取真实的、未归一化的效用值**
            # 这对于被 VecNormalize 封装的环境至关重要。
            original_reward = env.get_original_reward()
            
            info_dict = info[0]
            
            trajectory_data['age'].append(info_dict['age'])
            trajectory_data['wealth'].append(info_dict['absolute_wealth'])
            trajectory_data['consumption'].append(info_dict['absolute_consumption'])
            trajectory_data['income'].append(info_dict['absolute_income'])
            trajectory_data['alpha'].append(info_dict['risky_share'])
            
            # 使用真实的效用值来累加终身回报
            current_discounted_reward +=  (beta ** time_step) * original_reward[0] # 
            time_step += 1

        episode_rewards.append(current_discounted_reward)
        episode_df = pd.DataFrame(trajectory_data)
        all_episode_data.append(episode_df)
        
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
                self.best_mean_reward = mean_reward

                fig, axes = plt.subplots(2, 1, figsize=(12, 10))
                retirement_age = self.eval_env.get_attr("tr")[0]
                
                axes[0].plot(rl_profiles_df['age'], rl_profiles_df['wealth'], 'b--', label='Wealth (RL)')
                axes[0].plot(rl_profiles_df['age'], rl_profiles_df['consumption'], 'r-', label='Consumption (RL)')
                axes[0].plot(rl_profiles_df['age'], rl_profiles_df['income'], 'k:', label='Income (Mean)')
                axes[0].axvline(x=retirement_age, color='grey', linestyle='--', label='Retirement')
                axes[0].set_title('Mean Lifecycle Profiles')
                axes[0].set_ylabel('Value (Absolute Units)')
                axes[0].legend(); axes[0].grid(True)
                
                axes[1].plot(rl_profiles_df['age'], rl_profiles_df['alpha'], 'm-', label='Risky Share (RL)')
                axes[1].axvline(x=retirement_age, color='grey', linestyle='--', label='Retirement')
                axes[1].set_title('Mean Portfolio Choice'); axes[1].set_xlabel('Age')
                axes[1].set_ylabel('Share'); axes[1].set_ylim(-0.05, 1.05)
                axes[1].legend(); axes[1].grid(True)
                
                fig.suptitle(
                    f"Evaluation at Timestep {self.num_timesteps}\n"
                    f"Current RL Utility: {mean_reward:.4f} (Best: {self.best_mean_reward:.4f})",
                    fontsize=14
                )
                fig.tight_layout(rect=[0, 0.03, 1, 0.94])
                
                self.logger.record("eval/lifecycle_profiles", Figure(fig, close=True), exclude=("stdout", "log", "json", "csv"))
                plt.close(fig)

        return True