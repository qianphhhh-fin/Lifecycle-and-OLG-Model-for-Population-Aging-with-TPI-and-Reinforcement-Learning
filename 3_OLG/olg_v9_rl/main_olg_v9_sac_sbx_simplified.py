# --- 开始文件：main_olg_v9_sac_sbx_simplified.py (最终版) ---

"""
OLG Model V9 SAC Simplified Training Script - SBX (Stable Baselines Jax)
在一个固定的宏观环境下训练SAC智能体

[最终版]
- 评估回调 (SimplifiedEvalCallback) 调用在 _utils 中定义的、封装好的
  HHSimulation_olgm_rl_simplified 函数，实现代码复用和逻辑统一。
- 最终评估也调用同一个函数，确保与训练过程中的评估完全一致。
"""

import numpy as np
import jax
import pickle
import time
import os
from typing import Dict, Any, Tuple
from sbx import SAC
from stable_baselines3.common.callbacks import BaseCallback
import matplotlib.pyplot as plt

# 导入主工具库和简化版环境
from main_olg_v9_utils import OLG_V9_Utils, OLGEnvV9SACSimp

# --- 开始：为简化版训练脚本添加的专用评估回调 ---
class SimplifiedEvalCallback(BaseCallback):
    """
    [重构版] 专门用于简化版环境的评估回调。
    它调用封装好的 HHSimulation_olgm_rl_simplified 函数进行评估。
    """
    def __init__(self, cS_obj: Any, paramS_rl: Dict, M_fixed: Dict,
                 n_eval_episodes: int = 100, eval_freq: int = 1000, 
                 log_path: str = None, best_model_save_path: str = None,
                 verbose: int = 1, use_survival_prob_in_eval: bool = True):
        super(SimplifiedEvalCallback, self).__init__(verbose)
        self.cS_obj = cS_obj
        self.paramS_rl = paramS_rl
        self.M_fixed = M_fixed
        self.n_eval_episodes = n_eval_episodes
        self.eval_freq = eval_freq
        self.log_path = log_path
        self.best_model_save_path = best_model_save_path
        self.use_survival_prob_in_eval = use_survival_prob_in_eval
        
        self.best_mean_reward = -np.inf
        
        if self.log_path is not None:
            os.makedirs(os.path.dirname(self.log_path), exist_ok=True)
            self.evaluations_results = []
            self.evaluations_timesteps = []
            
        print("✅ [重构版] 专用简化评估回调已初始化。")

    def _on_step(self) -> bool:
        if self.n_calls > 0 and self.n_calls % self.eval_freq == 0:
            # 1. 调用封装好的模拟函数获取消费路径
            sim_results = OLG_V9_Utils.HHSimulation_olgm_rl_simplified(
                self.model, self.cS_obj, self.paramS_rl, self.M_fixed,
                n_sim=self.n_eval_episodes, random_seed=42
            )
            c_paths = sim_results['c_path_rl']
            
            # 2. 计算每个路径的效用
            episode_rewards = [
                self._calculate_lifetime_utility(c_paths[i,:], self.cS_obj, self.use_survival_prob_in_eval)
                for i in range(self.n_eval_episodes)
            ]

            mean_reward = np.mean(episode_rewards)
            std_reward = np.std(episode_rewards)

            if self.verbose > 0:
                print(f"Eval @ {self.n_calls} ts: mean_reward:{mean_reward:.4f} +/- {std_reward:.4f}")
            
            if self.log_path is not None:
                self.evaluations_timesteps.append(self.num_timesteps)
                self.evaluations_results.append(episode_rewards)
                np.savez(self.log_path, timesteps=self.evaluations_timesteps, results=self.evaluations_results)
            
            if mean_reward > self.best_mean_reward:
                self.best_mean_reward = mean_reward
                if self.verbose > 0: print("🏆 New best mean reward!")
                if self.best_model_save_path is not None:
                    self.model.save(os.path.join(self.best_model_save_path, "best_model"))
            
            self.logger.record("eval/mean_reward", self.best_mean_reward)

        return True

    def _calculate_lifetime_utility(self, c_path: np.ndarray, cS: Any, use_survival_prob: bool) -> float:
        beta = cS.beta
        s_transitionV = cS.s_1yr_transitionV.flatten()
        aD = len(c_path)
        utility_sum = 0.0
        cumulative_discount = 1.0
        for a_group in range(aD):
            _, u = OLG_V9_Utils.CES_utility(c_path[a_group], cS.sigma, cS)
            utility_sum += cumulative_discount * u
            if a_group < aD - 1:
                survival_factor = s_transitionV[a_group] if use_survival_prob else 1.0
                cumulative_discount *= (beta * survival_factor)
        return utility_sum

# --- 结束：专用评估回调 ---

def plot_training_stats(log_path: str, save_path: str):
    """绘制训练统计图"""
    log_file = os.path.join(log_path, 'evaluations.npz')
    if os.path.exists(log_file):
        data = np.load(log_file)
        timesteps = data['timesteps']
        results = data['results']
        mean_rewards = np.mean(results, axis=1)
        
        plt.figure(figsize=(8, 5))
        plt.plot(timesteps, mean_rewards, 'b-')
        plt.xlabel('Timesteps')
        plt.ylabel('Mean Reward (Discounted Utility)')
        plt.title('Simplified SAC Training Progress')
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(save_path, dpi=300)
        plt.show()
        print(f"训练统计图已保存到 {save_path}")

def main_simplified():
    """主训练函数 (简化版)"""
    print("="*60)
    print("=== OLG 模型 V9 - SAC 智能体训练 (简化版 - 固定宏观环境) ===")
    print(f"    (JAX设备: {jax.devices()})")
    print("="*60)
    
    # 1. 初始化参数
    cS = OLG_V9_Utils.ParameterValues_HuggettStyle()
    paramS_for_rl = {}
    (paramS_for_rl['leLogGridV'], paramS_for_rl['leTrProbM'], paramS_for_rl['leProb1V']) = OLG_V9_Utils.EarningProcess_olgm(cS)
    paramS_for_rl['leGridV'] = np.exp(paramS_for_rl['leLogGridV'])
    
    # 2. 定义固定的宏观经济环境
    M_FIXED = {
        'R_k_net_factor': 1.03, 'w_gross': 2.0, 'TR_total': 0.1,
        'b_payg_avg_retiree': 0.4, 'tau_l': 0.15, 'theta_payg_actual': 0.12
    }
    print("\n--- 训练和评估将在以下固定环境下进行 ---")
    for key, value in M_FIXED.items(): print(f"  - {key}: {value:.4f}")
    
    # 3. 创建简化版强化学习环境
    env = OLGEnvV9SACSimp(
        cS=cS, paramS_rl=paramS_for_rl, M_fixed=M_FIXED,
        training_mode=True, reward_shaping_scheme='exponential', reward_alpha=0.1
    )

    # 4. 创建SBX SAC Agent
    model_kwargs = {
        'policy': 'MlpPolicy', 'env': env, 'learning_rate': 3e-4, 'buffer_size': int(1e6),
        'batch_size': 256, 'tau': 5e-3, 'gamma': 0.97, 'ent_coef': 'auto',
        'policy_kwargs': {'net_arch': [256, 256]}, 'verbose': 1,
        'learning_starts': 5000, 'seed': 42,
        'tensorboard_log': './py/tensorboard_logs_sbx_simp/'
    }
    model = SAC(**model_kwargs)
    
    # 5. 设置评估回调
    eval_callback = SimplifiedEvalCallback(
        cS_obj=cS, paramS_rl=paramS_for_rl, M_fixed=M_FIXED,
        n_eval_episodes=100, eval_freq=5000,
        log_path='./py/logs_sbx_simp/evaluations.npz',
        best_model_save_path='./py/best_model_sbx_simp/',
        verbose=1, use_survival_prob_in_eval=True
    )
    
    # 6. 开始训练
    total_timesteps = 200_000
    print(f"\n--- 开始训练 (总步数: {total_timesteps}) ---")
    
    os.makedirs('./py/best_model_sbx_simp/', exist_ok=True)
    os.makedirs('./py/logs_sbx_simp/', exist_ok=True)
    
    model.learn(total_timesteps=total_timesteps, callback=eval_callback, progress_bar=True)
    
    # 7. 保存最终模型和配置
    model_save_path = './py/final_sac_agent_olg_sbx_simp'
    model.save(model_save_path)
    config = {'cS': cS, 'paramS_for_rl': paramS_for_rl, 'M_FIXED': M_FIXED, 'model_kwargs': model_kwargs}
    with open(f'{model_save_path}_config.pkl', 'wb') as f: pickle.dump(config, f)
    print(f"\n✅ 最终模型和配置已保存至: {model_save_path}.zip / .pkl")
    
    # 8. 最终评估
    print("\n--- 最终评估 ---")
    best_model_path = os.path.join('./py/best_model_sbx_simp/', "best_model.zip")
    final_model = SAC.load(best_model_path, env=env) if os.path.exists(best_model_path) else model

    final_sim_results = OLG_V9_Utils.HHSimulation_olgm_rl_simplified(
        final_model, cS, paramS_for_rl, M_FIXED, n_sim=500, random_seed=42)
    final_c_paths = final_sim_results['c_path_rl']
    final_rewards = [eval_callback._calculate_lifetime_utility(final_c_paths[i,:], cS, True) for i in range(500)]
    
    print("\n" + "="*60)
    print("📈 简化版模型最终评估结果")
    print("="*60)
    print(f"  - 评估指标: VFI等价折现总效用 (含生存概率)")
    print(f"  - 平均回报: {np.mean(final_rewards):.4f} ± {np.std(final_rewards):.4f}")
    print(f"  - 训练过程最佳回报: {eval_callback.best_mean_reward:.4f}")
    
    # 9. 绘制训练曲线
    plot_training_stats('./py/logs_sbx_simp/', './py/training_stats_sbx_simp.png')

if __name__ == "__main__":
    main_simplified()