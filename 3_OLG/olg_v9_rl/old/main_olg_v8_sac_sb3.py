"""
OLG Model V8 SAC Training Script - SB3 (Stable Baselines 3) Implementation
基于SB3 (PyTorch) 的SAC算法训练OLG模型智能体

SB3的优势:
- 纯PyTorch实现，生态系统成熟
- 广泛使用，社区支持良好
- 支持现代RL算法（SAC、TQC、PPO等）
"""

import numpy as np
import pickle
import time
import os
import json
from typing import Dict, Any, Tuple

# 切换到SB3
import torch as th
import onnx
from stable_baselines3 import SAC
from stable_baselines3.common.evaluation import evaluate_policy
from stable_baselines3.common.callbacks import EvalCallback, StopTrainingOnRewardThreshold
from stable_baselines3.common.monitor import Monitor

from olg_utils import OLGUtils
from olg_env_v8_sac import OLGEnvV8SAC

class OnnxablePolicy(th.nn.Module):
    """
    用于ONNX导出的包装类 (遵循SB3官方文档)
    它返回actor网络的确定性（未缩放）输出。
    """
    def __init__(self, actor: th.nn.Module):
        super().__init__()
        self.actor = actor

    def forward(self, observation: th.Tensor) -> th.Tensor:
        # NOTE: The deterministic action is the mean of the distribution
        # This is the raw output from the network, typically in [-1, 1]
        return self.actor(observation, deterministic=True)


def evaluate_policy_with_discount(model, env, n_eval_episodes=10, deterministic=True, 
                                 gamma=0.97, render=False, callback=None, 
                                 return_episode_rewards=False):
    """
    自定义评估函数，使用指定的折现因子计算episode总回报
    解决SB3的evaluate_policy默认无折现(γ=1.0)的问题
    """
    episode_rewards = []
    episode_lengths = []
    
    for i in range(n_eval_episodes):
        reset_result = env.reset()
        if isinstance(reset_result, tuple):
            obs, _ = reset_result
        else:
            obs = reset_result
            
        done = False
        episode_reward = 0.0
        episode_length = 0
        discount_factor = 1.0
        
        while not done:
            action, _ = model.predict(obs, deterministic=deterministic)
            
            step_result = env.step(action)
            if len(step_result) == 5:
                obs, reward, terminated, truncated, info = step_result
                done = terminated or truncated
            else:
                obs, reward, done, info = step_result
            
            episode_reward += discount_factor * reward
            discount_factor *= gamma
            
            episode_length += 1
            
            if render:
                env.render()
            
            if callback is not None:
                callback(locals(), globals())
        
        episode_rewards.append(episode_reward)
        episode_lengths.append(episode_length)
    
    mean_reward = np.mean(episode_rewards)
    std_reward = np.std(episode_rewards)
    
    if return_episode_rewards:
        return mean_reward, std_reward, episode_rewards
    else:
        return mean_reward, std_reward

def evaluate_policy_lifecycle_simulation(model, cS, paramS_for_rl, rng_M, 
                                        n_sim=50, deterministic=True, 
                                        gamma=0.97, random_seed=42, verbose=True):
    """
    使用生命周期模拟评估RL模型（与compare_rl_and_vfi.py和MATLAB版本完全一致的评估过程）
    """
    M_fixed = {
        'R_k_net_factor': 1.03, 'w_gross': 2.0, 'TR_total': 0.1,
        'b_payg_avg_retiree': 0.4, 'tau_l': 0.15, 'theta_payg_actual': 0.12
    }
    
    if verbose:
        print(f"\n🧬 生命周期模拟评估 (n_sim={n_sim}, seed={random_seed})")
    
    np.random.seed(random_seed)
    
    aD_new, aR_new, nw = int(cS['aD_new']), int(cS['aR_new']), int(cS['nw'])
    leGridV = np.array(paramS_for_rl['leGridV']).flatten()
    leTrProbM = np.array(paramS_for_rl['leTrProbM'])
    leProb1V = np.array(paramS_for_rl['leProb1V']).flatten()
    ageEffV_new = np.array(cS['ageEffV_new']).flatten()
    beta, sigma = float(cS['beta']), float(cS['sigma'])
    kMin, kppsMin, kppsMax = float(cS['kMin']), float(cS['kppsMin']), float(cS['kppsMax'])
    pps_active = bool(cS['pps_active'])
    
    env = OLGEnvV8SAC(cS, paramS_for_rl, rng_M)
    env.set_macro_parameters(M_fixed)
    
    lifetime_utility_rl = np.zeros(n_sim)
    k_path_rl = np.zeros((n_sim, aD_new))
    c_path_rl = np.zeros((n_sim, aD_new))
    cpps_path_rl = np.zeros((n_sim, aD_new))
    
    for i_sim in range(n_sim):
        k_current_rl, kpps_current_rl = kMin, kppsMin
        eps_idx_current = np.where(np.random.rand() <= np.cumsum(leProb1V))[0][0]
        utility_sum_rl = 0
        obs, _ = env.reset()
        env.set_macro_parameters(M_fixed)
        
        for age_idx in range(aD_new):
            if age_idx == aD_new - 1:
                # 最后一期特殊处理... (逻辑与之前版本相同)
                capital_income = k_current_rl * (M_fixed['R_k_net_factor'] - 1)
                age_efficiency = ageEffV_new[min(age_idx, len(ageEffV_new)-1)]
                labor_income = M_fixed['w_gross'] * age_efficiency * leGridV[eps_idx_current]
                transfer_income = M_fixed['TR_total']
                payg_income = M_fixed['b_payg_avg_retiree'] if age_idx >= aR_new else 0
                total_non_pps_resources = capital_income + labor_income + transfer_income + payg_income
                pps_withdrawal = 0
                if pps_active and age_idx >= aR_new:
                    pps_tax_rate_withdrawal = float(cS.get('pps_tax_rate_withdrawal', 0.03))
                    pps_withdrawal = kpps_current_rl * (1 - pps_tax_rate_withdrawal)
                total_available_resources = total_non_pps_resources + pps_withdrawal
                k_next_rl, cpps_rl = kMin, 0
                c_rl = max(float(cS.get('cFloor', 0.05)), total_available_resources / (1 + float(cS.get('tau_c', 0.10))))
                u_rl = (max(c_rl, 1e-10)**(1 - sigma)) / (1 - sigma) if abs(sigma - 1.0) > 1e-6 else np.log(max(c_rl, 1e-10))
            else:
                action, _ = model.predict(obs, deterministic=deterministic)
                next_obs, reward, terminated, truncated, info = env.step(action)
                c_rl, k_next_rl, cpps_rl = info.get('consumption', 0), info.get('k_prime', 0), info.get('c_pps', 0)
                u_rl, obs = reward, next_obs
            
            k_path_rl[i_sim, age_idx], c_path_rl[i_sim, age_idx], cpps_path_rl[i_sim, age_idx] = k_current_rl, c_rl, cpps_rl
            utility_sum_rl += (beta ** age_idx) * u_rl
            k_current_rl = k_next_rl
            
            if pps_active:
                pps_return_premium = float(cS.get('pps_return_rate_premium', 0))
                pps_return_factor = 1 + ((M_fixed['R_k_net_factor'] - 1) + pps_return_premium)
                kpps_current_rl = max(kppsMin, min(kppsMax, (kpps_current_rl + cpps_rl) * pps_return_factor))
            
            if age_idx < aD_new - 1:
                eps_idx_current = np.where(np.random.rand() <= np.cumsum(leTrProbM[eps_idx_current, :]))[0][0]
                if 'terminated' in locals() and (terminated or truncated): break
        
        lifetime_utility_rl[i_sim] = utility_sum_rl
    
    mean_utility_rl, std_utility_rl = np.mean(lifetime_utility_rl), np.std(lifetime_utility_rl)
    
    lifecycle_results = {
        'mean_utility_rl': mean_utility_rl, 'std_utility_rl': std_utility_rl,
        'lifetime_utility_rl': lifetime_utility_rl, 'k_path_rl': k_path_rl,
        'c_path_rl': c_path_rl, 'cpps_path_rl': cpps_path_rl, 'aD_new': aD_new
    }
    
    return mean_utility_rl, std_utility_rl, lifecycle_results

class EvalCallbackWithDiscount(EvalCallback):
    """
    自定义EvalCallback，使用生命周期模拟进行评估，并在找到最佳模型时导出为ONNX。
    """
    def __init__(self, eval_env, cS, paramS_for_rl, rng_M, gamma=0.97, **kwargs):
        super().__init__(eval_env, **kwargs)
        self.gamma = gamma
        self.cS = cS
        self.paramS_for_rl = paramS_for_rl
        self.rng_M = rng_M
        self.best_timestep = 0
        print(f"🔧 使用生命周期模拟评估回调 (γ={gamma})")

    def _export_best_model_to_onnx(self):
        """每当找到新的最佳模型时，将其导出为ONNX并保存参数。"""
        if self.best_model_save_path is None: return
        
        onnx_save_dir = self.best_model_save_path
        onnx_path = os.path.join(onnx_save_dir, "best_model.onnx")
        params_path = os.path.join(onnx_save_dir, "best_model_params.json")
        model_zip_path = os.path.join(onnx_save_dir, "best_model.zip")

        if self.verbose > 0:
            print(f"📦 检测到新最佳模型，导出到ONNX: {onnx_path}")

        try:
            model = SAC.load(model_zip_path, device='cpu')
            model_params = {}
            obs_space, act_space = model.observation_space, model.action_space
            model_params['input_dim'] = obs_space.shape[0]
            model_params['output_dim'] = act_space.shape[0]
            model_params['action_space_low'] = act_space.low.tolist()
            model_params['action_space_high'] = act_space.high.tolist()
            model_params['net_arch'] = model.policy.net_arch
            
            with open(params_path, 'w') as f:
                json.dump(model_params, f, indent=4)
            if self.verbose > 0: print(f"   - 参数已保存: {params_path}")

            onnxable_model = OnnxablePolicy(model.policy.actor)
            dummy_input = th.randn(1, *obs_space.shape, device='cpu')

            th.onnx.export(
                onnxable_model, dummy_input, onnx_path, opset_version=14,
                input_names=["observation"], output_names=["action"],
                dynamic_axes={"observation": {0: "batch_size"}, "action": {0: "batch_size"}}
            )
            onnx.checker.check_model(onnx.load(onnx_path))
            if self.verbose > 0: print(f"   - ✅ ONNX模型验证并保存成功。")
        except Exception as e:
            if self.verbose > 0: print(f"   - ❌ ONNX导出失败: {e}")

    def _on_step(self) -> bool:
        if self.eval_freq > 0 and self.n_calls % self.eval_freq == 0:
            eval_random_seed = 42 + (self.num_timesteps // 1000)
            mean_reward, std_reward, lifecycle_results = evaluate_policy_lifecycle_simulation(
                self.model, self.cS, self.paramS_for_rl, self.rng_M,
                n_sim=max(self.n_eval_episodes, 10), deterministic=self.deterministic,
                gamma=self.gamma, random_seed=eval_random_seed, verbose=False
            )
            
            if self.verbose >= 1:
                print(f"📊 Eval num_timesteps={self.num_timesteps}, Reward: {mean_reward:.4f} +/- {std_reward:.4f}")

            if mean_reward > self.best_mean_reward:
                if self.verbose >= 1: print("🏆 New best mean reward!")
                self.best_mean_reward = mean_reward
                self.best_timestep = self.num_timesteps
                if self.best_model_save_path is not None:
                    self.model.save(os.path.join(self.best_model_save_path, "best_model"))
                    self._export_best_model_to_onnx()
        return True

def main():
    """主训练函数"""
    print("=== OLG 模型 V8 - SAC 智能体训练 (SB3 - Stable Baselines 3 版本) ===")
    
    cS = OLGUtils.parameter_values_huggett_style()
    paramS_for_rl = {}
    (paramS_for_rl['leLogGridV'], paramS_for_rl['leTrProbM'], paramS_for_rl['leProb1V']) = OLGUtils.earning_process_olgm(cS)
    paramS_for_rl['leGridV'] = np.exp(paramS_for_rl['leLogGridV'])
    paramS_for_rl['ageEffV_new'] = cS['ageEffV_new']
    
    rng_M = {
        'R_k_net_factor': [1.01, 1.05], 'w_gross': [1.5, 2.5], 'TR_total': [0.0, 0.2],
        'b_payg_avg_retiree': [0.1, 0.8], 'tau_l': [0.05, 0.25], 'theta_payg_actual': [0.05, 0.20]
    }
    
    env = Monitor(OLGEnvV8SAC(cS, paramS_for_rl, rng_M))
    
    model_kwargs = {
        'policy': 'MlpPolicy', 'env': env, 'learning_rate': 1e-4, 'buffer_size': int(1e6),
        'batch_size': 512, 'tau': 5e-3, 'gamma': 0.97, 'ent_coef': 'auto',
        'gradient_steps': 1,
        'policy_kwargs': {'net_arch': [256, 256], 'activation_fn': th.nn.ReLU},
        'verbose': 1, 'learning_starts': 10000, 'seed': 42,
        'tensorboard_log': './py/tensorboard_logs_sb3/'
    }
    model = SAC(**model_kwargs)
    print("SB3 SAC Agent已创建。")
    
    total_timesteps = 150_000
    
    eval_env = Monitor(OLGEnvV8SAC(cS, paramS_for_rl, rng_M))
    
    eval_callback = EvalCallbackWithDiscount(
        eval_env, cS, paramS_for_rl, rng_M, gamma=0.97,
        best_model_save_path='./py/best_model_sb3/',
        log_path='./py/logs_sb3/', eval_freq=5_000, n_eval_episodes=100,
        deterministic=True
    )
    
    os.makedirs('./py/best_model_sb3/', exist_ok=True)
    os.makedirs('./py/logs_sb3/', exist_ok=True)
    os.makedirs('./py/tensorboard_logs_sb3/', exist_ok=True)
    
    start_time = time.time()
    model.learn(total_timesteps=total_timesteps, callback=eval_callback, progress_bar=True)
    training_time = time.time() - start_time
    print(f"训练用时: {training_time:.2f} 秒")
    
    model.save('./py/final_sac_agent_olg_sb3')
    config = {'cS': cS, 'paramS_for_rl': paramS_for_rl, 'rng_M': rng_M, 'model_kwargs': model_kwargs}
    with open('./py/training_config_sb3.pkl', 'wb') as f:
        pickle.dump(config, f)
    
    print("\n评估最终模型:")
    mean_reward_lifecycle, _, _ = evaluate_policy_lifecycle_simulation(
        model, cS, paramS_for_rl, rng_M, n_sim=100, deterministic=True, random_seed=42
    )
    print(f"生命周期模拟评估 (γ=0.97): {mean_reward_lifecycle:.4f}")
    print(f"训练过程最佳结果: {eval_callback.best_mean_reward:.4f} (timestep {eval_callback.best_timestep})")
    
if __name__ == "__main__":
    main() 