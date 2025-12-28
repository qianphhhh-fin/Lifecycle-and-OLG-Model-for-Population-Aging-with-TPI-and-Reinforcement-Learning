"""
OLG Model V8 SAC Training Script - SBX (Stable Baselines Jax) Implementation
基于SBX (SB3 + JAX) 的SAC算法训练OLG模型智能体

SBX的优势:
- 基于JAX，更高的计算性能
- 与SB3兼容的API，易于迁移
- 支持现代RL算法（SAC、TQC、PPO等）
"""

import numpy as np
import jax
import pickle
import time
import os
from typing import Dict, Any, Tuple
from sbx import SAC  # 使用SBX而不是SB3
from stable_baselines3.common.evaluation import evaluate_policy
from stable_baselines3.common.callbacks import EvalCallback, StopTrainingOnRewardThreshold
from stable_baselines3.common.monitor import Monitor
import matplotlib.pyplot as plt

from olg_utils import OLGUtils
from olg_env_v8_sac import OLGEnvV8SAC

def evaluate_policy_with_discount(model, env, n_eval_episodes=10, deterministic=True, 
                                 gamma=0.97, render=False, callback=None, 
                                 return_episode_rewards=False):
    """
    自定义评估函数，使用指定的折现因子计算episode总回报
    解决SBX/SB3的evaluate_policy默认无折现(γ=1.0)的问题
    
    Args:
        model: 要评估的模型
        env: 环境
        n_eval_episodes: 评估回合数
        deterministic: 是否使用确定性策略
        gamma: 折现因子（应与训练时的gamma保持一致）
        render: 是否渲染
        callback: 回调函数
        return_episode_rewards: 是否返回每个episode的奖励
        
    Returns:
        mean_reward: 平均折现回报
        std_reward: 标准差
        episode_rewards: (可选) 每个episode的折现回报列表
    """
    episode_rewards = []
    episode_lengths = []
    
    for i in range(n_eval_episodes):
        # 兼容不同的reset()返回值格式
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
            
            # 兼容不同的step()返回值格式
            step_result = env.step(action)
            if len(step_result) == 5:
                obs, reward, terminated, truncated, info = step_result
                done = terminated or truncated
            elif len(step_result) == 4:
                obs, reward, done, info = step_result
            else:
                raise ValueError(f"Unexpected step() return format: {len(step_result)} values")
            
            # 使用折现因子累积奖励
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
    
    🚨 重要更新：最后一期处理与MATLAB版本完全一致
    - 非最后一期：使用RL策略进行决策
    - 最后一期：强制消费所有资产，不储蓄，不缴费PPS（与main_olg_v8_utils.m一致）
    
    Args:
        model: 要评估的模型
        cS: 模型参数
        paramS_for_rl: RL专用参数
        rng_M: 宏观参数范围（用于环境初始化）
        n_sim: 模拟个体数量
        deterministic: 是否使用确定性策略
        gamma: 折现因子
        random_seed: 随机种子，确保结果可重现
        verbose: 是否显示详细信息
        
    Returns:
        mean_reward: 平均生命周期折现效用
        std_reward: 标准差
        lifecycle_results: 详细的生命周期模拟结果
    
    处理逻辑:
        1. 年龄0到aD_new-2：使用RL模型的策略决策
        2. 年龄aD_new-1（最后一期）：强制消费所有资产，不储蓄，不缴费PPS
           - k_next = kMin (不储蓄)
           - c_pps = 0 (不缴费PPS)
           - consumption = 所有可用资源 / (1 + tau_c)
           - utility = CRRA(consumption)
    """
    # 🎯 使用与compare_rl_and_vfi.py完全相同的固定测试参数
    M_fixed = {
        'R_k_net_factor': 1.03,
        'w_gross': 2.0,
        'TR_total': 0.1,
        'b_payg_avg_retiree': 0.4,
        'tau_l': 0.15,
        'theta_payg_actual': 0.12
    }
    
    if verbose:
        print(f"\n🧬 生命周期模拟评估 (n_sim={n_sim}, seed={random_seed})")
        print("🎯 使用与compare_rl_and_vfi.py完全相同的评估逻辑")
        print("📊 固定测试参数:")
        for key, value in M_fixed.items():
            print(f"  {key} = {value:.3f}")
    
    # 🎲 设置随机种子确保可重现性
    np.random.seed(random_seed)
    
    # 获取维度参数
    aD_new = int(cS['aD_new'])
    aR_new = int(cS['aR_new'])
    nk = int(cS['nk'])
    nkpps = int(cS['nkpps'])
    nw = int(cS['nw'])
    
    # 从参数中提取网格和转移概率
    kGridV = np.array(cS['kGridV']).flatten()
    kppsGridV = np.array(cS['kppsGridV']).flatten()
    leGridV = np.array(paramS_for_rl['leGridV']).flatten()
    leTrProbM = np.array(paramS_for_rl['leTrProbM'])
    leProb1V = np.array(paramS_for_rl['leProb1V']).flatten()
    ageEffV_new = np.array(cS['ageEffV_new']).flatten()
    beta = float(cS['beta'])
    sigma = float(cS['sigma'])
    kMin = float(cS['kMin'])
    kppsMin = float(cS['kppsMin'])
    kMax = float(cS['kMax'])
    kppsMax = float(cS['kppsMax'])
    pps_active = bool(cS['pps_active'])
    
    if verbose:
        print(f"📊 模型维度: aD_new={aD_new}, nk={nk}, nkpps={nkpps}, nw={nw}")
        print(f"📊 折现因子: β={beta:.3f} (应与γ={gamma:.3f}一致)")
    
    # 创建专用评估环境（与compare_rl_and_vfi.py中的设置完全相同）
    env = OLGEnvV8SAC(cS, paramS_for_rl, rng_M)
    env.set_macro_parameters(M_fixed)
    
    # 初始化结果存储
    lifetime_utility_rl = np.zeros(n_sim)
    
    # 生命周期轨迹存储
    k_path_rl = np.zeros((n_sim, aD_new))
    c_path_rl = np.zeros((n_sim, aD_new))
    cpps_path_rl = np.zeros((n_sim, aD_new))
    
    if verbose:
        print("🔄 开始生命周期模拟...")
    
    sim_start_time = time.time()
    
    for i_sim in range(n_sim):
        if verbose and (i_sim + 1) % 10 == 0:
            print(f"  进度: {i_sim + 1}/{n_sim}")
        
        # 初始状态
        k_current_rl = kMin
        kpps_current_rl = kppsMin
        
        # 初始效率冲击（与compare_rl_and_vfi.py完全相同的逻辑）
        eps_idx_current = np.where(np.random.rand() <= np.cumsum(leProb1V))[0]
        if len(eps_idx_current) > 0:
            eps_idx_current = eps_idx_current[0]
        else:
            eps_idx_current = 0
        
        utility_sum_rl = 0
        
        # RL环境重置
        obs, _ = env.reset()
        env.set_macro_parameters(M_fixed)
        
        for age_idx in range(aD_new):
            # 🚨 最后一期特殊处理：与MATLAB版本完全一致
            if age_idx == aD_new - 1:
                # 最后一期：强制消费所有资产，不储蓄，不缴费PPS
                if verbose and i_sim == 0:  # 只在第一个模拟时打印
                    print(f"    🏁 最后一期 (age_idx={age_idx}): 强制消费所有资产")
                
                # 计算可用资源（与MATLAB的HHIncome_Huggett逻辑一致）
                # 基础收入：资本收入 + 劳动收入 + 转移支付 + PAYG收益
                r_k_net = M_fixed['R_k_net_factor'] - 1  # 净利率
                capital_income = k_current_rl * r_k_net
                
                # 劳动收入（最后一期通常是退休期，劳动效率为0）
                age_efficiency = ageEffV_new[min(age_idx, len(ageEffV_new)-1)]
                epsilon_val = leGridV[eps_idx_current]
                labor_income = M_fixed['w_gross'] * age_efficiency * epsilon_val
                
                # 转移支付和PAYG收益
                transfer_income = M_fixed['TR_total']
                payg_income = M_fixed['b_payg_avg_retiree'] if age_idx >= aR_new else 0
                
                # 总的非PPS资源
                total_non_pps_resources = capital_income + labor_income + transfer_income + payg_income
                
                # PPS资源（如果激活且在提取期）
                pps_withdrawal = 0
                if pps_active and age_idx >= aR_new:  # 退休期可以提取PPS
                    pps_tax_rate_withdrawal = float(cS.get('pps_tax_rate_withdrawal', 0.03))
                    pps_withdrawal = kpps_current_rl * (1 - pps_tax_rate_withdrawal)
                
                # 总可用资源
                total_available_resources = total_non_pps_resources + pps_withdrawal
                
                # 强制设置决策变量（与MATLAB一致）
                k_next_rl = kMin  # 不储蓄
                cpps_rl = 0       # 不缴费PPS
                
                # 计算消费（扣除消费税）
                tau_c = float(cS.get('tau_c', 0.10))
                c_rl = max(float(cS.get('cFloor', 0.05)), total_available_resources / (1 + tau_c))
                
                # 计算效用（直接基于消费计算，与MATLAB的CES_utility一致）
                sigma_val = float(cS['sigma'])
                if abs(sigma_val - 1.0) < 1e-6:
                    u_rl = np.log(max(c_rl, 1e-10))  # log效用
                else:
                    u_rl = (max(c_rl, 1e-10)**(1 - sigma_val)) / (1 - sigma_val)  # CRRA效用
                
                if verbose and i_sim == 0:
                    print(f"      💰 总可用资源: {total_available_resources:.4f}")
                    print(f"      🍽️  最后一期消费: {c_rl:.4f}")
                    print(f"      📊 最后一期效用: {u_rl:.4f}")
                    print(f"      💾 储蓄: {k_next_rl:.4f} (强制为kMin)")
                    print(f"      🏦 PPS缴费: {cpps_rl:.4f} (强制为0)")
                
            else:
                # 非最后一期：使用RL策略
                action, _ = model.predict(obs, deterministic=deterministic)
                next_obs, reward, terminated, truncated, info = env.step(action)
                
                c_rl = info.get('consumption', 0)
                k_next_rl = info.get('k_prime', 0)
                cpps_rl = info.get('c_pps', 0)
                
                # RL的效用直接来自reward（与compare_rl_and_vfi.py一致）
                u_rl = reward
                
                obs = next_obs
            
            # 存储轨迹（适用于所有年龄）
            k_path_rl[i_sim, age_idx] = k_current_rl
            c_path_rl[i_sim, age_idx] = c_rl
            cpps_path_rl[i_sim, age_idx] = cpps_rl
            
            # 累积折现效用（与compare_rl_and_vfi.py完全相同的逻辑）
            discount_factor = beta ** age_idx
            utility_sum_rl += discount_factor * u_rl
            
            # 状态更新
            k_current_rl = k_next_rl
            
            # PPS资产演化（简化处理，与compare_rl_and_vfi.py一致）
            if pps_active:
                R_k_net_factor = M_fixed['R_k_net_factor']
                pps_return_premium = float(cS.get('pps_return_rate_premium', 0))
                pps_return_factor = 1 + ((R_k_net_factor - 1) + pps_return_premium)
                
                kpps_current_rl = (kpps_current_rl + cpps_rl) * pps_return_factor
                kpps_current_rl = max(kppsMin, min(kppsMax, kpps_current_rl))
            
            # 效率冲击演化（与compare_rl_and_vfi.py完全相同的逻辑）
            if age_idx < aD_new - 1:
                trans_probs = leTrProbM[eps_idx_current, :]
                eps_idx_current = np.where(np.random.rand() <= np.cumsum(trans_probs))[0]
                if len(eps_idx_current) > 0:
                    eps_idx_current = eps_idx_current[0]
                else:
                    eps_idx_current = min(eps_idx_current, nw-1)
            
            # 检查终止条件（仅对非最后一期）
            if age_idx < aD_new - 1 and (terminated or truncated):
                break
        
        lifetime_utility_rl[i_sim] = utility_sum_rl
    
    sim_time = time.time() - sim_start_time
    
    # 计算统计结果
    mean_utility_rl = np.mean(lifetime_utility_rl)
    std_utility_rl = np.std(lifetime_utility_rl)
    
    if verbose:
        print(f"✅ 生命周期模拟完成，耗时: {sim_time:.2f} 秒")
        print(f"📊 RL生命周期效用: {mean_utility_rl:.4f} ± {std_utility_rl:.4f}")
    
    # 构造详细结果
    lifecycle_results = {
        'success': True,
        'n_sim': n_sim,
        'random_seed': random_seed,
        'M_fixed': M_fixed,
        'simulation_time': sim_time,
        'mean_utility_rl': mean_utility_rl,
        'std_utility_rl': std_utility_rl,
        'lifetime_utility_rl': lifetime_utility_rl,
        'k_path_rl': k_path_rl,
        'c_path_rl': c_path_rl,
        'cpps_path_rl': cpps_path_rl,
        'aD_new': aD_new,
        'beta': beta,
        'gamma': gamma
    }
    
    return mean_utility_rl, std_utility_rl, lifecycle_results

class EvalCallbackWithDiscount(EvalCallback):
    """
    自定义EvalCallback，使用生命周期模拟进行评估
    🆕 使用evaluate_policy_lifecycle_simulation确保与最终评估完全一致
    """
    
    def __init__(self, eval_env, cS, paramS_for_rl, rng_M, gamma=0.97, 
                 use_lifecycle_simulation=True, **kwargs):
        """
        Args:
            eval_env: 评估环境
            cS: 模型参数
            paramS_for_rl: RL专用参数
            rng_M: 宏观参数范围
            gamma: 折现因子，应与训练时的gamma保持一致
            use_lifecycle_simulation: 是否使用生命周期模拟（默认True）
            **kwargs: 其他EvalCallback参数
        """
        super().__init__(eval_env, **kwargs)
        self.gamma = gamma
        self.cS = cS
        self.paramS_for_rl = paramS_for_rl
        self.rng_M = rng_M
        self.use_lifecycle_simulation = use_lifecycle_simulation
        
        if use_lifecycle_simulation:
            print(f"🔧 使用生命周期模拟评估回调 (γ={gamma})")
            print("🎯 确保与evaluate_policy_lifecycle_simulation完全一致")
        else:
            print(f"🔧 使用传统折现评估回调 (γ={gamma})")
        
        # 测试环境的返回格式
        try:
            reset_result = self.eval_env.reset()
            print(f"📋 环境reset()返回格式: {type(reset_result)} - {len(reset_result) if isinstance(reset_result, tuple) else '单值'}")
        except Exception as e:
            print(f"⚠️ 环境测试失败: {e}")
    
    def _on_step(self) -> bool:
        """
        重写_on_step方法，使用带折现的评估
        """
        continue_training = True

        if self.eval_freq > 0 and self.n_calls % self.eval_freq == 0:
            # 使用带折现的评估
            episode_rewards, episode_lengths = self._evaluate_with_discount()
            
            if self.log_path is not None:
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

            mean_reward, std_reward = np.mean(episode_rewards), np.std(episode_rewards)
            mean_ep_length, std_ep_length = np.mean(episode_lengths), np.std(episode_lengths)
            self.last_mean_reward = mean_reward

            if self.verbose >= 1:
                print(f"📊 Eval num_timesteps={self.num_timesteps}")
                print(f"   Current reward: {mean_reward:.2f} +/- {std_reward:.2f} (折现)")
                print(f"   Episode length: {mean_ep_length:.2f} +/- {std_ep_length:.2f}")
                print(f"🏆 Best reward: {self.best_mean_reward:.2f} (at timestep {getattr(self, 'best_timestep', 'N/A')})")
                
                # 显示改进情况
                if hasattr(self, 'best_timestep'):
                    improvement = mean_reward - self.best_mean_reward
                    if improvement > 0:
                        print(f"   📈 Improvement: +{improvement:.2f}")
                    elif improvement < 0:
                        print(f"   📉 Below best: {improvement:.2f}")
                    else:
                        print(f"   ➖ Same as best")

            # Add to current Logger
            self.logger.record("eval/mean_reward", float(mean_reward))
            self.logger.record("eval/mean_ep_length", mean_ep_length)

            if len(self._is_success_buffer) > 0:
                success_rate = np.mean(self._is_success_buffer)
                if self.verbose >= 1:
                    print(f"Success rate: {100 * success_rate:.2f}%")
                self.logger.record("eval/success_rate", success_rate)

            # Dump log so the evaluation results are printed with the correct timestep
            self.logger.record("time/total_timesteps", self.num_timesteps, exclude="tensorboard")
            self.logger.dump(self.num_timesteps)

            if mean_reward > self.best_mean_reward:
                if self.verbose >= 1:
                    print("🏆 New best mean reward!")
                if self.best_model_save_path is not None:
                    self.model.save(os.path.join(self.best_model_save_path, "best_model"))
                self.best_mean_reward = mean_reward
                self.best_timestep = self.num_timesteps  # 记录最佳结果的timestep

                # Trigger callback on new best model, if needed
                if self.callback_on_new_best is not None:
                    continue_training = self.callback_on_new_best.on_step()

        return continue_training
    
    def _evaluate_with_discount(self):
        """使用生命周期模拟或传统折现进行评估"""
        if self.use_lifecycle_simulation:
            return self._evaluate_with_lifecycle_simulation()
        else:
            return self._evaluate_with_traditional_discount()
    
    def _evaluate_with_lifecycle_simulation(self):
        """🆕 使用生命周期模拟进行评估"""
        # 生成基于timestep的随机种子，确保每次评估的可重现性
        eval_random_seed = 42 + (self.num_timesteps // 1000)  # 基于训练步数的种子
        
        # 调用统一的生命周期评估函数
        mean_reward, std_reward, lifecycle_results = evaluate_policy_lifecycle_simulation(
            self.model, self.cS, self.paramS_for_rl, self.rng_M,
            n_sim=max(self.n_eval_episodes,10), 
            deterministic=self.deterministic,
            gamma=self.gamma,
            random_seed=eval_random_seed,
            verbose=False  # 训练过程中不显示详细信息
        )
        
        # 生成与传统方法兼容的返回格式
        # 为了兼容原有的episode格式，我们创建n_eval_episodes个相同的结果
        episode_rewards = [mean_reward] * self.n_eval_episodes
        episode_lengths = [lifecycle_results['aD_new']] * self.n_eval_episodes  # 生命周期长度
        
        if self.verbose >= 1:
            print(f"   🧬 生命周期模拟 (n_sim={lifecycle_results['n_sim']}, seed={eval_random_seed})")
            print(f"   📊 模拟结果: {mean_reward:.4f} ± {std_reward:.4f}")
        
        return episode_rewards, episode_lengths
    
    def _evaluate_with_traditional_discount(self):
        """传统的折现评估方法（备用）"""
        episode_rewards = []
        episode_lengths = []
        
        for i in range(self.n_eval_episodes):
            # 兼容不同的reset()返回值格式
            reset_result = self.eval_env.reset()
            if isinstance(reset_result, tuple):
                obs, _ = reset_result
            else:
                obs = reset_result
                
            done = False
            episode_reward = 0.0
            episode_length = 0
            discount_factor = 1.0
            
            while not done:
                action, _ = self.model.predict(obs, deterministic=self.deterministic)
                
                # 兼容不同的step()返回值格式
                step_result = self.eval_env.step(action)
                if len(step_result) == 5:
                    obs, reward, terminated, truncated, info = step_result
                    done = terminated or truncated
                elif len(step_result) == 4:
                    obs, reward, done, info = step_result
                else:
                    raise ValueError(f"Unexpected step() return format: {len(step_result)} values")
                
                # 使用折现因子累积奖励
                episode_reward += discount_factor * reward
                discount_factor *= self.gamma
                
                episode_length += 1
                
                if self.render:
                    self.eval_env.render()
            
            episode_rewards.append(episode_reward)
            episode_lengths.append(episode_length)
        
        return episode_rewards, episode_lengths

def main():
    """主训练函数"""
    print("=== OLG 模型 V8 - SAC 智能体训练 (SBX - Stable Baselines Jax 版本) ===")
    print("    (在线RL，宏观变量M作为环境参数)")
    print("    (决策变量：PPS缴费比例, 非PPS储蓄比例)")
    print("    (网络架构：SBX JAX SAC实现)")
    print(f"    (JAX设备: {jax.devices()})")
    
    # 1. 初始化参数
    print("\n--- 1. 初始化参数 ---")
    cS = OLGUtils.parameter_values_huggett_style()
    
    # 计算RL相关参数
    paramS_for_rl = {}
    if ('leGridV' not in cS or 'leTrProbM' not in cS or 'leProb1V' not in cS):
        (paramS_for_rl['leLogGridV'], 
         paramS_for_rl['leTrProbM'], 
         paramS_for_rl['leProb1V']) = OLGUtils.earning_process_olgm(cS)
        paramS_for_rl['leGridV'] = np.exp(paramS_for_rl['leLogGridV'])
        cS['leGridV'] = paramS_for_rl['leGridV']
        cS['leTrProbM'] = paramS_for_rl['leTrProbM']
        cS['leProb1V'] = paramS_for_rl['leProb1V']
    else:
        paramS_for_rl['leGridV'] = cS['leGridV']
        paramS_for_rl['leTrProbM'] = cS['leTrProbM']
        paramS_for_rl['leProb1V'] = cS['leProb1V']
    
    paramS_for_rl['ageEffV_new'] = cS['ageEffV_new']
    
    # 2. 定义宏观状态M的采样范围
    print("\n--- 2. 定义宏观状态 M 的采样范围 ---")
    rng_M = {
        'R_k_net_factor': [1.01, 1.05],
        'w_gross': [1.5, 2.5],
        'TR_total': [0.0, 0.2],
        'b_payg_avg_retiree': [0.1, 0.8],
        'tau_l': [0.05, 0.25],
        'theta_payg_actual': [0.05, 0.20]
    }
    print("宏观参数采样范围已定义。")
    
    # 3. 创建强化学习环境
    print("\n--- 3. 创建强化学习环境 ---")
    env = OLGEnvV8SAC(cS, paramS_for_rl, rng_M)
    env = Monitor(env)  # 包装环境以记录统计信息
    
    print(f"观测空间: {env.observation_space}")
    print(f"动作空间: {env.action_space}")
    print("RL环境已创建。")
    
    # 4. 创建SBX SAC Agent
    print("\n--- 4. 创建 SBX SAC Agent ---")
    
    # SBX SAC超参数（JAX优化版本）
    model_kwargs = {
        'policy': 'MlpPolicy',
        'env': env,
        'learning_rate': 1e-4,           # SBX推荐的学习率
        'buffer_size': int(1e6),         # 经验回放缓冲区大小
        'batch_size': 512,               # SBX推荐的批量大小
        'tau': 5e-3,                     # 目标网络软更新系数
        'gamma': 0.97,                   # 折扣因子
        'ent_coef': 'auto',              # 自动调整熵系数
        'gradient_steps': 1,             # 每次更新的梯度步数
        'policy_kwargs': {
            'net_arch': [256, 256],      # 网络架构
            # SBX使用JAX，不需要指定activation_fn
        },
        'verbose': 1,
        'learning_starts': 10000,
        'seed': 42,
        'tensorboard_log': './py/tensorboard_logs_sbx/'  # 📊 TensorBoard日志路径
        # 注意：SBX不支持target_update_interval参数，使用默认值
    }
    
    print("创建SBX SAC模型...")
    model = SAC(**model_kwargs)
    print("SBX SAC Agent已创建。")
    print(f"使用JAX设备: {jax.devices()}")
    print(f"网络架构: {model_kwargs['policy_kwargs']}")
    
    # 5. 设置训练参数
    print("\n--- 5. 设置训练参数 ---")
    max_steps_per_episode = cS['aD_new']    # 每回合最大步数
    total_timesteps = 150_000               # 总训练步数（适中的数量用于测试）
    stop_training_value = -20               # 停止训练阈值
    eval_freq = 5_000                       # 评估频率
    n_eval_episodes = 100                    # 评估回合数
    
    print(f"每回合最大步数: {max_steps_per_episode}")
    print(f"总训练步数: {total_timesteps}")
    print(f"停止训练值: {stop_training_value}")
    
    # 6. 测试环境和智能体初始化
    print("\n--- 6. 测试环境和智能体初始化 ---")
    obs, _ = env.reset()
    print(f"环境重置成功，观察维度: {obs.shape}")
    
    # 测试智能体的初始动作
    action, _ = model.predict(obs, deterministic=False)
    print(f"智能体初始动作生成成功，动作维度: {action.shape}，动作值: [{action[0]:.4f}, {action[1]:.4f}]")
    
    # 7. 设置评估回调
    print("\n--- 7. 设置评估和回调 ---")
    
    # 创建评估环境（使用固定参数）
    eval_env_base = OLGEnvV8SAC(cS, paramS_for_rl, rng_M)
    M_eval = {
        'R_k_net_factor': 1.03,
        'w_gross': 2.0,
        'TR_total': 0.1,
        'b_payg_avg_retiree': 0.4,
        'tau_l': 0.15,
        'theta_payg_actual': 0.12
    }
    eval_env_base.set_macro_parameters(M_eval)
    eval_env = Monitor(eval_env_base)
    
    # 设置停止训练回调
    stop_callback = StopTrainingOnRewardThreshold(
        reward_threshold=stop_training_value, 
        verbose=1
    )
    
    # 设置评估回调
    eval_callback = EvalCallbackWithDiscount(
        eval_env,
        cS,                              # 🆕 传递模型参数
        paramS_for_rl,                   # 🆕 传递RL专用参数
        rng_M,                           # 🆕 传递宏观参数范围
        gamma=0.97,                      # 🔧 重要：使用与训练一致的折现因子
        use_lifecycle_simulation=True,   # 🆕 使用生命周期模拟
        best_model_save_path='./py/best_model_sbx/',
        log_path='./py/logs_sbx/',
        eval_freq=eval_freq,
        n_eval_episodes=n_eval_episodes,
        deterministic=True,
        render=False,
        callback_on_new_best=stop_callback,
        verbose=1
    )
    
    print("🧬 生命周期模拟评估回调已设置 (γ=0.97)。")
    print("🎯 训练过程评估将与最终评估使用完全相同的逻辑。")
    
    # 8. 开始训练
    print("\n--- 8. 开始训练 ---")
    print("使用SBX (Stable Baselines Jax) SAC算法训练...")
    
    # 创建保存目录
    os.makedirs('./py/best_model_sbx/', exist_ok=True)
    os.makedirs('./py/logs_sbx/', exist_ok=True)
    os.makedirs('./py/tensorboard_logs_sbx/', exist_ok=True)  # 📊 TensorBoard日志目录
    
    start_time = time.time()
    
    try:
        model.learn(
            total_timesteps=total_timesteps,
            callback=eval_callback,
            log_interval=100,
            progress_bar=True
        )
        print("训练完成。")
    except KeyboardInterrupt:
        print("训练被用户中断。")
    except Exception as e:
        print(f"训练过程中发生错误: {e}")
        import traceback
        traceback.print_exc()
    
    training_time = time.time() - start_time
    print(f"训练用时: {training_time:.2f} 秒")
    
    # 9. 保存最终模型
    print("\n--- 9. 保存最终模型 ---")
    model.save('./py/final_sac_agent_olg_sbx')

    # 保存参数和环境配置
    config = {
        'cS': cS,
        'paramS_for_rl': paramS_for_rl,
        'rng_M': rng_M,
        'M_eval': M_eval,
        'model_kwargs': model_kwargs,
        'training_time': training_time,
        'algorithm': 'SBX_SAC',
        'jax_devices': str(jax.devices())
    }

    # 导出Actor网络为ONNX格式
    print("\n--- 9.1 导出Actor网络为ONNX格式 ---")
    try:
        import torch as th
        import onnx
        import json
        
        class OnnxableSACPolicy(th.nn.Module):
            """用于ONNX导出的包装类"""
            def __init__(self, model):
                super().__init__()
                self.model = model
                
            def forward(self, observation: th.Tensor) -> th.Tensor:
                """只导出actor网络的前向传播"""
                # 使用确定性策略
                with th.no_grad():
                    # 这里我们只需要动作，不需要其他返回值
                    action = self.model.predict(observation.numpy(), deterministic=True)[0]
                    return th.as_tensor(action)
        
        # 提取模型参数
        model_params = {}
        
        # 获取观察空间和动作空间信息
        observation_space = model.observation_space
        action_space = model.action_space
        
        # 提取维度信息
        if hasattr(observation_space, 'shape'):
            model_params['input_dim'] = observation_space.shape[0]
        else:
            model_params['input_dim'] = 1
        
        if hasattr(action_space, 'shape'):
            model_params['output_dim'] = action_space.shape[0]
        else:
            model_params['output_dim'] = 1
        
        # 提取动作空间边界
        if hasattr(model.action_space, 'low') and hasattr(model.action_space, 'high'):
            model_params['action_space_low'] = model.action_space.low.tolist()
            model_params['action_space_high'] = model.action_space.high.tolist()
            print(f"提取动作空间边界: low={model_params['action_space_low']}, high={model_params['action_space_high']}")
        else:
            print("警告: 未能提取动作空间边界，MATLAB中可能需要手动设置。")
        
        # 尝试提取网络架构信息
        try:
            if hasattr(model, 'actor') and hasattr(model.actor, 'latent_pi_net'):
                # 提取actor网络架构
                actor_net = model.actor.latent_pi_net
                net_arch = []
                for module in actor_net.modules():
                    if isinstance(module, th.nn.Linear):
                        net_arch.append(module.out_features)
                model_params['net_arch'] = net_arch[:-1]  # 排除输出层
            else:
                # 默认网络架构
                model_params['net_arch'] = [256, 256]
        except Exception as e:
            print(f"提取网络架构失败: {e}")
            model_params['net_arch'] = [256, 256]  # 默认值
        
        # 创建ONNX导出包装器
        onnxable_model = OnnxableSACPolicy(model)
        
        # 创建示例输入
        observation_size = model.observation_space.shape
        dummy_input = th.zeros(1, *observation_size)
        
        # ONNX输出路径
        onnx_path = './py/final_sac_agent_olg_sbx.onnx'
        
        # 导出到ONNX
        print(f"导出模型到ONNX: {onnx_path}")
        th.onnx.export(
            onnxable_model,
            dummy_input,
            onnx_path,
            opset_version=14,
            input_names=["observation"],
            output_names=["action"],
            dynamic_axes={
                "observation": {0: "batch_size"},
                "action": {0: "batch_size"}
            }
        )
        
        # 验证ONNX模型
        onnx_model = onnx.load(onnx_path)
        onnx.checker.check_model(onnx_model)
        print("ONNX模型验证通过")
        
        # 从ONNX模型中提取额外信息
        for input_info in onnx_model.graph.input:
            if input_info.name == "observation":
                # 获取输入维度
                input_shape = []
                for dim in input_info.type.tensor_type.shape.dim:
                    if dim.dim_param:  # 动态维度
                        input_shape.append(-1)
                    else:
                        input_shape.append(dim.dim_value)
                input_shape = input_shape[1:]  # 去掉批次维度
                model_params['input_shape'] = input_shape
        
        for output_info in onnx_model.graph.output:
            if output_info.name == "action":
                # 获取输出维度
                output_shape = []
                for dim in output_info.type.tensor_type.shape.dim:
                    if dim.dim_param:  # 动态维度
                        output_shape.append(-1)
                    else:
                        output_shape.append(dim.dim_value)
                output_shape = output_shape[1:]  # 去掉批次维度
                model_params['output_shape'] = output_shape
        
        # 保存模型参数到JSON文件
        params_path = './py/final_sac_agent_olg_sbx_params.json'
        with open(params_path, 'w') as f:
            json.dump(model_params, f, indent=4)
        print(f"模型参数已保存到: {params_path}")
        
        # 同样为最佳模型导出ONNX和参数
        best_model_path = './py/best_model_sbx/best_model.onnx'
        best_model = SAC.load('./py/best_model_sbx/best_model')
        
        # 提取最佳模型参数（与上面相同的逻辑）
        best_model_params = {}
        best_model_params['input_dim'] = model_params['input_dim']  # 应该相同
        best_model_params['output_dim'] = model_params['output_dim']  # 应该相同
        best_model_params['net_arch'] = model_params['net_arch']  # 应该相同
        
        # 提取最佳模型的动作空间边界
        if hasattr(best_model.action_space, 'low') and hasattr(best_model.action_space, 'high'):
            best_model_params['action_space_low'] = best_model.action_space.low.tolist()
            best_model_params['action_space_high'] = best_model.action_space.high.tolist()
        
        onnxable_best_model = OnnxableSACPolicy(best_model)
        
        print(f"导出最佳模型到ONNX: {best_model_path}")
        th.onnx.export(
            onnxable_best_model,
            dummy_input,
            best_model_path,
            opset_version=14,
            input_names=["observation"],
            output_names=["action"],
            dynamic_axes={
                "observation": {0: "batch_size"},
                "action": {0: "batch_size"}
            }
        )
        
        # 验证最佳模型的ONNX
        best_onnx_model = onnx.load(best_model_path)
        onnx.checker.check_model(best_onnx_model)
        print("最佳模型ONNX验证通过")
        
        # 保存最佳模型参数
        best_params_path = './py/best_model_sbx/best_model_params.json'
        with open(best_params_path, 'w') as f:
            json.dump(best_model_params, f, indent=4)
        print(f"最佳模型参数已保存到: {best_params_path}")
        
        # 更新配置信息
        config['onnx_export'] = {
            'final_model': onnx_path,
            'final_model_params': params_path,
            'best_model': best_model_path,
            'best_model_params': best_params_path
        }
        
        # 打印模型参数摘要
        print("\n模型参数摘要:")
        print(f"输入维度: {model_params['input_dim']}")
        print(f"输出维度: {model_params['output_dim']}")
        print(f"网络架构: {model_params['net_arch']}")
        
    except ImportError as e:
        print(f"ONNX导出失败: {e}")
        print("请确保已安装PyTorch和ONNX包")
    except Exception as e:
        print(f"ONNX导出过程中出错: {e}")
        import traceback
        traceback.print_exc()

    # 保存参数和环境配置
    with open('./py/training_config_sbx.pkl', 'wb') as f:
        pickle.dump(config, f)
    
    print("最终模型和配置已保存。")
    
    # 10. 评估训练好的Agent
    print("\n--- 10. 评估训练好的 Agent ---")
    
    # 在固定参数下评估（参数已在前面设置）
    print("使用固定测试参数（与MATLAB版本一致）:")
    for key, value in M_eval.items():
        print(f"  {key} = {value:.3f}")
    
    # 对比不同折现率的评估结果
    print("\n🔧 折现率对比评估:")
    
    # 无折现评估（SB3默认）
    mean_reward_no_discount, std_reward_no_discount = evaluate_policy(
        model, eval_env, n_eval_episodes=100, deterministic=True
    )
    
    # 正确的折现评估（γ=0.97）
    mean_reward_discounted, std_reward_discounted = evaluate_policy_with_discount(
        model, eval_env, n_eval_episodes=100, deterministic=True, gamma=0.97
    )
    
    # 🆕 新增：生命周期模拟评估（与compare_rl_and_vfi.py完全相同）
    mean_reward_lifecycle, std_reward_lifecycle, lifecycle_results = evaluate_policy_lifecycle_simulation(
        model, cS, paramS_for_rl, rng_M, n_sim=100, deterministic=True, gamma=0.97, random_seed=42
    )
    
    print(f"❌ 无折现评估 (γ=1.0): {mean_reward_no_discount:.2f} ± {std_reward_no_discount:.2f}")
    print(f"✅ 正确折现评估，有随机宏观变量 (γ=0.97): {mean_reward_discounted:.2f} ± {std_reward_discounted:.2f}")
    print(f"🧬 生命周期模拟评估，无随机宏观变量 (γ=0.97): {mean_reward_lifecycle:.2f} ± {std_reward_lifecycle:.2f}")
    print(f"🏆 训练过程最佳结果: {eval_callback.best_mean_reward:.2f} (timestep {getattr(eval_callback, 'best_timestep', 'N/A')})")
    print(f"📊 最终 vs 最佳差异: {mean_reward_discounted - eval_callback.best_mean_reward:.2f}")
    print(f"💡 与VFI比较应使用生命周期模拟结果: {mean_reward_lifecycle:.2f}")
    
    print(f"\n📈 评估完成。训练过程中的评估已使用正确的折现因子 (γ=0.97)")
    
    # TensorBoard使用指南
    print(f"\n📊 TensorBoard 监控指南:")
    print(f"  启动命令: tensorboard --logdir=py/tensorboard_logs_sbx/ --port=6006")
    print(f"  访问地址: http://localhost:6006/")
    print(f"  重要指标:")
    print(f"    - train/actor_loss: Actor网络损失")
    print(f"    - train/critic_loss: Critic网络损失") 
    print(f"    - eval/mean_reward: 评估奖励（带折现）")
    print(f"    - eval/mean_ep_length: 平均回合长度")
    print(f"    - train/ent_coef: 熵系数变化")
    
    # 更新配置信息
    config['eval_gamma'] = 0.97
    config['mean_reward_no_discount'] = mean_reward_no_discount
    config['mean_reward_discounted'] = mean_reward_discounted
    config['mean_reward_lifecycle'] = mean_reward_lifecycle  # 🆕 新增生命周期评估结果
    config['std_reward_lifecycle'] = std_reward_lifecycle
    config['lifecycle_results'] = lifecycle_results  # 🆕 保存完整的生命周期结果
    config['best_training_reward'] = eval_callback.best_mean_reward
    config['best_training_timestep'] = getattr(eval_callback, 'best_timestep', None)
    config['tensorboard_log_path'] = './py/tensorboard_logs_sbx/'
    
    # 11. 绘制训练统计
    print("\n--- 11. 绘制训练统计 ---")
    try:
        plot_training_stats('./py/logs_sbx/', './py/training_stats_sbx.png')
    except Exception as e:
        print(f"绘图失败: {e}")
    
    print("SBX SAC Agent 训练和处理框架完成。")
    
    # 添加ONNX导出信息
    if 'onnx_export' in config:
        print("")
        print("🔄 ONNX模型导出信息:")
        print(f"  ✅ 最终模型ONNX: {config['onnx_export']['final_model']}")
        print(f"  ✅ 最终模型参数: {config['onnx_export']['final_model_params']}")
        print(f"  ✅ 最佳模型ONNX: {config['onnx_export']['best_model']}")
        print(f"  ✅ 最佳模型参数: {config['onnx_export']['best_model_params']}")
        print(f"  💡 ONNX模型和参数文件可直接在MATLAB中使用")
    
    print("="*60)
    
    return model, config, M_eval

def plot_training_stats(log_path: str, save_path: str):
    """绘制训练统计图"""
    
    # 读取训练日志
    log_file = os.path.join(log_path, 'evaluations.npz')
    if os.path.exists(log_file):
        data = np.load(log_file)
        timesteps = data['timesteps']
        results = data['results']
        
        plt.figure(figsize=(12, 4))
        
        # 评估奖励图
        plt.subplot(1, 2, 1)
        mean_rewards = np.mean(results, axis=1)
        std_rewards = np.std(results, axis=1)
        plt.plot(timesteps, mean_rewards, 'b-', label='Mean Reward', linewidth=2)
        plt.fill_between(timesteps, mean_rewards - std_rewards, 
                        mean_rewards + std_rewards, alpha=0.3)
        plt.xlabel('Timesteps')
        plt.ylabel('Reward')
        plt.title('SBX SAC Training Progress: Evaluation Rewards')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 奖励分布图
        plt.subplot(1, 2, 2)
        plt.boxplot(results.T, positions=timesteps[::max(1, len(timesteps)//10)])
        plt.xlabel('Timesteps')
        plt.ylabel('Reward')
        plt.title('Reward Distribution Over Training')
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"训练统计图已保存到 {save_path}")
    else:
        print(f"未找到训练日志文件: {log_file}")

if __name__ == "__main__":
    # 运行训练
    model, config, M_eval = main()
    
    # 输出最终结果摘要
    print("\n" + "="*60)
    print("🎯 SBX SAC 训练完成摘要")
    print("="*60)
    print(f"算法: {config['algorithm']}")
    print(f"JAX设备: {config['jax_devices']}")
    print(f"训练时间: {config['training_time']:.2f} 秒")
    print(f"网络架构: {config['model_kwargs']['policy_kwargs']['net_arch']}")
    print(f"训练折现因子: {config['model_kwargs']['gamma']}")
    print(f"评估折现因子: {config['eval_gamma']}")
    print(f"模型保存路径: ./py/final_sac_agent_olg_sbx.zip")
    print("")
    print("🏆 训练结果摘要:")
    print(f"  最佳训练结果: {config['best_training_reward']:.2f} (timestep {config['best_training_timestep']})")
    print(f"  最终评估结果: {config['mean_reward_discounted']:.2f}")
    print(f"  生命周期模拟结果: {config['mean_reward_lifecycle']:.2f} ± {config['std_reward_lifecycle']:.2f}")
    print(f"  性能退化: {config['mean_reward_discounted'] - config['best_training_reward']:.2f}")
    print("")
    print("🔧 折现率一致性确认:")
    print(f"  ✅ 训练时γ = {config['model_kwargs']['gamma']}")
    print(f"  ✅ 评估时γ = {config['eval_gamma']}")
    print(f"  ✅ 与VFI比较时应使用生命周期模拟结果: {config['mean_reward_lifecycle']:.2f}")
    print(f"  🎯 生命周期评估与compare_rl_and_vfi.py完全一致")
    
    # 添加ONNX导出信息
    if 'onnx_export' in config:
        print("")
        print("🔄 ONNX模型导出信息:")
        print(f"  ✅ 最终模型ONNX: {config['onnx_export']['final_model']}")
        print(f"  ✅ 最终模型参数: {config['onnx_export']['final_model_params']}")
        print(f"  ✅ 最佳模型ONNX: {config['onnx_export']['best_model']}")
        print(f"  ✅ 最佳模型参数: {config['onnx_export']['best_model_params']}")
        print(f"  💡 ONNX模型和参数文件可直接在MATLAB中使用")
    
    print("="*60) 