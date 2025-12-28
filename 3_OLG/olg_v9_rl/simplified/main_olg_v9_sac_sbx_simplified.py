"""
OLG Model V9 SAC Training Script - Simplified Version with Fixed Macro Parameters
简化版SAC训练脚本 - 使用固定宏观参数

🎯 核心简化：
- 宏观参数在环境初始化时固定
- RL状态变量与VFI基本相同：(k, k_pps, age, ε) - 4维
- 保持累积存活概率方法的理论等价性
- 更公平地比较RL和VFI性能

🎮 动作空间设计：
- 2维连续动作：[PPS缴费比例, 消费比例]
- 决策顺序：先PPS缴费 → 再消费 → 最后储蓄（自动）
- PPS缴费比例：[0, pps_max_contrib_frac]
- 消费比例：[0, 1.0]，作用于可用于消费的资源

💡 决策逻辑：
1. 根据PPS缴费比例确定PPS缴费
2. 根据消费比例确定消费支出（基础消费 + 比例 × 超额资源）
3. 剩余资源自动用于储蓄

主要变化：
1. 环境使用固定宏观参数，不再作为状态变量
2. 观测空间从10维降到4维
3. 动作空间语义改为[PPS缴费比例, 消费比例]
4. 其他训练逻辑保持不变
5. 保持累积存活概率方法
"""

import numpy as np
import jax
import pickle
import time
import os
from typing import Dict, Any, Tuple
from sbx import SAC  # 使用标准SBX SAC
from stable_baselines3.common.evaluation import evaluate_policy
from stable_baselines3.common.callbacks import EvalCallback, StopTrainingOnRewardThreshold
from stable_baselines3.common.monitor import Monitor
import matplotlib.pyplot as plt

# 导入简化版工具类
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from main_olg_v9_utils import OLG_V9_Utils
from simplified.main_olg_v9_utils_simplified import OLGEnvV9SACSimplified, OLGUtilsSimplified

def evaluate_policy_lifecycle_simulation_simplified(model, cS, paramS_for_rl, M_fixed, 
                                                   n_sim=50, deterministic=True, 
                                                   gamma=0.97, random_seed=42, verbose=True,
                                                   eIdxM_annual_input=None):
    """
    简化版生命周期模拟评估
    
    Args:
        model: RL模型
        cS: 模型参数
        paramS_for_rl: RL专用参数
        M_fixed: 固定宏观参数
        n_sim: 模拟个体数量
        deterministic: 是否使用确定性策略
        gamma: 折现因子（应该等于cS.beta）
        random_seed: 随机种子
        verbose: 是否显示详细信息
        eIdxM_annual_input: 可选的预生成效率冲击序列
        
    Returns:
        mean_reward: 平均生命周期折现效用
        std_reward: 标准差
        lifecycle_results: 详细的生命周期模拟结果
    """
    
    sim_start_time = time.time()

    if verbose:
        print(f"\n🧬 简化版生命周期模拟评估")
        print(f"   (n_sim={n_sim}, seed={random_seed})")
        print("📊 固定宏观参数:")
        for key, value in M_fixed.items():
            print(f"     {key} = {value:.3f}")

    # 设置随机种子并生成效率冲击路径
    np.random.seed(random_seed)
    aD_orig = int(cS.aD_orig)
    nw = int(cS.nw)
    leProb1V = np.array(paramS_for_rl['leProb1V']).flatten()
    leTrProbM = np.array(paramS_for_rl['leTrProbM'])
    
    if eIdxM_annual_input is None:
        if verbose: print("🔄 生成年度效率冲击路径...")
        random_numbers = np.random.rand(n_sim, aD_orig)
        eIdxM_annual = OLG_V9_Utils.MarkovChainSimulation(
            n_sim, aD_orig, leProb1V, leTrProbM, random_numbers
        )
    else:
        if verbose: print("🔄 使用预生成的年度效率冲击路径...")
        eIdxM_annual = eIdxM_annual_input
        if eIdxM_annual.shape[0] != n_sim:
            print(f"⚠️ 警告：预生成的效率冲击个体数 ({eIdxM_annual.shape[0]}) 与 n_sim ({n_sim}) 不符。将使用预生成路径的个体数。")
            n_sim = eIdxM_annual.shape[0]

    # 调用简化版核心模拟函数
    if verbose: print("🚀 调用简化版生命周期模拟...")
    k_path_rl_group, kpps_path_rl_group, c_path_rl_group, cpps_path_rl_group = OLGUtilsSimplified.HHSimulation_olgm_rl_simplified(
        model, cS, paramS_for_rl, M_fixed, eIdxM_annual
    )
    if verbose: print("✅ 核心模拟完成。")

    aD_new = int(cS.aD_new)
    
    # 初始化年度路径矩阵
    k_path_rl = np.zeros((n_sim, aD_orig))
    kpps_path_rl = np.zeros((n_sim, aD_orig))
    c_path_rl = np.zeros((n_sim, aD_orig))
    cpps_path_rl = np.zeros((n_sim, aD_orig))
    

    # 根据返回的消费路径计算生命周期效用
    if verbose: print("📊 计算生命周期效用...")
    lifetime_utility_rl = np.zeros(n_sim)
    beta = float(cS.beta)
    for i in range(n_sim):
        utility_sum = 0
        for a in range(aD_orig):
            c = c_path_rl[i, a]
            _, u = OLG_V9_Utils.CES_utility(c, cS.sigma, cS)
            utility_sum += (beta ** a) * u
        lifetime_utility_rl[i] = utility_sum

    mean_utility_rl = np.mean(lifetime_utility_rl)
    std_utility_rl = np.std(lifetime_utility_rl)
    
    sim_time = time.time() - sim_start_time

    if verbose:
        print(f"✅ 简化版生命周期评估完成，耗时: {sim_time:.2f} 秒")
        print(f"📊 RL生命周期效用: {mean_utility_rl:.4f} ± {std_utility_rl:.4f}")
    
    # 整理并返回结果
    ageToGroupMap = np.zeros(cS.aD_orig, dtype=int)
    for i_group, indices in enumerate(cS.physAgeMap):
        if indices:
            ageToGroupMap[np.array(indices)] = i_group

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
        'kpps_path_rl': kpps_path_rl,
        'c_path_rl': c_path_rl,
        'cpps_path_rl': cpps_path_rl,
        'eIdxM_annual': eIdxM_annual,
        'aD_orig': aD_orig,
        'aD_new': int(cS.aD_new),
        'ageToGroupMap': ageToGroupMap,
        'beta': beta,
        'gamma': gamma,
        'time_scale': 'annual_matlab_consistent',
        'matlab_compatible': True,
        'vfi_equivalent': True,
        'simplified_version': True,
        'state_space_dim': 4  # 简化版特征
    }
    
    return mean_utility_rl, std_utility_rl, lifecycle_results

class EvalCallbackSimplified(EvalCallback):
    """
    简化版评估回调
    """
    
    def __init__(self, eval_env, cS, paramS_for_rl, M_fixed, gamma=0.97, 
                 use_lifecycle_simulation=True, **kwargs):
        """
        Args:
            eval_env: 评估环境
            cS: 模型参数
            paramS_for_rl: RL专用参数
            M_fixed: 固定宏观参数
            gamma: 折现因子
            use_lifecycle_simulation: 是否使用生命周期模拟
            **kwargs: 其他EvalCallback参数
        """
        super().__init__(eval_env, **kwargs)
        self.gamma = gamma
        self.cS = cS
        self.paramS_for_rl = paramS_for_rl
        self.M_fixed = M_fixed
        self.use_lifecycle_simulation = use_lifecycle_simulation
        
        # 生成统一的效率冲击路径
        if use_lifecycle_simulation:
            self._generate_unified_efficiency_shocks()
        
        print(f"🔧 简化版评估回调 (状态空间4维, γ={gamma})")
        print("🎯 基于累积存活概率方法，与VFI理论完全等价")

    def _generate_unified_efficiency_shocks(self):
        """生成统一的效率冲击路径"""
        print("🎲 生成统一效率冲击序列（简化版）...")
        
        np.random.seed(42)
        n_sim_target = max(self.n_eval_episodes, 100)
        aD_orig = int(self.cS.aD_orig)
        
        leProb1V = np.array(self.paramS_for_rl['leProb1V']).flatten()
        leTrProbM = np.array(self.paramS_for_rl['leTrProbM'])
        
        random_numbers = np.random.rand(n_sim_target, aD_orig)
        # self.eIdxM_unified = OLG_V9_Utils.MarkovChainSimulation_(
        #     n_sim_target, aD_orig, leProb1V, leTrProbM, random_numbers
        # )

        self.eIdxM_unified = OLG_V9_Utils.LaborEndowSimulation_olgm_AgeGroup(self.cS, self.paramS_for_rl)
        
        print(f"✅ 统一效率冲击矩阵生成完成: {self.eIdxM_unified.shape}")

    def _on_step(self) -> bool:
        """重写_on_step方法"""
        continue_training = True

        if self.eval_freq > 0 and self.n_calls % self.eval_freq == 0:
            episode_rewards, episode_lengths = self._evaluate_with_discount()
            
            if self.log_path is not None:
                self.evaluations_timesteps.append(self.num_timesteps)
                self.evaluations_results.append(episode_rewards)
                self.evaluations_length.append(episode_lengths)

                kwargs = {}
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
                print(f"   Current reward: {mean_reward:.2f} +/- {std_reward:.2f} (简化版)")
                print(f"   Episode length: {mean_ep_length:.2f} +/- {std_ep_length:.2f}")
                print(f"🏆 Best reward: {self.best_mean_reward:.2f} (at timestep {getattr(self, 'best_timestep', 'N/A')})")

            self.logger.record("eval/mean_reward", float(mean_reward))
            self.logger.record("eval/mean_ep_length", mean_ep_length)

            if len(self._is_success_buffer) > 0:
                success_rate = np.mean(self._is_success_buffer)
                if self.verbose >= 1:
                    print(f"Success rate: {100 * success_rate:.2f}%")
                self.logger.record("eval/success_rate", success_rate)

            self.logger.record("time/total_timesteps", self.num_timesteps, exclude="tensorboard")
            self.logger.dump(self.num_timesteps)

            if mean_reward > self.best_mean_reward:
                if self.verbose >= 1:
                    print("🏆 New best mean reward!")
                if self.best_model_save_path is not None:
                    self.model.save(os.path.join(self.best_model_save_path, "best_model"))
                self.best_mean_reward = mean_reward
                self.best_timestep = self.num_timesteps

                if self.callback_on_new_best is not None:
                    continue_training = self.callback_on_new_best.on_step()

        return continue_training
    
    def _evaluate_with_discount(self):
        """使用生命周期模拟或标准折现进行评估"""
        if self.use_lifecycle_simulation:
            return self._evaluate_with_lifecycle_simulation()
        else:
            return self._evaluate_with_traditional_discount()
    
    def _evaluate_with_lifecycle_simulation(self):
        """使用生命周期模拟进行评估"""
        eval_random_seed = 42
        n_sim_actual = max(self.n_eval_episodes, 100)
        eIdxM_input = self.eIdxM_unified[:n_sim_actual, :] if hasattr(self, 'eIdxM_unified') else None
        
        mean_reward, std_reward, lifecycle_results = evaluate_policy_lifecycle_simulation_simplified(
            self.model, self.cS, self.paramS_for_rl, self.M_fixed,
            n_sim=n_sim_actual, 
            deterministic=self.deterministic,
            gamma=self.gamma,
            random_seed=eval_random_seed,
            verbose=False,
            eIdxM_annual_input=eIdxM_input
        )
        
        episode_rewards = [mean_reward] * self.n_eval_episodes
        episode_lengths = [lifecycle_results['aD_orig']] * self.n_eval_episodes
        
        if self.verbose >= 1:
            print(f"   🧬 简化版生命周期模拟 (n_sim={lifecycle_results['n_sim']}, seed={eval_random_seed})")
            print(f"   📊 模拟结果: {mean_reward:.4f} ± {std_reward:.4f}")
        
        return episode_rewards, episode_lengths
    
    def _evaluate_with_traditional_discount(self):
        """传统的折现评估方法（备用）"""
        episode_rewards = []
        episode_lengths = []
        
        for i in range(self.n_eval_episodes):
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
                
                step_result = self.eval_env.step(action)
                if len(step_result) == 5:
                    obs, reward, terminated, truncated, info = step_result
                    done = terminated or truncated
                elif len(step_result) == 4:
                    obs, reward, done, info = step_result
                else:
                    raise ValueError(f"Unexpected step() return format: {len(step_result)} values")
                
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
    print("=== OLG 模型 V9 - SAC 智能体训练 (简化版 - 固定宏观参数) ===")
    print("    🎯 状态空间：(k, k_pps, age, ε) - 4维，与VFI基本相同")
    print("    🔧 宏观参数固定，不作为状态变量")
    print("    💡 累积存活概率方法保持不变")
    print(f"    (JAX设备: {jax.devices()})")
    
    # 1. 初始化参数
    print("\n--- 1. 初始化参数 ---")
    cS = OLG_V9_Utils.ParameterValues_HuggettStyle()
    
    # 计算RL相关参数
    paramS_for_rl = {}
    if (not hasattr(cS, 'leGridV') or not hasattr(cS, 'leTrProbM') or not hasattr(cS, 'leProb1V')):
        (paramS_for_rl['leLogGridV'], 
         paramS_for_rl['leTrProbM'], 
         paramS_for_rl['leProb1V']) = OLG_V9_Utils.EarningProcess_olgm(cS)
        paramS_for_rl['leGridV'] = np.exp(paramS_for_rl['leLogGridV'])
        cS.leGridV = paramS_for_rl['leGridV']
        cS.leTrProbM = paramS_for_rl['leTrProbM']
        cS.leProb1V = paramS_for_rl['leProb1V']
    else:
        paramS_for_rl['leGridV'] = cS.leGridV
        paramS_for_rl['leTrProbM'] = cS.leTrProbM
        paramS_for_rl['leProb1V'] = cS.leProb1V
    
    paramS_for_rl['ageEffV_new'] = cS.ageEffV_new
    
    # 2. 定义固定宏观参数
    print("\n--- 2. 定义固定宏观参数 ---")
    M_fixed = {
        'R_k_net_factor': 1.03,
        'w_gross': 2.0,
        'TR_total': 0.1,
        'b_payg_avg_retiree': 0.4,
        'tau_l': 0.15,
        'theta_payg_actual': 0.12
    }
    
    print("🎯 固定宏观参数（与VFI测试参数一致）:")
    for key, value in M_fixed.items():
        print(f"  {key} = {value:.3f}")
    
    # 3. 创建强化学习环境
    print("\n--- 3. 创建简化版强化学习环境 ---")
    # 训练环境：使用训练模式
    env = OLGEnvV9SACSimplified(cS, paramS_for_rl, M_fixed, training_mode=True)
    env = Monitor(env)
    
    print(f"观测空间: {env.observation_space} (4维)")
    print(f"动作空间: {env.action_space} ([PPS缴费比例, 消费比例])")
    print("🎯 简化版训练环境已创建（固定宏观参数）")
    print("💡 动作空间：先PPS缴费，再消费，最后自动储蓄")
    
    # 4. 创建SBX SAC Agent
    print("\n--- 4. 创建 SBX SAC Agent ---")
    
    model_kwargs = {
        'policy': 'MlpPolicy',
        'env': env,
        'learning_rate': 3e-5,
        'buffer_size': int(1e6),
        'batch_size': 256,
        'tau': 5e-3,
        'gamma': 0.97,
        'ent_coef': 'auto',
        'gradient_steps': 1,
        'policy_kwargs': {
            'net_arch': [256, 256],  # 简化网络架构
        },
        'verbose': 1,
        'learning_starts': 5000,
        'seed': 42,
        'tensorboard_log': './simplified/tensorboard_logs_simplified/'
    }
    
    print("创建简化版SBX SAC模型...")
    model = SAC(**model_kwargs)
    print("简化版SBX SAC Agent已创建。")
    print(f"网络架构: {model_kwargs['policy_kwargs']}")
    
    # 5. 设置训练参数
    print("\n--- 5. 设置训练参数 ---")
    max_steps_per_episode = cS.aD_new
    total_timesteps = 300_000  # 简化版使用较少的训练步数
    stop_training_value = -20
    eval_freq = 1_000
    n_eval_episodes = 100
    
    print(f"每回合最大步数: {max_steps_per_episode}")
    print(f"总训练步数: {total_timesteps}")
    print(f"状态空间维度: 4 (简化版)")
    
    # 6. 测试环境和智能体初始化
    print("\n--- 6. 测试环境和智能体初始化 ---")
    obs, _ = env.reset()
    print(f"环境重置成功，观察维度: {obs.shape} (4维)")
    
    action, _ = model.predict(obs, deterministic=False)
    print(f"智能体初始动作生成成功，动作维度: {action.shape}")
    print(f"动作值: [PPS缴费比例={action[0]:.4f}, 消费比例={action[1]:.4f}]")
    
    # 7. 设置评估回调
    print("\n--- 7. 设置评估和回调 ---")
    
    # 创建评估环境（评估模式）
    eval_env_base = OLGEnvV9SACSimplified(cS, paramS_for_rl, M_fixed, training_mode=False)
    eval_env = Monitor(eval_env_base)
    
    stop_callback = StopTrainingOnRewardThreshold(
        reward_threshold=stop_training_value, 
        verbose=1
    )
    
    # 设置评估回调
    eval_callback = EvalCallbackSimplified(
        eval_env,
        cS,
        paramS_for_rl,
        M_fixed,
        gamma=0.97,
        use_lifecycle_simulation=True,
        best_model_save_path='./simplified/best_model_simplified/',
        log_path='./simplified/logs_simplified/',
        eval_freq=eval_freq,
        n_eval_episodes=n_eval_episodes,
        deterministic=True,
        render=False,
        callback_on_new_best=stop_callback,
        verbose=1
    )
    
    print("🧬 简化版生命周期模拟评估回调已设置。")
    
    # 8. 开始训练
    print("\n--- 8. 开始训练 ---")
    print("使用简化版SBX SAC算法训练...")
    
    # 创建保存目录
    os.makedirs('./simplified/best_model_simplified/', exist_ok=True)
    os.makedirs('./simplified/logs_simplified/', exist_ok=True)
    os.makedirs('./simplified/tensorboard_logs_simplified/', exist_ok=True)
    
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
    model.save('./simplified/final_sac_agent_olg_simplified')

    # 保存配置
    config = {
        'cS': cS,
        'paramS_for_rl': paramS_for_rl,
        'M_fixed': M_fixed,
        'model_kwargs': model_kwargs,
        'training_time': training_time,
        'algorithm': 'SBX_SAC_Simplified',
        'state_space_dim': 4,
        'simplified_version': True,
        'fixed_macro_params': True,
        'jax_devices': str(jax.devices())
    }

    with open('./simplified/training_config_simplified.pkl', 'wb') as f:
        pickle.dump(config, f)
    
    print("简化版模型和配置已保存。")
    
    # 10. 评估训练好的Agent
    print("\n--- 10. 评估训练好的 Agent ---")
    
    print("使用固定宏观参数:")
    for key, value in M_fixed.items():
        print(f"  {key} = {value:.3f}")
    
    # 1. 无折现评估
    mean_reward_no_discount, std_reward_no_discount = evaluate_policy(
        model, eval_env, n_eval_episodes=100, deterministic=True
    )
    
    # 2. 生命周期模拟评估（简化版）
    eIdxM_annual_unified = None
    if hasattr(eval_callback, 'eIdxM_annual_unified'):
        eIdxM_annual_unified = eval_callback.eIdxM_annual_unified[:100, :]
        print("🎯 使用训练过程中生成的统一效率冲击路径进行最终评估")
    
    mean_reward_lifecycle, std_reward_lifecycle, lifecycle_results = evaluate_policy_lifecycle_simulation_simplified(
        model, cS, paramS_for_rl, M_fixed, n_sim=100, deterministic=True, gamma=0.97, random_seed=42,
        eIdxM_annual_input=eIdxM_annual_unified
    )
    
    print(f"\n📊 简化版评估结果:")
    print(f"❌ 无折现 (γ=1.0): {mean_reward_no_discount:.2f} ± {std_reward_no_discount:.2f}")
    print(f"🧬 生命周期模拟 (简化版): {mean_reward_lifecycle:.2f} ± {std_reward_lifecycle:.2f}")
    print(f"🏆 训练过程最佳结果: {eval_callback.best_mean_reward:.2f} (timestep {getattr(eval_callback, 'best_timestep', 'N/A')})")
    
    print(f"\n🎯 简化版特性:")
    print(f"  📊 状态空间: (k, k_pps, age, ε) - 4维，与VFI基本相同")
    print(f"  🔧 固定宏观参数: 消除宏观不确定性")
    print(f"  🏋️ 训练时: reward = u(c) * ∏_{{i=1}}^{{t}} s(i) - 累积存活概率加权")
    print(f"  📊 评估时: reward = u(c) - 纯效用，与VFI可比")
    print(f"  💡 理论基础: 数学等价于VFI的Bellman方程")
    
    print(f"\n🔍 简化版优势:")
    print(f"  ✅ 状态空间与VFI更接近，更公平的比较")
    print(f"  ✅ 消除宏观参数不确定性的影响")
    print(f"  ✅ 保持累积存活概率方法的理论等价性")
    print(f"  💡 推荐评估方式: 生命周期模拟 ({mean_reward_lifecycle:.2f})")
    
    # 更新配置信息
    config['mean_reward_no_discount'] = mean_reward_no_discount
    config['mean_reward_lifecycle'] = mean_reward_lifecycle
    config['std_reward_no_discount'] = std_reward_no_discount
    config['std_reward_lifecycle'] = std_reward_lifecycle
    config['lifecycle_results'] = lifecycle_results
    config['best_training_reward'] = eval_callback.best_mean_reward
    config['best_training_timestep'] = getattr(eval_callback, 'best_timestep', None)
    
    # 重新保存配置
    with open('./simplified/training_config_simplified.pkl', 'wb') as f:
        pickle.dump(config, f)
    
    print("简化版SBX SAC Agent 训练完成。")
    print("="*60)
    
    return model, config, M_fixed

if __name__ == "__main__":
    # 运行训练
    model, config, M_fixed = main()
    
    # 输出最终结果摘要
    print("\n" + "="*60)
    print("🎯 简化版SBX SAC 训练完成摘要")
    print("="*60)
    print(f"算法: {config['algorithm']}")
    print(f"状态空间维度: {config['state_space_dim']} (简化版)")
    print(f"固定宏观参数: {config['fixed_macro_params']}")
    print(f"训练时间: {config['training_time']:.2f} 秒")
    print(f"模型保存路径: ./simplified/final_sac_agent_olg_simplified.zip")
    print("")
    print("🏆 训练结果摘要:")
    print(f"  最佳训练结果: {config['best_training_reward']:.2f} (timestep {config['best_training_timestep']})")
    print("")
    print("🎯 简化版评估:")
    print(f"  ❌ 无折现: {config['mean_reward_no_discount']:.2f} ± {config['std_reward_no_discount']:.2f}")
    print(f"  🧬 生命周期模拟: {config['mean_reward_lifecycle']:.2f} ± {config['std_reward_lifecycle']:.2f}")
    print("")
    print("🔍 简化版特性:")
    print(f"  ✅ 状态空间: (k, k_pps, age, ε) - 与VFI基本相同")
    print(f"  ✅ 固定宏观参数: 消除宏观不确定性")
    print(f"  ✅ 累积存活概率方法: 理论与VFI完全等价")
    print(f"  💡 更公平的RL vs VFI比较")
    print("="*60) 