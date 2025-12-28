"""
OLG Model V8 SAC Training Script - Python Implementation
转换自MATLAB版本的main_olg_v8_SAC.m

使用Stable Baselines 3的SAC算法训练OLG模型智能体
"""

import numpy as np
import torch
import pickle
import time
from typing import Dict, Any, Tuple
from stable_baselines3 import SAC
from stable_baselines3.common.evaluation import evaluate_policy
from stable_baselines3.common.callbacks import EvalCallback, StopTrainingOnRewardThreshold
from stable_baselines3.common.monitor import Monitor
from stable_baselines3.common.vec_env import DummyVecEnv
import matplotlib.pyplot as plt

from olg_utils import OLGUtils
from olg_env_v8_sac import OLGEnvV8SAC

def main():
    """主训练函数"""
    print("=== OLG 模型 V8 - SAC 智能体训练 (Stable Baselines 3版本) ===")
    print("    (在线RL，宏观变量M作为环境参数)")
    print("    (决策变量：PPS缴费比例, 非PPS储蓄比例)")
    print("    (网络架构：SB3标准SAC实现)")
    
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
    
    # 4. 创建SAC Agent
    print("\n--- 4. 创建 SAC Agent ---")
    
    # SAC超参数（对应MATLAB官方标准）
    model_kwargs = {
        'policy': 'MlpPolicy',
        'env': env,
        'learning_rate': 1e-3,           # 对应MATLAB官方推荐学习率
        'buffer_size': int(1e6),         # 对应MATLAB的ExperienceBufferLength
        'batch_size': 64,                # 对应MATLAB的MiniBatchSize
        'tau': 5e-3,                     # 对应MATLAB的TargetSmoothFactor
        'gamma': 0.97,                   # 对应MATLAB的DiscountFactor
        'ent_coef': 0.2,                 # 对应MATLAB的Alpha熵系数
        'target_update_interval': 1,
        'gradient_steps': 1,
        'optimize_memory_usage': False,  # 修复：避免与handle_timeout_termination冲突
        'policy_kwargs': {
            'net_arch': [256, 256],      # 对应MATLAB的网络架构
            'activation_fn': torch.nn.ReLU,
        },
        'verbose': 1,
        'seed': 42,
        'device': 'cuda' if torch.cuda.is_available() else 'cpu'
    }
    
    model = SAC(**model_kwargs)
    print("SAC Agent已创建。")
    print(f"使用设备: {model.device}")
    print(f"网络架构: {model.policy_kwargs}")
    
    # 5. 设置训练参数（对应MATLAB官方标准配置）
    print("\n--- 5. 设置训练参数 ---")
    max_steps_per_episode = cS['aD_new']    # 对应MATLAB的MaxStepsPerEpisode
    total_timesteps = 50000 # max_episodes * max_steps_per_episode
    stop_training_value = -20               # 对应MATLAB的StopTrainingValue
    eval_freq = 10000                        # 评估频率
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
    eval_env = OLGEnvV8SAC(cS, paramS_for_rl, rng_M)
    M_eval = {
        'R_k_net_factor': 1.03,
        'w_gross': 2.0,
        'TR_total': 0.1,
        'b_payg_avg_retiree': 0.4,
        'tau_l': 0.15,
        'theta_payg_actual': 0.12
    }
    eval_env.set_macro_parameters(M_eval)
    eval_env = Monitor(eval_env)
    
    # 设置停止训练回调
    stop_callback = StopTrainingOnRewardThreshold(
        reward_threshold=stop_training_value, 
        verbose=1
    )
    
    # 设置评估回调
    eval_callback = EvalCallback(
        eval_env,
        best_model_save_path='./py/best_model/',
        log_path='./py/logs/',
        eval_freq=eval_freq,
        n_eval_episodes=n_eval_episodes,
        deterministic=True,
        render=False,
        callback_on_new_best=stop_callback,
        verbose=1
    )
    
    print("评估回调已设置。")
    
    # 8. 开始训练
    print("\n--- 8. 开始训练 ---")
    print("使用Stable Baselines 3 SAC算法训练...")
    
    # 创建保存目录
    import os
    os.makedirs('./py/best_model/', exist_ok=True)
    os.makedirs('./py/logs/', exist_ok=True)
    
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
    model.save('./py/final_sac_agent_olg_sb3')
    
    # 保存参数和环境配置
    config = {
        'cS': cS,
        'paramS_for_rl': paramS_for_rl,
        'rng_M': rng_M,
        'M_eval': M_eval,
        'model_kwargs': model_kwargs,
        'training_time': training_time
    }
    
    with open('./py/training_config.pkl', 'wb') as f:
        pickle.dump(config, f)
    
    print("最终模型和配置已保存。")
    
    # 10. 评估训练好的Agent
    print("\n--- 10. 评估训练好的 Agent ---")
    
    # 在固定参数下评估（通过env属性访问原始环境）
    eval_env.env.set_macro_parameters(M_eval)
    print("使用固定测试参数（与MATLAB版本一致）:")
    for key, value in M_eval.items():
        print(f"  {key} = {value:.3f}")
    
    # 进行评估
    mean_reward, std_reward = evaluate_policy(
        model, eval_env, n_eval_episodes=100, deterministic=True
    )
    
    print(f"评估完成。在固定M下的平均回报 (来自 100 次模拟): {mean_reward:.2f} ± {std_reward:.2f}")
    
    # 11. 绘制训练统计
    print("\n--- 11. 绘制训练统计 ---")
    try:
        plot_training_stats()
    except Exception as e:
        print(f"绘图失败: {e}")
    
    print("SAC Agent 训练和处理框架完成 (SB3版本)。")

def plot_training_stats():
    """绘制训练统计图"""
    import pandas as pd
    import os
    
    # 读取训练日志
    log_path = './py/logs/evaluations.npz'
    if os.path.exists(log_path):
        data = np.load(log_path)
        timesteps = data['timesteps']
        results = data['results']
        
        plt.figure(figsize=(12, 4))
        
        # 评估奖励图
        plt.subplot(1, 2, 1)
        mean_rewards = np.mean(results, axis=1)
        std_rewards = np.std(results, axis=1)
        plt.plot(timesteps, mean_rewards, 'b-', label='Mean Reward')
        plt.fill_between(timesteps, mean_rewards - std_rewards, 
                        mean_rewards + std_rewards, alpha=0.3)
        plt.xlabel('Timesteps')
        plt.ylabel('Reward')
        plt.title('Training Progress: Evaluation Rewards')
        plt.legend()
        plt.grid(True)
        
        # 奖励分布图
        plt.subplot(1, 2, 2)
        plt.boxplot(results.T, positions=timesteps[::max(1, len(timesteps)//10)])
        plt.xlabel('Timesteps')
        plt.ylabel('Reward')
        plt.title('Reward Distribution Over Training')
        plt.grid(True)
        
        plt.tight_layout()
        plt.savefig('./py/training_stats.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        print("训练统计图已保存到 ./py/training_stats.png")
    else:
        print("未找到训练日志文件。")



if __name__ == "__main__":
    main() 