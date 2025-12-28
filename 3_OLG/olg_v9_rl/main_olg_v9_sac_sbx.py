# --- 开始文件：train_sac_agent_olg_v9_sbx.py (修正版) ---

"""
OLG Model V9 SAC Training Script - SBX (Stable Baselines Jax) Implementation

[修正版]
- 适配了重构后的 OLGEnvV9SAC (全功能环境) 和 HHSimulation_olgm_rl (全功能模拟器)。
- 训练环境现在是 OLGEnvV9SAC，它在每轮开始时会采样不同的宏观变量 M。
- 评估环境也是 OLGEnvV9SAC，但宏观变量被固定，以进行稳定评估。
- 生命周期模拟评估函数 evaluate_... 已更新，以调用新的 HHSimulation_olgm_rl，并正确传递参数。
- 确保了从训练到评估再到最终比较的逻辑一致性。
"""

import numpy as np
import jax
import pickle
import time
import os
from typing import Dict, Any, Tuple
from sbx import SAC
from stable_baselines3.common.evaluation import evaluate_policy
from stable_baselines3.common.callbacks import EvalCallback, StopTrainingOnRewardThreshold
from stable_baselines3.common.monitor import Monitor
import matplotlib.pyplot as plt

# [修正] 导入新的全功能环境和模拟器
from main_olg_v9_utils import OLG_V9_Utils, OLGEnvV9SAC, TempParamSHH


# --- [新增] 专门用于评估模式的绘图函数 ---
def plot_evaluation_results(results: Dict, cS: Any):
    """可视化评估模式的结果"""
    print("\n📈 生成评估结果图表...")
    
    fig, axes = plt.subplots(2, 2, figsize=(20, 16))
    fig.suptitle('预训练RL智能体性能评估', fontsize=20, y=0.97)
    axes = axes.flatten()
    
    # 图A: 终身效用分布
    mean_u = results['mean_utility_rl']
    std_u = results['std_utility_rl']
    axes[0].hist(results['lifetime_utility_rl'], bins=30, density=True, alpha=0.7, label='效用分布')
    axes[0].axvline(mean_u, color='r', linestyle='--', label=f'平均值: {mean_u:.2f}')
    axes[0].set_title(f'A. 终身效用分布 (μ={mean_u:.2f}, σ={std_u:.2f})', fontsize=16)
    axes[0].set_xlabel('终身期望效用', fontsize=12)
    axes[0].set_ylabel('概率密度', fontsize=12)
    axes[0].legend()
    axes[0].grid(True, linestyle='--', alpha=0.6)

    age_groups = np.arange(cS.aD_new)
    
    # 图B: 平均消费路径
    axes[1].plot(age_groups, np.mean(results['c_path_rl'], axis=0), 'o-', color='red', label='平均消费')
    axes[1].set_title('B. 平均消费生命周期路径', fontsize=16)
    axes[1].set_xlabel('年龄组', fontsize=12)
    axes[1].set_ylabel('消费 (c)', fontsize=12)
    axes[1].legend()
    axes[1].grid(True, linestyle='--', alpha=0.6)
    
    # 图C: 平均资产路径
    axes[2].plot(age_groups, np.mean(results['k_path_rl'], axis=0), 's--', color='blue', label='平均非PPS资产 (k)')
    axes[2].set_title('C. 平均资产生命周期路径', fontsize=16)
    axes[2].set_xlabel('年龄组', fontsize=12)
    axes[2].set_ylabel('资产 (k)', fontsize=12)
    axes[2].legend()
    axes[2].grid(True, linestyle='--', alpha=0.6)

    # 图D: 平均PPS缴费路径
    axes[3].plot(age_groups, np.mean(results['cpps_path_rl'], axis=0), '^:', color='green', label='平均PPS缴费')
    axes[3].axvline(x=cS.aR_new - 1, color='gray', linestyle='--', label='退休年龄')
    axes[3].set_title('D. 平均PPS缴费生命周期路径', fontsize=16)
    axes[3].set_xlabel('年龄组', fontsize=12)
    axes[3].set_ylabel('PPS缴费 (c_pps)', fontsize=12)
    axes[3].legend()
    axes[3].grid(True, linestyle='--', alpha=0.6)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    save_path = './py/rl_agent_evaluation.png'
    plt.savefig(save_path, dpi=300)
    print(f"✅ 评估图表已保存到: {save_path}")
    plt.show()


# --- [新增] 专门的评估模式主函数 ---
def run_evaluation_only():
    """仅加载并评估最佳模型"""
    print("\n" + "="*80)
    print("🚀 进入评估模式 (Evaluation-Only Mode)")
    print("="*80)

    # 1. 定义模型和配置路径
    best_model_path = './py/best_model_sbx_full/best_model.zip'
    config_path = best_model_path.replace('.zip', '_config.pkl')

    if not os.path.exists(best_model_path) or not os.path.exists(config_path):
        print(f"❌ 错误: 找不到模型 '{best_model_path}' 或配置文件 '{config_path}'。")
        print("请先运行训练模式以生成模型。")
        return

    # 2. 加载模型和配置
    print(f"📁 正在加载模型: {best_model_path}")
    model = SAC.load(best_model_path)
    print(f"📁 正在加载配置: {config_path}")
    with open(config_path, 'rb') as f:
        config = pickle.load(f)
    
    cS = config['cS']
    paramS_for_rl = config['paramS_for_rl']
    rng_M = config['rng_M']
    
    # 3. 运行大规模模拟评估
    # 使用与最终比较脚本(test_vfi_grid_search.py)完全一致的评估环境
    mean_reward, std_reward, lifecycle_results = evaluate_policy_lifecycle_simulation_age_group(
        model, cS, paramS_for_rl, rng_M, 
        n_sim=1000, # 使用更多的模拟个体以获得更平滑的路径
        random_seed=42, 
        verbose=True
    )
    
    print("\n" + "="*50)
    print("📈 最终评估结果")
    print(f"   平均终身效用: {mean_reward:.4f} ± {std_reward:.4f}")
    print("="*50)

    # 4. 可视化结果
    plot_evaluation_results(lifecycle_results, cS)

def evaluate_policy_lifecycle_simulation_age_group(
    model: Any, cS: Any, paramS_for_rl: Dict, rng_M: Dict, 
    n_sim: int = 50, deterministic: bool = True, 
    random_seed: int = 42, verbose: bool = True,
    eIdxM_group_input: np.ndarray = None) -> Tuple[float, float, Dict]:
    """
    [修正版] 使用生命周期模拟评估RL模型，适配全功能环境和模拟器。
    
    此函数现在调用新的 HHSimulation_olgm_rl，该函数内部处理与 OLGEnvV9SAC
    环境的交互，确保逻辑一致。
    """
    sim_start_time = time.time()

    # 1. 设置固定的宏观经济参数用于评估
    M_fixed = {
        'R_k_net_factor': 1.03, 'w_gross': 2.0, 'TR_total': 0.1,
        'b_payg_avg_retiree': 0.4, 'tau_l': 0.15, 'theta_payg_actual': 0.12
    }
    
    if verbose:
        print(f"\n🧬 年龄组生命周期模拟评估 (全功能环境版)")
        print(f"   (n_sim={n_sim}, seed={random_seed})")
        print("📊 固定测试参数:", M_fixed)

    # 2. 生成或使用预设的效率冲击路径
    np.random.seed(random_seed)
    aD_new = int(cS.aD_new)
    if eIdxM_group_input is None:
        if verbose: print("🔄 生成年龄组效率冲击路径...")
        eIdxM_group = OLG_V9_Utils.MarkovChainSimulation_AgeGroup(
            n_sim, cS, paramS_for_rl['leProb1V'], paramS_for_rl['leTrProbM'])
    else:
        if verbose: print("🔄 使用预生成的年龄组效率冲击路径...")
        eIdxM_group = eIdxM_group_input
    
    # [修正] eIdxM_group 需要是1-based索引
    eIdxM_group_1based = eIdxM_group + 1

    # 3. 准备模拟参数
    bV_payg_group = np.zeros(aD_new)
    if int(cS.aR_new) < aD_new:
        bV_payg_group[int(cS.aR_new):] = M_fixed['b_payg_avg_retiree']

    # [新] 模拟器需要一个包含所有训练时参数的配置字典
    rl_config = {'cS': cS, 'paramS_for_rl': paramS_for_rl, 'rng_M': rng_M}
    
    # [新] 模拟器还需要一个临时的 paramS 对象，包含当前迭代的税率等
    paramS_sim = TempParamSHH(
        M_fixed['tau_l'],
        M_fixed['theta_payg_actual'],
        cS.pps_active,
        cS.ageEffV_new
    )
    
    # 4. [核心修正] 调用与全功能环境适配的 HHSimulation_olgm_rl 模拟器
    if verbose: print("🚀 调用 HHSimulation_olgm_rl 进行核心模拟...")
    k_path_rl, kpps_path_rl, c_path_rl, cpps_path_rl = OLG_V9_Utils.HHSimulation_olgm_rl(
        model, rl_config, eIdxM_group_1based,
        M_fixed['R_k_net_factor'], M_fixed['w_gross'], M_fixed['TR_total'],
        bV_payg_group, paramS_sim, cS
    )
    if verbose: print("✅ 核心模拟完成。")

    # 5. 计算生命周期效用 (此部分逻辑不变)
    if verbose: print("📊 计算VFI等价的生命周期效用...")
    lifetime_utility_rl = np.zeros(n_sim)
    beta = float(cS.beta)
    s_transitionV = cS.s_1yr_transitionV.flatten()
    for i in range(n_sim):
        utility_sum = 0.0
        cumulative_discount = 1.0
        for a_group in range(aD_new):
            c = c_path_rl[i, a_group]
            _, u = OLG_V9_Utils.CES_utility(c, cS.sigma, cS)
            utility_sum += cumulative_discount * u
            if a_group < aD_new - 1:
                cumulative_discount *= (beta * s_transitionV[a_group])
        lifetime_utility_rl[i] = utility_sum

    mean_utility_rl = np.mean(lifetime_utility_rl)
    std_utility_rl = np.std(lifetime_utility_rl)
    sim_time = time.time() - sim_start_time

    if verbose:
        print(f"✅ 评估完成，耗时: {sim_time:.2f} 秒")
        print(f"📊 RL生命周期效用: {mean_utility_rl:.4f} ± {std_utility_rl:.4f}")
    
    # 6. 整理并返回结果
    lifecycle_results = {
        'mean_utility_rl': mean_utility_rl,
        'std_utility_rl': std_utility_rl,
        'lifetime_utility_rl': lifetime_utility_rl,
        'k_path_rl': k_path_rl,
        'c_path_rl': c_path_rl,
        'cpps_path_rl': cpps_path_rl,
        'eIdxM_group': eIdxM_group, # 返回0-based索引
    }
    return mean_utility_rl, std_utility_rl, lifecycle_results


class EvalCallbackWithDiscount(EvalCallback):
    """[修正版] 自定义EvalCallback，适配全功能环境评估"""
    
    def __init__(self, eval_env, cS, paramS_for_rl, rng_M, **kwargs):
        super().__init__(eval_env, **kwargs)
        self.cS = cS
        self.paramS_for_rl = paramS_for_rl
        self.rng_M = rng_M
        self._generate_unified_efficiency_shocks()
        print(f"🔧 使用生命周期模拟评估回调 (VFI等价)")

    def _generate_unified_efficiency_shocks(self):
        print("🎲 生成统一效率冲击序列（用于所有评估）...")
        np.random.seed(42) # 固定种子
        n_sim_target = max(self.n_eval_episodes, 100)
        self.eIdxM_group_unified = OLG_V9_Utils.MarkovChainSimulation_AgeGroup(
            n_sim_target, self.cS, self.paramS_for_rl['leProb1V'], self.paramS_for_rl['leTrProbM']
        )
        print(f"✅ 统一效率冲击矩阵生成完成: {self.eIdxM_group_unified.shape}")

    def _on_step(self) -> bool:
        # [简化] 只保留核心的评估和保存逻辑
        if self.eval_freq > 0 and self.n_calls % self.eval_freq == 0:
            
            # [核心修正] 调用适配新环境的评估函数
            mean_reward, std_reward, _ = evaluate_policy_lifecycle_simulation_age_group(
                self.model, self.cS, self.paramS_for_rl, self.rng_M,
                n_sim=self.n_eval_episodes, 
                deterministic=self.deterministic,
                random_seed=42,
                verbose=False,
                eIdxM_group_input=self.eIdxM_group_unified
            )
            
            if self.verbose > 0:
                print(f"Eval @ T={self.num_timesteps} - Mean reward: {mean_reward:.4f} +/- {std_reward:.4f}")
            
            self.logger.record("eval/mean_reward", mean_reward)
            self.logger.dump(self.num_timesteps)

            if mean_reward > self.best_mean_reward:
                self.best_mean_reward = mean_reward
                if self.verbose > 0:
                    print(f"🏆 New best mean reward: {self.best_mean_reward:.4f}")
                if self.best_model_save_path is not None:
                    self.model.save(os.path.join(self.best_model_save_path, "best_model"))
            
            # (可选) 可以在这里添加停止训练的回调逻辑
            if self.callback is not None:
                 return self.callback.on_step()
                 
        return True


def main(args):
    """主函数，根据命令行参数选择模式"""
    
    # [新增] 根据 --eval_only 标志选择执行路径
    if args.eval_only:
        run_evaluation_only()
        return # 评估后直接退出

    # --- 以下是原始的训练流程 ---
    print("=== OLG 模型 V9 - SAC 智能体训练 (全功能环境版) ===")
    
    # 1. 初始化参数
    print("\n--- 1. 初始化参数 ---")
    cS = OLG_V9_Utils.ParameterValues_HuggettStyle()
    leLogGridV, leTrProbM, leProb1V = OLG_V9_Utils.EarningProcess_olgm(cS)
    paramS_for_rl = {
        'leLogGridV': leLogGridV, 'leTrProbM': leTrProbM, 'leProb1V': leProb1V,
        'leGridV': np.exp(leLogGridV), 'ageEffV_new': cS.ageEffV_new
    }
    
    # 2. 定义宏观状态M的采样范围
    rng_M = {
        'R_k_net_factor': [1.01, 1.05], 'w_gross': [1.5, 2.5],
        'TR_total': [0.0, 0.2], 'b_payg_avg_retiree': [0.1, 0.8],
        'tau_l': [0.05, 0.25], 'theta_payg_actual': [0.05, 0.20]
    }
    
    # 3. 创建训练环境
    print("\n--- 3. 创建强化学习环境 ---")
    env = OLGEnvV9SAC(cS, paramS_for_rl, rng_M, training_mode=True)
    env = Monitor(env)
    
    # 4. 创建SAC Agent
    model_kwargs = {
        'policy': 'MlpPolicy', 
        'env': env, 
        'learning_rate': 1e-4,
        'buffer_size': int(1e6), 
        'batch_size': 512, 
        'tau': 5e-3,
        'gamma': 0.97, 
        'ent_coef': 'auto',
        'policy_kwargs': {'net_arch': [384, 384]},
        'verbose': 1, 
        'learning_starts': 5000, 
        'seed': 42,
        'tensorboard_log': './py/tensorboard_logs_sbx/'
    }
    model = SAC(**model_kwargs)
    
    # 5. 设置评估回调
    print("\n--- 5. 设置评估回调 ---")
    # [修正] 确保评估环境的宏观变量范围与训练时一致
    eval_env = OLGEnvV9SAC(cS, paramS_for_rl, rng_M, training_mode=False)
    eval_env = Monitor(eval_env)
    eval_callback = EvalCallbackWithDiscount(
        eval_env, cS, paramS_for_rl, rng_M,
        best_model_save_path='./py/best_model_sbx_full/',
        log_path='./py/logs_sbx_full/',
        eval_freq=5000, n_eval_episodes=100, deterministic=True, verbose=1
    )
    
    # 6. 开始训练
    print("\n--- 6. 开始训练 ---")
    os.makedirs('./py/best_model_sbx_full/', exist_ok=True)
    os.makedirs('./py/logs_sbx_full/', exist_ok=True)
    
    start_time = time.time()
    model.learn(total_timesteps=500_000, callback=eval_callback, progress_bar=True)
    training_time = time.time() - start_time
    print(f"训练用时: {training_time:.2f} 秒")
    
    # 7. 保存最终模型和配置
    print("\n--- 7. 保存最终模型和配置 ---")
    final_model_path = './py/final_sac_agent_olg_sbx_full.zip'
    model.save(final_model_path.replace('.zip', ''))
    config = {
        'cS': cS, 'paramS_for_rl': paramS_for_rl, 'rng_M': rng_M,
        'model_kwargs': {k: v for k, v in model_kwargs.items() if k != 'env'},
        'training_time': training_time, 'algorithm': 'SBX_SAC_Full_Env'
    }
    with open(final_model_path.replace('.zip', '_config.pkl'), 'wb') as f:
        pickle.dump(config, f)
    best_model_path = './py/best_model_sbx_full/best_model.zip'
    if os.path.exists(best_model_path):
        with open(best_model_path.replace('.zip', '_config.pkl'), 'wb') as f:
            pickle.dump(config, f)
            
    print(f"最终模型保存在: {final_model_path}")
    print("训练和配置保存完成。")
    
    # 8. 最终评估
    print("\n--- 8. 最终评估 ---")
    mean_reward, std_reward, _ = evaluate_policy_lifecycle_simulation_age_group(
        model, cS, paramS_for_rl, rng_M, n_sim=200, random_seed=42, verbose=True
    )
    print(f"\n最终模型在固定宏观环境下的生命周期效用: {mean_reward:.4f} +/- {std_reward:.4f}")


import argparse # [新增] 导入命令行参数模块

if __name__ == "__main__":
    # [新增] 解析命令行参数
    parser = argparse.ArgumentParser(description="训练或评估OLG模型的SAC智能体。")
    parser.add_argument(
        '--eval_only',
        action='store_true',
        help='如果设置此标志，则只运行评估模式，不进行训练。'
    )
    args = parser.parse_args()
    
    main(args)