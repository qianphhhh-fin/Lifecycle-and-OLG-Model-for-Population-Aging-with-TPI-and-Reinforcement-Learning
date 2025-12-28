# --- START OF FILE main_olg_v9_transition.py ---

"""
OLG 模型 V9 - 过渡动态分析 (使用RL Agent, 内生劳动税率)

目标:
- 从一个现实的初始人口分布开始，模拟经济随时间演化的路径。
- 采用适应性预期/有界理性的假设，使用预训练的RL智能体来模拟家庭行为。
- [核心功能] 在每个时间步 t，通过内层迭代求解使当期政府预算平衡的劳动税率 τ_l,t。
- 观察关键宏观变量（K, L, w, r, τ_l）如何从初始状态逐步收敛到新的伪稳态。
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
import pickle
from typing import Dict
from main_olg_v9_utils import OLG_V9_Utils, ModelParameters, TempParamSHH
from sbx import SAC as SBX_SAC
import sys

class TransitionalDynamicsAnalyzer:
    """
    一个用于运行和分析OLG模型过渡动态的类，具有内生政策求解功能。
    """
    def __init__(self, T=150, n_sim=5000):
        self.T = T
        self.n_sim = n_sim
        self.rl_model, self.rl_config = self._load_rl_model()
        
        # [核心修改] 存储宏观变量的训练范围
        if 'rng_M' not in self.rl_config:
            raise ValueError("RL配置文件 'rl_config' 中缺少 'rng_M'，无法进行范围检测。")
        self.macro_variable_ranges = self.rl_config['rng_M']
        
        self.cS = self.rl_config['cS']
        self.paramS = self._initialize_paramS()
        self.eIdxM = self._generate_shock_paths()


    def _load_rl_model(self):
        """加载训练好的RL模型和配置"""
        print("--- 正在加载RL模型 ---")
        best_model_path = './py/best_model_sbx_full/best_model.zip'
        if not os.path.exists(best_model_path):
            best_model_path = './py/final_sac_agent_olg_sbx_full.zip'
        
        if not os.path.exists(best_model_path):
            raise FileNotFoundError("未找到RL模型文件。请确保模型存在。")

        model = SBX_SAC.load(best_model_path)
        config_path = best_model_path.replace('.zip', '_config.pkl')
        with open(config_path, 'rb') as f:
            config = pickle.load(f)
        print("✅ RL模型和配置加载成功。")
        return model, config

    def _initialize_paramS(self):
        """初始化基础的paramS结构"""
        paramS = ModelParameters()
        paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V = OLG_V9_Utils.EarningProcess_olgm(self.cS)
        paramS.leGridV = np.exp(paramS.leLogGridV.flatten())
        paramS.ageEffV_new = self.cS.ageEffV_new
        return paramS

    def _generate_shock_paths(self):
        """生成统一的效率冲击路径"""
        print("--- 正在生成效率冲击路径 ---")
        self.cS.nSim = self.n_sim
        return OLG_V9_Utils.LaborEndowSimulation_olgm_AgeGroup(self.cS, self.paramS)

    def run_analysis(self):
        """执行完整的过渡动态分析"""
        # 1. 计算人口的完整过渡路径
        Z_path_norm = self._compute_population_path()

        # 2. 过渡动态主循环
        results = self._simulate_transition(Z_path_norm)
        
        # 3. 可视化结果
        self._plot_results(results)

    def _compute_population_path(self) -> np.ndarray:
        """从初始分布开始，计算人口随时间演化的路径"""
        print("\n--- 1. 计算人口过渡路径 ---")
        popS_initial = OLG_V9_Utils.initPopulation(self.cS)
        popS_transition = OLG_V9_Utils.populationDynamics(popS_initial, self.cS)
        
        Z_path = np.zeros((self.cS.aD_new, self.T + 1))
        num_pop_periods = popS_transition.Z.shape[1]
        
        for t in range(self.T + 1):
            if t < num_pop_periods:
                Z_path[:, t] = popS_transition.Z[:, t]
            else:
                Z_path[:, t] = popS_transition.Z[:, -1]
        
        Z_path_norm = Z_path / np.sum(Z_path, axis=0, keepdims=True)
        print("✅ 人口过渡路径计算完成。")
        return Z_path_norm

    def _simulate_consumption_at_t(self, M_t: Dict, b_t: float) -> float:
        """
        [内层函数] 在给定的宏观环境 M_t 下，使用RL模型模拟，并返回当期的总消费 C_t。
        """
        R_k_net_factor = 1 + M_t['r_net_t']
        bV_payg = np.zeros(self.cS.aD_new)
        bV_payg[self.cS.aR_new:] = b_t
        
        paramS_sim = TempParamSHH(M_t['tau_l_t'], 0.10, self.cS.pps_active, self.cS.ageEffV_new)
        
        _, _, c_hist, _ = OLG_V9_Utils.HHSimulation_olgm_rl(
            self.rl_model, self.rl_config, self.eIdxM,
            R_k_net_factor, M_t['w_t'], 0.0,
            bV_payg, paramS_sim, self.cS
        )
        
        C_t = np.dot(np.mean(c_hist, axis=0), M_t['Z_t_norm'])
        return C_t

    def _get_prices_and_policy_at_t(self, K_t: float, Z_t_norm: np.ndarray, rho_prime_target: float) -> Dict:
        """[v3] 求解 t 期的瞬时均衡，并进行范围检测"""
        paramS_t = self.paramS
        paramS_t.ageMassV = Z_t_norm
        paramS_t.popGrowthForDebt = 0.0
        _, L_t = OLG_V9_Utils.LaborSupply_Huggett(self.eIdxM, self.cS, paramS_t, Z_t_norm)
        paramS_t.L_per_capita = L_t
        R_mkt_factor, w_t = OLG_V9_Utils.HHPrices_Huggett(K_t, L_t, self.cS)
        r_mkt = R_mkt_factor - 1 - self.cS.ddk
        r_net_t = r_mkt * (1 - self.cS.tau_k)
        mass_workers_t = np.sum(Z_t_norm[:self.cS.aR_new])
        avg_wage_t = w_t * L_t / mass_workers_t if mass_workers_t > 0 else 0
        b_t = rho_prime_target * avg_wage_t
        Y_t = self.cS.A * (K_t**self.cS.alpha) * (L_t**(1-self.cS.alpha))
        G_t = self.cS.gov_exp_frac_Y * Y_t
        B_t = self.cS.gov_debt_frac_Y * Y_t

        # --- [核心修改] 内层迭代，集成进度条和实时范围检测 ---
        tau_l_guess = 0.15
        max_inner_iter, tol_inner, damp_inner = 20, 1e-4, 0.5
        
        progress_bar_started = False
        out_of_range_warnings = [] # 用于收集本轮迭代中的所有范围警告

        for inner_iter in range(max_inner_iter):
            # --- 1. 实时范围检测 ---
            #    将当前猜测的宏观变量与训练范围进行比较
            out_of_range_msg = ""
            current_macro_vars = {
                'R_k_net_factor': 1 + r_net_t,
                'w_gross': w_t,
                'tau_l': tau_l_guess, # 使用当前的猜测值
            }
            for var, value in current_macro_vars.items():
                if var in self.macro_variable_ranges:
                    low, high = self.macro_variable_ranges[var]
                    if not (low - 1e-5 <= value <= high + 1e-5):
                        # 如果超出范围，构建一个简短的警告信息
                        out_of_range_msg += f" 🚨{var}={value:.3f}!"

            # --- 2. 模拟与GBC检查 ---
            M_t_guess = {"r_net_t": r_net_t, "w_t": w_t, "tau_l_t": tau_l_guess, "Z_t_norm": Z_t_norm}
            C_t_model = self._simulate_consumption_at_t(M_t_guess, b_t)
            gbc_residual = OLG_V9_Utils.check_gbc_residual(K_t, C_t_model, Y_t, G_t, B_t, w_t, r_mkt, 0, tau_l_guess, b_t, 0, 0, self.cS, paramS_t)
            
            # --- 3. 更新动态进度条 (现在包含范围警告) ---
            progress_msg = f"    内循环 (t={self.current_t}): iter={inner_iter+1}/{max_inner_iter}, τ_l={tau_l_guess:.5f}, GBC残差={gbc_residual:.3e}{out_of_range_msg}"
            sys.stdout.write(f"\r{progress_msg:<100}") # 增加宽度以容纳警告
            sys.stdout.flush()
            progress_bar_started = True

            if out_of_range_msg and out_of_range_msg not in out_of_range_warnings:
                 out_of_range_warnings.append(out_of_range_msg.strip())

            # --- 4. 收敛判断与更新 ---
            if abs(gbc_residual) < tol_inner:
                break
            
            tau_l_update = -gbc_residual / (w_t * L_t + 1e-9)
            tau_l_guess += damp_inner * tau_l_update
            tau_l_guess = np.clip(tau_l_guess, self.cS.tau_l_min, self.cS.tau_l_max)
        
        # --- 5. 清理与总结 ---
        if progress_bar_started:
            sys.stdout.write(f"\r{' ' * 100}\r")
            sys.stdout.flush()

        # 打印在迭代过程中出现过的所有不重复的范围警告
        if out_of_range_warnings:
            unique_warnings = " | ".join(sorted(list(set(out_of_range_warnings))))
            print(f"  💡 提示: t={self.current_t}, τ_l求解中曾超出范围: {unique_warnings}")

        if inner_iter == max_inner_iter - 1 and abs(gbc_residual) >= tol_inner:
            print(f"  🚨 警告: t={self.current_t}, τ_l 内层迭代未在 {max_inner_iter} 次内收敛。最终GBC残差: {gbc_residual:.2e}")
        
        tau_l_t = tau_l_guess
        
        return {"r_net_t": r_net_t, "w_t": w_t, "tau_l_t": tau_l_t, "b_t": b_t, "L_t": L_t}

    def _simulate_one_period_and_get_next_K(self, M_t: Dict, Z_t_norm: np.ndarray) -> float:
        """使用RL模型，在给定的宏观环境 M_t 下，模拟并计算下一期的总资本 K_{t+1}。"""
        R_k_net_factor = 1 + M_t['r_net_t']
        bV_payg = np.zeros(self.cS.aD_new)
        bV_payg[self.cS.aR_new:] = M_t['b_t']
        
        paramS_sim = TempParamSHH(M_t['tau_l_t'], 0.10, self.cS.pps_active, self.cS.ageEffV_new)
        
        k_hist, _, _, _ = OLG_V9_Utils.HHSimulation_olgm_rl(
            self.rl_model, self.rl_config, self.eIdxM,
            R_k_net_factor, M_t['w_t'], 0.0,
            bV_payg, paramS_sim, self.cS
        )
        
        # K_{t+1} = E[k'] = sum over age a { mean(k'_{a,t}) * mass(a, t) }
        # k'_{a,t} 对应于 k_hist 的 k_{a+1}
        # 使用 t 期的人口分布 Z_t_norm 来加权 t 期的储蓄决策，以形成 t+1 期的总资本
        k_prime_by_age = np.mean(k_hist[:, 1:], axis=0) # mean k' for each age group > 1
        weights = Z_t_norm[:-1] # mass of people who will be in age group a+1 at t+1
        
        K_next = np.dot(k_prime_by_age, weights)
        return K_next

    def _simulate_transition(self, Z_path_norm: np.ndarray) -> Dict:
        """执行过渡动态的主向前模拟循环。"""
        print("\n--- 2. 开始RL过渡动态向前模拟 (内生τ_l) ---")
        paths = {key: np.zeros(self.T + 1) for key in ["K", "tau_l", "w", "r_net", "L"]}
        paths["K"][0] = 3.5
        rho_prime_target = 0.4
        
        print("\n" + "="*90)
        print(f"{'时期 (t)':>10s} | {'资本 (K_t)':>12s} | {'劳动 (L_t)':>12s} | {'工资 (w_t)':>12s} | {'利率 (r_net_t)':>14s} | {'税率 (τ_l,t)':>14s}")
        print("-"*90)

        for t in range(self.T):
            self.current_t = t # 用于在警告信息中打印
            M_t = self._get_prices_and_policy_at_t(paths["K"][t], Z_path_norm[:, t], rho_prime_target)
            
            for key in ["r_net", "w", "tau_l", "L"]:
                paths[key][t] = M_t[f'{key}_t']
            
            print(f"{t:>10d} | {paths['K'][t]:>12.4f} | {paths['L'][t]:>12.4f} | {paths['w'][t]:>12.4f} | {paths['r_net'][t]:>14.4%} | {paths['tau_l'][t]:>14.4%}")
            
            paths["K"][t+1] = self._simulate_one_period_and_get_next_K(M_t, Z_path_norm[:, t])
        
        final_M = self._get_prices_and_policy_at_t(paths["K"][self.T], Z_path_norm[:, self.T], rho_prime_target)
        for key in ["r_net", "w", "tau_l", "L"]:
            paths[key][self.T] = final_M[f'{key}_t']

        print(f"{self.T:>10d} | {paths['K'][self.T]:>12.4f} | {paths['L'][self.T]:>12.4f} | {paths['w'][self.T]:>12.4f} | {paths['r_net'][self.T]:>14.4%} | {paths['tau_l'][self.T]:>14.4%}")
        print("-"*90)
        
        print("✅ 过渡动态模拟完成。")
        return {f"{key}_path": val for key, val in paths.items()}
    
    def _plot_results(self, results: Dict):
        """可视化过渡动态的结果"""
        print("\n--- 3. 绘制过渡动态路径 ---")
        time_axis = np.arange(self.T + 1)
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle('OLG 模型过渡动态 (RL Agent, 内生τ_l)', fontsize=16)

        axes[0, 0].plot(time_axis, results['K_path'], 'b-', lw=2)
        axes[0, 0].set_title('资本存量路径 K_t'); axes[0, 0].set_xlabel('时期 (t)'); axes[0, 0].grid(True)
        
        axes[0, 1].plot(time_axis, results['w_path'], 'r-', lw=2)
        axes[0, 1].set_title('工资路径 w_t'); axes[0, 1].set_xlabel('时期 (t)'); axes[0, 1].grid(True)

        axes[1, 0].plot(time_axis, results['r_net_path'], 'g-', lw=2)
        axes[1, 0].set_title('净利率路径 r_t'); axes[1, 0].set_xlabel('时期 (t)'); axes[1, 0].grid(True)
        
        axes[1, 1].plot(time_axis, results['tau_l_path'], 'm-', lw=2)
        axes[1, 1].set_title('劳动税率路径 τ_l,t'); axes[1, 1].set_xlabel('时期 (t)'); axes[1, 1].grid(True)

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.show()

if __name__ == "__main__":
    analyzer = TransitionalDynamicsAnalyzer(T=150, n_sim=5000)
    analyzer.run_analysis()