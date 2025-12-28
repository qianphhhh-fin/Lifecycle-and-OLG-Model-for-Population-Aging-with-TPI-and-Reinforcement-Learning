# --- 开始文件：compare_rl_and_vfi_matlab.py (多环境评估版) ---

"""
在多个随机宏观环境下，比较VFI、全功能RL和朴素策略的优化结果

[多环境评估版]
- 脚本生成N个在预设范围内的随机宏观环境。
- 对每个环境：
  - VFI策略需要重新求解。
  - 同一个RL智能体被用于所有环境，测试其泛化能力。
  - 朴素策略规则不变，但行为随环境变化。
- 最终聚合所有环境下的评估结果，提供更鲁棒的性能比较。
- 所有模拟都使用相同的随机冲击路径，确保跨策略比较的公平性。
"""

import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pickle
import time
from typing import Dict, Any, Tuple, List
from scipy import stats

# MATLAB Engine导入
try:
    import matlab.engine
    MATLAB_AVAILABLE = True
except ImportError:
    MATLAB_AVAILABLE = False
    print("❌ MATLAB Engine 不可用，无法继续。")
    exit()

# 导入RL库和自定义模块
from sbx import SAC as SBX_SAC
from main_olg_v9_utils import OLG_V9_Utils, TempParamSHH

# 配置matplotlib中文字体支持
def setup_chinese_fonts():
    """设置matplotlib中文字体"""
    try:
        from matplotlib.font_manager import fontManager
        font_path = 'C:/Windows/Fonts/simhei.ttf' 
        if os.path.exists(font_path):
            fontManager.addfont(font_path)
            plt.rcParams['font.sans-serif'] = ['SimHei']
            print(f"✅ 使用中文字体: SimHei")
        else:
            chinese_fonts = ['SimHei', 'Microsoft YaHei', 'WenQuanYi Micro Hei', 'Heiti TC']
            available_fonts = [f.name for f in fontManager.ttflist]
            for font in chinese_fonts:
                if font in available_fonts:
                    plt.rcParams['font.sans-serif'] = [font]
                    print(f"✅ 使用中文字体: {font}")
                    return
            print("⚠️ 未找到指定中文字体，请确保已安装。")
    except Exception as e:
        print(f"⚠️ 设置中文字体时出错: {e}")
    plt.rcParams['axes.unicode_minus'] = False
setup_chinese_fonts()

class RLVFIComparatorMultiEnv:
    """[多环境版] 在一系列随机宏观环境下比较三种策略的性能。"""
    
    def __init__(self, num_test_envs: int = 10):

        self.vfi_grid_configs = [
            {'label': 'VFI_Med_hybrid (20x20)', 'nk': 20, 'nkpps': 20, 'nkprime': 20, 'npps': 20, 'solver': 'hybrid'},
        ]
                
        # [修改] 定义宏观参数的采样范围，而不是固定值
        self.rng_M = {
            'R_k_net_factor': [1.02, 1.04],
            'w_gross': [1.8, 2.2],
            'TR_total': [0.05, 0.15],
            'b_payg_avg_retiree': [0.3, 0.5],
            'tau_l': [0.10, 0.20],
            'theta_payg_actual': [0.08, 0.16]
        }
        self.num_test_envs = num_test_envs
        if num_test_envs == 1:
            self.rng_M = {
                'R_k_net_factor': [1.03,1.03],
                'w_gross': [1.8,1.8],
                'TR_total': [0.0,0.0],
                'b_payg_avg_retiree': [0.4,0.4],
                'tau_l': [0.15,0.15],
                'theta_payg_actual': [0.10,0.10]
            }

        self.test_environments = self._generate_test_environments()
        
        self.eng = matlab.engine.start_matlab()
        self.eng.addpath(os.getcwd(), nargout=0)
        print("✅ MATLAB Engine 已启动。")

    def _generate_test_environments(self) -> List[Dict[str, float]]:
        """生成N个用于测试的随机宏观环境"""
        np.random.seed(123) # 固定种子以复现环境
        environments = []
        print(f"🌍 正在生成 {self.num_test_envs} 个随机宏观测试环境...")
        for i in range(self.num_test_envs):
            env_params = {key: np.random.uniform(low, high) for key, (low, high) in self.rng_M.items()}
            environments.append(env_params)
            print(f"  环境 {i+1}: R_k={env_params['R_k_net_factor']:.3f}, w={env_params['w_gross']:.3f}, tau_l={env_params['tau_l']:.3f}")
        return environments

    def __del__(self):
        if hasattr(self, 'eng') and self.eng is not None:
            self.eng.quit()
            print("✅ MATLAB Engine已关闭。")

    def load_rl_model(self, use_best_model: bool = True) -> Tuple[Any, Dict]:
        # (此函数无需修改)
        model_path = './py/best_model_sbx_full/best_model.zip' if use_best_model else './py/final_sac_agent_olg_sbx_full.zip'
        if not os.path.exists(model_path): model_path = './py/final_sac_agent_olg_sbx_full.zip'
        config_path = model_path.replace('.zip', '_config.pkl')
        print(f"📁 正在加载全功能RL模型: {model_path}")
        model = SBX_SAC.load(model_path)
        with open(config_path, 'rb') as f: config = pickle.load(f)
        print("✅ 全功能模型和配置加载成功。")
        return model, config

    def run_matlab_vfi_for_env(self, M_env: Dict[str, float]) -> Dict[str, Any]:
        """[修改] 为指定的单个宏观环境求解VFI策略。"""
        # (逻辑与之前的 run_matlab_vfi 类似，但使用传入的 M_env)
        cS_python = OLG_V9_Utils.ParameterValues_HuggettStyle()
        (leLogGridV, leTrProbM, leProb1V) = OLG_V9_Utils.EarningProcess_olgm(cS_python)
        paramS_vfi_dict = {
            'leLogGridV': leLogGridV, 'leTrProbM': leTrProbM, 'leProb1V': leProb1V,
            'leGridV': np.exp(leLogGridV), 'ageEffV_new': cS_python.ageEffV_new,
            'tau_l': M_env['tau_l'],
            'theta_payg_actual_for_hh': M_env['theta_payg_actual'],
            'pps_tax_deferral_active': bool(cS_python.pps_active)
        }
        cS_matlab_dict = self._convert_dict_to_matlab_struct(cS_python.__dict__)
        paramS_vfi_matlab = self._convert_dict_to_matlab_struct(paramS_vfi_dict)
        bV_payg_vfi = np.zeros(cS_python.aD_new)
        if cS_python.aR_new < cS_python.aD_new:
            bV_payg_vfi[cS_python.aR_new:] = M_env['b_payg_avg_retiree']
        bV_payg_matlab = matlab.double(bV_payg_vfi.tolist())
        cS_matlab_dict['interpolation_method'] = 'spline'
        cPolM, kPolM, cPpsM, _ = self.eng.main_olg_v8_utils.HHSolution_VFI_Huggett(
            M_env['R_k_net_factor'], M_env['w_gross'], M_env['TR_total'],
            bV_payg_matlab, paramS_vfi_matlab, cS_matlab_dict, nargout=4)
        return {'kPolM': np.array(kPolM), 'cPpsPolM_choice': np.array(cPpsM),
                'M_test': M_env, 'cS_python_obj': cS_python, 'paramS_vfi_dict': paramS_vfi_dict}

    def _convert_dict_to_matlab_struct(self, py_dict: Dict) -> Dict:
        # (此函数无需修改)
        matlab_struct = {}
        for key, value in py_dict.items():
            if key in ['physAgeMap', 'interpolation_method']: continue
            if isinstance(value, np.ndarray): matlab_struct[key] = matlab.double(value.tolist())
            elif isinstance(value, list): matlab_struct[key] = matlab.double(value)
            elif isinstance(value, (int, float, bool)): matlab_struct[key] = float(value)
        return matlab_struct
        
    def _calculate_lifetime_utility(self, c_path: np.ndarray, cS: Any, use_survival_prob: bool) -> float:
        # (此函数无需修改)
        beta, aD = cS.beta, c_path.shape[0]
        s_transitionV = cS.s_1yr_transitionV.flatten()
        utility_sum, cumulative_discount = 0.0, 1.0
        for a_group in range(aD):
            c_val = c_path[a_group]
            _, u = OLG_V9_Utils.CES_utility(c_val, cS.sigma, cS)
            utility_sum += cumulative_discount * u
            if a_group < aD - 1:
                survival_factor = s_transitionV[a_group] if use_survival_prob else 1.0
                cumulative_discount *= (beta * survival_factor)
        return utility_sum

    def run_rule_of_thumb_simulation(self, M_env: Dict, cS_obj: Any, paramS_dict: Dict, 
                                     n_sim: int, eIdxM_group: np.ndarray) -> Dict:
        """[修改] 朴素策略模拟器，接收指定的宏观环境 M_env。"""
        SAVING_RATE, PPS_ALLOCATION_RATE = 0.20, 0.25
        leGridV, ageEffV_new = paramS_dict['leGridV'], paramS_dict['ageEffV_new']
        kMin, kppsMin, kMax, kppsMax = cS_obj.kMin, cS_obj.kppsMin, cS_obj.kMax, cS_obj.kppsMax
        aR_new, aD_new, tau_c, cFloor, pps_active = cS_obj.aR_new, cS_obj.aD_new, cS_obj.tau_c, cS_obj.cFloor, cS_obj.pps_active
        k_paths, kpps_paths, c_paths, cpps_paths = [], [], [], []
        paramS_hh = TempParamSHH(M_env['tau_l'], M_env['theta_payg_actual'], pps_active, ageEffV_new)
        bV_payg = np.zeros(aD_new)
        if aR_new < aD_new: bV_payg[aR_new:] = M_env['b_payg_avg_retiree']
        for i_sim in range(n_sim):
            k_current, kpps_current = kMin, kppsMin
            k_path, kpps_path, c_path, cpps_path = [], [], [], []
            for age_idx in range(aD_new):
                k_path.append(k_current); kpps_path.append(kpps_current)
                epsilon_val = leGridV[eIdxM_group[i_sim, age_idx] - 1]
                paramS_hh.current_pps_withdrawal = kpps_current * cS_obj.pps_withdrawal_rate if age_idx >= aR_new and pps_active else 0
                if age_idx < aR_new:
                    income, _, _ = OLG_V9_Utils.HHIncome_Huggett(k_current, M_env['R_k_net_factor'], M_env['w_gross'], M_env['TR_total'], bV_payg[age_idx], 0.0, age_idx, paramS_hh, cS_obj, epsilon_val)
                    total_savings = SAVING_RATE * income
                    c_pps = PPS_ALLOCATION_RATE * total_savings if pps_active else 0
                    k_prime = total_savings - c_pps
                    consumption_expenditure = income - total_savings
                else:
                    c_pps = 0
                    total_wealth = k_current + kpps_current * (1 - cS_obj.pps_tax_rate_withdrawal)
                    remaining_periods = aD_new - age_idx
                    consumption_expenditure = total_wealth / remaining_periods if remaining_periods > 0 else total_wealth
                    k_prime = 0
                current_c = max(cFloor, consumption_expenditure / (1 + tau_c))
                c_path.append(current_c); cpps_path.append(max(0, c_pps))
                k_current = max(kMin, min(kMax, k_prime))
                if pps_active:
                    pps_return_factor = 1 + ((M_env['R_k_net_factor'] - 1) + cS_obj.pps_return_rate_premium)
                    kpps_current = max(kppsMin, min(kppsMax, (kpps_current + c_pps - paramS_hh.current_pps_withdrawal) * pps_return_factor))
                else: kpps_current = kppsMin
            k_paths.append(k_path); kpps_paths.append(kpps_path); c_paths.append(c_path); cpps_paths.append(cpps_path)
        return {"k_path_rot": np.array(k_paths), "c_path_rot": np.array(c_paths), "cpps_path_rot": np.array(cpps_paths)}

# 在 RLVFIComparatorMultiEnv 类中，替换 run_comparison 函数

    def run_comparison(self, n_sim=500, random_seed=42, use_survival_prob_in_eval=True):
        """[重写] 主比较流程，在多个宏观环境下进行评估。"""
        print("\n" + "="*80)
        print(f"🔬 在 {self.num_test_envs} 个随机宏观环境下进行策略比较")
        print("=" * 80)

        # 1. 加载唯一的RL模型
        rl_model, rl_config = self.load_rl_model(use_best_model=True)
        
        # 2. 生成统一的随机冲击“剧本”
        temp_cs = OLG_V9_Utils.ParameterValues_HuggettStyle()
        _, tr_prob, p0 = OLG_V9_Utils.EarningProcess_olgm(temp_cs)
        eIdxM_group_0based = OLG_V9_Utils.MarkovChainSimulation_AgeGroup(n_sim, temp_cs, p0, tr_prob)
        eIdxM_group_for_sim = eIdxM_group_0based + 1
        
        # 3. 循环遍历每个测试环境
        all_results_by_env = []
        first_env_paths = {} # [新增] 用于存储第一个环境的路径

        for i, M_env in enumerate(self.test_environments):
            print(f"\n--- 正在评估环境 {i+1}/{self.num_test_envs} ---")
            print(f"    参数: R_k={M_env['R_k_net_factor']:.3f}, w={M_env['w_gross']:.3f}, tau_l={M_env['tau_l']:.3f}")

            # a. 为当前环境求解VFI策略
            # [修改] 确保返回 cPolM 以便进行最精确的模拟
            vfi_results = self.run_matlab_vfi_for_env(M_env)
            vfi_results['cPolM'] = self.eng.workspace['cPolM'] # 假设cPolM在MATLAB工作区
            
            cS_obj = vfi_results['cS_python_obj']
            
            # b. 在当前环境下模拟所有策略
            vfi_sim_paths = OLG_V9_Utils.HHSimulation_olgm_vfi_simplified(vfi_results, n_sim, eIdxM_group_for_sim)
            
            paramS_sim_for_rl = TempParamSHH(M_env['tau_l'], M_env['theta_payg_actual'], cS_obj.pps_active, cS_obj.ageEffV_new)
            bV_payg_for_rl = np.zeros(cS_obj.aD_new)
            if cS_obj.aR_new < cS_obj.aD_new: bV_payg_for_rl[cS_obj.aR_new:] = M_env['b_payg_avg_retiree']
            
            k_rl, kpps_rl, c_rl, cpps_rl = OLG_V9_Utils.HHSimulation_olgm_rl(
                rl_model, rl_config, eIdxM_group_for_sim,
                M_env['R_k_net_factor'], M_env['w_gross'], M_env['TR_total'],
                bV_payg_for_rl, paramS_sim_for_rl, cS_obj
            )
            
            rot_sim_paths = self.run_rule_of_thumb_simulation(M_env, cS_obj, vfi_results['paramS_vfi_dict'], n_sim, eIdxM_group_for_sim)

            # c. 计算当前环境下的效用
            utility_vfi = np.array([self._calculate_lifetime_utility(vfi_sim_paths['c_path_vfi'][j,:], cS_obj, use_survival_prob_in_eval) for j in range(n_sim)])
            utility_rl = np.array([self._calculate_lifetime_utility(c_rl[j,:], cS_obj, use_survival_prob_in_eval) for j in range(n_sim)])
            utility_rot = np.array([self._calculate_lifetime_utility(rot_sim_paths['c_path_rot'][j,:], cS_obj, use_survival_prob_in_eval) for j in range(n_sim)])
            
            all_results_by_env.append({
                'M_env': M_env,
                'mean_utility_vfi': np.mean(utility_vfi),
                'mean_utility_rl': np.mean(utility_rl),
                'mean_utility_rot': np.mean(utility_rot),
            })
            print(f"    环境 {i+1} 结果: VFI={np.mean(utility_vfi):.3f}, RL={np.mean(utility_rl):.3f}, RoT={np.mean(utility_rot):.3f}")

            # [新增] 如果是第一个环境，保存所有路径用于绘图
            if i == 0:
                first_env_paths['vfi'] = vfi_sim_paths
                first_env_paths['rl'] = {'k_path_rl': k_rl, 'kpps_path_rl': kpps_rl, 'c_path_rl': c_rl}
                first_env_paths['rot'] = rot_sim_paths
                first_env_paths['cS_obj'] = cS_obj # 保存该环境对应的cS对象

        # 4. 聚合所有环境的结果并分析
        self.analyze_and_plot(all_results_by_env, first_env_paths)

    def analyze_and_plot(self, all_results_by_env: List[Dict]):
        """[修改] 分析并绘制多环境评估的结果"""
        print("\n📈 分析与绘图 (多环境评估)...")
        
        # 提取每个策略在所有环境下的平均效用列表
        utilities_vfi = [res['mean_utility_vfi'] for res in all_results_by_env]
        utilities_rl = [res['mean_utility_rl'] for res in all_results_by_env]
        utilities_rot = [res['mean_utility_rot'] for res in all_results_by_env]
        
        # 计算跨环境的总体平均值和标准差
        mean_vfi = np.mean(utilities_vfi)
        std_vfi = np.std(utilities_vfi)
        mean_rl = np.mean(utilities_rl)
        std_rl = np.std(utilities_rl)
        mean_rot = np.mean(utilities_rot)
        std_rot = np.std(utilities_rot)
        
        # 使用配对t检验，因为每个策略都在相同的环境集上进行了评估
        t_rl_vfi, p_rl_vfi = stats.ttest_rel(utilities_rl, utilities_vfi)
        t_rl_rot, p_rl_rot = stats.ttest_rel(utilities_rl, utilities_rot)

        print("\n" + "=" * 80 + "\n📋 跨环境性能比较摘要\n" + "=" * 80)
        print(f" - VFI 平均效用 (跨环境):         {mean_vfi:.4f} ± {std_vfi:.4f}")
        print(f" - RL (全功能) 平均效用 (跨环境): {mean_rl:.4f} ± {std_rl:.4f}")
        print(f" - 朴素策略 平均效用 (跨环境):   {mean_rot:.4f} ± {std_rot:.4f}")
        print("-" * 40)
        print(f" - RL vs VFI:  效用差异 {mean_rl - mean_vfi:+.4f}, p-value={p_rl_vfi:.4f} ({'显著' if p_rl_vfi < 0.05 else '不显著'})")
        print(f" - RL vs 朴素: 效用差异 {mean_rl - mean_rot:+.4f}, p-value={p_rl_rot:.4f} ({'显著' if p_rl_rot < 0.05 else '不显著'})")
        
                # --- Part 2: 绘图 (扩展为 2x2 布局) ---
        fig, axes = plt.subplots(2, 2, figsize=(18, 14))
        fig.suptitle(f'跨 {self.num_test_envs} 个随机宏观环境的策略性能比较', fontsize=18, y=0.98)
        axes = axes.flatten() # 将2x2的axes数组展平为一维，方便索引

     # 图1: 效用箱线图 (与之前相同)
        ax = axes[0]
        ax.boxplot([utilities_vfi, utilities_rl, utilities_rot],
                labels=['VFI', 'RL (全功能)', '朴素策略'],
                patch_artist=True,
                boxprops=dict(facecolor='lightblue', color='blue'),
                medianprops=dict(color='red', linewidth=2))
        ax.set_ylabel('平均终身效用')
        ax.set_title('A. 策略性能的稳健性 (跨环境分布)', fontsize=14)
        ax.grid(True, linestyle='--', alpha=0.6)
        
        # 图2: 效用增益分布图 (与之前相同)
        ax = axes[1]
        rl_gain = (np.array(utilities_rl) - np.array(utilities_vfi)) / np.abs(utilities_vfi) * 100
        rot_gain = (np.array(utilities_rot) - np.array(utilities_vfi)) / np.abs(utilities_vfi) * 100
        ax.hist(rl_gain, bins=10, alpha=0.7, label=f'RL vs VFI (均值: {np.mean(rl_gain):.1f}%)', color='red')
        ax.hist(rot_gain, bins=10, alpha=0.7, label=f'朴素 vs VFI (均值: {np.mean(rot_gain):.1f}%)', color='green')
        ax.axvline(0, color='black', linestyle='--')
        ax.set_xlabel('效用增益 (%) [相对于VFI]')
        ax.set_ylabel('环境数量')
        ax.set_title('B. 效用增益分布', fontsize=14)
        ax.legend()
        ax.grid(True, linestyle='--', alpha=0.6)

        # --- [新增] Part 3: 绘制第一个环境的生命周期路径 ---
        if path_data:
            age_groups = np.arange(path_data['cS_obj'].aD_new)
            
            # 计算平均路径
            c_vfi_mean = np.mean(path_data['vfi']['c_path_vfi'], axis=0)
            k_vfi_mean = np.mean(path_data['vfi']['k_path_vfi'], axis=0)
            
            c_rl_mean = np.mean(path_data['rl']['c_path_rl'], axis=0)
            k_rl_mean = np.mean(path_data['rl']['k_path_rl'], axis=0)
            
            c_rot_mean = np.mean(path_data['rot']['c_path_rot'], axis=0)
            k_rot_mean = np.mean(path_data['rot']['k_path_rot'], axis=0)

            # 图3: 平均消费路径
            ax = axes[2]
            ax.plot(age_groups, c_vfi_mean, 'o-', label='VFI', color='blue')
            ax.plot(age_groups, c_rl_mean, 's--', label='RL (全功能)', color='red')
            ax.plot(age_groups, c_rot_mean, '^:', label='朴素策略', color='green')
            ax.set_xlabel('年龄组')
            ax.set_ylabel('平均消费')
            ax.set_title('C. 平均消费生命周期路径 (第一个环境)', fontsize=14)
            ax.legend()
            ax.grid(True, linestyle='--', alpha=0.6)

            # 图4: 平均资产路径
            ax = axes[3]
            ax.plot(age_groups, k_vfi_mean, 'o-', label='VFI', color='blue')
            ax.plot(age_groups, k_rl_mean, 's--', label='RL (全功能)', color='red')
            ax.plot(age_groups, k_rot_mean, '^:', label='朴素策略', color='green')
            ax.set_xlabel('年龄组')
            ax.set_ylabel('平均非PPS资产 (k)')
            ax.set_title('D. 平均资产生命周期路径 (第一个环境)', fontsize=14)
            ax.legend()
            ax.grid(True, linestyle='--', alpha=0.6)

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        save_path = './py/multi_env_comparison_with_paths.png'
        plt.savefig(save_path, dpi=300)
        print(f"📈 多环境比较图表已保存到: {save_path}")
        plt.show()

def main():
    if not MATLAB_AVAILABLE: return
    # 可以通过修改 num_test_envs 来控制评估环境的数量
    comparator = RLVFIComparatorMultiEnv(num_test_envs=1)
    comparator.run_comparison(n_sim=500, use_survival_prob_in_eval=True)

if __name__ == "__main__":
    main()