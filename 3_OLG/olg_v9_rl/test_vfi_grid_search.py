# --- 开始文件：test_vfi_grid_search.py (最终整合版) ---

"""
在固定的宏观环境下，对多种策略进行全面的性能和行为比较。

[最终整合版]
- 比较对象:
  1. 多种不同配置的VFI策略 (通过MATLAB求解)
  2. 一个预训练的全功能RL策略
  3. 一个基于规则的朴素策略 (Rule of Thumb)
- 评估维度:
  1. 最终性能: 平均终身效用
  2. 经济学行为: 平均生命周期路径 (资产、消费、PPS缴费)
- 核心特性:
  - 所有策略都在完全相同的固定宏观环境和随机冲击路径下进行评估，确保绝对公平。
  - 模块化设计，易于增删VFI配置或修改朴素策略规则。
  - 生成信息丰富的组合图表，一站式完成性能排序和行为分析。
"""

import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pickle
import time
from typing import Dict, Any, Tuple, List

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


class ComprehensiveStrategyComparator:
    """在固定宏观环境下，对VFI, RL, 朴素策略进行全面的性能和行为比较。"""
    
    def __init__(self):
        # 1. 定义要测试的VFI网格配置列表
        self.vfi_grid_configs = [
            # {'label': 'VFI_Hybrid (20x20)', 'nk': 5, 'nkpps': 5, 'nkprime': 5, 'npps': 5, 'solver': 'hybrid'},
            {'label': 'VFI_Vectorized (5x5)', 'nk': 5, 'nkpps': 5, 'nkprime': 5, 'npps': 5, 'solver': 'vectorized_grid'},
            {'label': 'VFI_Vectorized (20x20)', 'nk': 20, 'nkpps': 20, 'nkprime': 20, 'npps': 20, 'solver': 'vectorized_grid'},
            {'label': 'VFI_Vectorized (50x50)', 'nk': 50, 'nkpps': 50, 'nkprime': 50, 'npps': 50, 'solver': 'vectorized_grid'},
        ]
        
        # 2. 定义一个固定的宏观测试环境
        self.M_fixed = {
            'R_k_net_factor': 1.03, 'w_gross': 2.0, 'TR_total': 0.1,
            'b_payg_avg_retiree': 0.4, 'tau_l': 0.15, 'theta_payg_actual': 0.12
        }
        print("="*80)
        print("🌍 固定宏观测试环境已设置:")
        for key, val in self.M_fixed.items():
            print(f"  - {key}: {val}")
        print("="*80)

        # 3. 启动MATLAB引擎
        self.eng = matlab.engine.start_matlab()
        self.eng.addpath(os.getcwd(), nargout=0)
        print("✅ MATLAB Engine 已启动。")

    def __del__(self):
        if hasattr(self, 'eng') and self.eng is not None:
            self.eng.quit()
            print("✅ MATLAB Engine已关闭。")

    def load_rl_model(self, use_best_model: bool = True) -> Tuple[Any, Dict]:
        model_path = './py/best_model_sbx_full/best_model.zip' if use_best_model else './py/final_sac_agent_olg_sbx_full.zip'
        if not os.path.exists(model_path): model_path = './py/final_sac_agent_olg_sbx_full.zip'
        config_path = model_path.replace('.zip', '_config.pkl')
        print(f"📁 正在加载全功能RL模型: {model_path}")
        model = SBX_SAC.load(model_path)
        with open(config_path, 'rb') as f: config = pickle.load(f)
        print("✅ 全功能模型和配置加载成功。")
        return model, config

    def run_matlab_vfi_for_config(self, vfi_config: Dict[str, Any], base_cS: Any) -> Dict[str, Any]:
        """为指定的VFI网格配置求解策略，使用基准cS确保参数一致性。"""
        print(f"\n--- 求解VFI策略: {vfi_config['label']} (求解器: {vfi_config['solver']}) ---")
        start_time = time.time()
        
        import copy
        cS_python = copy.deepcopy(base_cS)
        cS_python.nk, cS_python.nkpps, cS_python.nkprime, cS_python.npps = vfi_config['nk'], vfi_config['nkpps'], vfi_config['nkprime'], vfi_config['npps']
        cS_python = OLG_V9_Utils.generateGrids(cS_python)
        
        (leLogGridV, leTrProbM, leProb1V) = OLG_V9_Utils.EarningProcess_olgm(cS_python)
        paramS_vfi_dict = {
            'leLogGridV': leLogGridV, 'leTrProbM': leTrProbM, 'leProb1V': leProb1V,
            'leGridV': np.exp(leLogGridV), 'ageEffV_new': cS_python.ageEffV_new,
            'tau_l': self.M_fixed['tau_l'],
            'theta_payg_actual_for_hh': self.M_fixed['theta_payg_actual'],
            'pps_tax_deferral_active': bool(cS_python.pps_active)
        }
        
        cS_matlab_dict = self._convert_dict_to_matlab_struct(cS_python.__dict__)
        paramS_vfi_matlab = self._convert_dict_to_matlab_struct(paramS_vfi_dict)
        
        bV_payg_vfi = np.zeros(cS_python.aD_new)
        if cS_python.aR_new < cS_python.aD_new:
            bV_payg_vfi[cS_python.aR_new:] = self.M_fixed['b_payg_avg_retiree']
        bV_payg_matlab = matlab.double(bV_payg_vfi.tolist())
        
        cPolM, kPolM, cPpsM, _ = self.eng.main_olg_v8_utils.HHSolution_VFI_Huggett(
            self.M_fixed['R_k_net_factor'], self.M_fixed['w_gross'], self.M_fixed['TR_total'],
            bV_payg_matlab, paramS_vfi_matlab, cS_matlab_dict, vfi_config['solver'], nargout=4)
        
        end_time = time.time()
        print(f"    ✅ VFI求解完成，耗时: {end_time - start_time:.2f} 秒")
        
        return {
            'cPolM': np.array(cPolM),
            'kPolM': np.array(kPolM), 
            'cPpsPolM_choice': np.array(cPpsM),
            'M_test': self.M_fixed,  # <<<<< 加上这一行
            'cS_python_obj': cS_python,
            'paramS_vfi_dict': paramS_vfi_dict
        }
    def run_rule_of_thumb_simulation(self, cS_obj: Any, paramS_dict: Dict, n_sim: int, eIdxM_group: np.ndarray) -> Dict:
        """朴素策略（固定储蓄率）模拟器。"""
        print("\n--- 模拟朴素策略 (Rule of Thumb) ---")
        SAVING_RATE, PPS_ALLOCATION_RATE = 0.20, 0.25
        leGridV = paramS_dict['leGridV']
        kMin, kppsMin, kMax, kppsMax = cS_obj.kMin, cS_obj.kppsMin, cS_obj.kMax, cS_obj.kppsMax
        aR_new, aD_new, tau_c, cFloor, pps_active = cS_obj.aR_new, cS_obj.aD_new, cS_obj.tau_c, cS_obj.cFloor, cS_obj.pps_active
        k_paths, kpps_paths, c_paths, cpps_paths = [], [], [], []
        
        paramS_hh = TempParamSHH(self.M_fixed['tau_l'], self.M_fixed['theta_payg_actual'], pps_active, cS_obj.ageEffV_new)
        bV_payg = np.zeros(aD_new)
        if aR_new < aD_new: bV_payg[aR_new:] = self.M_fixed['b_payg_avg_retiree']

        for i_sim in range(n_sim):
            k, kpps = kMin, kppsMin
            k_path, kpps_path, c_path, cpps_path = [], [], [], []
            for age_idx in range(aD_new):
                k_path.append(k); kpps_path.append(kpps)
                epsilon_val = leGridV[eIdxM_group[i_sim, age_idx] - 1]
                
                if age_idx < aR_new: # 工作期
                    income, _, _ = OLG_V9_Utils.HHIncome_Huggett(k, self.M_fixed['R_k_net_factor'], self.M_fixed['w_gross'], self.M_fixed['TR_total'], bV_payg[age_idx], 0.0, age_idx, paramS_hh, cS_obj, epsilon_val)
                    total_savings = SAVING_RATE * income
                    c_pps_choice = PPS_ALLOCATION_RATE * total_savings if pps_active else 0.0
                    k_prime = total_savings - c_pps_choice
                    c_expend = income - total_savings
                else: # 退休期
                    c_pps_choice = 0.0
                    paramS_hh.current_pps_withdrawal = kpps * cS_obj.pps_withdrawal_rate if pps_active else 0
                    income, _, _ = OLG_V9_Utils.HHIncome_Huggett(k, self.M_fixed['R_k_net_factor'], self.M_fixed['w_gross'], self.M_fixed['TR_total'], bV_payg[age_idx], 0.0, age_idx, paramS_hh, cS_obj, epsilon_val)
                    total_wealth = income + k # income already includes capital income
                    remaining_periods = aD_new - age_idx
                    c_expend = total_wealth / remaining_periods if remaining_periods > 0 else total_wealth
                    k_prime = 0.0

                c_val = max(cFloor, c_expend / (1 + tau_c))
                c_path.append(c_val); cpps_path.append(max(0, c_pps_choice))
                
                k = max(kMin, min(kMax, k_prime))
                if pps_active:
                    pps_return_factor = 1 + ((self.M_fixed['R_k_net_factor'] - 1) + cS_obj.pps_return_rate_premium)
                    kpps = max(kppsMin, min(kppsMax, (kpps + c_pps_choice - paramS_hh.current_pps_withdrawal) * pps_return_factor))
                else:
                    kpps = kppsMin
            
            k_paths.append(k_path); kpps_paths.append(kpps_path); c_paths.append(c_path); cpps_paths.append(cpps_path)
        
        print("    ✅ 朴素策略模拟完成。")
        return {"k_path": np.array(k_paths), "kpps_path": np.array(kpps_paths), "c_path": np.array(c_paths), "cpps_path": np.array(cpps_paths)}

    def _convert_dict_to_matlab_struct(self, py_dict: Dict) -> Dict:
        matlab_struct = {}
        for key, value in py_dict.items():
            if key in ['physAgeMap', 'interpolation_method', 'initial_pop', 'use_continuous_optimization']: continue
            if isinstance(value, np.ndarray): matlab_struct[key] = matlab.double(value.tolist())
            elif isinstance(value, list): matlab_struct[key] = matlab.double(value)
            elif isinstance(value, (int, float, bool)): matlab_struct[key] = float(value)
        return matlab_struct
        
    def _calculate_lifetime_utility(self, c_path: np.ndarray, cS: Any, use_survival_prob: bool) -> float:
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

    def run_comparison(self, n_sim=500, use_survival_prob_in_eval=True):
        """主比较流程：在一个脚本中完成所有策略的比较和分析。"""
        print("\n" + "="*80)
        print(f"🔬 开始全面的策略比较 (模拟个体数: {n_sim})")
        print("=" * 80)

        # 1. 加载RL模型，并将其配置确立为所有比较的唯一基准
        rl_model, rl_config = self.load_rl_model(use_best_model=True)
        baseline_cS = rl_config['cS']
        # baseline_cS.pps_return_rate_premium = -0.1 

        # 2. 生成统一的随机冲击“剧本”，使用基准cS确保随机过程一致
        _, tr_prob, p0 = OLG_V9_Utils.EarningProcess_olgm(baseline_cS)
        eIdxM_group_0based = OLG_V9_Utils.MarkovChainSimulation_AgeGroup(n_sim, baseline_cS, p0, tr_prob)
        eIdxM_group_for_sim = eIdxM_group_0based + 1
        
        all_results = []
        
        # 3. 循环遍历每个VFI配置进行评估
        for vfi_config in self.vfi_grid_configs:
            vfi_results = self.run_matlab_vfi_for_config(vfi_config, baseline_cS)
            sim_paths = OLG_V9_Utils.HHSimulation_olgm_vfi_simplified(vfi_results, n_sim, eIdxM_group_for_sim)
            utility = np.array([self._calculate_lifetime_utility(sim_paths['c_path_vfi'][j,:], baseline_cS, use_survival_prob_in_eval) for j in range(n_sim)])
            all_results.append({
                'label': vfi_config['label'], 'mean_utility': np.mean(utility),
                'paths': {'k': sim_paths['k_path_vfi'], 'c': sim_paths['c_path_vfi'], 'cpps': sim_paths['cpps_path_vfi']}
            })
            print(f"    📊 {vfi_config['label']} 评估结果: 平均终身效用 = {np.mean(utility):.4f}")

        # 4. 评估RL智能体
        print("\n--- 模拟RL智能体 (基准) ---")
        paramS_sim_for_rl = TempParamSHH(self.M_fixed['tau_l'], self.M_fixed['theta_payg_actual'], baseline_cS.pps_active, baseline_cS.ageEffV_new)
        bV_payg_for_rl = np.zeros(baseline_cS.aD_new)
        if baseline_cS.aR_new < baseline_cS.aD_new:
            bV_payg_for_rl[baseline_cS.aR_new:] = self.M_fixed['b_payg_avg_retiree']
        k_rl, _, c_rl, cpps_rl = OLG_V9_Utils.HHSimulation_olgm_rl(
            rl_model, rl_config, eIdxM_group_for_sim, self.M_fixed['R_k_net_factor'], 
            self.M_fixed['w_gross'], self.M_fixed['TR_total'], bV_payg_for_rl, paramS_sim_for_rl, baseline_cS
        )
        utility_rl = np.array([self._calculate_lifetime_utility(c_rl[j,:], baseline_cS, use_survival_prob_in_eval) for j in range(n_sim)])
        all_results.append({
            'label': 'RL (全功能)', 'mean_utility': np.mean(utility_rl),
            'paths': {'k': k_rl, 'c': c_rl, 'cpps': cpps_rl}
        })
        print(f"    📊 RL 评估结果: 平均终身效用 = {np.mean(utility_rl):.4f}")

        # 5. 评估朴素策略
        paramS_rot = {'leGridV': np.exp(OLG_V9_Utils.EarningProcess_olgm(baseline_cS)[0])}
        rot_paths = self.run_rule_of_thumb_simulation(baseline_cS, paramS_rot, n_sim, eIdxM_group_for_sim)
        utility_rot = np.array([self._calculate_lifetime_utility(rot_paths['c_path'][j,:], baseline_cS, use_survival_prob_in_eval) for j in range(n_sim)])
        all_results.append({
            'label': '朴素策略', 'mean_utility': np.mean(utility_rot),
            'paths': {'k': rot_paths['k_path'], 'c': rot_paths['c_path'], 'cpps': rot_paths['cpps_path']}
        })
        print(f"    📊 朴素策略 评估结果: 平均终身效用 = {np.mean(utility_rot):.4f}")
        
        # 6. 分析并绘图
        self.analyze_and_plot(all_results, baseline_cS)

# 在 ComprehensiveStrategyComparator 类中，替换此函数

    def analyze_and_plot(self, all_results: List[Dict], cS: Any):
        """分析并绘制所有策略的性能和行为路径。"""
        print("\n📈 分析与绘图 (综合比较)...")
        
        all_results.sort(key=lambda x: x['mean_utility'], reverse=True)
        
        # --- Part 1: 性能排序文本输出 ---
        print("\n" + "=" * 80 + "\n📋 策略性能排序\n" + "=" * 80)
        for i, res in enumerate(all_results):
            print(f"  {i+1}. {res['label']:<25}: {res['mean_utility']:.4f}")
        print("=" * 80)

        # --- Part 2: 绘图 (2x2 布局) ---
        fig, axes = plt.subplots(2, 2, figsize=(20, 16))
        fig.suptitle('固定宏观环境下的全面策略比较', fontsize=20, y=0.97)
        axes = axes.flatten()
        
        # [核心修复] 准备颜色和样式映射，使用简化键
        color_map = {'RL': 'red', '朴素策略': 'green'}
        style_map = {'RL': 's--', '朴素策略': '^:'}
        vfi_colors = plt.cm.Blues(np.linspace(0.9, 0.4, len(self.vfi_grid_configs) if self.vfi_grid_configs else 1))
        vfi_styles = ['o-', 'd-', 'v-', 'p-']
        
        # --- 图A: 平均终身效用条形图 ---
        labels = [res['label'] for res in all_results]
        utilities = [res['mean_utility'] for res in all_results]
        bar_colors = []
        vfi_idx = 0
        for label in labels:
            if 'RL' in label: bar_colors.append(color_map['RL'])
            elif '朴素' in label: bar_colors.append(color_map['朴素策略'])
            else:
                bar_colors.append(vfi_colors[vfi_idx % len(vfi_colors)])
                vfi_idx += 1
        
        bars = axes[0].bar(labels, utilities, color=bar_colors)
        axes[0].set_ylabel('平均终身效用', fontsize=14)
        axes[0].set_title('A. 最终性能排序', fontsize=16)
        axes[0].set_xticklabels(labels, rotation=35, ha='right', fontsize=12)
        axes[0].grid(True, axis='y', linestyle='--', alpha=0.7)
        for bar in bars:
            yval = bar.get_height()
            axes[0].text(bar.get_x() + bar.get_width()/2.0, yval, f'{yval:.3f}', va='bottom' if yval > 0 else 'top', ha='center', fontsize=10)

        # --- 准备路径数据 ---
        age_groups = np.arange(cS.aD_new)
        mean_paths = {}
        for res in all_results:
            mean_paths[res['label']] = {
                'k': np.mean(res['paths']['k'], axis=0),
                'c': np.mean(res['paths']['c'], axis=0),
                'cpps': np.mean(res['paths']['cpps'], axis=0)
            }
        
        # [核心修复] 循环绘图的逻辑
        def plot_paths(ax, path_key, ylabel, title):
            vfi_idx = 0
            for label, m_paths in mean_paths.items():
                path_to_plot = m_paths[path_key]
                if 'RL' in label:
                    ax.plot(age_groups, path_to_plot, style_map['RL'], label=label, color=color_map['RL'], lw=2.5)
                elif '朴素' in label:
                    ax.plot(age_groups, path_to_plot, style_map['朴素策略'], label=label, color=color_map['朴素策略'], lw=2.5)
                else:
                    ax.plot(age_groups, path_to_plot, vfi_styles[vfi_idx % len(vfi_styles)], label=label, color=vfi_colors[vfi_idx % len(vfi_colors)], lw=2)
                    vfi_idx += 1
            ax.set_xlabel('年龄组', fontsize=12)
            ax.set_ylabel(ylabel, fontsize=14)
            ax.set_title(title, fontsize=16)
            ax.legend(fontsize=11)
            ax.grid(True, linestyle='--', alpha=0.6)

        # --- 图B: 平均消费路径 ---
        plot_paths(axes[1], 'c', '平均消费 (c)', 'B. 平均消费生命周期路径')

        # --- 图C: 平均资产路径 ---
        plot_paths(axes[2], 'k', '平均非PPS资产 (k)', 'C. 平均资产生命周期路径')
        
        # --- 图D: 平均PPS缴费路径 ---
        plot_paths(axes[3], 'cpps', '平均PPS缴费 (c_pps)', 'D. 平均PPS缴费生命周期路径')
        axes[3].axvline(x=cS.aR_new - 1, color='gray', linestyle='--', label='退休年龄') # aR_new是1-based长度, 退休发生在age_idx = aR_new-1 之后
        axes[3].legend(fontsize=11) # 重新调用以显示退休线
        
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        save_path = './py/comprehensive_strategy_comparison.png'
        plt.savefig(save_path, dpi=300)
        print(f"\n📈 综合比较图表已保存到: {save_path}")
        plt.show()


def main():
    if not MATLAB_AVAILABLE:
        return
    comparator = ComprehensiveStrategyComparator()
    comparator.run_comparison(n_sim=1000, use_survival_prob_in_eval=True)

if __name__ == "__main__":
    main()