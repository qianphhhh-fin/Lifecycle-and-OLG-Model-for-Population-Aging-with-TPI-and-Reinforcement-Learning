# --- 开始文件：compare_rl_and_vfi_matlab_simplified.py (最终独立绘图版 - 已加入朴素策略) ---

"""
比较VFI、简化版RL和朴素策略的优化结果 - 固定宏观环境

[最终独立版 - 已加入朴素策略]
- 新增了一个基于“拇指法则”的朴素策略作为第三个比较基准。
- 朴素策略：工作期固定比例储蓄，退休期年金化消费。
- 所有三个策略的模拟都使用完全相同的宏观环境和随机冲击路径，确保公平。
- 分析和绘图功能已更新，以展示三个策略的对比。
"""

import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pickle
import json
import time
from typing import Dict, Any, Tuple
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
from main_olg_v9_utils import OLG_V9_Utils, TempParamSHH # 导入 TempParamSHH

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

class RLVFIComparatorSimplified:
    """[最终独立版] 简化版RL、VFI和朴素策略比较器。"""
    
    def __init__(self, use_sbx: bool = True):
        self.M_FIXED = {
            'R_k_net_factor': 1.03, 'w_gross': 2.0, 'TR_total': 0.1,
            'b_payg_avg_retiree': 0.4, 'tau_l': 0.15, 'theta_payg_actual': 0.12
        }
        self.eng = matlab.engine.start_matlab()
        self.eng.addpath(os.getcwd(), nargout=0)
        print("✅ MATLAB Engine 已启动。")
    
    def __del__(self):
        if hasattr(self, 'eng') and self.eng is not None:
            self.eng.quit()
            print("✅ MATLAB Engine已关闭。")

    def load_rl_model(self, use_best_model: bool = True) -> Tuple[Any, Dict]:
        # ... (此函数无需修改) ...
        model_path = './py/best_model_sbx_simp/best_model.zip' if use_best_model else './py/final_sac_agent_olg_sbx_simp.zip'
        if not os.path.exists(model_path):
            print(f"⚠️ 未找到最佳模型 '{model_path}'，回退到最终模型...")
            model_path = './py/final_sac_agent_olg_sbx_simp.zip'
        config_path = model_path.replace('.zip', '_config.pkl')
        if "best_model" in model_path and not os.path.exists(config_path):
             config_path = './py/final_sac_agent_olg_sbx_simp_config.pkl'
        print(f"📁 正在加载简化版RL模型: {model_path}")
        model = SBX_SAC.load(model_path)
        with open(config_path, 'rb') as f: config = pickle.load(f)
        print("✅ 模型和配置加载成功。")
        return model, config

    def run_matlab_vfi(self) -> Dict[str, Any]:
        # ... (此函数无需修改) ...
        print("\n🔧 VFI求解 (在固定的宏观环境下)...")
        cS_python = OLG_V9_Utils.ParameterValues_HuggettStyle()
        (leLogGridV, leTrProbM, leProb1V) = OLG_V9_Utils.EarningProcess_olgm(cS_python)
        leGridV = np.exp(leLogGridV)
        paramS_vfi_dict = {
            'leLogGridV': leLogGridV, 'leTrProbM': leTrProbM, 'leProb1V': leProb1V,
            'leGridV': leGridV, 'ageEffV_new': cS_python.ageEffV_new,
            'tau_l': self.M_FIXED['tau_l'],
            'theta_payg_actual_for_hh': self.M_FIXED['theta_payg_actual'],
            'pps_tax_deferral_active': bool(cS_python.pps_active)
        }
        cS_matlab_dict = self._convert_dict_to_matlab_struct(cS_python.__dict__)
        paramS_vfi_matlab = self._convert_dict_to_matlab_struct(paramS_vfi_dict)
        bV_payg_vfi = np.zeros(cS_python.aD_new)
        if cS_python.aR_new < cS_python.aD_new:
            bV_payg_vfi[cS_python.aR_new:] = self.M_FIXED['b_payg_avg_retiree']
        bV_payg_matlab = matlab.double(bV_payg_vfi.tolist())
        cS_matlab_dict['interpolation_method'] = 'spline'
        cPolM, kPolM, cPpsM, _ = self.eng.main_olg_v8_utils.HHSolution_VFI_Huggett(
            self.M_FIXED['R_k_net_factor'], self.M_FIXED['w_gross'], self.M_FIXED['TR_total'],
            bV_payg_matlab, paramS_vfi_matlab, cS_matlab_dict, nargout=4)
        print("✅ VFI策略求解完成。")
        return {'cPolM_q': np.array(cPolM), 'kPolM': np.array(kPolM), 'cPpsPolM_choice': np.array(cPpsM),
                'M_test': self.M_FIXED, 'bV_payg_eq': bV_payg_vfi,
                'cS_python_obj': cS_python, 'paramS_vfi_dict': paramS_vfi_dict}
    
    def _convert_dict_to_matlab_struct(self, py_dict: Dict) -> Dict:
        # ... (此函数无需修改) ...
        matlab_struct = {}
        for key, value in py_dict.items():
            if key in ['physAgeMap', 'interpolation_method']: continue
            if isinstance(value, np.ndarray): matlab_struct[key] = matlab.double(value.tolist())
            elif isinstance(value, list) and all(isinstance(i, (int, float, np.number)) for i in value): matlab_struct[key] = matlab.double(value)
            elif isinstance(value, (int, float, bool, np.number)): matlab_struct[key] = float(value)
        return matlab_struct
        
    def _calculate_lifetime_utility(self, c_path: np.ndarray, cS: Any, use_survival_prob: bool) -> float:
        # ... (此函数无需修改) ...
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

    # [新] 朴素策略模拟器
    def run_rule_of_thumb_simulation(self, cS_obj: Any, paramS_dict: Dict, n_sim: int, 
                                     eIdxM_group: np.ndarray) -> Dict:
        """
        在固定的宏观环境下，模拟一个基于“拇指法则”的朴素策略。
        
        Args:
            cS_obj: 模型参数对象。
            paramS_dict: 包含leGridV, ageEffV_new等的参数字典。
            n_sim: 模拟个体数。
            eIdxM_group: 预先生成的效率冲击路径 (n_sim x aD_new)，1-based索引。

        Returns:
            一个包含所有朴素策略模拟路径的字典。
        """
        print("👍 运行朴素策略模拟 (Rule-of-Thumb)...")

        # 1. 定义朴素策略的参数
        SAVING_RATE = 0.20  # 工作期储蓄率
        PPS_ALLOCATION_RATE = 0.25  # 储蓄中用于PPS的比例

        # 2. 提取必要的参数
        leGridV = paramS_dict['leGridV']
        ageEffV_new = paramS_dict['ageEffV_new']
        kMin, kppsMin = cS_obj.kMin, cS_obj.kppsMin
        kMax, kppsMax = cS_obj.kMax, cS_obj.kppsMax
        aR_new = cS_obj.aR_new
        aD_new = cS_obj.aD_new
        tau_c = cS_obj.tau_c
        cFloor = cS_obj.cFloor
        pps_active = cS_obj.pps_active
        
        # 3. 初始化结果存储
        k_paths, kpps_paths, c_paths, cpps_paths = [], [], [], []

        # 创建一个临时的paramS对象，用于调用HHIncome_Huggett
        paramS_hh = TempParamSHH(
            self.M_FIXED['tau_l'],
            self.M_FIXED['theta_payg_actual'],
            pps_active,
            ageEffV_new
        )
        
        bV_payg = np.zeros(aD_new)
        if aR_new < aD_new:
            bV_payg[aR_new:] = self.M_FIXED['b_payg_avg_retiree']

        # 4. 循环模拟每个个体的生命周期
        for i_sim in range(n_sim):
            k_current, kpps_current = kMin, kppsMin
            k_path, kpps_path, c_path, cpps_path = [], [], [], []

            for age_idx in range(aD_new): # 0-based
                k_path.append(k_current)
                kpps_path.append(kpps_current)
                
                eps_idx_1based = eIdxM_group[i_sim, age_idx]
                epsilon_val = leGridV[eps_idx_1based - 1]
                
                # 更新paramS_hh中的PPS提取部分
                if age_idx >= aR_new and pps_active:
                    paramS_hh.current_pps_withdrawal = kpps_current * cS_obj.pps_withdrawal_rate
                else:
                    paramS_hh.current_pps_withdrawal = 0

                if age_idx < aR_new: # 工作期
                    # 确定收入（此时还不知道要交多少PPS，所以先假设为0）
                    income, _, _ = OLG_V9_Utils.HHIncome_Huggett(
                        k_current, self.M_FIXED['R_k_net_factor'], self.M_FIXED['w_gross'],
                        self.M_FIXED['TR_total'], bV_payg[age_idx], 0.0,
                        age_idx, paramS_hh, cS_obj, epsilon_val
                    )
                    
                    total_savings = SAVING_RATE * income
                    
                    # 分配储蓄到PPS和非PPS
                    c_pps = PPS_ALLOCATION_RATE * total_savings if pps_active else 0
                    k_prime = total_savings - c_pps
                    
                    consumption_expenditure = income - total_savings
                    
                else: # 退休期
                    c_pps = 0 # 退休后不缴费
                    
                    # 退休期策略：年金化消费
                    # 计算总财富
                    total_wealth = (k_current + kpps_current * (1 - cS_obj.pps_tax_rate_withdrawal))
                    # 计算剩余生命年数
                    remaining_periods = aD_new - age_idx
                    # 年金化消费支出
                    consumption_expenditure = total_wealth / remaining_periods if remaining_periods > 0 else total_wealth
                    # 退休后消费所有财富，不再储蓄
                    k_prime = 0
                
                # 计算最终消费
                current_c = max(cFloor, consumption_expenditure / (1 + tau_c))
                c_path.append(current_c)
                cpps_path.append(max(0, c_pps))
                
                # 更新状态
                k_current = max(kMin, min(kMax, k_prime))
                if pps_active:
                    pps_return_factor = 1 + ((self.M_FIXED['R_k_net_factor'] - 1) + cS_obj.pps_return_rate_premium)
                    pps_withdrawal = paramS_hh.current_pps_withdrawal
                    next_k_pps = (kpps_current + c_pps - pps_withdrawal) * pps_return_factor
                    kpps_current = max(kppsMin, min(kppsMax, next_k_pps))
                else:
                    kpps_current = kppsMin

            k_paths.append(k_path)
            kpps_paths.append(kpps_path)
            c_paths.append(c_path)
            cpps_paths.append(cpps_path)
        
        print("✅ 朴素策略模拟完成。")
        return {
            "k_path_rot": np.array(k_paths),
            "c_path_rot": np.array(c_paths),
            "cpps_path_rot": np.array(cpps_paths),
        }

    def run_comparison(self, n_sim=500, random_seed=42, use_survival_prob_in_eval=True):
        """[修改] 主比较流程，加入朴素策略"""
        print("\n" + "="*80)
        print("🔬 VFI vs RL vs 朴素策略 优化方法比较 [最终独立版]")
        print("=" * 80)

        # 1. 加载/求解各种策略
        rl_model, rl_config = self.load_rl_model(use_best_model=True)
        vfi_results = self.run_matlab_vfi()
        
        cS_obj = vfi_results['cS_python_obj']
        paramS_dict = vfi_results['paramS_vfi_dict']
        
        print("\n--- 开始使用统一模拟器进行生命周期模拟 ---")
        
        # 2. 生成统一的随机冲击“剧本”
        eIdxM_group_0based = OLG_V9_Utils.MarkovChainSimulation_AgeGroup(n_sim, cS_obj, paramS_dict['leProb1V'], paramS_dict['leTrProbM'])
        eIdxM_group_for_sim = eIdxM_group_0based + 1
        
        # 3. 运行所有模拟
        vfi_sim_paths = OLG_V9_Utils.HHSimulation_olgm_vfi_simplified(vfi_results, n_sim, random_seed, eIdxM_group_for_sim)
        
        rl_sim_paths = OLG_V9_Utils.HHSimulation_olgm_rl_simplified(
            rl_model, cS_obj, rl_config['paramS_for_rl'], self.M_FIXED, n_sim, eIdxM_group_for_sim
        )
        
        rot_sim_paths = self.run_rule_of_thumb_simulation(
            cS_obj, paramS_dict, n_sim, eIdxM_group_for_sim
        )
        
        # 4. 计算所有策略的生涯效用
        print(f"\n📊 计算生涯效用 ({'包含' if use_survival_prob_in_eval else '不含'}生存概率)...")
        utility_vfi = np.array([self._calculate_lifetime_utility(vfi_sim_paths['c_path_vfi'][i,:], cS_obj, use_survival_prob_in_eval) for i in range(n_sim)])
        utility_rl = np.array([self._calculate_lifetime_utility(rl_sim_paths['c_path_rl'][i,:], cS_obj, use_survival_prob_in_eval) for i in range(n_sim)])
        utility_rot = np.array([self._calculate_lifetime_utility(rot_sim_paths['c_path_rot'][i,:], cS_obj, use_survival_prob_in_eval) for i in range(n_sim)])
        
        results_for_plot = {
            **vfi_sim_paths, **rl_sim_paths, **rot_sim_paths,
            'utility_vfi': utility_vfi, 'utility_rl': utility_rl, 'utility_rot': utility_rot,
            'n_sim': n_sim
        }
        self.analyze_and_plot(results_for_plot)

    def analyze_and_plot(self, results):
        """[修改] 分析并绘制三个策略的比较结果"""
        print("\n📈 分析与绘图...")
        mean_utility_vfi = np.mean(results['utility_vfi'])
        mean_utility_rl = np.mean(results['utility_rl'])
        mean_utility_rot = np.mean(results['utility_rot'])
        std_utility_vfi = np.std(results['utility_vfi'])
        std_utility_rl = np.std(results['utility_rl'])
        std_utility_rot = np.std(results['utility_rot'])
        
        # T-tests
        t_rl_vfi, p_rl_vfi = stats.ttest_ind(results['utility_rl'], results['utility_vfi'], equal_var=False)
        t_rl_rot, p_rl_rot = stats.ttest_ind(results['utility_rl'], results['utility_rot'], equal_var=False)

        print("\n" + "=" * 80)
        print("📋 VFI vs RL vs 朴素策略 比较分析摘要")
        print("=" * 80)
        print(f" - VFI 平均效用:           {mean_utility_vfi:.4f} ± {std_utility_vfi:.4f}")
        print(f" - RL (专家) 平均效用:     {mean_utility_rl:.4f} ± {std_utility_rl:.4f}")
        print(f" - 朴素策略 平均效用:    {mean_utility_rot:.4f} ± {std_utility_rot:.4f}")
        print("-" * 40)
        print(f" - RL vs VFI:  效用差异 {mean_utility_rl - mean_utility_vfi:+.4f}, p-value={p_rl_vfi:.4f} ({'显著' if p_rl_vfi < 0.05 else '不显著'})")
        print(f" - RL vs 朴素: 效用差异 {mean_utility_rl - mean_utility_rot:+.4f}, p-value={p_rl_rot:.4f} ({'显著' if p_rl_rot < 0.05 else '不显著'})")

        self.plot_comparison_results(results, save_path='./py/rl_vfi_rot_comparison.png')

    def plot_comparison_results(self, results, save_path):
        """[修改] 生成6个核心对比图表，包含三个策略"""
        fig, axes = plt.subplots(2, 3, figsize=(20, 12))
        fig.suptitle(f'VFI vs RL vs 朴素策略 在固定环境下的比较 (n={results["n_sim"]})', fontsize=18)
        
        # 1. 效用分布比较
        ax = axes[0, 0]
        ax.hist(results['utility_vfi'], bins=30, alpha=0.6, label='VFI', color='blue', density=True)
        ax.hist(results['utility_rl'], bins=30, alpha=0.6, label='RL', color='red', density=True)
        ax.hist(results['utility_rot'], bins=30, alpha=0.6, label='朴素策略', color='green', density=True)
        ax.axvline(np.mean(results['utility_vfi']), color='navy', linestyle='--', label=f'VFI均值: {np.mean(results["utility_vfi"]):.2f}')
        ax.axvline(np.mean(results['utility_rl']), color='darkred', linestyle='--', label=f'RL均值: {np.mean(results["utility_rl"]):.2f}')
        ax.axvline(np.mean(results['utility_rot']), color='darkgreen', linestyle='--', label=f'朴素策略均值: {np.mean(results["utility_rot"]):.2f}')
        ax.set_xlabel('生涯总效用', fontsize=12)
        ax.set_ylabel('密度', fontsize=12)
        ax.set_title('效用分布比较', fontsize=14)
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # 2. 平均资产路径
        ax = axes[0, 1]
        age_path = np.arange(1, results['k_path_vfi'].shape[1] + 1)
        ax.plot(age_path, np.mean(results['k_path_vfi'], axis=0), 'b-', label='VFI', lw=2)
        ax.plot(age_path, np.mean(results['k_path_rl'], axis=0), 'r--', label='RL', lw=2)
        ax.plot(age_path, np.mean(results['k_path_rot'], axis=0), 'g-.', label='朴素策略', lw=2)
        ax.set_xlabel('年龄组索引', fontsize=12)
        ax.set_ylabel('平均资产', fontsize=12)
        ax.set_title('平均资产路径', fontsize=14)
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # 3. 平均消费路径
        ax = axes[0, 2]
        ax.plot(age_path, np.mean(results['c_path_vfi'], axis=0), 'b-', label='VFI', lw=2)
        ax.plot(age_path, np.mean(results['c_path_rl'], axis=0), 'r--', label='RL', lw=2)
        ax.plot(age_path, np.mean(results['c_path_rot'], axis=0), 'g-.', label='朴素策略', lw=2)
        ax.set_xlabel('年龄组索引', fontsize=12)
        ax.set_ylabel('平均消费', fontsize=12)
        ax.set_title('平均消费路径', fontsize=14)
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # 4. RL vs VFI 个体效用对比
        ax = axes[1, 0]
        ax.scatter(results['utility_vfi'], results['utility_rl'], alpha=0.4, s=20, label='个体')
        min_val = min(np.min(results['utility_vfi']), np.min(results['utility_rl']))
        max_val = max(np.max(results['utility_vfi']), np.max(results['utility_rl']))
        ax.plot([min_val, max_val], [min_val, max_val], 'r--', label='45度线')
        ax.set_xlabel('VFI 生涯效用', fontsize=12)
        ax.set_ylabel('RL 生涯效用', fontsize=12)
        ax.set_title('个体效用对比: RL vs VFI', fontsize=14)
        ax.legend()
        ax.grid(True, alpha=0.3)

        # 5. RL vs 朴素策略 个体效用对比
        ax = axes[1, 1]
        ax.scatter(results['utility_rot'], results['utility_rl'], alpha=0.4, s=20, color='green', label='个体')
        min_val = min(np.min(results['utility_rot']), np.min(results['utility_rl']))
        max_val = max(np.max(results['utility_rot']), np.max(results['utility_rl']))
        ax.plot([min_val, max_val], [min_val, max_val], 'r--', label='45度线')
        ax.set_xlabel('朴素策略 生涯效用', fontsize=12)
        ax.set_ylabel('RL 生涯效用', fontsize=12)
        ax.set_title('个体效用对比: RL vs 朴素策略', fontsize=14)
        ax.legend()
        ax.grid(True, alpha=0.3)

        # 6. 平均PPS缴费路径
        ax = axes[1, 2]
        ax.plot(age_path, np.mean(results['cpps_path_vfi'], axis=0), 'b-', label='VFI', lw=2)
        ax.plot(age_path, np.mean(results['cpps_path_rl'], axis=0), 'r--', label='RL', lw=2)
        ax.plot(age_path, np.mean(results['cpps_path_rot'], axis=0), 'g-.', label='朴素策略', lw=2)
        ax.set_xlabel('年龄组索引', fontsize=12)
        ax.set_ylabel('平均PPS缴费', fontsize=12)
        ax.set_title('平均PPS缴费路径', fontsize=14)
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plt.savefig(save_path, dpi=300)
        print(f"📈 三策略比较图表已保存到: {save_path}")
        plt.show()

def main():
    if not MATLAB_AVAILABLE: return
    comparator = RLVFIComparatorSimplified()
    comparator.run_comparison(use_survival_prob_in_eval=True)

if __name__ == "__main__":
    main()