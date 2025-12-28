# --- 开始文件：test_vfi_grid_sensitivity.py (与 compare... 模拟逻辑完全一致版) ---

"""
VFI 网格密度敏感性测试脚本

[最终一致版]
目标：
- 研究 nk 和 nkpps 的密度对VFI策略性能的影响。
- [核心] 所有策略评估都统一使用在 compare_rl_and_vfi... 脚本中经过验证的
  HHSimulation_olgm_vfi_simplified 模拟器，以确保评估标准完全一致。

方法：
1. 定义要测试的 (nk, nkpps) 组合。
2. 对每种组合，调用MATLAB求解VFI策略。
3. 使用统一的、经过验证的 HHSimulation_olgm_vfi_simplified 函数模拟每种策略。
4. 可视化结果。
"""

import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import time
from typing import Dict, List, Tuple, Any

# MATLAB Engine导入
try:
    import matlab.engine
    MATLAB_AVAILABLE = True
except ImportError:
    MATLAB_AVAILABLE = False
    print("❌ MATLAB Engine 不可用，无法继续。")
    exit()

# 导入自定义模块
from main_olg_v9_utils import OLG_V9_Utils # 确保这个模块包含 HHSimulation_olgm_vfi_simplified

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

class VFIGridSensitivityTester:
    """
    一个专门用于测试VFI网格密度敏感性的类。
    """
    
    def __init__(self):
        self.M_FIXED = {
            'R_k_net_factor': 1.03, 'w_gross': 2.0, 'TR_total': 0.1,
            'b_payg_avg_retiree': 0.4, 'tau_l': 0.15, 'theta_payg_actual': 0.12
        }
        self.eng = matlab.engine.start_matlab()
        self.eng.addpath(os.getcwd(), nargout=0)
        print("✅ MATLAB Engine 已启动。")
    
    def __del__(self):
        if hasattr(self, 'eng') and self.eng is not None: self.eng.quit()

    def run_vfi_for_grid(self, nk: int, nkpps: int, interpolation_method: str = 'linear') -> Dict:
        """[修正] 为指定的网格密度运行VFI求解器，并返回与 `compare...` 脚本兼容的字典"""
        print(f"\n--- 求解VFI for (nk={nk}, nkpps={nkpps}, interp='{interpolation_method}') ---")
        
        # 1. 初始化并修改参数
        cS_python = OLG_V9_Utils.ParameterValues_HuggettStyle()
        cS_python.nk = nk
        cS_python.nkpps = nkpps
        
        # 2. 重新生成依赖于网格的参数
        power = 1.5
        cS_python.kGridV = cS_python.kMin + (cS_python.kMax - cS_python.kMin) * (np.linspace(0, 1, nk)**power)
        if nk > 0: cS_python.kGridV[0] = cS_python.kMin
        cS_python.kGridV = cS_python.kGridV.flatten()

        power_kpps = 1.5
        if nkpps > 1:
            kppsGridV_temp = cS_python.kppsMin + (cS_python.kppsMax - cS_python.kppsMin) * (np.linspace(0, 1, nkpps)**power_kpps)
            kppsGridV_temp[0] = cS_python.kppsMin
        elif nkpps == 1: kppsGridV_temp = np.array([cS_python.kppsMin])
        else: kppsGridV_temp = np.array([])
        cS_python.kppsGridV = kppsGridV_temp.flatten()

        (leLogGridV, leTrProbM, leProb1V) = OLG_V9_Utils.EarningProcess_olgm(cS_python)
        
        # [新] 添加 ageEffV_new 到 paramS_vfi_dict 中，以供模拟器使用
        paramS_vfi_dict = {
            'leLogGridV': leLogGridV, # 添加这个以备用
            'leGridV': np.exp(leLogGridV), 'leTrProbM': leTrProbM, 'leProb1V': leProb1V,
            'ageEffV_new': cS_python.ageEffV_new, # 添加 ageEffV_new
            'tau_l': self.M_FIXED['tau_l'],
            'theta_payg_actual_for_hh': self.M_FIXED['theta_payg_actual'],
            'pps_tax_deferral_active': bool(cS_python.pps_active)
        }
        
        # 3. 转换并调用MATLAB
        cS_python.interpolation_method = interpolation_method
        
        cS_matlab_dict = self._convert_dict_to_matlab_struct(cS_python.__dict__)
        paramS_vfi_matlab = self._convert_dict_to_matlab_struct(paramS_vfi_dict)
        
        bV_payg_vfi = np.zeros(cS_python.aD_new)
        if cS_python.aR_new <= cS_python.aD_new:
            bV_payg_vfi[cS_python.aR_new-1:] = self.M_FIXED['b_payg_avg_retiree']
        bV_payg_matlab = matlab.double(bV_payg_vfi.tolist())

        cPolM, kPolM, cPpsM, _ = self.eng.main_olg_v8_utils.HHSolution_VFI_Huggett(
            self.M_FIXED['R_k_net_factor'], self.M_FIXED['w_gross'], self.M_FIXED['TR_total'],
            bV_payg_matlab, paramS_vfi_matlab, cS_matlab_dict, nargout=4
        )
        
        print(f"✅ VFI求解完成 for (nk={nk}, nkpps={nkpps})")

        # [修正] 返回与 HHSimulation_olgm_vfi_simplified 兼容的字典
        return {
            'cPolM_q': np.array(cPolM), 
            'kPolM': np.array(kPolM), 
            'cPpsPolM_choice': np.array(cPpsM),
            'M_test': self.M_FIXED,
            'bV_payg_eq': bV_payg_vfi,
            'cS_python_obj': cS_python, 
            'paramS_vfi_dict': paramS_vfi_dict,
            'nk': nk, 'nkpps': nkpps # 保留网格信息
        }

    def _convert_dict_to_matlab_struct(self, py_dict: Dict) -> Dict:
        matlab_struct = {}
        for key, value in py_dict.items():
            if key in ['physAgeMap', 'interpolation_method']: continue
            if isinstance(value, np.ndarray):
                matlab_struct[key] = matlab.double(value.tolist())
            elif isinstance(value, list) and all(isinstance(i, (int, float, np.number)) for i in value):
                matlab_struct[key] = matlab.double(value)
            elif isinstance(value, (int, float, bool, np.number)):
                matlab_struct[key] = float(value)
        return matlab_struct

    def run_sensitivity_analysis(self, grid_configs: List[Tuple[int, int]], n_sim=500, random_seed=42):
        """对一系列网格配置进行敏感性分析"""
        all_results = []
        
        print("\n--- 生成统一的效率冲击路径用于所有评估 ---")
        temp_cs = OLG_V9_Utils.ParameterValues_HuggettStyle()
        (leLogGridV, temp_tr_prob, temp_p0) = OLG_V9_Utils.EarningProcess_olgm(temp_cs)
        
        # [修正] 生成0-based索引的路径，然后转换为1-based给模拟器
        eIdxM_group_0based = OLG_V9_Utils.MarkovChainSimulation_AgeGroup(n_sim, temp_cs, temp_p0, temp_tr_prob)
        eIdxM_group_global_1based = eIdxM_group_0based + 1
        print(f"✅ 全局效率路径已生成 (shape: {eIdxM_group_global_1based.shape})。")

        for nk, nkpps in grid_configs:
            for interp_method in ['linear', 'spline']:
                # 1. 为当前网格配置求解VFI
                vfi_results = self.run_vfi_for_grid(nk, nkpps, interpolation_method=interp_method)
                
                cS_obj = vfi_results['cS_python_obj']
                
                # 2. [核心修正] 使用与 compare... 脚本完全一致的模拟器进行评估
                #    这个函数需要一个特定的 vfi_results 字典和 eIdxM_group (1-based)
                sim_paths = OLG_V9_Utils.HHSimulation_olgm_vfi_simplified(
                    vfi_results,
                    n_sim,
                    eIdxM_group_global_1based
                )
                
                # 3. 计算终身效用
                utility_vfi = np.array([
                    self._calculate_lifetime_utility(sim_paths['c_path_vfi'][i,:], cS_obj, True)
                    for i in range(n_sim)
                ])
                
                # 4. 存储结果
                all_results.append({
                    'nk': nk, 'nkpps': nkpps, 'interp_method': interp_method,
                    'mean_utility': np.mean(utility_vfi), 'std_utility': np.std(utility_vfi),
                    'sim_paths': sim_paths
                })
                print(f"📈 结果 for (nk={nk}, nkpps={nkpps}, interp='{interp_method}'): 平均效用 = {np.mean(utility_vfi):.4f}")

        self.plot_sensitivity_results(all_results)
        return all_results

    def _calculate_lifetime_utility(self, c_path: np.ndarray, cS: Any, use_survival_prob: bool) -> float:
        beta, aD = cS.beta, c_path.shape[0]
        s_transitionV = cS.s_1yr_transitionV.flatten()
        utility_sum, cumulative_discount = 0.0, 1.0
        for a_group in range(aD):
            _, u = OLG_V9_Utils.CES_utility(c_path[a_group], cS.sigma, cS)
            utility_sum += cumulative_discount * u
            if a_group < aD - 1:
                survival_factor = s_transitionV[a_group] if use_survival_prob else 1.0
                cumulative_discount *= (beta * survival_factor)
        return utility_sum

    def plot_sensitivity_results(self, results: List[Dict]):
        print("\n--- 可视化网格敏感性结果 ---")
        
        linear_results = [r for r in results if r['interp_method'] == 'linear']
        spline_results = [r for r in results if r['interp_method'] == 'spline']
        
        fig, axes = plt.subplots(2, 2, figsize=(18, 12), constrained_layout=True)
        fig.suptitle('VFI性能对网格密度和插值方法的敏感性分析', fontsize=16)

        # 1. 效用 vs 网格点总数
        ax = axes[0, 0]
        if linear_results:
            total_points = [r['nk'] * r['nkpps'] for r in linear_results]
            mean_utilities = [r['mean_utility'] for r in linear_results]
            ax.plot(total_points, mean_utilities, 'o-', markersize=8, label='interp = linear')
        if spline_results:
            total_points = [r['nk'] * r['nkpps'] for r in spline_results]
            mean_utilities = [r['mean_utility'] for r in spline_results]
            ax.plot(total_points, mean_utilities, 's--', markersize=8, label='interp = spline')
        
        ax.set_xlabel('总状态点数 (nk * nkpps)')
        ax.set_ylabel('平均终身效用')
        ax.set_title('VFI性能 vs. 网格总点数')
        ax.legend()
        ax.grid(True)

        # 2. 平均消费路径对比 (linear)
        ax = axes[1, 0]
        colors = plt.cm.viridis(np.linspace(0, 1, len(linear_results)))
        for i, r in enumerate(linear_results):
            mean_c_path = np.mean(r['sim_paths']['c_path_vfi'], axis=0)
            ax.plot(mean_c_path, label=f"({r['nk']},{r['nkpps']})", color=colors[i])
        ax.set_title("平均消费路径 (interp='linear')")
        ax.set_xlabel('年龄组索引')
        ax.set_ylabel('平均消费')
        ax.legend(fontsize='small')
        ax.grid(True)
        
        # 3. 平均消费路径对比 (spline)
        ax = axes[1, 1]
        colors = plt.cm.plasma(np.linspace(0, 1, len(spline_results)))
        for i, r in enumerate(spline_results):
            mean_c_path = np.mean(r['sim_paths']['c_path_vfi'], axis=0)
            ax.plot(mean_c_path, label=f"({r['nk']},{r['nkpps']})", color=colors[i])
        ax.set_title("平均消费路径 (interp='spline')")
        ax.set_xlabel('年龄组索引')
        ax.set_ylabel('平均消费')
        ax.legend(fontsize='small')
        ax.grid(True)
        
        # [新] 4. 平均资产路径对比 (spline, 作为示例)
        ax = axes[0, 1]
        colors = plt.cm.plasma(np.linspace(0, 1, len(spline_results)))
        for i, r in enumerate(spline_results):
            mean_k_path = np.mean(r['sim_paths']['k_path_vfi'], axis=0)
            ax.plot(mean_k_path, label=f"({r['nk']},{r['nkpps']})", color=colors[i])
        ax.set_title("平均资产路径 (interp='spline')")
        ax.set_xlabel('年龄组索引')
        ax.set_ylabel('平均资产')
        ax.legend(fontsize='small')
        ax.grid(True)

        plt.savefig('./py/vfi_grid_sensitivity_analysis.png', dpi=300)
        print("📈 敏感性分析图表已保存。")
        plt.show()

def main():
    if not MATLAB_AVAILABLE: return

    print("⚠️ 注意: 请确保您的 main_olg_v8_utils.m 已被修改，")
    print("   以便 HHSolutionByAge... 函数能从 cS 结构体中读取 'interpolation_method'。")

    tester = VFIGridSensitivityTester()

    grid_configurations = [
        (5, 5), (10, 10), (20, 20)
    ]

    tester.run_sensitivity_analysis(grid_configurations)

if __name__ == "__main__":
    main()