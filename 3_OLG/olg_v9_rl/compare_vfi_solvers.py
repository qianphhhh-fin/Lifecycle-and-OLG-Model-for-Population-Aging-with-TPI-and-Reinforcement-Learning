# --- 开始文件：compare_vfi_solvers.py ---

"""
[Python版] VFI求解器敏感性与性能比较脚本

目标:
- 在Python环境中，调用MATLAB引擎来运行不同配置的VFI求解器。
- 系统性地比较 'grid' 和 'hybrid' 两种求解方法在不同网格密度下的性能。
- 所有模拟和评估都在Python端完成，与RL的评估框架完全对齐。
- 这有助于诊断和理解MATLAB端不同求解器策略的最终效果。
"""

import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

import numpy as np
import matplotlib.pyplot as plt
import time
from typing import Dict, Any, List

# MATLAB Engine导入
try:
    import matlab.engine
    MATLAB_AVAILABLE = True
except ImportError:
    MATLAB_AVAILABLE = False
    print("❌ MATLAB Engine 不可用，无法继续。")
    exit()

# 导入与 `compare_rl_and_vfi_matlab.py` 一致的自定义模块
from main_olg_v9_utils import OLG_V9_Utils, TempParamSHH

# 配置matplotlib中文字体支持 (与另一个脚本一致)
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
            # 备用字体逻辑
            pass
    except Exception as e:
        print(f"⚠️ 设置中文字体时出错: {e}")
    plt.rcParams['axes.unicode_minus'] = False

setup_chinese_fonts()


class VFISolverComparator:
    """在Python中比较MATLAB VFI求解器的性能"""
    
    def __init__(self):
        self.eng = matlab.engine.start_matlab()
        self.eng.addpath(os.getcwd(), nargout=0)
        print("✅ MATLAB Engine 已启动。")

    def __del__(self):
        if hasattr(self, 'eng') and self.eng is not None:
            self.eng.quit()
            print("✅ MATLAB Engine已关闭。")

    def _convert_dict_to_matlab_struct(self, py_dict: Dict) -> Dict:
        """将Python字典转换为MATLAB结构体 (与RL比较脚本一致)"""
        matlab_struct = {}
        for key, value in py_dict.items():
            if key in ['physAgeMap', 'interpolation_method']: continue
            if isinstance(value, np.ndarray):
                # 确保数据类型和形状正确
                if value.ndim == 1:
                    matlab_struct[key] = matlab.double(value.tolist())
                else:
                    matlab_struct[key] = matlab.double(value.tolist())
            elif isinstance(value, list):
                matlab_struct[key] = matlab.double(value)
            elif isinstance(value, (int, float, bool)):
                matlab_struct[key] = float(value)
        return matlab_struct

    def run_matlab_vfi(self, M_env: Dict, cS_config: Dict, solver_method: str) -> Dict[str, Any]:
        """
        为指定的宏观环境和cS配置，调用指定的MATLAB VFI求解器。
        
        Args:
            M_env (Dict): 宏观经济环境参数。
            cS_config (Dict): 控制模型结构和网格的参数 (如 nk, nkpps)。
            solver_method (str): 要使用的求解器 ('grid' 或 'hybrid')。

        Returns:
            Dict: 包含策略矩阵和相关参数的结果字典。
        """
        # 1. 创建基础cS和paramS对象
        cS_python = OLG_V9_Utils.ParameterValues_HuggettStyle()
        
        # 2. 应用自定义配置
        for key, value in cS_config.items():
            setattr(cS_python, key, value)
        
        # 3. 重新生成依赖网格的参数 (在Python端完成)
        cS_python = OLG_V9_Utils.generateGrids(cS_python)

        # 4. 创建 paramS
        (leLogGridV, leTrProbM, leProb1V) = OLG_V9_Utils.EarningProcess_olgm(cS_python)
        paramS_vfi_dict = {
            'leLogGridV': leLogGridV, 'leTrProbM': leTrProbM, 'leProb1V': leProb1V,
            'leGridV': np.exp(leLogGridV),
            'ageEffV_new': cS_python.ageEffV_new,
            'tau_l': M_env['tau_l'],
            'theta_payg_actual_for_hh': 0.10, # 假设值，可以从M_env传入
            'pps_tax_deferral_active': bool(cS_python.pps_active),
        }

        # 5. 转换为MATLAB格式
        cS_matlab_dict = self._convert_dict_to_matlab_struct(cS_python.__dict__)
        paramS_vfi_matlab = self._convert_dict_to_matlab_struct(paramS_vfi_dict)
        
        bV_payg_vfi = np.zeros(cS_python.aD_new)
        if cS_python.aR_new < cS_python.aD_new:
            bV_payg_vfi[cS_python.aR_new:] = M_env['b_payg_avg_retiree']
        bV_payg_matlab = matlab.double(bV_payg_vfi.tolist())
        
        # 6. [核心] 调用MATLAB VFI求解器，并传入solver_method
        print(f"  调用MATLAB: HHSolution_VFI_Huggett (solver='{solver_method}')")
        cPolM, kPolM, cPpsM, _ = self.eng.main_olg_v8_utils.HHSolution_VFI_Huggett(
            M_env['R_k_net_factor'], M_env['w_gross'], M_env['TR_total'],
            bV_payg_matlab, paramS_vfi_matlab, cS_matlab_dict, solver_method, nargout=4)
        
        return {
            'cPolM': np.array(cPolM),
            'kPolM': np.array(kPolM),
            'cPpsPolM_choice': np.array(cPpsM),
            'cS_python_obj': cS_python,
            'paramS_vfi_dict': paramS_vfi_dict, # <<<< 添加这一行
            'M_test': M_env, # <<<< 添加这一行，M_env 就是我们需要的宏观环境
        }

    def _calculate_lifetime_utility(self, c_path: np.ndarray, cS: Any, use_survival_prob: bool) -> float:
        """计算终身效用 (与RL比较脚本完全一致)"""
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

# --- 在 compare_vfi_solvers.py 的 VFISolverComparator class 内部 ---
# ... (所有其他方法保持不变) ...

    def run_comparison(self, n_sim=5000, random_seed=42):
        """主比较流程"""
        print("\n" + "="*80)
        print("🔬 Python端 VFI 求解器性能比较")
        print("=" * 80)
        
        # 1. 定义测试配置
        solver_methods_to_test = ['grid', 'hybrid']
        grid_configurations = [
            {'nk': 20, 'nkpps': 10},
            {'nk': 5, 'nkpps': 5},
        ]

        # 2. 定义固定的宏观环境 (与MATLAB脚本对齐)
        M_FIXED = {
            'R_k_net_factor': 1.03,
            'w_gross': 1.8,
            'TR_total': 0.0,
            'b_payg_avg_retiree': 0.4,
            'tau_l': 0.15,
            'theta_payg_actual': 0.10, # <<<< 添加这一行，可以使用一个合理的默认值
        }
        
        # 3. 生成统一的随机冲击“剧本” (1-based for simulators)
        # 注意：这里需要确保 MarkovChainSimulation_AgeGroup 存在并能正确工作
        temp_cs = OLG_V9_Utils.ParameterValues_HuggettStyle()
        temp_cs.nSim = n_sim
        _, tr_prob, p0 = OLG_V9_Utils.EarningProcess_olgm(temp_cs)
        eIdxM_group_0based = OLG_V9_Utils.MarkovChainSimulation_AgeGroup(n_sim, temp_cs, p0, tr_prob)
        eIdxM_group_for_sim = eIdxM_group_0based + 1

        # 4. 循环遍历所有配置并收集结果
        all_results = []
        for cS_config in grid_configurations:
            for solver_method in solver_methods_to_test:
                print(f"\n--- 正在评估配置: {cS_config}, solver='{solver_method}' ---")
                
                # a. 调用MATLAB求解VFI (这一步返回的 vfi_results 包含列向量 kGridV)
                start_time = time.time()
                vfi_results = self.run_matlab_vfi(M_FIXED, cS_config, solver_method)
                solve_time = time.time() - start_time
                print(f"  VFI求解耗时: {solve_time:.2f} 秒。")

                # ===============================
                # === 在这里添加核心修正代码 ===
                # ===============================
                # 在调用模拟器之前，手动修正 vfi_results 内部的 cS 对象，
                # 将网格向量从 (n, 1) 转换为 (n,)，以满足 scipy 的要求。
                # 这样做可以避免修改 main_olg_v9_utils.py 文件。
                
                print("  预处理数据：将网格向量转换为1D格式以适配模拟器...")
                cS_obj_for_sim = vfi_results['cS_python_obj']
                
                if hasattr(cS_obj_for_sim, 'kGridV') and cS_obj_for_sim.kGridV.ndim > 1:
                    cS_obj_for_sim.kGridV = cS_obj_for_sim.kGridV.flatten()
                
                if hasattr(cS_obj_for_sim, 'kppsGridV') and cS_obj_for_sim.kppsGridV.ndim > 1:
                    cS_obj_for_sim.kppsGridV = cS_obj_for_sim.kppsGridV.flatten()
                
                # 将修正后的 cS 对象放回 vfi_results 字典中
                vfi_results['cS_python_obj'] = cS_obj_for_sim
                # ===============================
                # === 修正结束 ===
                # ===============================

                # b. 在Python端进行模拟 (现在它会接收到格式正确的网格)
                start_time = time.time()
                sim_paths = OLG_V9_Utils.HHSimulation_olgm_vfi_simplified(
                    vfi_results, n_sim, eIdxM_group_for_sim
                )
                sim_time = time.time() - start_time
                print(f"  Python端模拟耗时: {sim_time:.2f} 秒。")
                
                # c. 计算效用
                c_path_vfi = sim_paths['c_path_vfi']
                utility_vfi = np.array([self._calculate_lifetime_utility(c_path_vfi[j,:], cS_obj_for_sim, True) for j in range(n_sim)])
                
                # d. 存储结果
                all_results.append({
                    'nk': cS_config['nk'],
                    'nkpps': cS_config['nkpps'],
                    'solver_method': solver_method,
                    'mean_utility': np.mean(utility_vfi),
                    'std_utility': np.std(utility_vfi),
                })
                print(f"📈 结果: 平均效用 = {np.mean(utility_vfi):.4f} (标准差 = {np.std(utility_vfi):.4f})")
        
        # 5. 分析和绘图
        self.analyze_and_plot(all_results, grid_configurations)


    def analyze_and_plot(self, all_results: List[Dict], grid_configs: List[Dict]):
        """分析并绘制结果 (与MATLAB脚本的绘图逻辑对齐)"""
        fig, ax = plt.subplots(1, 1, figsize=(10, 7))
        
        grid_labels = [f"nk={c['nk']}, nkpps={c['nkpps']}" for c in grid_configs]
        n_grids = len(grid_configs)
        bar_width = 0.35
        index = np.arange(n_grids)
        
        grid_solver_means = [r['mean_utility'] for r in all_results if r['solver_method'] == 'grid']
        hybrid_solver_means = [r['mean_utility'] for r in all_results if r['solver_method'] == 'hybrid']
        
        bar1 = ax.bar(index - bar_width/2, grid_solver_means, bar_width, label='Grid Solver', color='royalblue')
        bar2 = ax.bar(index + bar_width/2, hybrid_solver_means, bar_width, label='Hybrid Solver', color='seagreen')

        ax.set_ylabel('平均终身效用')
        ax.set_title('VFI 求解器性能比较 (Python评估框架)')
        ax.set_xticks(index)
        ax.set_xticklabels(grid_labels)
        ax.legend()
        ax.grid(True, linestyle='--', alpha=0.6)

        # 在每个bar上显示数值
        for bar in [bar1, bar2]:
            for rect in bar:
                height = rect.get_height()
                ax.annotate(f'{height:.2f}',
                            xy=(rect.get_x() + rect.get_width() / 2, height),
                            xytext=(0, 3),  # 3 points vertical offset
                            textcoords="offset points",
                            ha='center', va='bottom')

        plt.tight_layout()
        save_path = './py/vfi_solver_comparison_python.png'
        plt.savefig(save_path, dpi=300)
        print(f"\n📈 比较图表已保存到: {save_path}")
        plt.show()


def main():
    if not MATLAB_AVAILABLE: return
    comparator = VFISolverComparator()
    comparator.run_comparison()

if __name__ == "__main__":
    main()