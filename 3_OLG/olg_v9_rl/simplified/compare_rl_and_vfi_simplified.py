"""
简化版RL和VFI比较脚本 - 固定宏观参数版本

🎯 核心目标：
- 使用固定宏观参数，RL状态空间与VFI基本相同
- 公平比较RL和VFI的性能
- 验证RL是否能在相同条件下达到VFI的结果

🎮 RL动作空间：
- 2维连续动作：[PPS缴费比例, 消费比例]
- 决策顺序：先PPS缴费 → 再消费 → 最后储蓄（自动）
- 与VFI的离散网格决策形成对比

📚 使用方法：
1. 完整VFI vs RL比较（默认）:
   python compare_rl_and_vfi_simplified.py
   
2. 仅RL模型评估（快速模式）:
   python compare_rl_and_vfi_simplified.py --rl-only
   
3. 仅VFI模型评估（调试模式）:
   python compare_rl_and_vfi_simplified.py --vfi-only
   
🚀 eva_rl_only模式优势：
- 跳过MATLAB VFI求解，启动更快
- 专注于RL模型性能评估
- 适用于模型调试和快速验证
- 详细的生命周期路径分析

🔧 eva_vfi_only模式优势：
- 跳过RL模型加载，专注VFI分析
- 快速验证VFI求解结果
- 详细的VFI策略函数分析
- 适用于VFI参数调试和验证
"""

import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

import numpy as np
import matplotlib.pyplot as plt
import pickle
import time
from typing import Dict, Any, Tuple

# MATLAB Engine导入
try:
    import matlab.engine
    MATLAB_AVAILABLE = True
except ImportError:
    MATLAB_AVAILABLE = False

# 导入必要的模块
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sbx import SAC as SBX_SAC
from main_olg_v9_utils import OLG_V9_Utils
from simplified.main_olg_v9_utils_simplified import OLGUtilsSimplified
from simplified.main_olg_v9_sac_sbx_simplified import evaluate_policy_lifecycle_simulation_simplified

import matplotlib
import matplotlib.font_manager as fm
# 配置matplotlib中文字体支持
def setup_chinese_fonts():
    """设置matplotlib中文字体"""
    chinese_fonts = ['SimHei', 'Microsoft YaHei', 'WenQuanYi Micro Hei', 'DejaVu Sans']
    available_fonts = [f.name for f in fm.fontManager.ttflist]
    
    for font in chinese_fonts:
        if font in available_fonts:
            matplotlib.rcParams['font.sans-serif'] = [font] + matplotlib.rcParams['font.sans-serif']
            print(f"✅ 使用中文字体: {font}")
            break
    else:
        print("⚠️ 未找到中文字体，可能显示为方框")
        matplotlib.rcParams['font.sans-serif'] = ['DejaVu Sans']
    
    matplotlib.rcParams['axes.unicode_minus'] = False  # 正确显示负号

# 初始化中文字体
setup_chinese_fonts()

class RLVFIComparatorSimplified:
    """简化版RL和VFI比较器"""
    
    def __init__(self, use_matlab=True):
        """
        初始化比较器
        
        Args:
            use_matlab: 是否使用MATLAB Engine（RL单独评估时可设为False）
        """
        # 固定的测试参数（与VFI保持一致）
        self.M_fixed = {
            'R_k_net_factor': 1.03,
            'w_gross': 2.0,
            'TR_total': 0.1,
            'b_payg_avg_retiree': 0.4,
            'tau_l': 0.15,
            'theta_payg_actual': 0.12
        }
        
        print("🎯 简化版RL评估器")
        print("📊 固定宏观参数:")
        for key, value in self.M_fixed.items():
            print(f"  {key} = {value:.3f}")
        
        # 初始化MATLAB Engine（可选）
        self.eng = None
        if use_matlab:
            if not MATLAB_AVAILABLE:
                raise ImportError("MATLAB Engine不可用")
            self.eng = matlab.engine.start_matlab()
            self.eng.addpath(os.getcwd(), nargout=0)
            print("✅ MATLAB Engine启动成功")
        else:
            print("⚠️ 跳过MATLAB Engine初始化（仅RL评估模式）")
    
    def load_rl_model_simplified(self, model_path: str = None):
        """加载简化版RL模型"""
        if model_path is None:
            # 检查最佳模型
            best_model_path = './simplified/best_model_simplified/best_model.zip'
            final_model_path = './simplified/final_sac_agent_olg_simplified.zip'
            
            if os.path.exists(best_model_path):
                model_path = best_model_path
                print(f"🏆 使用最佳模型: {model_path}")
            elif os.path.exists(final_model_path):
                model_path = final_model_path
                print(f"📁 使用最终模型: {model_path}")
            else:
                raise FileNotFoundError("未找到简化版模型文件")
        
        model = SBX_SAC.load(model_path)
        
        # 加载配置
        config_path = './simplified/training_config_simplified.pkl'
        config = {}
        if os.path.exists(config_path):
            with open(config_path, 'rb') as f:
                config = pickle.load(f)
        
        return model, config
    
    def run_vfi_simplified(self):
        """运行VFI（使用固定宏观参数）"""
        print("\n🔧 运行VFI（固定宏观参数）...")
        
        start_time = time.time()
        
        # 使用Python初始化参数
        cS = OLG_V9_Utils.ParameterValues_HuggettStyle()
        leLogGridV, leTrProbM, leProb1V = OLG_V9_Utils.EarningProcess_olgm(cS)
        leGridV = np.exp(leLogGridV)
        
        # 构建MATLAB兼容的参数
        paramS_vfi_matlab = {
            'leLogGridV': matlab.double(leLogGridV.tolist()),
            'leTrProbM': matlab.double(leTrProbM.tolist()),
            'leProb1V': matlab.double(leProb1V.tolist()),
            'leGridV': matlab.double(leGridV.tolist()),
            'ageEffV_new': matlab.double(cS.ageEffV_new.tolist()),
            'tau_l': self.M_fixed['tau_l'],
            'theta_payg_actual_for_hh': self.M_fixed['theta_payg_actual'],
            'pps_tax_deferral_active': bool(cS.pps_active)
        }
        
        # 构建PAYG福利向量
        bV_payg = np.zeros(cS.aD_new)
        if cS.aR_new < cS.aD_new:
            bV_payg[cS.aR_new:] = self.M_fixed['b_payg_avg_retiree']
        bV_payg_matlab = matlab.double(bV_payg.tolist())
        
        # 构建MATLAB兼容的cS
        cS_matlab_dict = {}
        for attr_name in dir(cS):
            if not attr_name.startswith('_'):
                attr_value = getattr(cS, attr_name)
                try:
                    if isinstance(attr_value, (int, float, bool, np.integer, np.floating)):
                        cS_matlab_dict[attr_name] = float(attr_value)
                    elif isinstance(attr_value, np.ndarray):
                        if attr_value.ndim <= 2:
                            cS_matlab_dict[attr_name] = matlab.double(attr_value.tolist())
                    elif attr_name == 'physAgeMap':
                        matlab_cell_list = []
                        for sublist in attr_value:
                            if hasattr(sublist, '__iter__'):
                                matlab_indices = [idx + 1 for idx in sublist]
                                matlab_cell_list.append(matlab.double(matlab_indices))
                        cS_matlab_dict[attr_name] = matlab_cell_list
                except (ValueError, TypeError):
                    continue
        
        # 调用MATLAB VFI
        print("🔧 调用MATLAB VFI求解...")
        vfi_start = time.time()
        
        cPolM_vfi, kPolM_vfi, cPpsPolM_vfi, VPolM_vfi = self.eng.main_olg_v8_utils.HHSolution_VFI_Huggett(
            self.M_fixed['R_k_net_factor'],
            self.M_fixed['w_gross'],
            self.M_fixed['TR_total'],
            bV_payg_matlab,
            paramS_vfi_matlab,
            cS_matlab_dict,
            nargout=4
        )
        
        # 转换为numpy数组
        # 转换MATLAB结果为numpy数组
        cPolM_vfi = np.array(cPolM_vfi)
        kPolM_vfi = np.array(kPolM_vfi)
        cPpsPolM_vfi = np.array(cPpsPolM_vfi)
        VPolM_vfi = np.array(VPolM_vfi)
    
        
        vfi_time = time.time() - vfi_start
        total_time = time.time() - start_time
        
        print(f"✅ VFI求解完成，耗时: {vfi_time:.2f} 秒")
        
        # 准备返回结果
        result_dict = {
            'cPolM': cPolM_vfi,
            'kPolM': kPolM_vfi,
            'cPpsPolM': cPpsPolM_vfi,
            'cPpsPolM_choice': cPpsPolM_vfi,
            'cPolM_q': cPolM_vfi,
            'VPolM': VPolM_vfi,
            'M_fixed': self.M_fixed,
            'cS': cS,
            'paramS_vfi': {
                'leLogGridV': leLogGridV,
                'leTrProbM': leTrProbM,
                'leProb1V': leProb1V,
                'leGridV': leGridV,
                'ageEffV_new': cS.ageEffV_new,
                'tau_l': self.M_fixed['tau_l'],
                'theta_payg_actual_for_hh': self.M_fixed['theta_payg_actual']
            },
            'vfi_time': vfi_time,
            'total_time': total_time,
            'bV_payg': bV_payg
        }
        
        return result_dict
    
    def simulate_lifecycle_comparison_simplified(self, rl_model, vfi_results, rl_config, 
                                               n_sim=100, random_seed=42):
        """简化版生命周期比较"""
        print(f"\n🔄 简化版生命周期比较 (n_sim={n_sim}, seed={random_seed})...")
        
        np.random.seed(random_seed)
        
        cS = vfi_results['cS']
        paramS_vfi = vfi_results['paramS_vfi']
        
        # VFI策略
        kPolM_vfi = vfi_results['kPolM']
        cPpsPolM_vfi = vfi_results['cPpsPolM_choice']
        cPolM_vfi = vfi_results['cPolM_q']
        
        # 生成效率冲击序列
        aD_orig = int(cS.aD_orig)
        random_numbers = np.random.rand(n_sim, aD_orig)
        eIdxM_annual = OLG_V9_Utils.MarkovChainSimulation(
            n_sim, aD_orig, paramS_vfi['leProb1V'], paramS_vfi['leTrProbM'], random_numbers
        )
        
        # VFI模拟 - 创建参数对象
        print("🎯 VFI模拟...")
        
        # 创建参数对象（与main_olg_v9_utils.py中的TempParamSHH兼容）
        class TempParamSVFI:
            def __init__(self, params_dict):
                for key, value in params_dict.items():
                    setattr(self, key, value)
        
        paramS_vfi_obj = TempParamSVFI(paramS_vfi)
        
        kHistM_vfi, kPpsHistM_vfi, cHistM_vfi, cppsHistM_vfi = OLG_V9_Utils.HHSimulation_olgm(
            kPolM_vfi, cPpsPolM_vfi, cPolM_vfi, eIdxM_annual,
            self.M_fixed['R_k_net_factor'], self.M_fixed['w_gross'], self.M_fixed['TR_total'],
            vfi_results['bV_payg'], paramS_vfi_obj, cS
        )
        
        # RL模拟（简化版）
        print("🤖 RL模拟（简化版）...")
        paramS_for_rl = {
            'leLogGridV': paramS_vfi['leLogGridV'],
            'leTrProbM': paramS_vfi['leTrProbM'],
            'leProb1V': paramS_vfi['leProb1V'],
            'leGridV': paramS_vfi['leGridV'],
            'ageEffV_new': paramS_vfi['ageEffV_new']
        }
        
        kHistM_rl, kPpsHistM_rl, cHistM_rl, cppsHistM_rl = OLGUtilsSimplified.HHSimulation_olgm_rl_simplified(
            rl_model, cS, paramS_for_rl, self.M_fixed, eIdxM_annual
        )
        
        # 计算生涯效用
        print("📊 计算生涯效用...")
        lifetime_utility_vfi = self._calculate_lifetime_utility(cHistM_vfi, cS)
        lifetime_utility_rl = self._calculate_lifetime_utility(cHistM_rl, cS)
        
        # 统计结果
        mean_utility_vfi = np.mean(lifetime_utility_vfi)
        mean_utility_rl = np.mean(lifetime_utility_rl)
        std_utility_vfi = np.std(lifetime_utility_vfi)
        std_utility_rl = np.std(lifetime_utility_rl)
        utility_diff = mean_utility_rl - mean_utility_vfi
        
        print(f"\n📈 简化版比较结果:")
        print(f"   VFI平均效用: {mean_utility_vfi:.4f} ± {std_utility_vfi:.4f}")
        print(f"   RL平均效用: {mean_utility_rl:.4f} ± {std_utility_rl:.4f}")
        print(f"   效用差异(RL-VFI): {utility_diff:.6f} ({utility_diff/abs(mean_utility_vfi)*100:.3f}%)")
        
        # 返回结果
        return {
            'n_sim': n_sim,
            'random_seed': random_seed,
            'k_path_vfi': kHistM_vfi,
            'k_path_rl': kHistM_rl,
            'c_path_vfi': cHistM_vfi,
            'c_path_rl': cHistM_rl,
            'cpps_path_vfi': cppsHistM_vfi,
            'cpps_path_rl': cppsHistM_rl,
            'lifetime_utility_vfi': lifetime_utility_vfi,
            'lifetime_utility_rl': lifetime_utility_rl,
            'mean_utility_vfi': mean_utility_vfi,
            'mean_utility_rl': mean_utility_rl,
            'std_utility_vfi': std_utility_vfi,
            'std_utility_rl': std_utility_rl,
            'utility_diff': utility_diff,
            'M_fixed': self.M_fixed,
            'simplified_version': True,
            'state_space_dim': 4
        }
    
    def _calculate_lifetime_utility(self, cHistM, cS):
        """计算生涯效用"""
        n_sim, aD_orig = cHistM.shape
        lifetime_utility = np.zeros(n_sim)
        
        for i_sim in range(n_sim):
            utility_sum = 0
            for age_idx in range(aD_orig):
                c_val = cHistM[i_sim, age_idx]
                _, u_val = OLG_V9_Utils.CES_utility(c_val, cS.sigma, cS)
                discount_factor = (cS.beta ** age_idx)
                utility_sum += discount_factor * u_val
            lifetime_utility[i_sim] = utility_sum
        
        return lifetime_utility
    
    def evaluate_vfi_only_simplified(self, n_sim=100, random_seed=42):
        """
        仅评估VFI模型性能（不与RL比较）
        
        Args:
            n_sim: 模拟个体数
            random_seed: 随机种子
            
        Returns:
            dict: VFI评估结果
        """
        print(f"\n🔧 VFI模型单独评估 (n_sim={n_sim}, seed={random_seed})...")
        
        # 运行VFI求解
        vfi_results = self.run_vfi_simplified()
        
        np.random.seed(random_seed)
        
        # 初始化参数
        cS = vfi_results['cS']
        paramS_vfi = vfi_results['paramS_vfi']
        
        # VFI策略
        kPolM_vfi = vfi_results['kPolM']
        cPpsPolM_vfi = vfi_results['cPpsPolM_choice']
        cPolM_vfi = vfi_results['cPolM_q']
        
        # 生成效率冲击序列
        aD_orig = int(cS.aD_orig)
        random_numbers = np.random.rand(n_sim, aD_orig)
        eIdxM_annual = OLG_V9_Utils.MarkovChainSimulation(
            n_sim, aD_orig, paramS_vfi['leProb1V'], paramS_vfi['leTrProbM'], random_numbers
        )
        
        # VFI模拟
        print("🔧 运行VFI生命周期模拟...")
        
        # 创建参数对象（与main_olg_v9_utils.py中的TempParamSHH兼容）
        class TempParamSVFI:
            def __init__(self, params_dict):
                for key, value in params_dict.items():
                    setattr(self, key, value)
        
        paramS_vfi_obj = TempParamSVFI(paramS_vfi)
        
        kHistM_vfi, kPpsHistM_vfi, cHistM_vfi, cppsHistM_vfi = OLG_V9_Utils.HHSimulation_olgm(
            kPolM_vfi, cPpsPolM_vfi, cPolM_vfi, eIdxM_annual,
            self.M_fixed['R_k_net_factor'], self.M_fixed['w_gross'], self.M_fixed['TR_total'],
            vfi_results['bV_payg'], paramS_vfi_obj, cS
        )

        # 打印VFI决策变量的描述性统计
        print("\n📊 VFI决策变量路径描述性统计:")
        print(f"消费决策矩阵形状: {cHistM_vfi.shape}")
        print(f"消费决策 (cHistM_vfi): 均值={np.mean(cHistM_vfi):.4f}, 中位数={np.median(cHistM_vfi):.4f}, 最小值={np.min(cHistM_vfi):.4f}, 最大值={np.max(cHistM_vfi):.4f}")
        print(f"资产决策 (kHistM_vfi): 均值={np.mean(kHistM_vfi):.4f}, 中位数={np.median(kHistM_vfi):.4f}, 最小值={np.min(kHistM_vfi):.4f}, 最大值={np.max(kHistM_vfi):.4f}")
        print(f"PPS缴费决策 (cppsHistM_vfi): 均值={np.mean(cppsHistM_vfi):.4f}, 中位数={np.median(cppsHistM_vfi):.4f}, 最小值={np.min(cppsHistM_vfi):.4f}, 最大值={np.max(cppsHistM_vfi):.4f}")

        
        # 计算生涯效用
        print("📊 计算VFI生涯效用...")
        lifetime_utility_vfi = self._calculate_lifetime_utility(cHistM_vfi, cS)
        
        # 统计结果
        mean_utility_vfi = np.mean(lifetime_utility_vfi)
        std_utility_vfi = np.std(lifetime_utility_vfi)
        median_utility_vfi = np.median(lifetime_utility_vfi)
        min_utility_vfi = np.min(lifetime_utility_vfi)
        max_utility_vfi = np.max(lifetime_utility_vfi)
        
        # 计算其他统计指标
        mean_consumption = np.mean(cHistM_vfi)
        mean_savings = np.mean(kHistM_vfi)
        mean_pps_savings = np.mean(kPpsHistM_vfi)
        mean_pps_contrib = np.mean(cppsHistM_vfi)
        
        # 计算策略函数统计
        mean_consumption_policy = np.mean(cPolM_vfi)
        mean_savings_policy = np.mean(kPolM_vfi)
        mean_pps_policy = np.mean(cPpsPolM_vfi)
        
        print(f"\n📈 VFI模型评估结果:")
        print(f"   平均生涯效用: {mean_utility_vfi:.4f} ± {std_utility_vfi:.4f}")
        print(f"   中位数效用: {median_utility_vfi:.4f}")
        print(f"   效用范围: [{min_utility_vfi:.4f}, {max_utility_vfi:.4f}]")
        print(f"   平均消费: {mean_consumption:.4f}")
        print(f"   平均储蓄: {mean_savings:.4f}")
        print(f"   平均PPS储蓄: {mean_pps_savings:.4f}")
        print(f"   平均PPS缴费: {mean_pps_contrib:.4f}")
        print(f"   VFI求解时间: {vfi_results['vfi_time']:.2f}秒")
        
        # 返回结果
        return {
            'n_sim': n_sim,
            'random_seed': random_seed,
            'k_path_vfi': kHistM_vfi,
            'kpps_path_vfi': kPpsHistM_vfi,
            'c_path_vfi': cHistM_vfi,
            'cpps_path_vfi': cppsHistM_vfi,
            'lifetime_utility_vfi': lifetime_utility_vfi,
            'mean_utility_vfi': mean_utility_vfi,
            'std_utility_vfi': std_utility_vfi,
            'median_utility_vfi': median_utility_vfi,
            'min_utility_vfi': min_utility_vfi,
            'max_utility_vfi': max_utility_vfi,
            'mean_consumption': mean_consumption,
            'mean_savings': mean_savings,
            'mean_pps_savings': mean_pps_savings,
            'mean_pps_contrib': mean_pps_contrib,
            'mean_consumption_policy': mean_consumption_policy,
            'mean_savings_policy': mean_savings_policy,
            'mean_pps_policy': mean_pps_policy,
            'k_policy_vfi': kPolM_vfi,
            'c_policy_vfi': cPolM_vfi,
            'cpps_policy_vfi': cPpsPolM_vfi,
            'V_policy_vfi': vfi_results['VPolM'],
            'M_fixed': self.M_fixed,
            'vfi_time': vfi_results['vfi_time'],
            'cS': cS,
            'paramS_vfi': paramS_vfi,
            'simplified_version': True,
            'vfi_only_evaluation': True
        }
    
    def evaluate_rl_only_simplified(self, rl_model, rl_config=None, n_sim=100, random_seed=42):
        """
        仅评估RL模型性能（不与VFI比较）
        
        Args:
            rl_model: RL模型
            rl_config: RL配置（可选）
            n_sim: 模拟个体数
            random_seed: 随机种子
            
        Returns:
            dict: RL评估结果
        """
        print(f"\n🤖 RL模型单独评估 (n_sim={n_sim}, seed={random_seed})...")
        
        np.random.seed(random_seed)
        
        # 初始化参数
        cS = OLG_V9_Utils.ParameterValues_HuggettStyle()
        leLogGridV, leTrProbM, leProb1V = OLG_V9_Utils.EarningProcess_olgm(cS)
        leGridV = np.exp(leLogGridV)
        
        # 生成效率冲击序列
        aD_orig = int(cS.aD_orig)
        random_numbers = np.random.rand(n_sim, aD_orig)
        eIdxM_annual = OLG_V9_Utils.MarkovChainSimulation(
            n_sim, aD_orig, leProb1V, leTrProbM, random_numbers
        )
        
        # RL模拟（简化版）
        print("🚀 运行RL生命周期模拟...")
        paramS_for_rl = {
            'leLogGridV': leLogGridV,
            'leTrProbM': leTrProbM,
            'leProb1V': leProb1V,
            'leGridV': leGridV,
            'ageEffV_new': cS.ageEffV_new
        }
        
        kHistM_rl, kPpsHistM_rl, cHistM_rl, cppsHistM_rl = OLGUtilsSimplified.HHSimulation_olgm_rl_simplified(
            rl_model, cS, paramS_for_rl, self.M_fixed, eIdxM_annual
        )

        # 打印RL决策变量的描述性统计
        print("\n📊 RL决策变量路径描述性统计:")
        print(f"消费决策矩阵形状: {cHistM_rl.shape}")
        print(f"消费决策 (cHistM_rl): 均值={np.mean(cHistM_rl):.4f}, 中位数={np.median(cHistM_rl):.4f}, 最小值={np.min(cHistM_rl):.4f}, 最大值={np.max(cHistM_rl):.4f}")
        print(f"资产决策 (kHistM_rl): 均值={np.mean(kHistM_rl):.4f}, 中位数={np.median(kHistM_rl):.4f}, 最小值={np.min(kHistM_rl):.4f}, 最大值={np.max(kHistM_rl):.4f}")
        print(f"PPS缴费决策 (cppsHistM_rl): 均值={np.mean(cppsHistM_rl):.4f}, 中位数={np.median(cppsHistM_rl):.4f}, 最小值={np.min(cppsHistM_rl):.4f}, 最大值={np.max(cppsHistM_rl):.4f}")
        

        
        # 计算生涯效用
        print("📊 计算RL生涯效用...")
        lifetime_utility_rl = self._calculate_lifetime_utility(cHistM_rl, cS)
        
        # 统计结果
        mean_utility_rl = np.mean(lifetime_utility_rl)
        std_utility_rl = np.std(lifetime_utility_rl)
        median_utility_rl = np.median(lifetime_utility_rl)
        min_utility_rl = np.min(lifetime_utility_rl)
        max_utility_rl = np.max(lifetime_utility_rl)
        
        # 计算其他统计指标
        mean_consumption = np.mean(cHistM_rl)
        mean_savings = np.mean(kHistM_rl)
        mean_pps_contrib = np.mean(cppsHistM_rl)
        
        print(f"\n📈 RL模型评估结果:")
        print(f"   平均生涯效用: {mean_utility_rl:.4f} ± {std_utility_rl:.4f}")
        print(f"   中位数效用: {median_utility_rl:.4f}")
        print(f"   效用范围: [{min_utility_rl:.4f}, {max_utility_rl:.4f}]")
        print(f"   平均消费: {mean_consumption:.4f}")
        print(f"   平均储蓄: {mean_savings:.4f}")
        print(f"   平均PPS缴费: {mean_pps_contrib:.4f}")
        
        # 返回结果
        return {
            'n_sim': n_sim,
            'random_seed': random_seed,
            'k_path_rl': kHistM_rl,
            'kpps_path_rl': kPpsHistM_rl,
            'c_path_rl': cHistM_rl,
            'cpps_path_rl': cppsHistM_rl,
            'lifetime_utility_rl': lifetime_utility_rl,
            'mean_utility_rl': mean_utility_rl,
            'std_utility_rl': std_utility_rl,
            'median_utility_rl': median_utility_rl,
            'min_utility_rl': min_utility_rl,
            'max_utility_rl': max_utility_rl,
            'mean_consumption': mean_consumption,
            'mean_savings': mean_savings,
            'mean_pps_contrib': mean_pps_contrib,
            'M_fixed': self.M_fixed,
            'simplified_version': True,
            'state_space_dim': 4,
            'rl_only_evaluation': True
        }
    
    def plot_rl_only_results(self, results, save_path='./simplified/rl_only_evaluation.png'):
        """绘制RL单独评估结果"""
        print("\n📊 绘制RL评估结果...")
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle(f'简化版RL模型评估 (n={results["n_sim"]})', fontsize=14)
        
        # 1. 生涯效用分布
        axes[0, 0].hist(results['lifetime_utility_rl'], bins=20, alpha=0.7, color='red', edgecolor='black')
        axes[0, 0].axvline(results['mean_utility_rl'], color='darkred', linestyle='-', linewidth=2, label=f'均值: {results["mean_utility_rl"]:.4f}')
        axes[0, 0].axvline(results['median_utility_rl'], color='orange', linestyle='--', linewidth=2, label=f'中位数: {results["median_utility_rl"]:.4f}')
        axes[0, 0].set_xlabel('生涯效用')
        axes[0, 0].set_ylabel('频数')
        axes[0, 0].set_title('生涯效用分布')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        # 2. 平均生命周期消费路径
        mean_c_rl = np.mean(results['c_path_rl'], axis=0)
        std_c_rl = np.std(results['c_path_rl'], axis=0)
        ages = np.arange(20, 20 + len(mean_c_rl))
        
        axes[0, 1].plot(ages, mean_c_rl, 'r-', linewidth=2, label='平均消费')
        axes[0, 1].fill_between(ages, mean_c_rl - std_c_rl, mean_c_rl + std_c_rl, alpha=0.3, color='red', label='±1标准差')
        axes[0, 1].set_xlabel('年龄')
        axes[0, 1].set_ylabel('消费')
        axes[0, 1].set_title('生命周期消费路径')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
        
        # 3. 平均储蓄路径
        mean_k_rl = np.mean(results['k_path_rl'], axis=0)
        std_k_rl = np.std(results['k_path_rl'], axis=0)
        
        axes[1, 0].plot(ages, mean_k_rl, 'b-', linewidth=2, label='平均储蓄')
        axes[1, 0].fill_between(ages, mean_k_rl - std_k_rl, mean_k_rl + std_k_rl, alpha=0.3, color='blue', label='±1标准差')
        axes[1, 0].set_xlabel('年龄')
        axes[1, 0].set_ylabel('储蓄')
        axes[1, 0].set_title('生命周期储蓄路径')
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)
        
        # 4. 统计摘要箱线图
        utility_stats = [results['lifetime_utility_rl']]
        axes[1, 1].boxplot(utility_stats, labels=['生涯效用'])
        axes[1, 1].set_ylabel('效用值')
        axes[1, 1].set_title('效用分布箱线图')
        axes[1, 1].grid(True, alpha=0.3)
        
        # 添加统计信息文本
        stats_text = f"""统计摘要:
        平均值: {results['mean_utility_rl']:.4f}
        标准差: {results['std_utility_rl']:.4f}
        中位数: {results['median_utility_rl']:.4f}
        最小值: {results['min_utility_rl']:.4f}
        最大值: {results['max_utility_rl']:.4f}"""
        
        axes[1, 1].text(0.02, 0.98, stats_text, transform=axes[1, 1].transAxes, 
                        verticalalignment='top', fontsize=8, 
                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"📈 RL评估图表已保存到: {save_path}")
    
    def plot_vfi_only_results(self, results, save_path='./simplified/vfi_only_evaluation.png'):
        """绘制VFI单独评估结果"""
        print("\n📊 绘制VFI评估结果...")
        
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        fig.suptitle(f'简化版VFI模型评估 (n={results["n_sim"]})', fontsize=14)
        
        # 1. 生涯效用分布
        axes[0, 0].hist(results['lifetime_utility_vfi'], bins=20, alpha=0.7, color='blue', edgecolor='black')
        axes[0, 0].axvline(results['mean_utility_vfi'], color='darkblue', linestyle='-', linewidth=2, label=f'均值: {results["mean_utility_vfi"]:.4f}')
        axes[0, 0].axvline(results['median_utility_vfi'], color='orange', linestyle='--', linewidth=2, label=f'中位数: {results["median_utility_vfi"]:.4f}')
        axes[0, 0].set_xlabel('生涯效用')
        axes[0, 0].set_ylabel('频数')
        axes[0, 0].set_title('生涯效用分布')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        # 2. 平均生命周期消费路径
        mean_c_vfi = np.mean(results['c_path_vfi'], axis=0)
        std_c_vfi = np.std(results['c_path_vfi'], axis=0)
        ages = np.arange(20, 20 + len(mean_c_vfi))
        
        axes[0, 1].plot(ages, mean_c_vfi, 'b-', linewidth=2, label='平均消费')
        axes[0, 1].fill_between(ages, mean_c_vfi - std_c_vfi, mean_c_vfi + std_c_vfi, alpha=0.3, color='blue', label='±1标准差')
        axes[0, 1].set_xlabel('年龄')
        axes[0, 1].set_ylabel('消费')
        axes[0, 1].set_title('生命周期消费路径')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
        
        # 3. 平均储蓄路径
        mean_k_vfi = np.mean(results['k_path_vfi'], axis=0)
        std_k_vfi = np.std(results['k_path_vfi'], axis=0)
        
        axes[0, 2].plot(ages, mean_k_vfi, 'g-', linewidth=2, label='平均储蓄')
        axes[0, 2].fill_between(ages, mean_k_vfi - std_k_vfi, mean_k_vfi + std_k_vfi, alpha=0.3, color='green', label='±1标准差')
        axes[0, 2].set_xlabel('年龄')
        axes[0, 2].set_ylabel('储蓄')
        axes[0, 2].set_title('生命周期储蓄路径')
        axes[0, 2].legend()
        axes[0, 2].grid(True, alpha=0.3)
        
        # 4. PPS储蓄路径
        mean_kpps_vfi = np.mean(results['kpps_path_vfi'], axis=0)
        std_kpps_vfi = np.std(results['kpps_path_vfi'], axis=0)
        
        axes[1, 0].plot(ages, mean_kpps_vfi, 'm-', linewidth=2, label='平均PPS储蓄')
        axes[1, 0].fill_between(ages, mean_kpps_vfi - std_kpps_vfi, mean_kpps_vfi + std_kpps_vfi, alpha=0.3, color='magenta', label='±1标准差')
        axes[1, 0].set_xlabel('年龄')
        axes[1, 0].set_ylabel('PPS储蓄')
        axes[1, 0].set_title('生命周期PPS储蓄路径')
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)
        
        # 5. PPS缴费路径
        mean_cpps_vfi = np.mean(results['cpps_path_vfi'], axis=0)
        std_cpps_vfi = np.std(results['cpps_path_vfi'], axis=0)
        
        axes[1, 1].plot(ages, mean_cpps_vfi, 'c-', linewidth=2, label='平均PPS缴费')
        axes[1, 1].fill_between(ages, mean_cpps_vfi - std_cpps_vfi, mean_cpps_vfi + std_cpps_vfi, alpha=0.3, color='cyan', label='±1标准差')
        axes[1, 1].set_xlabel('年龄')
        axes[1, 1].set_ylabel('PPS缴费')
        axes[1, 1].set_title('生命周期PPS缴费路径')
        axes[1, 1].legend()
        axes[1, 1].grid(True, alpha=0.3)
        
        # 6. 统计摘要和策略函数信息
        axes[1, 2].axis('off')
        
        # 添加统计信息文本
        stats_text = f"""生涯效用统计:
        平均值: {results['mean_utility_vfi']:.4f}
        标准差: {results['std_utility_vfi']:.4f}
        中位数: {results['median_utility_vfi']:.4f}
        最小值: {results['min_utility_vfi']:.4f}
        最大值: {results['max_utility_vfi']:.4f}
        
        平均路径统计:
        平均消费: {results['mean_consumption']:.4f}
        平均储蓄: {results['mean_savings']:.4f}
        平均PPS储蓄: {results['mean_pps_savings']:.4f}
        平均PPS缴费: {results['mean_pps_contrib']:.4f}
        
        策略函数统计:
        策略消费: {results['mean_consumption_policy']:.4f}
        策略储蓄: {results['mean_savings_policy']:.4f}
        策略PPS: {results['mean_pps_policy']:.4f}
        
        VFI求解时间: {results['vfi_time']:.2f}秒"""
        
        axes[1, 2].text(0.05, 0.95, stats_text, transform=axes[1, 2].transAxes, 
                        verticalalignment='top', fontsize=9, 
                        bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
        axes[1, 2].set_title('VFI评估统计摘要')
        
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"📈 VFI评估图表已保存到: {save_path}")
    
    def plot_comparison_simplified(self, results, save_path='./simplified/comparison_simplified.png'):
        """绘制简化版比较结果"""
        print("\n📊 绘制简化版比较结果...")
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle(f'简化版VFI vs RL比较 (固定宏观参数, n={results["n_sim"]})', fontsize=14)
        
        # 1. 效用分布比较
        axes[0, 0].hist(results['lifetime_utility_vfi'], bins=20, alpha=0.7, label='VFI', color='blue')
        axes[0, 0].hist(results['lifetime_utility_rl'], bins=20, alpha=0.7, label='RL', color='red')
        axes[0, 0].axvline(results['mean_utility_vfi'], color='blue', linestyle='--')
        axes[0, 0].axvline(results['mean_utility_rl'], color='red', linestyle='--')
        axes[0, 0].set_xlabel('生涯效用')
        axes[0, 0].set_ylabel('频数')
        axes[0, 0].set_title('效用分布比较')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        # 2. 效用散点图
        axes[0, 1].scatter(results['lifetime_utility_vfi'], results['lifetime_utility_rl'], alpha=0.6)
        min_val = min(np.min(results['lifetime_utility_vfi']), np.min(results['lifetime_utility_rl']))
        max_val = max(np.max(results['lifetime_utility_vfi']), np.max(results['lifetime_utility_rl']))
        axes[0, 1].plot([min_val, max_val], [min_val, max_val], 'r--', linewidth=2)
        axes[0, 1].set_xlabel('VFI生涯效用')
        axes[0, 1].set_ylabel('RL生涯效用')
        axes[0, 1].set_title('个体效用对比')
        axes[0, 1].grid(True, alpha=0.3)
        
        # 3. 平均消费路径
        mean_c_vfi = np.mean(results['c_path_vfi'], axis=0)
        mean_c_rl = np.mean(results['c_path_rl'], axis=0)
        ages = np.arange(20, 20 + len(mean_c_vfi))
        
        axes[1, 0].plot(ages, mean_c_vfi, 'b-', linewidth=2, label='VFI')
        axes[1, 0].plot(ages, mean_c_rl, 'r--', linewidth=2, label='RL')
        axes[1, 0].set_xlabel('年龄')
        axes[1, 0].set_ylabel('平均消费')
        axes[1, 0].set_title('平均消费路径')
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)
        
        # 4. 效用差异分布
        utility_diff_individual = results['lifetime_utility_rl'] - results['lifetime_utility_vfi']
        axes[1, 1].hist(utility_diff_individual, bins=20, alpha=0.7, color='green')
        axes[1, 1].axvline(np.mean(utility_diff_individual), color='darkgreen', linewidth=2)
        axes[1, 1].axvline(0, color='black', linestyle='--', alpha=0.7)
        axes[1, 1].set_xlabel('效用差异 (RL - VFI)')
        axes[1, 1].set_ylabel('频数')
        axes[1, 1].set_title('个体效用差异分布')
        axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"📈 简化版比较图表已保存到: {save_path}")
    
    def __del__(self):
        """清理MATLAB Engine"""
        if hasattr(self, 'eng') and self.eng is not None:
            try:
                self.eng.quit()
            except:
                pass  # 忽略清理时的错误

def main(eva_rl_only=False, eva_vfi_only=False):
    """
    主函数
    
    Args:
        eva_rl_only: 如果为True，只评估RL模型，不运行VFI比较
        eva_vfi_only: 如果为True，只评估VFI模型，不运行RL比较
    """
    print("=" * 60)
    if eva_rl_only:
        print("🤖 简化版RL模型单独评估")
        print("🎯 只评估RL模型性能，跳过VFI比较")
        print("🎮 动作空间：[PPS缴费比例, 消费比例]")
    elif eva_vfi_only:
        print("🔧 简化版VFI模型单独评估")
        print("🎯 只评估VFI模型性能，跳过RL比较")
        print("🎮 策略空间：离散网格搜索")
    else:
        print("🔬 简化版VFI vs RL比较分析")
        print("🎯 固定宏观参数，状态空间与VFI基本相同")
        print("🎮 RL动作空间：[PPS缴费比例, 消费比例] vs VFI离散网格")
    print("=" * 60)
    
    if eva_rl_only:
        # RL单独评估模式
        print("\n1️⃣ 初始化RL评估器（跳过MATLAB）...")
        comparator = RLVFIComparatorSimplified(use_matlab=False)
        
        print("\n2️⃣ 加载简化版RL模型...")
        rl_model, rl_config = comparator.load_rl_model_simplified()
        print("✅ 简化版RL模型加载成功")
        
        print("\n3️⃣ 运行RL单独评估...")
        rl_results = comparator.evaluate_rl_only_simplified(
            rl_model, rl_config, n_sim=100, random_seed=42
        )
        print("✅ RL评估完成")
        
        print("\n4️⃣ 绘制RL评估结果...")
        comparator.plot_rl_only_results(rl_results)
        
        print("\n5️⃣ 保存RL评估结果...")
        with open('./simplified/rl_only_results_simplified.pkl', 'wb') as f:
            pickle.dump(rl_results, f)
        print("✅ 结果保存成功")
        
        # 输出RL评估摘要
        print("\n" + "=" * 60)
        print("📋 简化版RL模型评估摘要")
        print("=" * 60)
        print(f"🎯 状态空间: 4维 (k, k_pps, age, ε)")
        print(f"🔧 宏观参数: 固定，不作为状态变量")
        print(f"📊 模拟个体数: {rl_results['n_sim']}")
        print(f"🎲 随机种子: {rl_results['random_seed']}")
        print(f"📊 平均生涯效用: {rl_results['mean_utility_rl']:.4f} ± {rl_results['std_utility_rl']:.4f}")
        print(f"📊 中位数效用: {rl_results['median_utility_rl']:.4f}")
        print(f"📊 效用范围: [{rl_results['min_utility_rl']:.4f}, {rl_results['max_utility_rl']:.4f}]")
        print(f"💰 平均消费: {rl_results['mean_consumption']:.4f}")
        print(f"💼 平均储蓄: {rl_results['mean_savings']:.4f}")
        print(f"🏦 平均PPS缴费: {rl_results['mean_pps_contrib']:.4f}")
        
        print("\n💡 RL单独评估特性:")
        print("  ✅ 快速评估，无需VFI求解")
        print("  ✅ 详细的RL模型性能分析")
        print("  ✅ 生命周期路径可视化")
        print("  🎮 动作空间：[PPS缴费比例, 消费比例]")
        print("  🚀 适用于RL模型调试和验证")
        
        print("=" * 60)
        
        return rl_results
        
    elif eva_vfi_only:
        # VFI单独评估模式
        print("\n1️⃣ 初始化VFI评估器...")
        comparator = RLVFIComparatorSimplified(use_matlab=True)
        
        print("\n2️⃣ 运行VFI单独评估...")
        vfi_results = comparator.evaluate_vfi_only_simplified(n_sim=100, random_seed=42)
        print("✅ VFI评估完成")
        
        print("\n3️⃣ 绘制VFI评估结果...")
        comparator.plot_vfi_only_results(vfi_results)
        
        print("\n4️⃣ 保存VFI评估结果...")
        with open('./simplified/vfi_only_results_simplified.pkl', 'wb') as f:
            pickle.dump(vfi_results, f)
        print("✅ 结果保存成功")
        
        # 输出VFI评估摘要
        print("\n" + "=" * 60)
        print("📋 简化版VFI模型评估摘要")
        print("=" * 60)
        print(f"🎯 策略空间: 离散网格 (k×k_pps×age×ε)")
        print(f"🔧 宏观参数: 固定")
        print(f"📊 模拟个体数: {vfi_results['n_sim']}")
        print(f"🎲 随机种子: {vfi_results['random_seed']}")
        print(f"📊 平均生涯效用: {vfi_results['mean_utility_vfi']:.4f} ± {vfi_results['std_utility_vfi']:.4f}")
        print(f"📊 中位数效用: {vfi_results['median_utility_vfi']:.4f}")
        print(f"📊 效用范围: [{vfi_results['min_utility_vfi']:.4f}, {vfi_results['max_utility_vfi']:.4f}]")
        print(f"💰 平均消费: {vfi_results['mean_consumption']:.4f}")
        print(f"💼 平均储蓄: {vfi_results['mean_savings']:.4f}")
        print(f"🏦 平均PPS储蓄: {vfi_results['mean_pps_savings']:.4f}")
        print(f"💳 平均PPS缴费: {vfi_results['mean_pps_contrib']:.4f}")
        print(f"⏱️ VFI求解时间: {vfi_results['vfi_time']:.2f}秒")
        
        print("\n💡 VFI单独评估特性:")
        print("  ✅ 快速验证VFI求解结果")
        print("  ✅ 详细的VFI策略函数分析")
        print("  ✅ 生命周期路径可视化")
        print("  🎮 策略空间：离散网格搜索")
        print("  🚀 适用于VFI模型调试和验证")
        print("  🔧 无需RL模型，专注VFI分析")
        
        print("=" * 60)
        
        return vfi_results
        
    else:
        # 完整的VFI vs RL比较模式
        print("\n1️⃣ 初始化简化版比较器...")
        comparator = RLVFIComparatorSimplified(use_matlab=True)
        
        print("\n2️⃣ 加载简化版RL模型...")
        rl_model, rl_config = comparator.load_rl_model_simplified()
        print("✅ 简化版RL模型加载成功")
        
        print("\n3️⃣ 运行VFI（固定宏观参数）...")
        vfi_results = comparator.run_vfi_simplified()
        print("✅ VFI运行成功")
        
        print("\n4️⃣ 运行简化版比较模拟...")
        comparison_results = comparator.simulate_lifecycle_comparison_simplified(
            rl_model, vfi_results, rl_config, n_sim=100, random_seed=42
        )
        print("✅ 简化版比较模拟完成")
        
        print("\n5️⃣ 绘制简化版比较结果...")
        comparator.plot_comparison_simplified(comparison_results)
        
        print("\n6️⃣ 保存简化版比较结果...")
        with open('./simplified/comparison_results_simplified.pkl', 'wb') as f:
            pickle.dump(comparison_results, f)
        print("✅ 结果保存成功")
        
        # 输出VFI vs RL比较摘要
        print("\n" + "=" * 60)
        print("📋 简化版VFI vs RL比较摘要")
        print("=" * 60)
        print(f"🎯 状态空间: 4维 (k, k_pps, age, ε)")
        print(f"🔧 宏观参数: 固定，不作为状态变量")
        print(f"📊 模拟个体数: {comparison_results['n_sim']}")
        print(f"🎲 随机种子: {comparison_results['random_seed']}")
        print(f"📊 VFI平均效用: {comparison_results['mean_utility_vfi']:.4f} ± {comparison_results['std_utility_vfi']:.4f}")
        print(f"📊 RL平均效用: {comparison_results['mean_utility_rl']:.4f} ± {comparison_results['std_utility_rl']:.4f}")
        print(f"🔍 效用差异(RL-VFI): {comparison_results['utility_diff']:.6f}")
        print(f"📈 相对差异: {comparison_results['utility_diff']/abs(comparison_results['mean_utility_vfi'])*100:.3f}%")
        
        utility_diff = comparison_results['utility_diff']
        if abs(utility_diff) < 0.001:
            print("🎉 RL与VFI结果非常接近！")
        elif utility_diff > 0:
            print("🏆 RL表现优于VFI")
        else:
            print("🏆 VFI表现优于RL")
        
        print("\n💡 简化版特性:")
        print("  ✅ 状态空间与VFI基本相同，更公平的比较")
        print("  ✅ 固定宏观参数，消除宏观不确定性影响")
        print("  ✅ 保持累积存活概率方法的理论等价性")
        print("  🎮 RL动作空间：[PPS缴费比例, 消费比例] vs VFI离散网格")
        print("  🎯 验证RL在相同条件下是否能达到VFI的结果")
        
        print("=" * 60)
        
        return comparison_results
    
    # 清理
    if 'comparator' in locals():
        del comparator

if __name__ == "__main__":
    # 可以通过命令行参数或直接修改这里来控制评估模式
    import sys
    
    # 检查命令行参数
    eva_rl_only = True
    eva_vfi_only = False
    
    if len(sys.argv) > 1:
        if sys.argv[1] == '--rl-only':
            eva_rl_only = True
        elif sys.argv[1] == '--vfi-only':
            eva_vfi_only = True
    
    results = main(eva_rl_only=eva_rl_only, eva_vfi_only=eva_vfi_only)

# 便捷调用函数
def evaluate_rl_only():
    """便捷函数：仅评估RL模型"""
    return main(eva_rl_only=True, eva_vfi_only=False)

def evaluate_vfi_only():
    """便捷函数：仅评估VFI模型"""
    return main(eva_rl_only=False, eva_vfi_only=True)

def compare_rl_vfi():
    """便捷函数：完整VFI vs RL比较"""
    return main(eva_rl_only=False, eva_vfi_only=False) 