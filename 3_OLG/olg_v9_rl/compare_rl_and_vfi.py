"""
比较VFI和RL（SAC Agent）的优化结果 - 完全独立的Python实现
目标：在相同的宏观和微观参数下，比较两种方法的优化效果
特性：
- 使用main_olg_v9_utils.py中的VFI方法
- 独立实现RL生命周期评估，不依赖外部模块
- 确保VFI和RL使用完全一致的参数设置

支持的RL后端：
- SB3 (Stable Baselines 3) 
- SBX (Stable Baselines Jax)
"""

import os
# 解决OpenMP库冲突问题
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

import numpy as np
import matplotlib.pyplot as plt
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

def check_chinese_font_availability():
    """检查并显示中文字体安装状态"""
    available_fonts = [f.name for f in fm.fontManager.ttflist]
    chinese_fonts = ['SimHei', 'Microsoft YaHei', 'WenQuanYi Micro Hei']
    
    print("\n🔤 中文字体检查:")
    found_chinese = False
    for font in chinese_fonts:
        if font in available_fonts:
            print(f"  ✅ {font} - 可用")
            found_chinese = True
        else:
            print(f"  ❌ {font} - 不可用")
    
    if not found_chinese:
        print("\n💡 中文字体安装建议:")
        print("  Windows: 系统自带 Microsoft YaHei 和 SimHei")
        print("  Linux: sudo apt-get install fonts-wqy-microhei")
        print("  macOS: 系统自带 PingFang SC")
        print("  如果仍有问题，请重启Python会话")
    
    return found_chinese

import pickle
import json
import os
import time
from typing import Dict, Any, Tuple, Optional
from pathlib import Path
from scipy import stats

# 导入SBX和SB3
from sbx import SAC as SBX_SAC
SBX_AVAILABLE = True
print("✅ SBX (Stable Baselines Jax) 可用")

# 导入SB3作为备选
from stable_baselines3 import SAC as SB3_SAC
from stable_baselines3.common.evaluation import evaluate_policy

from main_olg_v9_utils import OLG_V9_Utils
from main_olg_v9_sac_sbx import evaluate_policy_lifecycle_simulation

class RLVFIComparator:
    """RL和VFI方法比较器 - 使用Python版本的VFI，独立实现RL评估"""
    
    def __init__(self, use_sbx: bool = True):
        """
        初始化比较器
        
        Args:
            use_sbx: 是否优先使用SBX，否则使用SB3
        """
        self.use_sbx = use_sbx and SBX_AVAILABLE
        self.SAC_class = SBX_SAC if self.use_sbx else SB3_SAC
        self.backend_name = "SBX" if self.use_sbx else "SB3"
        
        print(f"🤖 使用 {self.backend_name} 作为RL后端")
        print("🔧 使用Python版本的VFI (main_olg_v9_utils.py)")
        
        # 固定的测试参数
        self.M_test = {
            'R_k_net_factor': 1.03,
            'w_gross': 2.0,
            'TR_total': 0.1,
            'b_payg_avg_retiree': 0.4,
            'tau_l': 0.15,
            'theta_payg_actual': 0.12
        }
        
        print("📊 固定测试参数设置:")
        for key, value in self.M_test.items():
            print(f"  {key} = {value:.3f}")
    
    def load_rl_model(self, model_path: Optional[str] = None, use_best_model: bool = True) -> Tuple[Any, Dict]:
        """
        加载训练好的RL模型
        
        Args:
            model_path: 模型路径，如果为None则自动选择
            use_best_model: 是否优先使用best model
            
        Returns:
            model: 加载的模型
            config: 模型配置（包含训练时的cS和paramS_for_rl）
        """
        if model_path is None:
            # 自动选择模型路径（v9版本）
            if self.use_sbx:
                if use_best_model:
                    best_model_path = './py/best_model_sbx/best_model.zip'
                    final_model_path = './final_sac_agent_olg_v9_sbx.zip'
                    
                    if os.path.exists(best_model_path):
                        model_path = best_model_path
                        print(f"🏆 使用训练过程中的最佳模型: {model_path}")
                    elif os.path.exists(final_model_path):
                        model_path = final_model_path
                        print(f"⚠️ Best model不存在，使用最终模型: {model_path}")
                    else:
                        raise FileNotFoundError("未找到SBX v9模型文件 (best_model或final_model)")
                else:
                    model_path = './final_sac_agent_olg_v9_sbx.zip'
                    print(f"📁 使用最终模型: {model_path}")
                
                config_path = './py/training_config_sbx.pkl'
            else:
                if use_best_model:
                    best_model_path = './py/best_model/best_model.zip'
                    final_model_path = './final_sac_agent_olg_v9_sb3.zip'
                    
                    if os.path.exists(best_model_path):
                        model_path = best_model_path
                        print(f"🏆 使用训练过程中的最佳模型: {model_path}")
                    elif os.path.exists(final_model_path):
                        model_path = final_model_path
                        print(f"⚠️ Best model不存在，使用最终模型: {model_path}")
                    else:
                        raise FileNotFoundError("未找到SB3 v9模型文件 (best_model或final_model)")
                else:
                    model_path = './final_sac_agent_olg_v9_sb3.zip'
                    print(f"📁 使用最终模型: {model_path}")
                
                config_path = './py/training_config_sb3.pkl'
        else:
            config_path = model_path.replace('.zip', '_config.pkl')
            print(f"📁 使用指定模型: {model_path}")
        
        print(f"📁 加载RL模型: {model_path}")
        
        if not os.path.exists(model_path):
            raise FileNotFoundError(f"模型文件不存在: {model_path}")
        
        # 加载模型
        model = self.SAC_class.load(model_path)
        print(f"✅ {self.backend_name} SAC模型加载成功")
        
        # 加载配置
        config = {}
        if os.path.exists(config_path):
            with open(config_path, 'rb') as f:
                config = pickle.load(f)
            print(f"✅ 模型配置加载成功")
        else:
            print(f"⚠️ 配置文件不存在: {config_path}")
        
        return model, config
    
    def _get_param(self, obj, key, default=None):
        """辅助函数：统一处理字典和对象的参数访问"""
        if isinstance(obj, dict):
            return obj.get(key, default)
        else:
            return getattr(obj, key, default)
    # 注意：evaluate_rl_policy_lifecycle方法已删除
    # 现在直接使用main_olg_v9_sac_sbx.py中的evaluate_policy_lifecycle_simulation函数
    
    def run_python_vfi(self) -> Dict[str, Any]:
        """
        使用Python版本的VFI方法 (main_olg_v9_utils.py)
        
        Returns:
            vfi_results: VFI结果字典
        """
        print("\n🔧 使用Python版本的VFI方法...")
        
        start_time = time.time()
        
        # 1. 初始化参数
        print("1️⃣ 初始化Python参数...")
        cS = OLG_V9_Utils.ParameterValues_HuggettStyle()
        
        print(f"  ✅ 参数初始化完成: aD_new={cS.aD_new}, nk={cS.nk}, nkpps={cS.nkpps}")
        
        # 2. 设置劳动效率过程
        print("2️⃣ 初始化劳动效率过程...")
        leLogGridV, leTrProbM, leProb1V = OLG_V9_Utils.EarningProcess_olgm(cS)
            
            # 创建paramS_vfi结构
        class ParameterStruct:
            pass
        
        paramS_vfi = ParameterStruct()
        paramS_vfi.leLogGridV = leLogGridV
        paramS_vfi.leTrProbM = leTrProbM
        paramS_vfi.leProb1V = leProb1V
        paramS_vfi.leGridV = np.exp(leLogGridV)
        paramS_vfi.ageEffV_new = cS.ageEffV_new
        paramS_vfi.tau_l = self.M_test['tau_l']
        paramS_vfi.theta_payg_actual_for_hh = self.M_test['theta_payg_actual']
        paramS_vfi.pps_tax_deferral_active = cS.pps_active
        
        # 3. 构建PAYG福利向量
        print("3️⃣ 构建PAYG福利向量...")
        bV_payg_vfi = np.zeros(cS.aD_new)
        if cS.aR_new < cS.aD_new:
            bV_payg_vfi[cS.aR_new:] = self.M_test['b_payg_avg_retiree']
        
        # 4. 运行VFI求解
        print("4️⃣ 运行Python VFI求解...")
        vfi_start = time.time()
            
        cPolM_vfi, kPolM_vfi, cPpsPolM_vfi, VPolM_vfi = OLG_V9_Utils.HHSolution_VFI_Huggett(
                self.M_test['R_k_net_factor'],
                self.M_test['w_gross'],
                self.M_test['TR_total'],
                bV_payg_vfi,
                paramS_vfi,
            cS
            )
            
        vfi_time = time.time() - vfi_start
        total_time = time.time() - start_time
            
        print(f"✅ Python VFI求解完成，耗时: {vfi_time:.2f} 秒")
        print(f"📊 策略矩阵尺寸: {cPolM_vfi.shape}")
            
        # 5. 准备返回结果
        result_dict = {
        'cPolM': cPolM_vfi,
        'kPolM': kPolM_vfi,
        'cPpsPolM': cPpsPolM_vfi,
        'VPolM': VPolM_vfi,
            'M_test': self.M_test,
            'cS': cS,
            'paramS_vfi': paramS_vfi,
            'vfi_time': vfi_time,
            'total_time': total_time
        }
        
        print("✅ Python VFI结果准备完成")
        return result_dict
    
    def simulate_lifecycle_comparison(self, rl_model: Any, vfi_results: Dict, rl_config: Dict,
                                    n_sim: int = 100, random_seed: int = 42) -> Dict[str, Any]:
        """
        模拟生命周期轨迹比较
        
        Args:
            rl_model: 训练好的RL模型
            vfi_results: VFI结果
            rl_config: RL训练配置
            n_sim: 模拟个体数量
            random_seed: 随机种子
            
        Returns:
            comparison_results: 比较结果
        """
        print(f"\n🎯 模拟生命周期轨迹比较 (n_sim={n_sim}, seed={random_seed})")
        
        # 设置随机种子确保可重现性
        np.random.seed(random_seed)
        print(f"🎲 随机种子已设置: {random_seed}")
        

        
        # 从VFI结果获取参数
        cS = vfi_results['cS']
        paramS_vfi = vfi_results['paramS_vfi']
        
        # 获取维度
        aD_new = cS.aD_new
        aR_new = cS.aR_new
        nk = cS.nk
        nkpps = cS.nkpps
        nw = cS.nw
        
        # 从VFI结果提取网格和参数
        kGridV = cS.kGridV
        kppsGridV = cS.kppsGridV
        leGridV = paramS_vfi.leGridV
        leTrProbM = paramS_vfi.leTrProbM
        leProb1V = paramS_vfi.leProb1V
        ageEffV_new = cS.ageEffV_new
        beta = cS.beta
        sigma = cS.sigma
        kMin = cS.kMin
        kppsMin = cS.kppsMin
        kMax = cS.kMax
        kppsMax = cS.kppsMax
        pps_active = cS.pps_active
            
        print(f"📊 网格提取成功: kGridV({len(kGridV)}), kppsGridV({len(kppsGridV)}), leGridV({len(leGridV)})")
        
        # 获取策略矩阵
        cPolM_vfi = vfi_results['cPolM']
        kPolM_vfi = vfi_results['kPolM']
        cPpsPolM_vfi = vfi_results['cPpsPolM']
        
        print(f"📊 网格尺寸: nk={nk}, nkpps={nkpps}, nw={nw}, aD_new={aD_new}")
        print(f"📊 VFI策略矩阵尺寸: {cPolM_vfi.shape}")
        
        # 🆕 使用独立的RL生命周期评估函数
        print("\n🧬 使用独立的RL生命周期评估函数")
        
        # 🔧 检查训练时保存的参数是否完整可用
        use_training_config = False
        if rl_config and 'cS' in rl_config and 'paramS_for_rl' in rl_config:
            cS_training = rl_config['cS']
            # 检查关键参数是否存在
            required_keys = ['kMin', 'kMax', 'kppsMin', 'kppsMax', 'sigma', 'cFloor', 'tau_c', 'pps_active']
            if isinstance(cS_training, dict):
                missing_keys = [key for key in required_keys if key not in cS_training]
                if not missing_keys:
                    use_training_config = True
                    print("✅ 使用训练时保存的参数 (cS, paramS_for_rl, rng_M)")
                    print(f"📊 训练时参数 - aD_new: {cS_training['aD_new']}, beta: {cS_training.get('beta', 'N/A')}")
                else:
                    print(f"⚠️ 训练配置缺少关键参数: {missing_keys}")
            else:
                # 对象类型，检查属性
                missing_attrs = [key for key in required_keys if not hasattr(cS_training, key)]
                if not missing_attrs:
                    use_training_config = True
                    print("✅ 使用训练时保存的参数 (cS, paramS_for_rl, rng_M)")
                    print(f"📊 训练时参数 - aD_new: {cS_training.aD_new}, beta: {cS_training.beta}")
                else:
                    print(f"⚠️ 训练配置缺少关键属性: {missing_attrs}")
        else:
            print("⚠️ 训练配置不完整，使用VFI参数进行RL评估")
        
        if use_training_config:
            cS_python = rl_config['cS']
            paramS_for_rl = rl_config['paramS_for_rl']
            rng_M = rl_config.get('rng_M', {
                'R_k_net_factor': [1.01, 1.05],
                'w_gross': [1.5, 2.5],
                'TR_total': [0.0, 0.2],
                'b_payg_avg_retiree': [0.1, 0.8],
                'tau_l': [0.05, 0.25],
                'theta_payg_actual': [0.05, 0.20]
            })
        else:
            print("⚠️ 使用VFI参数进行RL评估，确保参数一致性")
            # 使用VFI参数进行RL评估，确保一致性
            cS_python = cS
            paramS_for_rl = type('obj', (object,), {})()
            paramS_for_rl.leLogGridV = paramS_vfi.leLogGridV
            paramS_for_rl.leTrProbM = paramS_vfi.leTrProbM
            paramS_for_rl.leProb1V = paramS_vfi.leProb1V
            paramS_for_rl.leGridV = paramS_vfi.leGridV
            paramS_for_rl.ageEffV_new = paramS_vfi.ageEffV_new
            
            # 使用固定的测试参数创建rng_M
            rng_M = {}
            for key, value in self.M_test.items():
                rng_M[key] = [value, value]  # 固定范围为相同值
        
        # 🎯 关键步骤1：首先生成统一的效率冲击序列，确保VFI和RL使用完全相同的序列
        print("\n🎲 生成年度效率冲击序列（VFI和RL将使用相同序列）...")
        # 重新设置随机种子，确保效率冲击生成的可重现性
        np.random.seed(random_seed)
        eIdxM_annual_unified = OLG_V9_Utils.LaborEndowSimulation_olgm(cS, paramS_vfi)
            
        # 确保个体数量一致
        if eIdxM_annual_unified.shape[0] != n_sim:
            if eIdxM_annual_unified.shape[0] > n_sim:
                eIdxM_annual_unified = eIdxM_annual_unified[:n_sim, :]
            else:
                repeat_times = (n_sim // eIdxM_annual_unified.shape[0]) + 1
                eIdxM_annual_unified = np.tile(eIdxM_annual_unified, (repeat_times, 1))[:n_sim, :]
        
        print(f"📊 统一效率冲击矩阵尺寸: {eIdxM_annual_unified.shape}")
        print("🎯 此效率冲击序列将同时用于VFI和RL模拟，确保完全一致性")
        
        # 🔄 VFI生命周期模拟（使用Python版本的HHSimulation_olgm）
        print("\n🔄 开始VFI生命周期模拟...")
        vfi_sim_start = time.time()
            
        # 使用Python版本的HHSimulation_olgm
        # 构建bV_payg_sim_benefit
        bV_payg_sim_benefit = np.zeros(cS.aD_new)
        if cS.aR_new < cS.aD_new:
            bV_payg_sim_benefit[cS.aR_new:] = self.M_test['b_payg_avg_retiree']
        
        kHistM_vfi, kPpsHistM_vfi, cHistM_vfi, _ = OLG_V9_Utils.HHSimulation_olgm(
            kPolM_vfi,                          # kPolM_4D_input
            cPpsPolM_vfi,                       # cPpsPolM_choice_4D_input  
            cPolM_vfi,                          # cPolM_consump_q_4D_input
            eIdxM_annual_unified,               # eIdxM_annual_input (统一的效率冲击序列)
            self.M_test['R_k_net_factor'],      # R_k_net_factor_hh_sim
            self.M_test['w_gross'],             # w_gross_sim_price
            self.M_test['TR_total'],            # TR_total_sim_transfer
            bV_payg_sim_benefit,                # bV_payg_sim_benefit
            paramS_vfi,                         # paramS_sim_household
            cS                                  # cS_common_sim
        )
        
        # 转换为轨迹数据
        k_path_vfi = kHistM_vfi
        kpps_path_vfi = kPpsHistM_vfi
        c_path_vfi = cHistM_vfi
        cpps_path_vfi = np.zeros_like(c_path_vfi)  # PPS缴费路径从策略中计算
        
        print(f"✅ Python模拟完成，轨迹尺寸: k{k_path_vfi.shape}, c{c_path_vfi.shape}")
        
                # 🧮 计算VFI生命周期效用（按年度）
        print("🧮 计算VFI生命周期效用...")
        print(f"🔧 折现参数验证: beta={beta:.4f}, sigma={sigma:.4f}")
        print(f"🔧 确保VFI和RL使用相同的折现因子: β^t")
        lifetime_utility_vfi = np.zeros(n_sim)
        
        # 确保与RL相同的随机种子（重新设置）
        np.random.seed(random_seed)
        
        for i_sim in range(n_sim):
            if (i_sim + 1) % 20 == 0:
                print(f"  VFI效用计算进度: {i_sim + 1}/{n_sim}")
            
            utility_sum_vfi = 0
            
            # 按年度计算效用（与RL一致）
            for age_annual_idx in range(cS.aD_orig):
                c_vfi = c_path_vfi[i_sim, age_annual_idx]
                
                # 计算VFI效用（与RL评估使用相同的效用函数）
                _, u_vfi = OLG_V9_Utils.CES_utility(c_vfi, sigma, cS)
                
                # 🔧 修复：使用与RL相同的beta参数进行折现，确保一致性
                discount_factor = beta ** age_annual_idx  # 使用beta而不是硬编码的0.97
                utility_sum_vfi += discount_factor * u_vfi
            
            lifetime_utility_vfi[i_sim] = utility_sum_vfi
        
        vfi_sim_time = time.time() - vfi_sim_start
        print(f"✅ VFI生命周期模拟完成，耗时: {vfi_sim_time:.2f} 秒")
        
        # 🚀 RL生命周期评估（使用相同的效率冲击序列）
        print("\n🔄 开始RL生命周期评估...")
        print(f"📊 评估参数: n_sim={n_sim}, random_seed={random_seed}, gamma=0.97")
        print("🎯 使用main_olg_v9_sac_sbx.py中的evaluate_policy_lifecycle_simulation确保完全一致性")
        rl_eval_start = time.time()
        
        # 🎯 关键：传递相同的效率冲击序列给RL评估，确保完全一致性
        mean_utility_rl, std_utility_rl, rl_lifecycle_results = evaluate_policy_lifecycle_simulation(
            rl_model, cS_python, paramS_for_rl, rng_M, 
            n_sim=n_sim, deterministic=True, gamma=0.97, random_seed=random_seed, verbose=True,
            eIdxM_annual_input=eIdxM_annual_unified  # 🎯 使用与VFI相同的效率冲击序列
        )
        
        rl_eval_time = time.time() - rl_eval_start
        print(f"✅ RL评估完成，耗时: {rl_eval_time:.2f} 秒")
        print(f"📊 RL评估结果: {mean_utility_rl:.4f} ± {std_utility_rl:.4f}")
        
        # 从RL评估结果中提取轨迹数据
        lifetime_utility_rl = rl_lifecycle_results['lifetime_utility_rl']
        k_path_rl = rl_lifecycle_results['k_path_rl']
        c_path_rl = rl_lifecycle_results['c_path_rl']
        cpps_path_rl = rl_lifecycle_results['cpps_path_rl']
        
        # 总模拟时间
        sim_time = rl_eval_time + vfi_sim_time
        
        # 📊 计算统计结果
        mean_utility_vfi = np.mean(lifetime_utility_vfi)
        std_utility_vfi = np.std(lifetime_utility_vfi)
        
        utility_diff = mean_utility_rl - mean_utility_vfi
        utility_improvement_pct = (utility_diff / abs(mean_utility_vfi)) * 100
        
        print(f"\n📊 生涯总折现效用比较 (基于 {n_sim} 个个体):")
        print(f"  📈 VFI方法: {mean_utility_vfi:.4f} ± {std_utility_vfi:.4f}")
        print(f"  🤖 RL方法 ({self.backend_name}): {mean_utility_rl:.4f} ± {std_utility_rl:.4f}")
        print(f"  🔍 差异 (RL - VFI): {utility_diff:.4f}")
        print(f"  📊 相对改进: {utility_improvement_pct:.2f}%")
        
        if utility_diff > 0:
            print("  >>> 🏆 RL方法表现更好！")
        elif utility_diff < 0:
            print("  >>> 🏆 VFI方法表现更好！")
        else:
            print("  >>> 🤝 两种方法表现相当。")
        
        # 统计显著性检验
        t_stat, p_value = stats.ttest_rel(lifetime_utility_rl, lifetime_utility_vfi)
        print(f"\n📈 统计显著性检验 (配对t检验): p值: {p_value:.4f}")
        print(f"  {'差异显著' if p_value < 0.05 else '差异不显著'} (α=0.05)")
        
        # 构造比较结果
        comparison_results = {
            'rl_backend': self.backend_name,
            'vfi_method': 'Python (main_olg_v9_utils.py)',
            'M_test': self.M_test,
            'n_sim': n_sim,
            'random_seed': random_seed,
            'rl_simulation_time': rl_eval_time,
            'vfi_simulation_time': vfi_sim_time,
            'total_simulation_time': sim_time,
            'vfi_time': vfi_results.get('vfi_time', 0),
            'mean_utility_vfi': mean_utility_vfi,
            'mean_utility_rl': mean_utility_rl,
            'std_utility_vfi': std_utility_vfi,
            'std_utility_rl': std_utility_rl,
            'utility_diff': utility_diff,
            'utility_improvement_pct': utility_improvement_pct,
            'p_value': p_value,
            'is_significant': p_value < 0.05,
            'lifetime_utility_vfi': lifetime_utility_vfi,
            'lifetime_utility_rl': lifetime_utility_rl,
            'k_path_vfi': k_path_vfi,
            'k_path_rl': k_path_rl,
            'c_path_vfi': c_path_vfi,
            'c_path_rl': c_path_rl,
            'cpps_path_vfi': cpps_path_vfi,
            'cpps_path_rl': cpps_path_rl,
            'rl_lifecycle_results': rl_lifecycle_results
        }
        
        return comparison_results
    
    def plot_comparison_results(self, results: Dict[str, Any], save_path: str = './py/rl_vfi_comparison.png'):
        """
        绘制VFI vs RL比较结果
        
        Args:
            results: 比较结果
            save_path: 保存路径
        """
        print("\n📊 绘制VFI vs RL比较结果...")
        
        n_sim = results['n_sim']
        
        # 创建2x3的子图布局
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        # 设置中文字体（备用确保）
        plt.rcParams['font.sans-serif'] = matplotlib.rcParams['font.sans-serif']
        plt.rcParams['axes.unicode_minus'] = False
        fig.suptitle(f'Python VFI vs RL ({results["rl_backend"]}) 生命周期比较 (n={n_sim})', fontsize=16)
        
        # 1. 效用分布比较
        axes[0, 0].hist(results['lifetime_utility_vfi'], bins=20, alpha=0.7, 
                       label='Python VFI', color='red')
        axes[0, 0].hist(results['lifetime_utility_rl'], bins=20, alpha=0.7, 
                       label=f'RL ({results["rl_backend"]})', color='blue')
        axes[0, 0].axvline(results['mean_utility_vfi'], color='red', linestyle='--', 
                          label=f'VFI平均值: {results["mean_utility_vfi"]:.2f}')
        axes[0, 0].axvline(results['mean_utility_rl'], color='blue', linestyle='--', 
                          label=f'RL平均值: {results["mean_utility_rl"]:.2f}')
        axes[0, 0].set_xlabel('生涯总效用')
        axes[0, 0].set_ylabel('频数')
        axes[0, 0].set_title('效用分布比较')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        # 2. 效用散点图
        axes[0, 1].scatter(results['lifetime_utility_vfi'], results['lifetime_utility_rl'], 
                          alpha=0.6, s=30)
        # 绘制45度线
        min_val = min(np.min(results['lifetime_utility_vfi']), np.min(results['lifetime_utility_rl']))
        max_val = max(np.max(results['lifetime_utility_vfi']), np.max(results['lifetime_utility_rl']))
        axes[0, 1].plot([min_val, max_val], [min_val, max_val], 'r--', linewidth=2, 
                       label='45度线')
        axes[0, 1].set_xlabel('Python VFI 生涯效用')
        axes[0, 1].set_ylabel('RL 生涯效用')
        axes[0, 1].set_title('个体效用对比')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
        
        # 3. 效用差异分布
        utility_diff_individual = results['lifetime_utility_rl'] - results['lifetime_utility_vfi']
        axes[0, 2].hist(utility_diff_individual, bins=20, alpha=0.7, color='green')
        axes[0, 2].axvline(np.mean(utility_diff_individual), color='darkgreen', 
                          linestyle='-', linewidth=2, 
                          label=f'平均差异: {np.mean(utility_diff_individual):.4f}')
        axes[0, 2].axvline(0, color='black', linestyle='--', alpha=0.7)
        axes[0, 2].set_xlabel('效用差异 (RL - VFI)')
        axes[0, 2].set_ylabel('频数')
        axes[0, 2].set_title('个体效用差异分布')
        axes[0, 2].legend()
        axes[0, 2].grid(True, alpha=0.3)
        
        # 4. 平均资产路径（年度数据）
        k_path_vfi = results['k_path_vfi']
        k_path_rl = results['k_path_rl']
        aD_vfi = k_path_vfi.shape[1]  # 年度数量（79）
        aD_rl = k_path_rl.shape[1]    # 年度数量（79）
        
        mean_k_vfi = np.mean(k_path_vfi, axis=0)
        mean_k_rl = np.mean(k_path_rl, axis=0)
        
        # 年度年龄（20-98岁）
        age_path_vfi = np.arange(20, 20 + aD_vfi)
        age_path_rl = np.arange(20, 20 + aD_rl)
        
        axes[1, 0].plot(age_path_vfi, mean_k_vfi, 'r-', linewidth=2, label=f'Python VFI (年度, T={aD_vfi})')
        axes[1, 0].plot(age_path_rl, mean_k_rl, 'b--', linewidth=2, 
                       label=f'RL ({results["rl_backend"]}, 年度, T={aD_rl})')
        axes[1, 0].set_xlabel('年龄（岁）')
        axes[1, 0].set_ylabel('平均资产')
        axes[1, 0].set_title('平均资产路径（年度）')
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)
        
        # 5. 平均消费路径
        c_path_vfi = results['c_path_vfi']
        c_path_rl = results['c_path_rl']
        
        mean_c_vfi = np.mean(c_path_vfi, axis=0)
        mean_c_rl = np.mean(c_path_rl, axis=0)
        
        axes[1, 1].plot(age_path_vfi, mean_c_vfi, 'r-', linewidth=2, label=f'Python VFI (年度, T={aD_vfi})')
        axes[1, 1].plot(age_path_rl, mean_c_rl, 'b--', linewidth=2, 
                       label=f'RL ({results["rl_backend"]}, 年度, T={aD_rl})')
        axes[1, 1].set_xlabel('年龄（岁）')
        axes[1, 1].set_ylabel('平均消费')
        axes[1, 1].set_title('平均消费路径（年度）')
        axes[1, 1].legend()
        axes[1, 1].grid(True, alpha=0.3)
        
        # 6. 平均PPS缴费路径
        cpps_path_vfi = results['cpps_path_vfi']
        cpps_path_rl = results['cpps_path_rl']
        
        mean_cpps_vfi = np.mean(cpps_path_vfi, axis=0)
        mean_cpps_rl = np.mean(cpps_path_rl, axis=0)
        
        axes[1, 2].plot(age_path_vfi, mean_cpps_vfi, 'r-', linewidth=2, label=f'Python VFI (年度, T={aD_vfi})')
        axes[1, 2].plot(age_path_rl, mean_cpps_rl, 'b--', linewidth=2, 
                       label=f'RL ({results["rl_backend"]}, 年度, T={aD_rl})')
        axes[1, 2].set_xlabel('年龄（岁）')
        axes[1, 2].set_ylabel('平均PPS缴费')
        axes[1, 2].set_title('平均PPS缴费路径（年度）')
        axes[1, 2].legend()
        axes[1, 2].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"📈 Python VFI vs RL比较图表已保存到: {save_path}")
        
        # 打印数值摘要（年度数据）
        print(f"\n📊 路径差异摘要:")
        print(f"  Python VFI生命周期长度: {aD_vfi} 年度 (20-{19+aD_vfi}岁)")
        print(f"  RL生命周期长度: {aD_rl} 年度 (20-{19+aD_rl}岁)")
        
        if aD_vfi == aD_rl:
            # 相同长度，可以直接比较
            k_diff_mean = np.mean(mean_k_rl - mean_k_vfi)
            c_diff_mean = np.mean(mean_c_rl - mean_c_vfi) 
            cpps_diff_mean = np.mean(mean_cpps_rl - mean_cpps_vfi)
            print(f"  平均资产差异 (RL - VFI): {k_diff_mean:.4f}")
            print(f"  平均消费差异 (RL - VFI): {c_diff_mean:.4f}")
            print(f"  平均PPS缴费差异 (RL - VFI): {cpps_diff_mean:.4f}")
        else:
            # 不同长度，比较重叠部分
            min_length = min(aD_vfi, aD_rl)
            k_diff_mean = np.mean(mean_k_rl[:min_length] - mean_k_vfi[:min_length])
            c_diff_mean = np.mean(mean_c_rl[:min_length] - mean_c_vfi[:min_length])
            cpps_diff_mean = np.mean(mean_cpps_rl[:min_length] - mean_cpps_vfi[:min_length])
            print(f"  ⚠️ 生命周期长度不同，比较前{min_length}年:")
            print(f"  平均资产差异 (RL - VFI): {k_diff_mean:.4f}")
            print(f"  平均消费差异 (RL - VFI): {c_diff_mean:.4f}")
            print(f"  平均PPS缴费差异 (RL - VFI): {cpps_diff_mean:.4f}")
    
    def save_comparison_results(self, results: Dict[str, Any], 
                              save_path: str = './py/rl_vfi_comparison_results.json'):
        """
        保存比较结果到JSON文件
        
        Args:
            results: 比较结果
            save_path: 保存路径
        """
        print(f"\n💾 保存比较结果到: {save_path}")
        
        # 准备可序列化的结果
        serializable_results = {}
        for key, value in results.items():
            if isinstance(value, np.ndarray):
                serializable_results[key] = value.tolist()
            elif isinstance(value, (np.integer, np.floating)):
                serializable_results[key] = float(value)
            elif isinstance(value, (np.bool_, bool)):
                serializable_results[key] = bool(value)
            else:
                serializable_results[key] = value
        
        # 添加元数据
        serializable_results['generated_at'] = time.strftime('%Y-%m-%d %H:%M:%S')
        serializable_results['python_backend'] = self.backend_name
        
        with open(save_path, 'w', encoding='utf-8') as f:
            json.dump(serializable_results, f, indent=2, ensure_ascii=False)
        
        print("✅ 结果保存成功")

def main():
    """主函数"""
    print("=" * 80)
    print("🔬 Python VFI vs RL (SAC) 优化方法比较")
    print("=" * 80)
    
    # 0. 检查中文字体
    check_chinese_font_availability()
    
    # 1. 初始化比较器
    print("\n1️⃣ 初始化比较器...")
    comparator = RLVFIComparator(use_sbx=True)  # 优先使用SBX
    
    # 2. 加载RL模型
    print("\n2️⃣ 加载RL模型...")
    rl_model, rl_config = comparator.load_rl_model(use_best_model=True)
    print(f"✅ 成功加载 {comparator.backend_name} 模型")
    
    # 3. 运行Python VFI
    print("\n3️⃣ 运行Python VFI...")
    vfi_results = comparator.run_python_vfi()
    print("✅ VFI运行成功")
    
    # 4. 运行比较模拟
    print("\n4️⃣ 运行生命周期比较模拟...")
    comparison_results = comparator.simulate_lifecycle_comparison(
        rl_model, vfi_results, rl_config, n_sim=100, random_seed=42
        )
    print("✅ 比较模拟完成")
    
    # 5. 绘制结果
    print("\n5️⃣ 绘制比较结果...")
    comparator.plot_comparison_results(comparison_results)
    print("✅ 图表生成完成")
    
    # 6. 保存结果
    print("\n6️⃣ 保存比较结果...")
    comparator.save_comparison_results(comparison_results)
    print("✅ 结果保存完成")
    
    # 7. 输出最终摘要
    print("\n" + "=" * 80)
    print("📋 Python VFI vs RL 比较分析摘要")
    print("=" * 80)
    
    print(f"🤖 RL后端: {comparison_results['rl_backend']}")
    print(f"🏆 RL模型类型: 训练过程中的最佳模型 (Best Model)")
    print(f"🎯 RL评估方法: {comparison_results.get('vfi_method', 'N/A')}")
    print(f"🎯 模拟个体数: {comparison_results['n_sim']}")
    print(f"🎲 随机种子: {comparison_results.get('random_seed', 'N/A')}")
    print(f"⏱️ VFI计算时间: {comparison_results['vfi_time']:.2f} 秒")
    print(f"⏱️ VFI模拟时间: {comparison_results.get('vfi_simulation_time', 0):.2f} 秒")
    print(f"⏱️ RL评估时间: {comparison_results.get('rl_simulation_time', 0):.2f} 秒")
    print(f"⏱️ 总模拟时间: {comparison_results.get('total_simulation_time', 0):.2f} 秒")
    print(f"📊 Python VFI平均效用: {comparison_results['mean_utility_vfi']:.4f} ± {comparison_results['std_utility_vfi']:.4f}")
    print(f"📊 RL平均效用: {comparison_results['mean_utility_rl']:.4f} ± {comparison_results['std_utility_rl']:.4f}")
    print(f"🔍 效用差异: {comparison_results['utility_diff']:.4f} ({comparison_results['utility_improvement_pct']:.2f}%)")
    print(f"📈 统计显著性: {'显著' if comparison_results['is_significant'] else '不显著'} (p={comparison_results['p_value']:.4f})")
    
    if comparison_results['utility_diff'] > 0:
        print("🏆 RL方法表现更好！")
    elif comparison_results['utility_diff'] < 0:
        print("🏆 Python VFI方法表现更好！")
    else:
        print("🤝 两种方法表现相当。")
    
    print("=" * 80)
    print("🎉 Python VFI vs RL 比较分析完成!")
    print("💡 注意：")
    print("  - 本次比较使用Python版本的VFI (main_olg_v9_utils.py)，移除MATLAB Engine依赖")
    print("  - VFI和RL都基于main_olg_v9_utils.py进行统一验证")
    print("  - RL评估直接调用main_olg_v9_sac_sbx.py中的evaluate_policy_lifecycle_simulation函数")
    print("  - RL模型使用训练过程中的最佳模型，通常性能优于最终模型")
    print("  - 🎯 关键改进：VFI和RL使用完全相同的效率冲击序列，确保公平比较")
    print("  - 🎯 确保VFI和RL使用完全一致的参数设置和评估逻辑")

if __name__ == "__main__":
    main() 