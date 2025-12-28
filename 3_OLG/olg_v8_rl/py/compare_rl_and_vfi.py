"""
比较VFI和RL（SAC Agent）的优化结果 - Python实现
目标：在相同的宏观和微观参数下，比较两种方法的优化效果
特性：使用MATLAB Engine调用MATLAB的VFI方法，与Python的SBX SAC进行比较

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
    # 尝试不同的中文字体
    chinese_fonts = ['SimHei', 'Microsoft YaHei', 'WenQuanYi Micro Hei', 'DejaVu Sans']
    available_fonts = [f.name for f in fm.fontManager.ttflist]
    
    # 找到可用的中文字体
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

# 导入MATLAB Engine
try:
    import matlab.engine
    MATLAB_AVAILABLE = True
    print("✅ MATLAB Engine 可用")
except ImportError:
    MATLAB_AVAILABLE = False
    print("❌ MATLAB Engine 不可用，请安装: pip install matlabengine")

# 尝试导入SBX，如果失败则使用SB3
try:
    from sbx import SAC as SBX_SAC
    SBX_AVAILABLE = True
    print("✅ SBX (Stable Baselines Jax) 可用")
except ImportError:
    SBX_AVAILABLE = False
    print("⚠️ SBX 不可用，将使用 SB3")

# 导入SB3作为备选
from stable_baselines3 import SAC as SB3_SAC
from stable_baselines3.common.evaluation import evaluate_policy

# 导入生命周期评估函数
from main_olg_v8_sac_sbx import evaluate_policy_lifecycle_simulation

from olg_utils import OLGUtils
from olg_env_v8_sac import OLGEnvV8SAC

class RLVFIComparator:
    """RL和VFI方法比较器 - 使用MATLAB Engine"""
    
    def __init__(self, use_sbx: bool = True):
        """
        初始化比较器
        
        Args:
            use_sbx: 是否优先使用SBX，否则使用SB3
        """
        if not MATLAB_AVAILABLE:
            raise RuntimeError("MATLAB Engine不可用，请安装: pip install matlabengine")
        
        self.use_sbx = use_sbx and SBX_AVAILABLE
        self.SAC_class = SBX_SAC if self.use_sbx else SB3_SAC
        self.backend_name = "SBX" if self.use_sbx else "SB3"
        
        print(f"🤖 使用 {self.backend_name} 作为RL后端")
        
        # 启动MATLAB Engine
        print("🔧 启动MATLAB Engine...")
        self.eng = matlab.engine.start_matlab()
        print("✅ MATLAB Engine启动成功")
        
        # 固定的测试参数（与MATLAB版本一致）
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
    
    def __del__(self):
        """析构函数，关闭MATLAB Engine"""
        if hasattr(self, 'eng') and self.eng is not None:
            self.eng.quit()
            print("🔧 MATLAB Engine已关闭")
    
    def load_rl_model(self, model_path: Optional[str] = None, use_best_model: bool = True) -> Tuple[Any, Dict]:
        """
        加载训练好的RL模型
        
        Args:
            model_path: 模型路径，如果为None则自动选择
            use_best_model: 是否优先使用best model（默认True）
                          - True: 优先使用训练过程中评估性能最佳的模型
                          - False: 使用训练结束时的最终模型
            
        Returns:
            model: 加载的模型
            config: 模型配置（包含训练时的cS和paramS_for_rl）
        
        Note:
            Best model通常在评估性能上表现更好，因为它避免了过拟合
            Final model是训练结束时的状态，可能存在性能退化
            🎯 关键：返回的config包含训练时的参数，确保评估一致性
        """
        if model_path is None:
            # 自动选择模型路径
            if self.use_sbx:
                if use_best_model:
                    best_model_path = './py/best_model_sbx/best_model.zip'
                    final_model_path = './py/final_sac_agent_olg_sbx.zip'
                    
                    # 优先使用best model，如果不存在则使用final model
                    if os.path.exists(best_model_path):
                        model_path = best_model_path
                        print(f"🏆 使用训练过程中的最佳模型: {model_path}")
                    elif os.path.exists(final_model_path):
                        model_path = final_model_path
                        print(f"⚠️ Best model不存在，使用最终模型: {model_path}")
                    else:
                        raise FileNotFoundError("未找到SBX模型文件 (best_model或final_model)")
                else:
                    model_path = './py/final_sac_agent_olg_sbx.zip'
                    print(f"📁 使用最终模型: {model_path}")
                
                config_path = './py/training_config_sbx.pkl'
            else:
                if use_best_model:
                    best_model_path = './py/best_model/best_model.zip'
                    final_model_path = './py/final_sac_agent_olg_sb3.zip'
                    
                    # 优先使用best model，如果不存在则使用final model
                    if os.path.exists(best_model_path):
                        model_path = best_model_path
                        print(f"🏆 使用训练过程中的最佳模型: {model_path}")
                    elif os.path.exists(final_model_path):
                        model_path = final_model_path
                        print(f"⚠️ Best model不存在，使用最终模型: {model_path}")
                    else:
                        raise FileNotFoundError("未找到SB3模型文件 (best_model或final_model)")
                else:
                    model_path = './py/final_sac_agent_olg_sb3.zip'
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
    
    def call_matlab_vfi(self) -> Dict[str, Any]:
        """
        使用MATLAB Engine调用VFI方法
        
        Returns:
            vfi_results: VFI结果字典
        """
        print("\n🔧 使用MATLAB Engine调用VFI方法...")
        
        try:
            start_time = time.time()
            
            # 1. 初始化参数
            print("1️⃣ 初始化MATLAB参数...")
            cS = self.eng.main_olg_v8_utils.ParameterValues_HuggettStyle()
            
            # 检查网格是否正确生成
            if 'kppsGridV' not in cS or len(cS['kppsGridV']) == 0:
                print("  重新生成网格...")
                cS = self.eng.main_olg_v8_utils.generateGrids(cS)
            
            # 2. 设置宏观参数
            M_test_matlab = matlab.double([
                [self.M_test['R_k_net_factor']],
                [self.M_test['w_gross']],
                [self.M_test['TR_total']],
                [self.M_test['b_payg_avg_retiree']],
                [self.M_test['tau_l']],
                [self.M_test['theta_payg_actual']]
            ])
            
            # 3. 初始化VFI参数
            print("2️⃣ 初始化VFI参数...")
            [leLogGridV, leTrProbM, leProb1V] = self.eng.main_olg_v8_utils.EarningProcess_olgm(cS, nargout=3)
            
            # 创建paramS_vfi结构
            paramS_vfi = {}
            paramS_vfi['leLogGridV'] = leLogGridV
            paramS_vfi['leTrProbM'] = leTrProbM
            paramS_vfi['leProb1V'] = leProb1V
            paramS_vfi['leGridV'] = self.eng.exp(leLogGridV)
            paramS_vfi['ageEffV_new'] = cS['ageEffV_new']
            paramS_vfi['tau_l'] = self.M_test['tau_l']
            paramS_vfi['theta_payg_actual_for_hh'] = self.M_test['theta_payg_actual']
            paramS_vfi['pps_tax_deferral_active'] = cS['pps_active']
            
            # 4. 构建PAYG福利向量
            print("3️⃣ 构建PAYG福利向量...")
            aD_new = int(cS['aD_new'])
            aR_new = int(cS['aR_new'])
            bV_payg_vfi = matlab.double(np.zeros((1, aD_new)))
            
            if aR_new < aD_new:
                # 创建PAYG福利向量
                payg_benefits = np.zeros(aD_new)
                payg_benefits[aR_new:] = self.M_test['b_payg_avg_retiree']
                bV_payg_vfi = matlab.double([payg_benefits])
            
            # 5. 运行VFI求解
            print("4️⃣ 运行VFI求解...")
            vfi_start = time.time()
            
            [cPolM_vfi, kPolM_vfi, cPpsPolM_vfi, VPolM_vfi] = self.eng.main_olg_v8_utils.HHSolution_VFI_Huggett(
                self.M_test['R_k_net_factor'],
                self.M_test['w_gross'],
                self.M_test['TR_total'],
                bV_payg_vfi,
                paramS_vfi,
                cS,
                nargout=4
            )
            
            vfi_time = time.time() - vfi_start
            total_time = time.time() - start_time
            
            print(f"✅ VFI求解完成，耗时: {vfi_time:.2f} 秒")
            print(f"📊 策略矩阵尺寸: {np.array(cPolM_vfi).shape}")
            
            # 6. 准备返回结果
            result_dict = {
                'success': True,
                'cPolM': np.array(cPolM_vfi),
                'kPolM': np.array(kPolM_vfi),
                'cPpsPolM': np.array(cPpsPolM_vfi),
                'VPolM': np.array(VPolM_vfi),
                'M_test': self.M_test,
                'cS': cS,
                'paramS_vfi': paramS_vfi,
                'vfi_time': vfi_time,
                'total_time': total_time
            }
            
            print("✅ VFI结果准备完成")
            return result_dict
            
        except Exception as e:
            print(f"❌ MATLAB VFI调用失败: {e}")
            return {
                'success': False,
                'error': str(e)
            }
    
    def simulate_lifecycle_comparison(self, rl_model: Any, vfi_results: Dict, rl_config: Dict,
                                    n_sim: int = 100, random_seed: int = 42) -> Dict[str, Any]:
        """
        模拟生命周期轨迹比较（期望化死亡率版本，与MATLAB版本一致）
        
        Args:
            rl_model: 训练好的RL模型
            vfi_results: VFI结果
            rl_config: RL训练配置（包含训练时的cS和paramS_for_rl）
            n_sim: 模拟个体数量
            random_seed: 随机种子，确保结果可重现
            
        Returns:
            comparison_results: 比较结果
        """
        print(f"\n🎯 模拟生命周期轨迹比较 (n_sim={n_sim}, seed={random_seed})")
        print("注意：使用期望化死亡率，确保与MATLAB版本一致")
        
        # 🎲 设置随机种子确保可重现性
        np.random.seed(random_seed)
        print(f"🎲 随机种子已设置: {random_seed}")
        
        if not vfi_results.get('success', False):
            print("❌ VFI结果不可用，无法进行比较")
            return {'success': False, 'error': 'VFI结果不可用'}
        
        # 从VFI结果获取参数
        cS = vfi_results['cS']
        paramS_vfi = vfi_results['paramS_vfi']
        
        # 获取维度
        aD_new = int(cS['aD_new'])
        aR_new = int(cS['aR_new'])
        nk = int(cS['nk'])
        nkpps = int(cS['nkpps'])
        nw = int(cS['nw'])
        
        # 从MATLAB结果提取网格
        kGridV = np.array(cS['kGridV']).flatten()
        kppsGridV = np.array(cS['kppsGridV']).flatten()
        leGridV = np.array(paramS_vfi['leGridV']).flatten()
        leTrProbM = np.array(paramS_vfi['leTrProbM'])
        leProb1V = np.array(paramS_vfi['leProb1V']).flatten()
        ageEffV_new = np.array(cS['ageEffV_new']).flatten()
        beta = float(cS['beta'])
        sigma = float(cS['sigma'])
        kMin = float(cS['kMin'])
        kppsMin = float(cS['kppsMin'])
        kMax = float(cS['kMax'])
        kppsMax = float(cS['kppsMax'])
        pps_active = bool(cS['pps_active'])
        
        # 获取策略矩阵
        cPolM_vfi = vfi_results['cPolM']
        kPolM_vfi = vfi_results['kPolM']
        cPpsPolM_vfi = vfi_results['cPpsPolM']
        
        print(f"📊 网格尺寸: nk={nk}, nkpps={nkpps}, nw={nw}, aD_new={aD_new}")
        print(f"📊 VFI策略矩阵尺寸: {cPolM_vfi.shape}")
        
        # 🆕 使用统一的生命周期评估函数进行RL评估
        print("\n🧬 使用统一生命周期评估函数 (evaluate_policy_lifecycle_simulation)")
        print("🎯 确保与main_olg_v8_sac_sbx.py完全一致的评估逻辑")
        
        # 🔧 使用训练时保存的参数，确保评估一致性
        if rl_config and 'cS' in rl_config and 'paramS_for_rl' in rl_config:
            print("✅ 使用训练时保存的参数 (cS, paramS_for_rl, rng_M)")
            cS_python = rl_config['cS']
            paramS_for_rl = rl_config['paramS_for_rl']
            rng_M = rl_config['rng_M']
            print(f"📊 训练时参数 - aD_new: {cS_python['aD_new']}, beta: {cS_python['beta']}")
        else:
            print("⚠️ 训练配置不完整，回退到重新创建参数")
            # 回退方案：重新创建参数（与原版本一致）
        cS_python = OLGUtils.parameter_values_huggett_style()
        paramS_for_rl = {}
        (paramS_for_rl['leLogGridV'], 
         paramS_for_rl['leTrProbM'], 
         paramS_for_rl['leProb1V']) = OLGUtils.earning_process_olgm(cS_python)
        paramS_for_rl['leGridV'] = np.exp(paramS_for_rl['leLogGridV'])
        paramS_for_rl['ageEffV_new'] = cS_python['ageEffV_new']
        
        rng_M = {
            'R_k_net_factor': [1.01, 1.05],
            'w_gross': [1.5, 2.5],
            'TR_total': [0.0, 0.2],
            'b_payg_avg_retiree': [0.1, 0.8],
            'tau_l': [0.05, 0.25],
            'theta_payg_actual': [0.05, 0.20]
        }
        
        # 🚀 调用统一的生命周期评估函数进行RL评估
        print("🔄 开始RL生命周期评估...")
        print(f"📊 评估参数: n_sim={n_sim}, random_seed={random_seed}, gamma=0.97")
        rl_eval_start = time.time()
        
        mean_utility_rl, std_utility_rl, rl_lifecycle_results = evaluate_policy_lifecycle_simulation(
            rl_model, cS_python, paramS_for_rl, rng_M, 
            n_sim=n_sim, deterministic=True, gamma=0.97, random_seed=random_seed, verbose=True
        )
        
        rl_eval_time = time.time() - rl_eval_start
        print(f"✅ RL评估完成，耗时: {rl_eval_time:.2f} 秒")
        print(f"📊 RL评估结果: {mean_utility_rl:.4f} ± {std_utility_rl:.4f}")
        
        # 🔍 检查参数一致性
        using_training_params = rl_config and 'cS' in rl_config and 'paramS_for_rl' in rl_config
        print(f"🎯 使用训练时参数: {'是' if using_training_params else '否'}")
        
        # 从RL评估结果中提取轨迹数据
        lifetime_utility_rl = rl_lifecycle_results['lifetime_utility_rl']
        k_path_rl = rl_lifecycle_results['k_path_rl']
        c_path_rl = rl_lifecycle_results['c_path_rl']
        cpps_path_rl = rl_lifecycle_results['cpps_path_rl']
        
        # 🔄 VFI生命周期模拟（保持原有逻辑）
        print("\n🔄 开始VFI生命周期模拟...")
        vfi_sim_start = time.time()
        
        # 初始化VFI结果存储
        lifetime_utility_vfi = np.zeros(n_sim)
        
        # VFI生命周期轨迹存储
        k_path_vfi = np.zeros((n_sim, aD_new))
        c_path_vfi = np.zeros((n_sim, aD_new))
        cpps_path_vfi = np.zeros((n_sim, aD_new))
        
        # 确保与RL相同的随机种子（重新设置）
        np.random.seed(random_seed)
        
        for i_sim in range(n_sim):
            if (i_sim + 1) % 20 == 0:
                print(f"  VFI进度: {i_sim + 1}/{n_sim}")
            
            # 初始状态
            k_current_vfi = kMin
            kpps_current_vfi = kppsMin
            
            # 初始效率冲击（与RL评估使用相同逻辑）
            eps_idx_current = np.where(np.random.rand() <= np.cumsum(leProb1V))[0]
            if len(eps_idx_current) > 0:
                eps_idx_current = eps_idx_current[0]
            else:
                eps_idx_current = 0
            
            utility_sum_vfi = 0
            
            for age_idx in range(aD_new):
                # VFI策略
                k_idx_vfi = np.argmin(np.abs(kGridV - k_current_vfi))
                kpps_idx_vfi = np.argmin(np.abs(kppsGridV - kpps_current_vfi))
                
                c_vfi = cPolM_vfi[k_idx_vfi, kpps_idx_vfi, eps_idx_current, age_idx]
                k_next_vfi = kPolM_vfi[k_idx_vfi, kpps_idx_vfi, eps_idx_current, age_idx]
                cpps_vfi = cPpsPolM_vfi[k_idx_vfi, kpps_idx_vfi, eps_idx_current, age_idx]
                
                # 存储VFI轨迹
                k_path_vfi[i_sim, age_idx] = k_current_vfi
                c_path_vfi[i_sim, age_idx] = c_vfi
                cpps_path_vfi[i_sim, age_idx] = cpps_vfi
                
                # 计算VFI效用（使用MATLAB的CES效用函数）
                try:
                    u_vfi = self.eng.main_olg_v8_utils.CES_utility(float(c_vfi), float(sigma), cS, nargout=2)[1]
                    u_vfi = float(u_vfi)
                except:
                    # 备选：简单CES效用函数
                    if sigma == 1:
                        u_vfi = np.log(max(c_vfi, 1e-8))
                    else:
                        u_vfi = (max(c_vfi, 1e-8)**(1-sigma) - 1) / (1-sigma)
                
                # 累积折现效用（确保使用相同的折现率0.97）
                discount_factor = beta ** age_idx
                utility_sum_vfi += discount_factor * u_vfi
                
                # VFI状态更新
                k_current_vfi = k_next_vfi
                
                # PPS资产演化（简化处理）
                if pps_active:
                    R_k_net_factor = self.M_test['R_k_net_factor']
                    pps_return_premium = float(cS.get('pps_return_rate_premium', 0))
                    pps_return_factor = 1 + ((R_k_net_factor - 1) + pps_return_premium)
                    
                    kpps_current_vfi = (kpps_current_vfi + cpps_vfi) * pps_return_factor
                    kpps_current_vfi = max(kppsMin, min(kppsMax, kpps_current_vfi))
                
                # 效率冲击演化
                if age_idx < aD_new - 1:
                    trans_probs = leTrProbM[eps_idx_current, :]
                    eps_idx_current = np.where(np.random.rand() <= np.cumsum(trans_probs))[0]
                    if len(eps_idx_current) > 0:
                        eps_idx_current = eps_idx_current[0]
                    else:
                        eps_idx_current = min(eps_idx_current, nw-1)
            
            lifetime_utility_vfi[i_sim] = utility_sum_vfi
        
        vfi_sim_time = time.time() - vfi_sim_start
        print(f"✅ VFI生命周期模拟完成，耗时: {vfi_sim_time:.2f} 秒")
        
        # 总模拟时间
        sim_time = rl_eval_time + vfi_sim_time
        
        # 📊 计算统计结果（RL结果来自统一评估函数）
        mean_utility_vfi = np.mean(lifetime_utility_vfi)
        std_utility_vfi = np.std(lifetime_utility_vfi)
        # RL结果已从统一函数获得
        # mean_utility_rl, std_utility_rl 已在上面定义
        
        utility_diff = mean_utility_rl - mean_utility_vfi
        utility_improvement_pct = (utility_diff / abs(mean_utility_vfi)) * 100
        
        print(f"\n📊 生涯总折现效用比较 (基于 {n_sim} 个个体):")
        print(f"  📈 VFI方法: {mean_utility_vfi:.4f} ± {std_utility_vfi:.4f}")
        print(f"  🤖 RL方法 ({self.backend_name}): {mean_utility_rl:.4f} ± {std_utility_rl:.4f}")
        print(f"  🔍 差异 (RL - VFI): {utility_diff:.4f}")
        print(f"  📊 相对改进: {utility_improvement_pct:.2f}%")
        print(f"  🎯 RL评估使用统一函数: evaluate_policy_lifecycle_simulation")
        
        if utility_diff > 0:
            print("  >>> 🏆 RL方法表现更好！")
        elif utility_diff < 0:
            print("  >>> 🏆 VFI方法表现更好！")
        else:
            print("  >>> 🤝 两种方法表现相当。")
        
        # 统计显著性检验
        from scipy import stats
        t_stat, p_value = stats.ttest_rel(lifetime_utility_rl, lifetime_utility_vfi)
        print(f"\n📈 统计显著性检验 (配对t检验): p值: {p_value:.4f}")
        print(f"  {'差异显著' if p_value < 0.05 else '差异不显著'} (α=0.05)")
        
        # 🆕 构造比较结果（整合统一评估函数的结果）
        comparison_results = {
            'success': True,
            'rl_backend': self.backend_name,
            'rl_evaluation_method': 'evaluate_policy_lifecycle_simulation',  # 🆕 标记使用统一函数
            'using_training_params': using_training_params,  # 🆕 标记是否使用训练时参数
            'M_test': self.M_test,
            'n_sim': n_sim,
            'random_seed': random_seed,
            'rl_simulation_time': rl_eval_time,  # 🆕 仅RL评估时间
            'vfi_simulation_time': vfi_sim_time,  # 🆕 仅VFI模拟时间
            'total_simulation_time': sim_time,   # 🆕 总时间
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
            'rl_lifecycle_results': rl_lifecycle_results  # 🆕 保存完整的RL评估结果
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
        
        if not results.get('success', False):
            print("❌ 结果数据不可用，无法绘图")
            return
        
        n_sim = results['n_sim']
        aD_new = results['k_path_rl'].shape[1]
        
        # 创建2x3的子图布局
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        # 设置中文字体（备用确保）
        plt.rcParams['font.sans-serif'] = matplotlib.rcParams['font.sans-serif']
        plt.rcParams['axes.unicode_minus'] = False
        fig.suptitle(f'VFI vs RL ({results["rl_backend"]}) 生命周期比较 (n={n_sim})', fontsize=16)
        
        # 1. 效用分布比较
        axes[0, 0].hist(results['lifetime_utility_vfi'], bins=20, alpha=0.7, 
                       label='VFI', color='red')
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
        axes[0, 1].set_xlabel('VFI 生涯效用')
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
        
        # 4. 平均资产路径
        age_path = np.arange(1, aD_new + 1)
        mean_k_vfi = np.mean(results['k_path_vfi'], axis=0)
        mean_k_rl = np.mean(results['k_path_rl'], axis=0)
        axes[1, 0].plot(age_path, mean_k_vfi, 'r-', linewidth=2, label='VFI')
        axes[1, 0].plot(age_path, mean_k_rl, 'b--', linewidth=2, 
                       label=f'RL ({results["rl_backend"]})')
        axes[1, 0].set_xlabel('年龄组')
        axes[1, 0].set_ylabel('平均资产')
        axes[1, 0].set_title('平均资产路径')
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)
        
        # 5. 平均消费路径
        mean_c_vfi = np.mean(results['c_path_vfi'], axis=0)
        mean_c_rl = np.mean(results['c_path_rl'], axis=0)
        axes[1, 1].plot(age_path, mean_c_vfi, 'r-', linewidth=2, label='VFI')
        axes[1, 1].plot(age_path, mean_c_rl, 'b--', linewidth=2, 
                       label=f'RL ({results["rl_backend"]})')
        axes[1, 1].set_xlabel('年龄组')
        axes[1, 1].set_ylabel('平均消费')
        axes[1, 1].set_title('平均消费路径')
        axes[1, 1].legend()
        axes[1, 1].grid(True, alpha=0.3)
        
        # 6. 平均PPS缴费路径
        mean_cpps_vfi = np.mean(results['cpps_path_vfi'], axis=0)
        mean_cpps_rl = np.mean(results['cpps_path_rl'], axis=0)
        axes[1, 2].plot(age_path, mean_cpps_vfi, 'r-', linewidth=2, label='VFI')
        axes[1, 2].plot(age_path, mean_cpps_rl, 'b--', linewidth=2, 
                       label=f'RL ({results["rl_backend"]})')
        axes[1, 2].set_xlabel('年龄组')
        axes[1, 2].set_ylabel('平均PPS缴费')
        axes[1, 2].set_title('平均PPS缴费路径')
        axes[1, 2].legend()
        axes[1, 2].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"📈 VFI vs RL比较图表已保存到: {save_path}")
        
        # 打印数值摘要
        print(f"\n📊 路径差异摘要:")
        k_diff_mean = np.mean(mean_k_rl - mean_k_vfi)
        c_diff_mean = np.mean(mean_c_rl - mean_c_vfi) 
        cpps_diff_mean = np.mean(mean_cpps_rl - mean_cpps_vfi)
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
    print("🔬 VFI vs RL (SAC) 优化方法比较")
    print("=" * 80)
    
    # 0. 检查中文字体
    check_chinese_font_availability()
    
    # 1. 初始化比较器
    print("\n1️⃣ 初始化比较器...")
    comparator = RLVFIComparator(use_sbx=True)  # 优先使用SBX
    
    # 2. 加载RL模型
    print("\n2️⃣ 加载RL模型...")
    try:
        rl_model, rl_config = comparator.load_rl_model(use_best_model=True)
        print(f"✅ 成功加载 {comparator.backend_name} 模型")
    except Exception as e:
        print(f"❌ 加载RL模型失败: {e}")
        print("💡 请先运行训练脚本生成模型")
        return
    
    # 3. 运行MATLAB VFI
    print("\n3️⃣ 运行MATLAB VFI...")
    try:
        vfi_results = comparator.call_matlab_vfi()
        if vfi_results.get('success', False):
            print("✅ VFI运行成功")
        else:
            print(f"❌ VFI运行失败: {vfi_results.get('error', '未知错误')}")
    except Exception as e:
        print(f"❌ VFI调用失败: {e}")
        vfi_results = {'success': False, 'error': str(e)}
    
    # 4. 运行比较模拟
    print("\n4️⃣ 运行生命周期比较模拟...")
    try:
        comparison_results = comparator.simulate_lifecycle_comparison(
            rl_model, vfi_results, rl_config, n_sim=500, random_seed=42
        )
        print("✅ 比较模拟完成")
    except Exception as e:
        print(f"❌ 比较模拟失败: {e}")
        return
    
    # 5. 绘制结果
    print("\n5️⃣ 绘制比较结果...")
    try:
        comparator.plot_comparison_results(comparison_results)
        print("✅ 图表生成完成")
    except Exception as e:
        print(f"❌ 绘图失败: {e}")
    
    # 6. 保存结果
    print("\n6️⃣ 保存比较结果...")
    try:
        comparator.save_comparison_results(comparison_results)
        print("✅ 结果保存完成")
    except Exception as e:
        print(f"❌ 保存失败: {e}")
    
    # 7. 输出最终摘要
    print("\n" + "=" * 80)
    print("📋 VFI vs RL 比较分析摘要")
    print("=" * 80)
    
    if comparison_results.get('success', False):
        print(f"🤖 RL后端: {comparison_results['rl_backend']}")
        print(f"🏆 RL模型类型: 训练过程中的最佳模型 (Best Model)")
        print(f"🎯 RL评估方法: {comparison_results.get('rl_evaluation_method', 'N/A')}")
        print(f"🔧 使用训练时参数: {'是' if comparison_results.get('using_training_params', False) else '否'}")
        print(f"🎯 模拟个体数: {comparison_results['n_sim']}")
        print(f"🎲 随机种子: {comparison_results.get('random_seed', 'N/A')}")
        print(f"⏱️ VFI计算时间: {comparison_results['vfi_time']:.2f} 秒")
        print(f"⏱️ VFI模拟时间: {comparison_results.get('vfi_simulation_time', 0):.2f} 秒")
        print(f"⏱️ RL评估时间: {comparison_results.get('rl_simulation_time', 0):.2f} 秒")
        print(f"⏱️ 总模拟时间: {comparison_results.get('total_simulation_time', 0):.2f} 秒")
        print(f"📊 VFI平均效用: {comparison_results['mean_utility_vfi']:.4f} ± {comparison_results['std_utility_vfi']:.4f}")
        print(f"📊 RL平均效用: {comparison_results['mean_utility_rl']:.4f} ± {comparison_results['std_utility_rl']:.4f}")
        print(f"🔍 效用差异: {comparison_results['utility_diff']:.4f} ({comparison_results['utility_improvement_pct']:.2f}%)")
        print(f"📈 统计显著性: {'显著' if comparison_results['is_significant'] else '不显著'} (p={comparison_results['p_value']:.4f})")
        
        if comparison_results['utility_diff'] > 0:
            print("🏆 RL方法表现更好！")
        elif comparison_results['utility_diff'] < 0:
            print("🏆 VFI方法表现更好！")
        else:
            print("🤝 两种方法表现相当。")
    else:
        print("❌ 比较分析失败")
    
    print("=" * 80)
    print("🎉 VFI vs RL 比较分析完成!")
    print("💡 注意：")
    print("  - 本次比较使用期望化死亡率，确保与MATLAB VFI版本一致")
    print("  - RL模型使用训练过程中的最佳模型，通常性能优于最终模型")
    print("  - 🆕 RL评估使用统一函数evaluate_policy_lifecycle_simulation")
    print("  - 🎯 确保与main_olg_v8_sac_sbx.py完全一致的评估逻辑")

if __name__ == "__main__":
    main() 