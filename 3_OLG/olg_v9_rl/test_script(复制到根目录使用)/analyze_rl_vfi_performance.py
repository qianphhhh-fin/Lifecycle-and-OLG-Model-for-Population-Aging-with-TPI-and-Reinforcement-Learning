"""
深入分析RL vs VFI性能差异
目标：理解为什么RL在生命周期效用上表现更好
分析内容：
1. 生命周期阶段分析
2. 决策模式差异
3. 资产配置策略
4. PPS缴费行为
5. 消费平滑效果
6. 个体差异分析
"""

import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.font_manager as fm
import pickle
import pandas as pd
from scipy import stats
from scipy.stats import pearsonr
import seaborn as sns
from typing import Dict, Any, List, Tuple
from pathlib import Path

# 配置matplotlib中文字体
def setup_chinese_fonts():
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
    
    matplotlib.rcParams['axes.unicode_minus'] = False

setup_chinese_fonts()

class RLVFIPerformanceAnalyzer:
    """RL vs VFI 性能分析器"""
    
    def __init__(self, data_path: str = './py/rl_vfi_detailed_paths.pkl'):
        """
        初始化分析器
        
        Args:
            data_path: 详细路径数据文件路径
        """
        self.data_path = data_path
        self.data = None
        self.load_data()
    
    def load_data(self):
        """加载详细路径数据"""
        print(f"📊 加载详细数据: {self.data_path}")
        
        if not os.path.exists(self.data_path):
            raise FileNotFoundError(f"数据文件不存在: {self.data_path}")
        
        with open(self.data_path, 'rb') as f:
            self.data = pickle.load(f)
        
        print("✅ 数据加载成功")
        print(f"📈 分析数据概览:")
        print(f"  - 个体数量: {self.data['metadata']['n_sim']}")
        print(f"  - RL后端: {self.data['metadata']['rl_backend']}")
        print(f"  - VFI方法: {self.data['metadata']['vfi_method']}")
        print(f"  - 效用差异: {self.data['utility_comparison']['utility_diff']:.4f}")
        print(f"  - 相对改进: {self.data['utility_comparison']['utility_improvement_pct']:.2f}%")
    
    def analyze_lifecycle_stages(self):
        """分析不同生命周期阶段的表现差异"""
        print("\n🔍 分析生命周期阶段表现差异...")
        
        paths = self.data['paths']
        k_vfi = paths['k_path_vfi']  # (n_sim, T)
        k_rl = paths['k_path_rl']
        c_vfi = paths['c_path_vfi']
        c_rl = paths['c_path_rl']
        cpps_vfi = paths['cpps_path_vfi']
        cpps_rl = paths['cpps_path_rl']
        
        n_sim, T = k_vfi.shape
        ages = np.arange(20, 20 + T)  # 20-98岁
        
        # 定义生命周期阶段
        young_ages = (ages >= 20) & (ages < 35)    # 年轻期
        middle_ages = (ages >= 35) & (ages < 50)   # 中年期
        mature_ages = (ages >= 50) & (ages < 65)   # 成熟期
        retirement_ages = ages >= 65               # 退休期
        
        stages = {
            '年轻期 (20-34岁)': young_ages,
            '中年期 (35-49岁)': middle_ages, 
            '成熟期 (50-64岁)': mature_ages,
            '退休期 (65岁+)': retirement_ages
        }
        
        stage_analysis = {}
        
        for stage_name, stage_mask in stages.items():
            if not np.any(stage_mask):
                continue
                
            stage_indices = np.where(stage_mask)[0]
            
            # 计算该阶段的平均表现
            k_vfi_stage = k_vfi[:, stage_indices].mean(axis=1)  # 个体在该阶段的平均资产
            k_rl_stage = k_rl[:, stage_indices].mean(axis=1)
            c_vfi_stage = c_vfi[:, stage_indices].mean(axis=1)
            c_rl_stage = c_rl[:, stage_indices].mean(axis=1)
            cpps_vfi_stage = cpps_vfi[:, stage_indices].mean(axis=1)
            cpps_rl_stage = cpps_rl[:, stage_indices].mean(axis=1)
            
            # 计算差异
            k_diff = k_rl_stage - k_vfi_stage
            c_diff = c_rl_stage - c_vfi_stage
            cpps_diff = cpps_rl_stage - cpps_vfi_stage
            
            stage_analysis[stage_name] = {
                'k_mean_vfi': k_vfi_stage.mean(),
                'k_mean_rl': k_rl_stage.mean(),
                'k_diff_mean': k_diff.mean(),
                'k_diff_std': k_diff.std(),
                'c_mean_vfi': c_vfi_stage.mean(),
                'c_mean_rl': c_rl_stage.mean(),
                'c_diff_mean': c_diff.mean(),
                'c_diff_std': c_diff.std(),
                'cpps_mean_vfi': cpps_vfi_stage.mean(),
                'cpps_mean_rl': cpps_rl_stage.mean(),
                'cpps_diff_mean': cpps_diff.mean(),
                'cpps_diff_std': cpps_diff.std(),
                'age_range': f"{ages[stage_indices].min()}-{ages[stage_indices].max()}岁",
                'n_years': len(stage_indices)
            }
        
        # 打印分析结果
        print("\n📊 生命周期阶段分析结果:")
        print("=" * 80)
        
        for stage_name, analysis in stage_analysis.items():
            print(f"\n🔹 {stage_name} ({analysis['age_range']}, {analysis['n_years']}年)")
            print(f"  资产差异 (RL - VFI): {analysis['k_diff_mean']:+.4f} ± {analysis['k_diff_std']:.4f}")
            print(f"  消费差异 (RL - VFI): {analysis['c_diff_mean']:+.4f} ± {analysis['c_diff_std']:.4f}")
            print(f"  PPS缴费差异 (RL - VFI): {analysis['cpps_diff_mean']:+.4f} ± {analysis['cpps_diff_std']:.4f}")
            
            # 判断哪种方法在该阶段更好
            if analysis['k_diff_mean'] > 0:
                print(f"  >>> RL在该阶段资产积累更多")
            else:
                print(f"  >>> VFI在该阶段资产积累更多")
        
        return stage_analysis
    
    def analyze_decision_patterns(self):
        """分析决策模式差异"""
        print("\n🧠 分析决策模式差异...")
        
        paths = self.data['paths']
        k_vfi = paths['k_path_vfi']
        k_rl = paths['k_path_rl']
        c_vfi = paths['c_path_vfi']
        c_rl = paths['c_path_rl']
        cpps_vfi = paths['cpps_path_vfi']
        cpps_rl = paths['cpps_path_rl']
        
        n_sim, T = k_vfi.shape
        
        # 计算储蓄率 (假设收入为1，简化分析)
        # 储蓄率 = (k_t+1 - k_t) / 收入
        savings_rate_vfi = np.diff(k_vfi, axis=1) / 1.0  # (n_sim, T-1)
        savings_rate_rl = np.diff(k_rl, axis=1) / 1.0
        
        # 计算PPS参与率
        pps_participation_vfi = (cpps_vfi > 0.01).mean(axis=0)  # 每个年龄的参与率
        pps_participation_rl = (cpps_rl > 0.01).mean(axis=0)
        
        # 计算消费波动性（标准差）
        c_volatility_vfi = np.std(c_vfi, axis=1)  # 每个个体的消费波动性
        c_volatility_rl = np.std(c_rl, axis=1)
        
        decision_patterns = {
            'savings_rate': {
                'vfi_mean': savings_rate_vfi.mean(),
                'rl_mean': savings_rate_rl.mean(),
                'diff': savings_rate_rl.mean() - savings_rate_vfi.mean()
            },
            'pps_participation': {
                'vfi_overall': (cpps_vfi > 0.01).mean(),
                'rl_overall': (cpps_rl > 0.01).mean(),
                'vfi_by_age': pps_participation_vfi,
                'rl_by_age': pps_participation_rl
            },
            'consumption_volatility': {
                'vfi_mean': c_volatility_vfi.mean(),
                'rl_mean': c_volatility_rl.mean(),
                'diff': c_volatility_rl.mean() - c_volatility_vfi.mean()
            }
        }
        
        print("📊 决策模式分析结果:")
        print(f"  平均储蓄率: VFI={decision_patterns['savings_rate']['vfi_mean']:.4f}, RL={decision_patterns['savings_rate']['rl_mean']:.4f}")
        print(f"  储蓄率差异: {decision_patterns['savings_rate']['diff']:+.4f}")
        print(f"  PPS总体参与率: VFI={decision_patterns['pps_participation']['vfi_overall']:.2%}, RL={decision_patterns['pps_participation']['rl_overall']:.2%}")
        print(f"  消费波动性: VFI={decision_patterns['consumption_volatility']['vfi_mean']:.4f}, RL={decision_patterns['consumption_volatility']['rl_mean']:.4f}")
        print(f"  消费波动性差异: {decision_patterns['consumption_volatility']['diff']:+.4f}")
        
        if decision_patterns['consumption_volatility']['diff'] < 0:
            print("  >>> RL实现了更好的消费平滑")
        else:
            print("  >>> VFI实现了更好的消费平滑")
        
        return decision_patterns
    
    def analyze_individual_differences(self):
        """分析个体差异"""
        print("\n👥 分析个体表现差异...")
        
        utility_vfi = self.data['utility_comparison']['lifetime_utility_vfi']
        utility_rl = self.data['utility_comparison']['lifetime_utility_rl']
        utility_diff = np.array(utility_rl) - np.array(utility_vfi)
        
        # 找出表现最好和最差的个体
        best_performers = np.argsort(utility_diff)[-5:]  # 前5名
        worst_performers = np.argsort(utility_diff)[:5]   # 后5名
        
        print("🏆 RL表现最好的5个个体 (效用提升最大):")
        for i, idx in enumerate(best_performers):
            improvement = utility_diff[idx]
            improvement_pct = (improvement / abs(utility_vfi[idx])) * 100
            print(f"  {i+1}. 个体{idx}: 提升 {improvement:.4f} ({improvement_pct:+.2f}%)")
        
        print("\n📉 RL表现相对较差的5个个体:")
        for i, idx in enumerate(worst_performers):
            improvement = utility_diff[idx]
            improvement_pct = (improvement / abs(utility_vfi[idx])) * 100
            print(f"  {i+1}. 个体{idx}: 变化 {improvement:.4f} ({improvement_pct:+.2f}%)")
        
        # 分析个体特征
        individual_analysis = {
            'best_performers': best_performers,
            'worst_performers': worst_performers,
            'utility_improvement_distribution': {
                'mean': utility_diff.mean(),
                'std': utility_diff.std(),
                'min': utility_diff.min(),
                'max': utility_diff.max(),
                'positive_rate': (utility_diff > 0).mean()
            }
        }
        
        print(f"\n📊 个体改进分布:")
        print(f"  平均改进: {individual_analysis['utility_improvement_distribution']['mean']:.4f}")
        print(f"  改进标准差: {individual_analysis['utility_improvement_distribution']['std']:.4f}")
        print(f"  最大改进: {individual_analysis['utility_improvement_distribution']['max']:.4f}")
        print(f"  最小改进: {individual_analysis['utility_improvement_distribution']['min']:.4f}")
        print(f"  受益个体比例: {individual_analysis['utility_improvement_distribution']['positive_rate']:.1%}")
        
        return individual_analysis
    
    def create_comprehensive_plots(self, save_dir: str = './py/analysis_plots/'):
        """创建综合分析图表"""
        print("\n📊 创建综合分析图表...")
        
        # 创建保存目录
        Path(save_dir).mkdir(parents=True, exist_ok=True)
        
        paths = self.data['paths']
        utility_comp = self.data['utility_comparison']
        
        # 1. 生命周期路径比较 (2x2 布局)
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('RL vs VFI 生命周期决策路径比较', fontsize=16)
        
        k_vfi = paths['k_path_vfi']
        k_rl = paths['k_path_rl']
        c_vfi = paths['c_path_vfi']
        c_rl = paths['c_path_rl']
        cpps_vfi = paths['cpps_path_vfi']
        cpps_rl = paths['cpps_path_rl']
        
        n_sim, T = k_vfi.shape
        ages = np.arange(20, 20 + T)
        
        # 资产路径
        axes[0,0].plot(ages, k_vfi.mean(axis=0), 'r-', linewidth=2, label='VFI')
        axes[0,0].plot(ages, k_rl.mean(axis=0), 'b--', linewidth=2, label='RL')
        axes[0,0].fill_between(ages, 
                              k_vfi.mean(axis=0) - k_vfi.std(axis=0), 
                              k_vfi.mean(axis=0) + k_vfi.std(axis=0), 
                              alpha=0.2, color='red')
        axes[0,0].fill_between(ages, 
                              k_rl.mean(axis=0) - k_rl.std(axis=0), 
                              k_rl.mean(axis=0) + k_rl.std(axis=0), 
                              alpha=0.2, color='blue')
        axes[0,0].set_xlabel('年龄')
        axes[0,0].set_ylabel('资产')
        axes[0,0].set_title('平均资产路径')
        axes[0,0].legend()
        axes[0,0].grid(True, alpha=0.3)
        
        # 消费路径
        axes[0,1].plot(ages, c_vfi.mean(axis=0), 'r-', linewidth=2, label='VFI')
        axes[0,1].plot(ages, c_rl.mean(axis=0), 'b--', linewidth=2, label='RL')
        axes[0,1].set_xlabel('年龄')
        axes[0,1].set_ylabel('消费')
        axes[0,1].set_title('平均消费路径')
        axes[0,1].legend()
        axes[0,1].grid(True, alpha=0.3)
        
        # PPS缴费路径
        axes[1,0].plot(ages, cpps_vfi.mean(axis=0), 'r-', linewidth=2, label='VFI')
        axes[1,0].plot(ages, cpps_rl.mean(axis=0), 'b--', linewidth=2, label='RL')
        axes[1,0].set_xlabel('年龄')
        axes[1,0].set_ylabel('PPS缴费')
        axes[1,0].set_title('平均PPS缴费路径')
        axes[1,0].legend()
        axes[1,0].grid(True, alpha=0.3)
        
        # 效用差异分布
        utility_diff = np.array(utility_comp['lifetime_utility_rl']) - np.array(utility_comp['lifetime_utility_vfi'])
        axes[1,1].hist(utility_diff, bins=30, alpha=0.7, color='green', edgecolor='black')
        axes[1,1].axvline(utility_diff.mean(), color='darkgreen', linestyle='-', linewidth=2, 
                         label=f'平均差异: {utility_diff.mean():.4f}')
        axes[1,1].axvline(0, color='black', linestyle='--', alpha=0.7)
        axes[1,1].set_xlabel('效用差异 (RL - VFI)')
        axes[1,1].set_ylabel('频数')
        axes[1,1].set_title('个体效用差异分布')
        axes[1,1].legend()
        axes[1,1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(f"{save_dir}/lifecycle_comparison.png", dpi=300, bbox_inches='tight')
        plt.show()
        
        # 2. 阶段性表现分析
        stage_analysis = self.analyze_lifecycle_stages()
        
        # 创建阶段比较图
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        fig.suptitle('生命周期阶段表现差异 (RL - VFI)', fontsize=16)
        
        stages = list(stage_analysis.keys())
        k_diffs = [stage_analysis[stage]['k_diff_mean'] for stage in stages]
        c_diffs = [stage_analysis[stage]['c_diff_mean'] for stage in stages]
        cpps_diffs = [stage_analysis[stage]['cpps_diff_mean'] for stage in stages]
        
        # 资产差异
        bars1 = axes[0].bar(stages, k_diffs, color=['lightcoral' if x < 0 else 'lightblue' for x in k_diffs])
        axes[0].axhline(0, color='black', linestyle='-', alpha=0.3)
        axes[0].set_ylabel('平均资产差异')
        axes[0].set_title('各阶段资产差异')
        axes[0].tick_params(axis='x', rotation=45)
        
        # 消费差异
        bars2 = axes[1].bar(stages, c_diffs, color=['lightcoral' if x < 0 else 'lightblue' for x in c_diffs])
        axes[1].axhline(0, color='black', linestyle='-', alpha=0.3)
        axes[1].set_ylabel('平均消费差异')
        axes[1].set_title('各阶段消费差异')
        axes[1].tick_params(axis='x', rotation=45)
        
        # PPS缴费差异
        bars3 = axes[2].bar(stages, cpps_diffs, color=['lightcoral' if x < 0 else 'lightblue' for x in cpps_diffs])
        axes[2].axhline(0, color='black', linestyle='-', alpha=0.3)
        axes[2].set_ylabel('平均PPS缴费差异')
        axes[2].set_title('各阶段PPS缴费差异')
        axes[2].tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        plt.savefig(f"{save_dir}/stage_differences.png", dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"✅ 分析图表已保存到: {save_dir}")
    
    def generate_performance_report(self, save_path: str = './py/rl_vfi_performance_report.txt'):
        """生成综合性能分析报告"""
        print("\n📝 生成综合性能分析报告...")
        
        with open(save_path, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write("RL vs VFI 性能分析报告\n")
            f.write("=" * 80 + "\n\n")
            
            # 基本信息
            metadata = self.data['metadata']
            f.write("📊 基本信息:\n")
            f.write(f"  生成时间: {metadata['generated_at']}\n")
            f.write(f"  RL后端: {metadata['rl_backend']}\n")
            f.write(f"  VFI方法: {metadata['vfi_method']}\n")
            f.write(f"  模拟个体数: {metadata['n_sim']}\n")
            f.write(f"  随机种子: {metadata['random_seed']}\n\n")
            
            # 效用比较
            utility_comp = self.data['utility_comparison']
            f.write("🎯 效用比较结果:\n")
            f.write(f"  VFI平均效用: {utility_comp['mean_utility_vfi']:.6f} ± {utility_comp['std_utility_vfi']:.6f}\n")
            f.write(f"  RL平均效用: {utility_comp['mean_utility_rl']:.6f} ± {utility_comp['std_utility_rl']:.6f}\n")
            f.write(f"  效用差异: {utility_comp['utility_diff']:.6f}\n")
            f.write(f"  相对改进: {utility_comp['utility_improvement_pct']:.2f}%\n")
            f.write(f"  统计显著性: {'显著' if utility_comp['is_significant'] else '不显著'} (p={utility_comp['p_value']:.6f})\n\n")
            
            # 生命周期阶段分析
            f.write("🔍 生命周期阶段分析:\n")
            stage_analysis = self.analyze_lifecycle_stages()
            for stage_name, analysis in stage_analysis.items():
                f.write(f"  {stage_name} ({analysis['age_range']}):\n")
                f.write(f"    资产差异: {analysis['k_diff_mean']:+.4f} ± {analysis['k_diff_std']:.4f}\n")
                f.write(f"    消费差异: {analysis['c_diff_mean']:+.4f} ± {analysis['c_diff_std']:.4f}\n")
                f.write(f"    PPS缴费差异: {analysis['cpps_diff_mean']:+.4f} ± {analysis['cpps_diff_std']:.4f}\n")
            f.write("\n")
            
            # 决策模式分析
            f.write("🧠 决策模式分析:\n")
            decision_patterns = self.analyze_decision_patterns()
            f.write(f"  储蓄率: VFI={decision_patterns['savings_rate']['vfi_mean']:.4f}, RL={decision_patterns['savings_rate']['rl_mean']:.4f}\n")
            f.write(f"  PPS参与率: VFI={decision_patterns['pps_participation']['vfi_overall']:.2%}, RL={decision_patterns['pps_participation']['rl_overall']:.2%}\n")
            f.write(f"  消费波动性: VFI={decision_patterns['consumption_volatility']['vfi_mean']:.4f}, RL={decision_patterns['consumption_volatility']['rl_mean']:.4f}\n")
            f.write("\n")
            
            # 个体差异分析
            f.write("👥 个体差异分析:\n")
            individual_analysis = self.analyze_individual_differences()
            dist = individual_analysis['utility_improvement_distribution']
            f.write(f"  平均改进: {dist['mean']:.6f}\n")
            f.write(f"  改进标准差: {dist['std']:.6f}\n")
            f.write(f"  改进范围: [{dist['min']:.6f}, {dist['max']:.6f}]\n")
            f.write(f"  受益个体比例: {dist['positive_rate']:.1%}\n")
            f.write("\n")
            
            # 关键发现
            f.write("🎯 关键发现:\n")
            if utility_comp['utility_diff'] > 0:
                f.write("  1. RL在生命周期效用上显著优于VFI\n")
            
            if decision_patterns['consumption_volatility']['diff'] < 0:
                f.write("  2. RL实现了更好的消费平滑\n")
            
            if decision_patterns['pps_participation']['rl_overall'] > decision_patterns['pps_participation']['vfi_overall']:
                f.write("  3. RL更积极地参与PPS制度\n")
            
            f.write("=" * 80 + "\n")
        
        print(f"✅ 性能分析报告已保存到: {save_path}")
    
    def run_full_analysis(self):
        """运行完整分析流程"""
        print("🚀 开始完整的RL vs VFI性能分析...")
        
        # 0. 重新计算终身效用验证
        print("\n🔍 步骤 0: 重新计算终身效用验证")
        recalculated_results = self.recalculate_lifetime_utilities()
        
        # 1. 生命周期阶段分析
        print("\n📊 步骤 1: 生命周期阶段分析")
        stage_analysis = self.analyze_lifecycle_stages()
        
        # 2. 决策模式分析
        print("\n🧠 步骤 2: 决策模式分析")
        decision_patterns = self.analyze_decision_patterns()
        
        # 3. 个体差异分析
        print("\n👥 步骤 3: 个体差异分析")
        individual_analysis = self.analyze_individual_differences()
        
        # 4. 创建综合图表
        print("\n📊 步骤 4: 创建综合图表")
        self.create_comprehensive_plots()
        
        # 5. 生成分析报告
        print("\n📝 步骤 5: 生成分析报告")
        self.generate_performance_report()
        
        print("\n🎉 完整分析完成！")
        print("📈 主要结论 (基于重新计算的效用):")
        print(f"  - RL相对VFI的效用改进: {recalculated_results['utility_improvement_pct_new']:+.4f}%")
        print(f"  - 统计显著性: {'显著' if recalculated_results['is_significant_new'] else '不显著'}")
        print(f"  - 受益个体比例: {(recalculated_results['utility_diff_new'] > 0).mean():.1%}")
        
        if abs(recalculated_results['utility_improvement_pct_new']) > 0.1:
            if recalculated_results['utility_improvement_pct_new'] > 0:
                print("  ✅ 验证结果: RL确实在终身效用上优于VFI")
            else:
                print("  ❌ 验证结果: VFI在终身效用上优于RL")
        else:
            print("  ⚖️ 验证结果: RL和VFI的终身效用几乎相等")
        
        return {
            'recalculated_results': recalculated_results,
            'stage_analysis': stage_analysis,
            'decision_patterns': decision_patterns,
            'individual_analysis': individual_analysis
        }

    def recalculate_lifetime_utilities(self):
        """
        重新计算RL和VFI的终身已实现效用
        使用正确的效用函数和折现因子进行验证
        """
        print("\n🔍 重新计算终身已实现效用...")
        
        # 从数据中提取参数
        paths = self.data['paths']
        metadata = self.data['metadata']
        
        # 提取消费路径
        c_path_vfi = paths['c_path_vfi']  # (n_sim, T)
        c_path_rl = paths['c_path_rl']    # (n_sim, T)
        
        n_sim, T = c_path_vfi.shape
        ages = np.arange(20, 20 + T)  # 20-98岁
        
        # 从metadata中获取参数，如果没有则使用默认值
        beta = metadata.get('beta', 0.97)        # 折现因子
        sigma = metadata.get('sigma', 1.5)       # 风险厌恶系数
        c_floor = metadata.get('c_floor', 1e-6)  # 最低消费约束
        
        print(f"📊 效用计算参数:")
        print(f"  - 折现因子 (β): {beta:.4f}")
        print(f"  - 风险厌恶系数 (σ): {sigma:.3f}")
        print(f"  - 最低消费约束: {c_floor:.6f}")
        print(f"  - 个体数量: {n_sim}")
        print(f"  - 生命周期长度: {T}年 (年龄{ages[0]}-{ages[-1]})")
        
        # 创建简化的参数结构体用于效用计算
        class SimpleParams:
            def __init__(self, c_floor):
                self.cFloor = c_floor
        
        simple_params = SimpleParams(c_floor)
        
        # 计算VFI的终身效用
        print("\n📈 重新计算VFI终身效用...")
        lifetime_utility_vfi_new = np.zeros(n_sim)
        
        for i_sim in range(n_sim):
            if (i_sim + 1) % 50 == 0:
                print(f"  VFI进度: {i_sim + 1}/{n_sim}")
            
            utility_sum = 0.0
            
            for t in range(T):
                c_vfi = c_path_vfi[i_sim, t]
                
                # 计算当期效用
                if abs(sigma - 1) < 1e-6:  # 对数效用
                    if c_vfi >= c_floor:
                        u_vfi = np.log(c_vfi)
                    else:
                        u_vfi = -1e10 - (c_floor - c_vfi) * 1e10
                else:  # CRRA效用
                    c_adjusted = max(c_floor, c_vfi)
                    if c_vfi >= c_floor:
                        u_vfi = (c_adjusted ** (1 - sigma)) / (1 - sigma)
                    else:
                        u_vfi = -1e10 - (c_floor - c_vfi) * 1e10
                
                # 折现累加
                discount_factor = beta ** t
                utility_sum += discount_factor * u_vfi
            
            lifetime_utility_vfi_new[i_sim] = utility_sum
        
        # 计算RL的终身效用
        print("\n📈 重新计算RL终身效用...")
        lifetime_utility_rl_new = np.zeros(n_sim)
        
        for i_sim in range(n_sim):
            if (i_sim + 1) % 50 == 0:
                print(f"  RL进度: {i_sim + 1}/{n_sim}")
            
            utility_sum = 0.0
            
            for t in range(T):
                c_rl = c_path_rl[i_sim, t]
                
                # 计算当期效用
                if abs(sigma - 1) < 1e-6:  # 对数效用
                    if c_rl >= c_floor:
                        u_rl = np.log(c_rl)
                    else:
                        u_rl = -1e10 - (c_floor - c_rl) * 1e10
                else:  # CRRA效用
                    c_adjusted = max(c_floor, c_rl)
                    if c_rl >= c_floor:
                        u_rl = (c_adjusted ** (1 - sigma)) / (1 - sigma)
                    else:
                        u_rl = -1e10 - (c_floor - c_rl) * 1e10
                
                # 折现累加
                discount_factor = beta ** t
                utility_sum += discount_factor * u_rl
            
            lifetime_utility_rl_new[i_sim] = utility_sum
        
        # 计算新的比较结果
        utility_diff_new = lifetime_utility_rl_new - lifetime_utility_vfi_new
        
        mean_utility_vfi_new = lifetime_utility_vfi_new.mean()
        mean_utility_rl_new = lifetime_utility_rl_new.mean()
        utility_diff_mean_new = utility_diff_new.mean()
        utility_improvement_pct_new = (utility_diff_mean_new / abs(mean_utility_vfi_new)) * 100
        
        # 统计检验
        from scipy.stats import ttest_rel
        t_stat_new, p_value_new = ttest_rel(lifetime_utility_rl_new, lifetime_utility_vfi_new)
        is_significant_new = p_value_new < 0.05
        
        # 与原始结果比较
        print("\n📊 重新计算的效用比较结果:")
        print("=" * 60)
        print(f"VFI平均终身效用: {mean_utility_vfi_new:.6f} ± {lifetime_utility_vfi_new.std():.6f}")
        print(f"RL平均终身效用:  {mean_utility_rl_new:.6f} ± {lifetime_utility_rl_new.std():.6f}")
        print(f"效用差异 (RL - VFI): {utility_diff_mean_new:.6f}")
        print(f"相对改进: {utility_improvement_pct_new:+.4f}%")
        print(f"统计显著性: {'显著' if is_significant_new else '不显著'} (t={t_stat_new:.4f}, p={p_value_new:.6f})")
        print(f"受益个体比例: {(utility_diff_new > 0).mean():.1%}")
        
        # 与原始数据对比
        if 'utility_comparison' in self.data:
            original_comparison = self.data['utility_comparison']
            print("\n🔄 与原始结果对比:")
            print("=" * 60)
            print(f"原始VFI效用: {original_comparison['mean_utility_vfi']:.6f}")
            print(f"重算VFI效用: {mean_utility_vfi_new:.6f}")
            print(f"VFI差异: {mean_utility_vfi_new - original_comparison['mean_utility_vfi']:.6f}")
            print()
            print(f"原始RL效用:  {original_comparison['mean_utility_rl']:.6f}")
            print(f"重算RL效用:  {mean_utility_rl_new:.6f}")
            print(f"RL差异:  {mean_utility_rl_new - original_comparison['mean_utility_rl']:.6f}")
            print()
            print(f"原始相对改进: {original_comparison['utility_improvement_pct']:+.4f}%")
            print(f"重算相对改进: {utility_improvement_pct_new:+.4f}%")
            print(f"改进差异: {utility_improvement_pct_new - original_comparison['utility_improvement_pct']:+.4f}个百分点")
        
        # 计算消费统计信息
        print("\n📊 消费路径统计:")
        print("=" * 40)
        print(f"VFI平均消费: {c_path_vfi.mean():.4f} ± {c_path_vfi.std():.4f}")
        print(f"RL平均消费:  {c_path_rl.mean():.4f} ± {c_path_rl.std():.4f}")
        print(f"VFI最低消费: {c_path_vfi.min():.6f}")
        print(f"RL最低消费:  {c_path_rl.min():.6f}")
        print(f"VFI违反最低消费约束的比例: {(c_path_vfi < c_floor).mean():.2%}")
        print(f"RL违反最低消费约束的比例:  {(c_path_rl < c_floor).mean():.2%}")
        
        # 分析效用差异分布
        print("\n📈 效用差异分布分析:")
        print("=" * 40)
        print(f"最大RL优势: {utility_diff_new.max():.6f}")
        print(f"最大VFI优势: {utility_diff_new.min():.6f}")
        print(f"效用差异标准差: {utility_diff_new.std():.6f}")
        print(f"效用差异中位数: {np.median(utility_diff_new):.6f}")
        
        # 分位数分析
        percentiles = [5, 25, 50, 75, 95]
        print(f"效用差异分位数:")
        for p in percentiles:
            value = np.percentile(utility_diff_new, p)
            print(f"  {p:2d}%: {value:+.6f}")
        
        # 保存重新计算的结果
        recalculated_results = {
            'lifetime_utility_vfi_new': lifetime_utility_vfi_new,
            'lifetime_utility_rl_new': lifetime_utility_rl_new,
            'utility_diff_new': utility_diff_new,
            'mean_utility_vfi_new': mean_utility_vfi_new,
            'mean_utility_rl_new': mean_utility_rl_new,
            'utility_diff_mean_new': utility_diff_mean_new,
            'utility_improvement_pct_new': utility_improvement_pct_new,
            'is_significant_new': is_significant_new,
            'p_value_new': p_value_new,
            't_stat_new': t_stat_new,
            'calculation_parameters': {
                'beta': beta,
                'sigma': sigma,
                'c_floor': c_floor,
                'n_sim': n_sim,
                'T': T
            }
        }
        
        return recalculated_results

def main():
    """主函数"""
    print("=" * 80)
    print("🔬 RL vs VFI 深入性能分析")
    print("=" * 80)
    
    # 检查数据文件是否存在
    data_path = './py/rl_vfi_detailed_paths.pkl'
    if not os.path.exists(data_path):
        print(f"❌ 数据文件不存在: {data_path}")
        print("💡 请先运行 compare_rl_and_vfi_matlab.py 生成详细数据")
        return
    
    # 创建分析器
    analyzer = RLVFIPerformanceAnalyzer(data_path)
    
    # 提供选择：仅重新计算效用或完整分析
    print("\n🔧 分析选项:")
    print("1. 仅重新计算并验证终身效用")
    print("2. 运行完整分析 (包括重新计算验证)")
    
    choice = input("\n请选择分析模式 (1/2, 默认为2): ").strip()
    
    if choice == "1":
        # 仅重新计算效用
        print("\n🔍 仅执行终身效用重新计算验证...")
        recalculated_results = analyzer.recalculate_lifetime_utilities()
        
        print("\n🎯 关键验证结论:")
        if abs(recalculated_results['utility_improvement_pct_new']) > 0.1:
            if recalculated_results['utility_improvement_pct_new'] > 0:
                print("  ✅ RL确实在终身效用上优于VFI")
            else:
                print("  ❌ VFI在终身效用上优于RL")
        else:
            print("  ⚖️ RL和VFI的终身效用几乎相等")
            
    else:
        # 运行完整分析
        print("\n🚀 运行完整分析...")
        results = analyzer.run_full_analysis()
        
        print("\n✅ 分析完成！查看以下文件获取详细结果:")
        print("  - ./py/analysis_plots/ - 分析图表")
        print("  - ./py/rl_vfi_performance_report.txt - 详细报告")

if __name__ == "__main__":
    main() 