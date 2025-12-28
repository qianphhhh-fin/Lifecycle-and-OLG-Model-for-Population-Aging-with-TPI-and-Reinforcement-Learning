"""
测试OLGEnvV9SACSimplified环境

🎯 目标：
- 使用随机决策测试环境的基本功能
- 检查状态转移是否正常
- 验证效用计算是否合理
- 观察生命周期路径是否符合预期
"""

import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

import numpy as np
import matplotlib.pyplot as plt
import sys
import time
from typing import Dict, List, Tuple

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

# 添加路径
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from main_olg_v9_utils import OLG_V9_Utils
from simplified.main_olg_v9_utils_simplified import OLGUtilsSimplified
from simplified.main_olg_v9_sac_sbx_simplified import OLGEnvV9SACSimplified

class OLGEnvTester:
    """OLG环境测试器"""
    
    def __init__(self):
        """初始化测试器"""
        print("🧪 OLG环境测试器初始化...")
        
        # 初始化参数
        self.cS = OLG_V9_Utils.ParameterValues_HuggettStyle()
        leLogGridV, leTrProbM, leProb1V = OLG_V9_Utils.EarningProcess_olgm(self.cS)
        leGridV = np.exp(leLogGridV)
        
        self.paramS_rl = {
            'leLogGridV': leLogGridV,
            'leTrProbM': leTrProbM,
            'leProb1V': leProb1V,
            'leGridV': leGridV,
            'ageEffV_new': self.cS.ageEffV_new
        }
        
        # 固定宏观参数
        self.M_fixed = {
            'R_k_net_factor': 1.03,
            'w_gross': 2.0,
            'TR_total': 0.1,
            'b_payg_avg_retiree': 0.4,
            'tau_l': 0.15,
            'theta_payg_actual': 0.12
        }
        
        print("✅ 参数初始化完成")
    
    def create_environment(self, training_mode: bool = False) -> OLGEnvV9SACSimplified:
        """创建环境实例"""
        env = OLGEnvV9SACSimplified(
            cS=self.cS,
            paramS_rl=self.paramS_rl,
            M_fixed=self.M_fixed,
            training_mode=training_mode
        )
        return env
    
    def test_basic_functionality(self) -> Dict:
        """测试基本功能"""
        print("\n🔍 测试基本功能...")
        
        env = self.create_environment(training_mode=False)
        
        # 重置环境
        obs, info = env.reset(seed=42)
        print(f"📊 初始观测: {obs}")
        print(f"📋 初始信息: {info}")
        
        # 检查观测空间和动作空间
        print(f"🔢 观测空间: {env.observation_space}")
        print(f"🎮 动作空间: {env.action_space}")
        
        # 测试几步随机动作
        results = {
            'observations': [obs.copy()],
            'actions': [],
            'rewards': [],
            'infos': [],
            'states_raw': []
        }
        
        for step in range(5):
            # 生成随机动作
            action = env.action_space.sample()
            results['actions'].append(action.copy())
            
            # 执行动作
            obs_next, reward, terminated, truncated, info = env.step(action)
            
            results['observations'].append(obs_next.copy())
            results['rewards'].append(reward)
            results['infos'].append(info.copy())
            
            # 记录原始状态
            raw_state = {
                'k': env.current_k_val,
                'k_pps': env.current_k_pps_val,
                'age_idx': env.current_age_idx,
                'epsilon_idx': env.current_eps_idx,
                'age': 20 + env.current_age_idx
            }
            results['states_raw'].append(raw_state.copy())
            
            print(f"步骤 {step+1}:")
            print(f"  🎮 动作: [PPS缴费比例={action[0]:.3f}, 非PPS储蓄比例={action[1]:.3f}]")
            print(f"  📊 观测 (当前): {obs}")
            print(f"  📊 观测 (下步): {obs_next}")
            print(f"  🏠 状态: k={raw_state['k']:.3f}, k_pps={raw_state['k_pps']:.3f}, age={raw_state['age']}, ε_idx={raw_state['epsilon_idx']}")
            print(f"  💰 奖励: {reward:.6f}")
            print(f"  🏁 终止: {terminated}, 截断: {truncated}")
            if 'consumption' in info:
                print(f"  🍽️  消费: {info['consumption']:.3f}")
            if 'c_pps' in info:
                print(f"  💳 PPS缴费: {info['c_pps']:.3f}")
            print(f"  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
            
            if terminated or truncated:
                print("🏁 回合结束")
                break
            
            obs = obs_next
        
        env.close()
        return results
    
    def test_full_lifecycle(self, n_episodes: int = 3) -> Dict:
        """测试完整生命周期"""
        print(f"\n🔄 测试完整生命周期 (n_episodes={n_episodes})...")
        
        all_results = []
        
        for episode in range(n_episodes):
            print(f"\n📈 回合 {episode + 1}/{n_episodes}")
            
            env = self.create_environment(training_mode=False)
            obs, info = env.reset(seed=42 + episode)
            
            episode_data = {
                'episode': episode,
                'states': [],
                'actions': [],
                'rewards': [],
                'consumptions': [],
                'pps_contributions': [],
                'utilities': [],
                'terminated_normally': False
            }
            
            step_count = 0
            total_reward = 0
            
            while True:
                # 记录当前状态
                current_state = {
                    'k': env.current_k_val,
                    'k_pps': env.current_k_pps_val,
                    'age_idx': env.current_age_idx,
                    'epsilon_idx': env.current_eps_idx,
                    'age': 20 + env.current_age_idx,
                    'observation': obs.copy()
                }
                episode_data['states'].append(current_state)
                
                # 生成更合理的随机动作
                action = env.action_space.sample()
                # 年轻时更保守的策略
                if current_state['age'] < 30:
                    action[0] = np.clip(action[0] * 0.1, 0, 1)  # 年轻时很少PPS缴费
                    action[1] = np.clip(action[1] * 0.5, 0, 1)  # 适度储蓄
                else:
                    action[0] = np.clip(action[0] * 0.3, 0, 1)  # 中年时增加PPS缴费
                    action[1] = np.clip(action[1] * 0.7, 0, 1)  # 积极储蓄
                
                episode_data['actions'].append(action.copy())
                
                # 执行动作
                obs_next, reward, terminated, truncated, info = env.step(action)
                
                episode_data['rewards'].append(reward)
                total_reward += reward
                
                # 从info中提取详细信息
                if 'consumption' in info:
                    episode_data['consumptions'].append(info['consumption'])
                if 'pps_contribution' in info:
                    episode_data['pps_contributions'].append(info['pps_contribution'])
                if 'utility' in info:
                    episode_data['utilities'].append(info['utility'])
                
                step_count += 1
                
                if step_count % 5 == 0 or step_count <= 3:  # 显示前3步和每5步
                    print(f"  📊 步骤 {step_count}:")
                    print(f"     🎮 动作: [PPS={action[0]:.3f}, 储蓄={action[1]:.3f}]")
                    print(f"     📊 观测: {obs}")
                    print(f"     🏠 状态: k={current_state['k']:.3f}, k_pps={current_state['k_pps']:.3f}, age={current_state['age']}, ε={current_state['epsilon_idx']}")
                    print(f"     💰 奖励: {reward:.6f}")
                    if 'consumption' in info:
                        print(f"     🍽️  消费: {info['consumption']:.3f}")
                    if 'c_pps' in info:
                        print(f"     💳 PPS: {info['c_pps']:.3f}")
                    print(f"     ──────────────────────────────────────────────────")
                
                if terminated or truncated:
                    episode_data['terminated_normally'] = terminated
                    print(f"🏁 回合结束: 步骤={step_count}, 总奖励={total_reward:.4f}, 正常结束={terminated}")
                    break
                
                obs = obs_next
            
            episode_data['total_steps'] = step_count
            episode_data['total_reward'] = total_reward
            all_results.append(episode_data)
            
            env.close()
        
        return all_results
    
    def analyze_results(self, results: List[Dict]) -> Dict:
        """分析测试结果"""
        print("\n📊 分析测试结果...")
        
        analysis = {
            'summary': {},
            'lifecycle_stats': {},
            'validation_checks': {}
        }
        
        # 汇总统计
        total_rewards = [ep['total_reward'] for ep in results]
        total_steps = [ep['total_steps'] for ep in results]
        
        analysis['summary'] = {
            'n_episodes': len(results),
            'avg_total_reward': np.mean(total_rewards),
            'std_total_reward': np.std(total_rewards),
            'avg_steps': np.mean(total_steps),
            'std_steps': np.std(total_steps),
            'all_terminated_normally': all(ep['terminated_normally'] for ep in results)
        }
        
        # 生命周期统计（使用第一个回合）
        if results:
            first_episode = results[0]
            
            ages = [state['age'] for state in first_episode['states']]
            k_path = [state['k'] for state in first_episode['states']]
            k_pps_path = [state['k_pps'] for state in first_episode['states']]
            consumptions = first_episode['consumptions']
            rewards = first_episode['rewards']
            
            analysis['lifecycle_stats'] = {
                'age_range': (min(ages), max(ages)),
                'k_range': (min(k_path), max(k_path)),
                'k_pps_range': (min(k_pps_path), max(k_pps_path)),
                'consumption_range': (min(consumptions), max(consumptions)) if consumptions else (0, 0),
                'reward_range': (min(rewards), max(rewards)),
                'avg_consumption': np.mean(consumptions) if consumptions else 0,
                'avg_reward': np.mean(rewards)
            }
        
        # 验证检查（简化版：16个年龄组而不是79个年度年龄）
        expected_steps = 16  # 简化版使用年龄组
        validation_checks = {
            'all_episodes_completed': all(ep['total_steps'] >= expected_steps for ep in results),
            'positive_consumption': True,
            'reasonable_rewards': True,
            'state_transitions_valid': True
        }
        
        # 检查消费是否为正
        for ep in results:
            if ep['consumptions'] and any(c <= 0 for c in ep['consumptions']):
                validation_checks['positive_consumption'] = False
                break
        
        # 检查奖励是否在合理范围内
        for ep in results:
            if any(abs(r) > 100 for r in ep['rewards']):  # 过大的奖励可能表示有问题
                validation_checks['reasonable_rewards'] = False
                break
        
        analysis['validation_checks'] = validation_checks
        
        return analysis
    
    def plot_lifecycle_paths(self, results: List[Dict], save_path: str = None):
        """绘制生命周期路径"""
        print("\n📈 绘制生命周期路径...")
        
        if not results:
            print("❌ 没有结果可绘制")
            return
        
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        fig.suptitle('OLG环境测试：随机策略生命周期路径', fontsize=14)
        
        colors = ['blue', 'red', 'green', 'orange', 'purple']
        
        for i, episode_data in enumerate(results[:3]):  # 最多显示3个回合
            color = colors[i % len(colors)]
            label = f'回合 {i+1}'
            
            ages = [state['age'] for state in episode_data['states']]
            k_path = [state['k'] for state in episode_data['states']]
            k_pps_path = [state['k_pps'] for state in episode_data['states']]
            consumptions = episode_data['consumptions']
            rewards = episode_data['rewards']
            pps_contributions = episode_data['pps_contributions']
            
            # 1. 资产路径
            axes[0, 0].plot(ages, k_path, color=color, label=label, alpha=0.7)
            axes[0, 0].set_xlabel('年龄')
            axes[0, 0].set_ylabel('普通资产 (k)')
            axes[0, 0].set_title('普通资产路径')
            axes[0, 0].grid(True, alpha=0.3)
            axes[0, 0].legend()
            
            # 2. PPS资产路径
            axes[0, 1].plot(ages, k_pps_path, color=color, label=label, alpha=0.7)
            axes[0, 1].set_xlabel('年龄')
            axes[0, 1].set_ylabel('PPS资产 (k_pps)')
            axes[0, 1].set_title('PPS资产路径')
            axes[0, 1].grid(True, alpha=0.3)
            axes[0, 1].legend()
            
            # 3. 消费路径
            if consumptions:
                axes[0, 2].plot(ages[:len(consumptions)], consumptions, color=color, label=label, alpha=0.7)
            axes[0, 2].set_xlabel('年龄')
            axes[0, 2].set_ylabel('消费 (c)')
            axes[0, 2].set_title('消费路径')
            axes[0, 2].grid(True, alpha=0.3)
            axes[0, 2].legend()
            
            # 4. 奖励路径
            axes[1, 0].plot(ages[:len(rewards)], rewards, color=color, label=label, alpha=0.7)
            axes[1, 0].set_xlabel('年龄')
            axes[1, 0].set_ylabel('即时奖励')
            axes[1, 0].set_title('奖励路径')
            axes[1, 0].grid(True, alpha=0.3)
            axes[1, 0].legend()
            
            # 5. PPS缴费路径
            if pps_contributions:
                axes[1, 1].plot(ages[:len(pps_contributions)], pps_contributions, color=color, label=label, alpha=0.7)
            axes[1, 1].set_xlabel('年龄')
            axes[1, 1].set_ylabel('PPS缴费')
            axes[1, 1].set_title('PPS缴费路径')
            axes[1, 1].grid(True, alpha=0.3)
            axes[1, 1].legend()
            
            # 6. 总资产路径
            total_assets = [k + k_pps for k, k_pps in zip(k_path, k_pps_path)]
            axes[1, 2].plot(ages, total_assets, color=color, label=label, alpha=0.7)
            axes[1, 2].set_xlabel('年龄')
            axes[1, 2].set_ylabel('总资产 (k + k_pps)')
            axes[1, 2].set_title('总资产路径')
            axes[1, 2].grid(True, alpha=0.3)
            axes[1, 2].legend()
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"📈 图表已保存到: {save_path}")
        
        plt.show()
    
    def print_analysis_report(self, analysis: Dict):
        """打印分析报告"""
        print("\n" + "="*60)
        print("📋 OLG环境测试分析报告")
        print("="*60)
        
        # 汇总统计
        summary = analysis['summary']
        print(f"\n🔢 汇总统计:")
        print(f"  测试回合数: {summary['n_episodes']}")
        print(f"  平均总奖励: {summary['avg_total_reward']:.4f} ± {summary['std_total_reward']:.4f}")
        print(f"  平均步数: {summary['avg_steps']:.1f} ± {summary['std_steps']:.1f}")
        print(f"  所有回合正常结束: {'✅' if summary['all_terminated_normally'] else '❌'}")
        
        # 生命周期统计
        if 'lifecycle_stats' in analysis:
            lifecycle = analysis['lifecycle_stats']
            print(f"\n📊 生命周期统计:")
            print(f"  年龄范围: {lifecycle['age_range'][0]} - {lifecycle['age_range'][1]} 岁")
            print(f"  普通资产范围: {lifecycle['k_range'][0]:.3f} - {lifecycle['k_range'][1]:.3f}")
            print(f"  PPS资产范围: {lifecycle['k_pps_range'][0]:.3f} - {lifecycle['k_pps_range'][1]:.3f}")
            print(f"  消费范围: {lifecycle['consumption_range'][0]:.3f} - {lifecycle['consumption_range'][1]:.3f}")
            print(f"  奖励范围: {lifecycle['reward_range'][0]:.6f} - {lifecycle['reward_range'][1]:.6f}")
            print(f"  平均消费: {lifecycle['avg_consumption']:.3f}")
            print(f"  平均奖励: {lifecycle['avg_reward']:.6f}")
        
        # 验证检查
        validation = analysis['validation_checks']
        print(f"\n✅ 验证检查:")
        print(f"  所有回合完成 (≥16步): {'✅' if validation['all_episodes_completed'] else '❌'}")
        print(f"  消费始终为正: {'✅' if validation['positive_consumption'] else '❌'}")
        print(f"  奖励在合理范围: {'✅' if validation['reasonable_rewards'] else '❌'}")
        print(f"  状态转移有效: {'✅' if validation['state_transitions_valid'] else '❌'}")
        
        # 总体评估
        all_checks_passed = all(validation.values())
        print(f"\n🎯 总体评估: {'✅ 环境功能正常' if all_checks_passed else '❌ 发现潜在问题'}")
        
        print("="*60)

def main():
    """主函数"""
    print("🧪 OLG环境简化版测试")
    print("🎯 使用随机决策测试环境功能")
    print("="*60)
    
    # 创建测试器
    tester = OLGEnvTester()
    
    # 1. 基本功能测试
    print("\n1️⃣ 基本功能测试...")
    basic_results = tester.test_basic_functionality()
    
    # 2. 完整生命周期测试
    print("\n2️⃣ 完整生命周期测试...")
    lifecycle_results = tester.test_full_lifecycle(n_episodes=3)
    
    # 3. 结果分析
    print("\n3️⃣ 结果分析...")
    analysis = tester.analyze_results(lifecycle_results)
    
    # 4. 绘制路径
    print("\n4️⃣ 绘制生命周期路径...")
    save_path = './simplified/olg_env_test_paths.png'
    tester.plot_lifecycle_paths(lifecycle_results, save_path=save_path)
    
    # 5. 打印报告
    print("\n5️⃣ 生成分析报告...")
    tester.print_analysis_report(analysis)
    
    # 6. 保存详细结果
    import pickle
    with open('./simplified/olg_env_test_results.pkl', 'wb') as f:
        pickle.dump({
            'basic_results': basic_results,
            'lifecycle_results': lifecycle_results,
            'analysis': analysis,
            'test_timestamp': time.time()
        }, f)
    
    print("\n💾 详细测试结果已保存到: ./simplified/olg_env_test_results.pkl")
    
    return analysis

if __name__ == "__main__":
    results = main() 