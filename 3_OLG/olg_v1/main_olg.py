import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

class FullOLGModel:
    def __init__(self):
        # 模型参数初始化
        self.periods = 30  # 总模拟时期数
        self.age_groups = 16  # 16个五年期年龄组：22-101岁
        
        # 人口参数
        self.x = -0.02  # 新生代人口增长率调整因子
        
        # 存活概率：从一个年龄组到下一个年龄组的存活概率
        # 进一步提高高龄组存活率以增加老龄人口
        self.beta = np.array([
            0.995, 0.99, 0.985, 0.98, 0.975, 0.97, 0.965, 0.96, 
            0.95, 0.94, 0.92, 0.89, 0.85, 0.80, 0.75  # 显著提高高龄存活率
        ])
        
        # 养老金系统参数
        self.tau_g = 0.16  # 进一步降低公共养老金缴费率，从0.15降至0.10
        self.lambda_pension = 0.45 # 进一步提高养老金替代率，从0.65提高到0.75
        
        # 个人养老金参数
        self.tau_p = 0.03  # 个人养老金提取税率
        self.phi = 0.005  # 费率优惠因子
        
        # 生产函数参数
        self.alpha = 0.36  # 资本产出弹性
        self.delta = 0.08  # 资本折旧率
        self.A = 1.0  # 初始技术水平
        self.growth_rate = 0.008  # 进一步降低技术进步率，从0.01降至0.008
        
        # 劳动效率参数 - 进一步降低劳动效率差异
        # 0-7组(22-61岁)为工作人口，8-15组(62-101岁)为退休人口
        self.h = np.array([
            2.2, 2.5, 2.7, 2.9, 3.0, 3.1, 3.1, 3.0,  # 进一步缩小劳动效率差距
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0   # 退休年龄组 (62-101岁)
        ])
        
        # 效用函数参数
        self.gamma = 2.0  # 相对风险厌恶系数
        self.theta = 0.04  # 时间偏好率
        
        # 市场参数
        self.Rf = 1.02  # 无风险收益率
        self.Rp_mean = 1.04  # 风险资产收益率期望
        self.Rp_std = 0.18  # 风险资产收益率标准差
        
        # 政府参数
        self.tau_c = 0.12  # 消费税率
        self.G_ratio = 0.18  # 政府消费占GDP比例
        
        # 退休参数
        self.retirement_age = 8  # 退休年龄组索引 (对应62岁开始)
    
    def initialize_state(self):
        # 初始化状态变量 - 使用2023年中国人口结构数据
        initial_pop = np.array([
            [76.20905535, 86.45596319, 113.8702459, 98.60198303, 86.64117824, 102.7890433, 112.0217783, 99.04620047,
             64.05142331, 66.93157492, 44.16815149, 25.40848066, 14.97325553, 6.872421945, 1.743059943, 0.216184341]
        ])
        
        self.Z = initial_pop  # 人口分布，行为时期，列为年龄组
        self.K = np.array([1200.0])  # 初始资本存量
        
        # 个人养老金资产分布 - 16个年龄组，降低初始积累
        initial_pension = np.zeros(16)
        for i in range(8):  # 工作年龄组0-7
            initial_pension[i] = 5.0 * (i + 1)  # 降低初始积累
        
        self.P = np.array([initial_pension])  # 个人养老金资产
        self.D = np.array([5.0])  # 初始公共养老金累计赤字
        
        # 存储模拟结果
        self.r_history = []  # 利率历史
        self.w_history = []  # 工资率历史
        self.Y_history = []  # 产出历史
        self.C_history = []  # 消费历史
        self.T_history = []  # 养老金收入历史
        self.B_history = []  # 养老金支出历史
        self.deficit_history = []  # 养老金赤字历史
        self.replacement_rate_history = []  # 替代率历史
        self.dependency_ratio_history = []  # 抚养比历史
    
    def update_population(self, t):
        # 更新人口分布 - 16个年龄组
        Z_new = np.zeros(16)
        
        # 最年轻组 (22-26岁) - 新生人口
        # 随时间变化的增长率：前期小负增长，后期大负增长
        growth_rate = self.x
        if t < 6:
            # 前6期较缓和的负增长
            growth_rate = -0.01 - 0.003 * t  # 从-1%逐渐增大负增长率
        elif t >= 6:
            # 后期强化负增长
            growth_rate = -0.03 - 0.004 * min(t-6, 10)  # 最终达到-7%的负增长率
        
        Z_new[0] = self.Z[t, 0] * (1 + growth_rate)
        
        # 其他年龄组 - 从上一组存活下来的人口
        for a in range(1, 16):
            Z_new[a] = self.Z[t, a-1] * self.beta[a-1]
        
        self.Z = np.vstack((self.Z, Z_new))
    
    def production_function(self, K, L):
        # Cobb-Douglas 生产函数
        return K**self.alpha * (self.A * L)**(1-self.alpha)
    
    def calculate_factor_prices(self, t):
        # 计算有效劳动力 - 仅考虑工作年龄组 (0-7)
        L = 0
        for a in range(8):  # 0-7组为工作人口
            L += self.h[a] * self.Z[t, a]
        
        # 计算总产出
        Y = self.production_function(self.K[t], L)
        
        # 要素价格：工资和利率
        w = (1 - self.alpha) * Y / L
        r = self.alpha * Y / self.K[t] - self.delta
        
        return Y, w, r, L
    
    def personal_pension_account(self, t, w, r):
        # 简化的个人养老金决策，降低缴费率
        # 工作人口各年龄组的缴费比例 - 降低为原来的70%
        q_rates = np.array([0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.06])
        
        # 风险资产配置比例（降低为20%）
        xi = 0.20
        
        # 模拟风险资产收益
        Rp = self.Rp_mean + np.random.normal(0, self.Rp_std)
        
        # 个人养老金账户收益率
        R_pension = xi * (Rp + self.phi) + (1 - xi) * self.Rf
        
        # 初始化新的养老金账户余额和退休金收益
        P_new = np.zeros(16)
        personal_pension_benefits = np.zeros(16)
        personal_pension_taxes = np.zeros(16)
        
        # 更新工作人口的养老金账户 (0-7组)
        for a in range(8):
            P_new[a] = self.P[t, a] * R_pension + q_rates[a] * w * self.h[a]
        
        # 更新退休人口的养老金账户 (8-15组)，加快提取速度
        withdrawal_rate = 0.15  # 提高提取率，从0.125提高到0.15
        for a in range(8, 16):
            # 退休者从账户中提取资金
            withdrawal = self.P[t, a] * withdrawal_rate
            personal_pension_benefits[a] = withdrawal * (1 - self.tau_p)
            personal_pension_taxes[a] = withdrawal * self.tau_p
            P_new[a] = max(0, self.P[t, a] * (1 - withdrawal_rate))
        
        # 计算总个人养老金收益和税收
        total_benefit = sum(personal_pension_benefits * self.Z[t, :])
        total_tax = sum(personal_pension_taxes * self.Z[t, :])
        
        return P_new, total_benefit, total_tax, personal_pension_benefits
    
    def consumption_decision(self, t, w, r, pension_benefit, personal_pension_benefits):
        # 简化的消费决策 (各年龄组消费为收入的固定比例)
        # 提高消费比例以减少储蓄
        
        # 初始化各年龄组消费
        C_by_age = np.zeros(16)
        
        # 工作人口消费 (0-7组)，提高消费比例
        for a in range(8):
            income = w * self.h[a]
            C_by_age[a] = (0.85 - 0.01 * a) * income  # 提高消费比例
        
        # 退休人口消费 (8-15组)
        for a in range(8, 16):
            income = pension_benefit + personal_pension_benefits[a]
            C_by_age[a] = 0.95 * income  # 提高退休者消费比例
        
        # 总消费
        C_total = sum(C_by_age * self.Z[t, :])
        
        return C_total
    
    def pension_system(self, t, w):
        # 公共养老金收入 - 来自工作人口 (0-7组)
        T_g = 0
        for a in range(8):  # 工作年龄组
            T_g += self.tau_g * w * self.h[a] * self.Z[t, a]
        
        # 退休时的工资率（简化为当期工资）
        retirement_wage = w
        
        # 公共养老金支出 - 支付给退休人口 (8-15组)
        pension_benefit = self.lambda_pension * retirement_wage
        B_g = 0
        for a in range(8, 16):  # 退休年龄组
            B_g += pension_benefit * self.Z[t, a]
        
        # 养老金赤字
        deficit = max(0, B_g - T_g)
        
        return T_g, B_g, deficit, pension_benefit
    
    def capital_market_clearing(self, t, Y, C, deficit, personal_pension_tax):
        # 政府消费
        G = self.G_ratio * Y
        
        # 投资
        I = Y - C - G
        
        # 更新资本存量
        K_next = (1 - self.delta) * self.K[t] + I
        
        # 更新养老金累计赤字
        D_next = self.D[t] * (1 + self.Rf - 1) + deficit
        
        return K_next, D_next
    
    def calculate_dependency_ratio(self, t):
        # 计算抚养比：退休人口/工作人口
        retired_population = sum(self.Z[t, 8:])
        working_population = sum(self.Z[t, :8])
        
        if working_population > 0:
            return retired_population / working_population
        else:
            return float('inf')  # 避免除零错误
    
    def run_simulation(self):
        # 初始化状态
        self.initialize_state()
        
        # 运行模拟
        for t in range(self.periods):
            # 计算当期的生产和要素价格
            Y, w, r, L = self.calculate_factor_prices(t)
            
            # 更新技术水平
            self.A *= (1 + self.growth_rate)
            
            # 公共养老金系统
            T_g, B_g, deficit, pension_benefit = self.pension_system(t, w)
            
            # 个人养老金账户更新
            P_new, personal_pension_benefit, personal_pension_tax, personal_pension_benefits = self.personal_pension_account(t, w, r)
            
            # 消费决策
            C = self.consumption_decision(t, w, r, pension_benefit, personal_pension_benefits)
            
            # 资本市场出清
            K_next, D_next = self.capital_market_clearing(t, Y, C, deficit, personal_pension_tax)
            
            # 更新状态变量
            self.K = np.append(self.K, K_next)
            self.D = np.append(self.D, D_next)
            self.P = np.vstack((self.P, P_new))
            
            # 更新人口
            if t < self.periods - 1:
                self.update_population(t)
            
            # 计算抚养比
            dependency_ratio = self.calculate_dependency_ratio(t)
            
            # 记录历史数据
            self.r_history.append(r)
            self.w_history.append(w)
            self.Y_history.append(Y)
            self.C_history.append(C)
            self.T_history.append(T_g)
            self.B_history.append(B_g)
            self.deficit_history.append(deficit)
            self.dependency_ratio_history.append(dependency_ratio)
            
            # 计算总替代率（公共+个人）
            # 使用人均养老金收益计算替代率
            avg_personal_pension = personal_pension_benefit / max(1, sum(self.Z[t, 8:]))
            total_replacement_rate = (pension_benefit + avg_personal_pension) / w
            self.replacement_rate_history.append(total_replacement_rate)
            
            print(f"Period {t}: Y={Y:.2f}, w={w:.2f}, r={r:.4f}, deficit={deficit:.2f}, dependency_ratio={dependency_ratio:.2f}")
    
    def plot_results(self):
        # 绘制模拟结果
        t = np.arange(self.periods)
        
        # 设置中文显示支持
        plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
        plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号
        
        # 计算各人口组数量和比例
        young_workers = np.sum(self.Z[:self.periods, 0:4], axis=1)  # 22-41岁
        older_workers = np.sum(self.Z[:self.periods, 4:8], axis=1)  # 42-61岁
        retirees = np.sum(self.Z[:self.periods, 8:], axis=1)  # 62-101岁
        total_population = young_workers + older_workers + retirees
        
        young_ratio = young_workers / total_population * 100
        older_ratio = older_workers / total_population * 100
        retirees_ratio = retirees / total_population * 100
        
        # 创建图表
        fig = plt.figure(figsize=(16, 12))
        
        # 1. 养老金赤字主图
        plt.subplot(2, 2, 1)
        plt.plot(t, self.T_history, label='养老金收入', color='blue', linewidth=2.5)
        plt.plot(t, self.B_history, label='养老金支出', color='orange', linewidth=2.5)
        plt.fill_between(t, self.B_history, self.T_history, 
                        where=(np.array(self.B_history) > np.array(self.T_history)), 
                        color='red', alpha=0.4, interpolate=True, label='养老金赤字')
        plt.title('养老金收支与赤字', fontsize=14)
        plt.xlabel('时期', fontsize=12)
        plt.ylabel('金额', fontsize=12)
        plt.legend(fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.3)
        
        # 2. 人口结构比例变化
        plt.subplot(2, 2, 2)
        plt.stackplot(t, young_ratio, older_ratio, retirees_ratio,
                     labels=['年轻工作者(22-41岁)', '老年工作者(42-61岁)', '退休者(62-101岁)'],
                     colors=['#3498db', '#e67e22', '#2ecc71'])
        plt.title('人口结构比例变化 (%)', fontsize=14)
        plt.ylabel('人口比例 (%)', fontsize=12)
        plt.xlabel('时期', fontsize=12)
        plt.legend(loc='upper right', fontsize=10)
        plt.ylim(0, 100)
        plt.grid(True, linestyle='--', alpha=0.3)
        
        # 3. 退休人口比例和抚养比
        plt.subplot(2, 2, 3)
        fig.ax1 = plt.gca()
        fig.ax2 = fig.ax1.twinx()
        
        fig.ax1.plot(t, retirees_ratio, 'g-', linewidth=2, label='退休人口比例')
        fig.ax1.set_ylabel('退休人口比例 (%)', color='g', fontsize=12)
        fig.ax1.tick_params(axis='y', labelcolor='g')
        
        fig.ax2.plot(t, self.dependency_ratio_history, 'r-', linewidth=2, label='抚养比')
        fig.ax2.set_ylabel('抚养比', color='r', fontsize=12)
        fig.ax2.tick_params(axis='y', labelcolor='r')
        
        plt.title('退休人口比例和抚养比', fontsize=14)
        plt.xlabel('时期', fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.3)
        
        # 添加两个图例
        lines1, labels1 = fig.ax1.get_legend_handles_labels()
        lines2, labels2 = fig.ax2.get_legend_handles_labels()
        plt.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=10)
        
        # 4. 养老金替代率和累计赤字占GDP比
        plt.subplot(2, 2, 4)
        
        # 计算累计赤字
        cumulative_deficit = np.zeros(self.periods)
        for i in range(self.periods):
            if i == 0:
                cumulative_deficit[i] = self.deficit_history[i]
            else:
                cumulative_deficit[i] = cumulative_deficit[i-1] * 1.02 + self.deficit_history[i]
        
        # 计算累计赤字占GDP比
        deficit_to_gdp = np.array(cumulative_deficit) / np.array(self.Y_history) * 100
        
        fig.ax3 = plt.gca()
        fig.ax4 = fig.ax3.twinx()
        
        fig.ax3.plot(t, self.replacement_rate_history, 'b-', linewidth=2, label='养老金替代率')
        fig.ax3.set_ylabel('替代率', color='b', fontsize=12)
        fig.ax3.tick_params(axis='y', labelcolor='b')
        
        fig.ax4.plot(t, deficit_to_gdp, 'm-', linewidth=2, label='累计赤字占GDP比')
        fig.ax4.set_ylabel('赤字占GDP比 (%)', color='m', fontsize=12)
        fig.ax4.tick_params(axis='y', labelcolor='m')
        
        plt.title('养老金替代率和累计赤字占GDP比', fontsize=14)
        plt.xlabel('时期', fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.3)
        
        # 添加两个图例
        lines3, labels3 = fig.ax3.get_legend_handles_labels()
        lines4, labels4 = fig.ax4.get_legend_handles_labels()
        plt.legend(lines3 + lines4, labels3 + labels4, loc='upper left', fontsize=10)
        
        plt.tight_layout()
        plt.savefig('fig/olg_results_with_deficit.png')
        
        # 单独创建养老金赤字图
        plt.figure(figsize=(12, 6))
        plt.plot(t, self.deficit_history, 'r-', linewidth=2.5, label='养老金赤字')
        plt.fill_between(t, 0, self.deficit_history, color='red', alpha=0.3)
        
        # 添加赤字起始点标记
        deficit_start = next((i for i, x in enumerate(self.deficit_history) if x > 0), None)
        if deficit_start is not None:
            plt.axvline(x=deficit_start, color='black', linestyle='--')
            plt.text(deficit_start+0.5, max(self.deficit_history)*0.8, 
                  f'赤字开始于第{deficit_start}期', fontsize=12)
        
        plt.title('养老金赤字演变', fontsize=16)
        plt.xlabel('时期', fontsize=14)
        plt.ylabel('赤字金额', fontsize=14)
        plt.grid(True, linestyle='--', alpha=0.3)
        plt.legend(fontsize=12)
        plt.savefig('fig/pension_deficit_only.png')
        
        plt.show()

if __name__ == "__main__":
    # 设置随机数种子以保证结果可重复
    np.random.seed(42)
    
    # 初始化并运行模型
    model = FullOLGModel()
    model.run_simulation()
    model.plot_results()
    
    print("模拟完成。结果已保存到 fig/olg_results_with_deficit.png")
