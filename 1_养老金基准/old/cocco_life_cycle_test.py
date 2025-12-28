import numpy as np
from scipy import interpolate
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import os
from scipy.stats import norm
from scipy.interpolate import interp1d

# 基础参数设置
tb = 20  # 开始工作年龄
tr = 61  # 退休年龄
td = 100  # 死亡年龄 
tn = td - tb + 1  # 总期数

# 资产配置和消费的网格设置
ncash = 51  # 现金网格点数
n = 5  # 收入冲击的网格点数

# 收入函数系数
aa = (-2.170042 + 2.700381)
b1 = 0.16818
b2 = -0.0323371/10
b3 = 0.0019704/100

# 其他参数
ret_fac = 0.6827  # 退休后固定支付的工资
pension_pct = 0  # 工资扣除用于养老保险的比例

# 标准差参数
smay = np.sqrt(0.01)  # 暂时性shock的标准差
smav = np.sqrt(0.01)  # 永久性shock的标准差
corr_z_epsilon = 0.0  # 永久性收入增长率与股票回报率的相关系数
corr_u_epsilon = 0.0  # 暂时性AR(1)收入与股票回报率的相关系数

# 效用函数参数
gamma = 3.84  # Epstein-Zin的风险厌恶系数
beta = 0.95  # Epstein-Zin的贴现因子
psi = 0.15  # Epstein-Zin的跨期替代弹性

# 资产回报参数
rf = 1.02  # 无风险利率
mu = 0.04  # 风险溢价
sigr = 0.27  # 风险资产收益率的标准差

# 初始化数组
survprob = np.zeros(tn-1)  # 生存概率
grid = np.zeros(n)  # 收入冲击网格
weig = np.zeros(n)  # 权重
gret = np.zeros(n)  # 回报率网格
yp = np.zeros((n,n))  # 永久性收入冲击
yh = np.zeros((n,n))  # 暂时性收入冲击

# 生存概率数据(这里只展示部分数据作为示例)
survprob = [
0.99975, 0.99974, 0.99973, 0.99972, 0.99971, 0.99969, 0.99968, 0.99966,
0.99963, 0.99961, 0.99958, 0.99954, 0.99951, 0.99947, 0.99943, 0.99938,
0.99934, 0.99928, 0.99922, 0.99916, 0.99908, 0.999, 0.99891, 0.99881,
0.9987, 0.99857, 0.99843, 0.99828, 0.99812, 0.99794, 0.99776, 0.99756,
0.99735, 0.99713, 0.9969, 0.99666, 0.9964, 0.99613, 0.99583, 0.99551,
0.99515, 0.99476, 0.99432, 0.99383, 0.9933, 0.9927, 0.99205, 0.99133,
0.99053, 0.98961, 0.98852, 0.98718, 0.98553, 0.98346, 0.98089, 0.97772,
0.97391, 0.96943, 0.96429, 0.95854, 0.95221, 0.94537, 0.93805, 0.93027,
0.92202, 0.91327, 0.90393, 0.89389, 0.88304, 0.87126, 0.85846, 0.84452,
0.82935, 0.81282, 0.79485, 0.77543, 0.75458, 0.7324, 0.70893, 0.68424]


def discretize_AR1_Tauchen(mu, rho, sigma, n, m=2):
    """
    Tauchen方法离散化AR(1)过程
    参数:
    mu: 均值
    rho: 自回归系数
    sigma: 标准差
    n: 网格点数
    m: 标准差的倍数范围
    """
    sigma_y = sigma / np.sqrt(1 - rho**2)
    y_max = m * sigma_y
    y_min = -y_max
    step = (y_max - y_min) / (n - 1)
    
    # 生成网格点
    grid = np.linspace(y_min, y_max, n)
    
    # 计算转移概率

    
    trans_mat = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if j == 0:
                trans_mat[i,j] = norm.cdf((grid[j] - mu - rho*grid[i] + step/2) / sigma)
            elif j == n-1:
                trans_mat[i,j] = 1 - norm.cdf((grid[j] - mu - rho*grid[i] - step/2) / sigma)
            else:
                trans_mat[i,j] = (norm.cdf((grid[j] - mu - rho*grid[i] + step/2) / sigma) - 
                                norm.cdf((grid[j] - mu - rho*grid[i] - step/2) / sigma))
    
    return grid, trans_mat

# 生成收入过程
# def generate_income_process():
    # 工作期收入增长率
f_y = np.zeros(tr-tb+1)
for i in range(tb, tr+1):
    f_y[i-tb] = np.exp(aa + b1*i + b2*i**2 + b3*i**3)

# 计算增长率
gy = np.zeros(tr-tb)
for i in range(tb, tr):
    gy[i-tb] = f_y[i-tb+1]/f_y[i-tb] - 1.0

# 生成永久性和暂时性收入冲击
gamma00 = 0  # AR自相关系数
mew = 0 # AR的常数项
sigma = 1 # 白噪声的标准差
grid, weig = discretize_AR1_Tauchen(mew, gamma00, sigma, n)
weig = np.diag(weig)

# 生成收入冲击矩阵
yh = np.zeros((n, n))
yp = np.zeros((n, n))
for i in range(n):
    grid2 = grid[i]*corr_u_epsilon + grid*np.sqrt(1-corr_u_epsilon**2)
    yh[:,i] = np.exp(grid2 * smay)
    
    grid2 = grid[i]*corr_z_epsilon + grid*np.sqrt(1-corr_z_epsilon**2)
    yp[:,i] = grid2 * smav

# 生成归一化后的收入增长率
gyp = np.zeros((n, n, tn-1))

# 工作期
for t in range(tr-tb):
    for i in range(n):
        gyp[:,i,t] = np.exp(gy[t] + yp[:,i])

# 退休期
for t in range(tr-tb, tn-1):
    gyp[:,:,t] = 1.0
    
    
# 现金网格生成

maxcash = 100
mincash = 0.25
l_maxcash = np.log(maxcash)
l_mincash = np.log(mincash)
stepcash = (l_maxcash - l_mincash)/(ncash-1)

lgcash = np.array([l_mincash + i*stepcash for i in range(ncash)])
gcash = np.exp(lgcash)
    


# 生成风险资产回报率
gret = np.zeros(n)
for i in range(n):
    gret[i] = rf + mu + grid[i]*sigr

# 生成权重矩阵

nweig1 = np.zeros((n,n,n))
for i in range(n):
    for j in range(n):
        for k in range(n):
            nweig1[i,j,k] = weig[i]*weig[j]*weig[k]



def retired_value_function(policy, cash, next_period_value, gret, rf, ret_fac,
                         gamma, beta, gcash, survprob, weig):
    """
    退休期的值函数
    参考 cocco_fun_valuefunc_retired.m 实现
    """
    c, alpha = policy
    if c <= 0 or c >= cash or alpha < 0 or alpha > 1:
        return 1e10
        
    # 计算储蓄
    # sav = cash - c
    sav = cash * (1-c)
    
    # 计算效用
    # u = (c**(1-gamma))/(1-gamma)
    u = ((cash*c)**(1-gamma))/(1-gamma)
    
    # 计算下期现金
    cash_next = (rf*(1-alpha) + gret*alpha)*sav + ret_fac
    
    # 插值计算下期值函数
    v_next = np.interp(cash_next, gcash, next_period_value)
    
    # 计算期望值
    auxVV = np.sum(weig * survprob * v_next)
    
    if auxVV == 0:
        v = -(u)
    else:
        v = -(u + beta * auxVV)
        
    return v

def working_value_function(policy, cash, gyp, next_period_value, yh, gret, rf,
                         gamma, beta, gcash, survprob, n, nweig1):
    """
    工作期的值函数
    参考 cocco_fun_valuefunc_work.m 实现
    """
    c, alpha = policy
    if c <= 0 or c >= cash or alpha < 0 or alpha > 1:
        return 1e10
    
    # 计算效用
    u = ((cash*c)**(1-gamma))/(1-gamma)

    # 计算期望值
    auxVV = 0
    for i in range(n): # z_t, 来自P_t=P_t-1+z_t
        for j in range(n): # epsilon_t, 即风险收益随机项
            for k in range(n): # u_t, 来自logY_t = f_t + P_t + u_t
                # 计算下期的cash-on-hand
                # sav = (cash - c)/gyp[i,j]
                sav = cash * (1-c)/gyp[i,j]
                cash_next = (rf*(1-alpha) + gret[j]*alpha)*sav + yh[k,j]
                cash_next = np.clip(cash_next, gcash[0], gcash[-1])               
                # 插值计算下期值函数
                v_next = np.interp(cash_next, gcash, next_period_value)
                
                # 累加计算期望值
                auxVV += nweig1[i,j,k] * survprob * (v_next * gyp[i,j])
    
    if auxVV == 0:
        v = -(u)
    else:
        v = -(u + beta * auxVV)
        
    return v

def solve_lifecycle_problem(gcash, V, C, A, tr, tb, tn, survprob):
    """
    求解生命周期问题
    """
    # 终末期处理
    # C[:, tn-1] = gcash  # 最后一期的消费等于现金
    C[:, tn-1] = 1 # 最后一期消费比例为1
    A[:, tn-1] = 0.0    # 最后一期不进行投资
    # V[:, tn-1] = (C[:, tn-1]**(1-gamma))/(1-gamma)  # 最后一期的值函数
    V[:, tn-1] = ((gcash**(1-gamma))/(1-gamma))  # 最后一期的值函数
    # 退休期求解
    print("Solving retirement periods...")
    for t in range(tn-2, tr-tb-1, -1):
        print(f"Period {t}")
        for i_cash in range(len(gcash)):
            cash = gcash[i_cash]
            
            # 初始猜测值
            x0 = [0.2, 0.2]
            
            # 约束条件
            # if i_cash == 0:
            #     bounds = [(0, 0.999*cash), (0, 1)]
            # else:
            #     bounds = [(C[i_cash-1, t], 0.999*cash), (0, 1)]
            if i_cash == 0:
                bounds = [(0,1), (0,1)]
            else:
                bounds = [(C[i_cash-1, t]*gcash[i_cash-1]/gcash[i_cash], 1), (0, 1)]
            
            # 求解最优化问题
            # 调用MATLAB的fmincon函数
            result = minimize(retired_value_function, x0,
                args=(cash, V[:, t+1], gret, rf, ret_fac, gamma, beta,
                        gcash, survprob[t], weig),
                method='SLSQP',
                bounds=bounds,
                options={'disp': False})
                
            # import matlab.engine
            # eng = matlab.engine.start_matlab()
            
            # # 将Python参数转换为MATLAB格式
            # x0_mat = matlab.double(x0)
            # cash_mat = matlab.double([cash])
            # V_next_mat = matlab.double(V[:, t+1].tolist())
            # gret_mat = matlab.double(gret.tolist())
            # rf_mat = matlab.double([rf])
            # ret_fac_mat = matlab.double([ret_fac])
            # gamma_mat = matlab.double([gamma])
            # beta_mat = matlab.double([beta])
            # gcash_mat = matlab.double(gcash.tolist())
            # survprob_mat = matlab.double([survprob[t]])
            # weig_mat = matlab.double(weig.tolist())
            
            # # 设置fmincon的选项
            # options = {'Display': 'off'}
            # lb = matlab.double([bounds[0][0], bounds[1][0]])
            # ub = matlab.double([bounds[0][1], bounds[1][1]])
            
            # # 调用fmincon
            # result = eng.fmincon('retired_value_function_mat', x0_mat, [], [], [], [], 
            #                    lb, ub, [], options,
            #                    cash_mat, V_next_mat, gret_mat, rf_mat, ret_fac_mat,
            #                    gamma_mat, beta_mat, gcash_mat, survprob_mat, weig_mat)
            
            # # 将结果转换回Python格式
            # result = {'x': np.array(result).flatten(), 'success': True}
            # eng.quit()
 
            
            if result.success:
                C[i_cash, t] = result.x[0]
                A[i_cash, t] = result.x[1]
                V[i_cash, t] = -result.fun
            else:
                print(f"Optimization failed at t={t}, cash={cash}")
    
    # 工作期求解
    print("Solving working periods...")
    
    for t in range(tr-tb-1, -1, -1):
        print(f"Period {t}")
        for i_cash in range(len(gcash)):
            cash = gcash[i_cash]
            
            # 初始猜测值
            # x0 = [0.1*cash, 0.1]
            
            # 约束条件
            # if i_cash == 0:
            #     bounds = [(0, 0.999*cash), (0, 1)]
            # else:
            #     bounds = [(C[i_cash-1, t], 0.999*cash), (0, 1)]
            x0 = [0.2, 0.2]
            if i_cash == 0:
                bounds = [(0,1), (0,1)]
            else:
                bounds = [(C[i_cash-1, t]*gcash[i_cash-1]/gcash[i_cash], 1), (0, 1)]    

                     
            # 求解最优化问题
            result = minimize(working_value_function, x0,
                            args=(cash, gyp[:,:,t], V[:, t+1], yh, gret, rf, gamma, 
                                 beta, gcash, survprob[t], n, nweig1),
                            method='SLSQP',
                            bounds=bounds)

            if result.success:
                C[i_cash, t] = result.x[0]
                A[i_cash, t] = result.x[1]
                V[i_cash, t] = -result.fun
            else:
                print(f"Optimization failed at t={t}, cash={cash}")
    
    return V, C, A

def run_lifecycle_model():
    """
    运行生命周期模型
    """
    # 初始化数组
    V = np.zeros((ncash, tn))
    C = np.zeros((ncash, tn))
    A = np.ones((ncash, tn))
    
    # 求解生命周期问题
    V, C, A = solve_lifecycle_problem(gcash, V, C, A, tr, tb, tn, survprob)
    
    return V, C, A

def simulate_lifecycle(C, A, gcash, nsim=10000):
    """
    生命周期模型的模拟
    参数:
    C: 最优消费策略
    A: 最优投资组合策略
    gcash: 现金网格
    nsim: 模拟路径数量
    """
    # 初始化模拟数组
    simPY = np.zeros((tn, nsim))  # 永久收入
    simGPY = np.zeros((tn, nsim))  # 收入增长
    simY = np.zeros((tn, nsim))    # 收入
    simC = np.zeros((tn, nsim))    # 消费
    simW = np.zeros((tn, nsim))    # 财富
    simA = np.zeros((tn, nsim))    # 投资组合份额
    simS = np.zeros((tn, nsim))    # 风险资产投资
    simB = np.zeros((tn, nsim))    # 无风险资产投资
    simW_Y = np.zeros((tn, nsim))  # 财富收入比
    simR = np.zeros((tn, nsim))    # 投资回报率

    # 初始化财富
    simW[:,:] = 0.2

    # 1. 模拟收入过程
    for i in range(nsim//2):
        # 工作期第一期
        eps_y = np.random.normal(0, 1)
        simPY[0, i] = eps_y * smav
        simPY[0, nsim//2 + i] = -eps_y * smav
        simGPY[0, i] = 1.0
        simGPY[0, nsim//2 + i] = 1.0
        
        simTY = np.random.normal(0, 1)
        simY[0, i] = np.exp(simTY * smay)
        simY[0, nsim//2 + i] = np.exp(-simTY * smay)

        # 工作期后续各期
        for t in range(1, tr-tb):
            eps_y = np.random.normal(0, 1)
            simPY[t, i] = eps_y * smav + simPY[t-1, i]
            simPY[t, nsim//2 + i] = -eps_y * smav + simPY[t-1, nsim//2 + i]
            
            simGPY[t, i] = np.exp(gy[t-1]) * np.exp(simPY[t, i]) / np.exp(simPY[t-1, i])
            simGPY[t, nsim//2 + i] = np.exp(gy[t-1]) * np.exp(simPY[t, nsim//2 + i]) / np.exp(simPY[t-1, nsim//2 + i])
            
            simTY = np.random.normal(0, 1)
            simY[t, i] = np.exp(simTY * smay)
            simY[t, nsim//2 + i] = np.exp(-simTY * smay)

    # 退休期收入
    simY[tr-tb:, :] = ret_fac
    simGPY[tr-tb:, :] = 1.0

    # 2. 模拟投资回报率
    for t in range(tn):
        for i in range(nsim//2):
            eps_r = np.random.normal(0, 1)
            simR[t, i] = mu + rf + sigr * eps_r
            simR[t, nsim//2 + i] = mu + rf - sigr * eps_r

    # 3. 模拟最优决策
    for t in range(tn):
        for i in range(nsim):
            # 计算当期状态
            simW_Y[t, i] = simW[t, i] / simY[t, i]
            cash = simW[t, i] + simY[t, i]

            # 插值得到最优决策
            simC[t, i] = np.interp(cash, gcash, C[:, t])
            simA[t, i] = np.interp(cash, gcash, A[:, t])

            # 确保决策在合理范围内
            simC[t, i] = np.clip(simC[t, i], 0, 1)
            simC[t, i] = simC[t, i] * cash # 实际消费
            simA[t, i] = np.clip(simA[t, i], 0, 1)

            # 计算储蓄和投资
            sav = cash - simC[t, i]
            simS[t, i] = simA[t, i] * sav
            simS[t, i] = np.minimum(simS[t, i], sav)
            simB[t, i] = sav - simS[t, i]

            # 更新下期财富
            if t < tn-1:
                simW[t+1, i] = (simB[t, i] * rf + simS[t, i] * simR[t, i]) / simGPY[t+1, i]

    return simY, simC, simW, simA, simS, simB, simW_Y

def calculate_statistics(simY, simC, simW, simA, simS, simB, simW_Y):
    """
    计算模拟结果的统计量
    """
    meanC = np.mean(simC, axis=1)
    meanY = np.mean(simY, axis=1)
    meanW = np.mean(simW, axis=1)
    meanS = np.mean(simS, axis=1)
    meanB = np.mean(simB, axis=1)
    meanWY = np.mean(simW_Y, axis=1)
    meanalpha = np.mean(simA, axis=1)

    return {
        'meanC': meanC,
        'meanY': meanY,
        'meanW': meanW,
        'meanS': meanS,
        'meanB': meanB,
        'meanWY': meanWY,
        'meanalpha': meanalpha,
    }

def plot_results(stats):
    """
    绘制模拟结果
    """
    # 创建图形
    # 设置中文字体
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
    plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号
    
    fig, ax = plt.subplots(figsize=(6, 5))
    
    # 绘制消费和投资组合份额
    ax.plot(stats['meanC'], label='消费')
    ax.plot(stats['meanalpha'], label='投资组合份额')
    ax.set_title('平均消费和投资组合份额')
    ax.legend()

    # 保存图片
    current_dir = os.path.dirname(os.path.abspath(__file__))
    plt.savefig(os.path.join(current_dir, 'result/lifecycle_plot.png'))
    
    plt.tight_layout()
    plt.show()

def save_results(stats, filepath):
    """
    保存模拟结果
    """
    np.savez(filepath,
             meanC=stats['meanC'],
             meanY=stats['meanY'],
             meanW=stats['meanW'],
             meanS=stats['meanS'],
             meanB=stats['meanB'],
             meanWY=stats['meanWY'],
             meanalpha=stats['meanalpha'])

# 运行模拟
def run_simulation(C, A, gcash, nsim=10000):
    """
    运行完整的模拟过程
    """
    # 执行模拟
    print("开始模拟...")
    sim_results = simulate_lifecycle(C, A, gcash, nsim)
    
    # 计算统计量
    print("计算统计量...")
    stats = calculate_statistics(*sim_results)
    
    # 绘制结果
    print("绘制结果...")
    plot_results(stats)
    
    # 保存结果
    print("保存结果...")
    current_dir = os.path.dirname(os.path.abspath(__file__))
    save_results(stats, os.path.join(current_dir, 'result/baseline_nopension.npz'))
    
    return sim_results, stats

# 使用示例
if __name__ == "__main__":
# 运行模型
    V, C, A = run_lifecycle_model()
    # 保存V,C,A为xlsx
    import pandas as pd
    
    # 创建DataFrame
    df_V = pd.DataFrame(V)
    df_C = pd.DataFrame(C) 
    df_A = pd.DataFrame(A)
    
    # 获取当前文件所在目录
    current_dir = os.path.dirname(os.path.abspath(__file__))
    
    # 保存到xlsx文件
    with pd.ExcelWriter(os.path.join(current_dir, 'result/VCA_results_2.xlsx')) as writer:
        df_V.to_excel(writer, sheet_name='Value_Function')
        df_C.to_excel(writer, sheet_name='Consumption')
        df_A.to_excel(writer, sheet_name='Portfolio_Share')
    
    print("V,C,A已保存到xlsx文件")



    # 假设已经求解得到最优策略C和A
    sim_results, stats = run_simulation(C, A, gcash)

