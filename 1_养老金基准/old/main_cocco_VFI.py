import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import minimize
from scipy.interpolate import interp1d
import pandas as pd
import os
from tqdm import tqdm
import argparse

# 辅助函数
def f_randn(size):
    """生成标准正态分布随机数"""
    return np.random.randn(size)

def f_ntoil_modified(value, grid, n):
    """寻找值在网格中最接近的位置"""
    if value <= grid[0]:
        return 0
    if value >= grid[n-1]:
        return n-1
    
    for i in range(n-1):
        if value >= grid[i] and value < grid[i+1]:
            return i
    
    return n-1

# AR(1)过程离散化的Tauchen方法
def discretizeAR1_Tauchen(mew, rho, sigma, znum, Tauchen_q, tauchenoptions=None):
    """
    离散化AR(1)过程为马尔可夫链
    
    参数:
    mew - 常数项
    rho - 自相关系数
    sigma - 创新的标准差
    znum - 离散化后的状态数（必须是奇数）
    Tauchen_q - 定义网格点的超参数
    tauchenoptions - 额外选项
    
    返回:
    z_grid - 离散的状态值
    P - 转移矩阵
    """
    if tauchenoptions is None:
        tauchenoptions = {'parallel': 0}
    
    if znum == 1:
        z_grid = np.array([mew / (1 - rho)])
        P = np.array([[1]])
        return z_grid, P
    
    # 标准方法
    zstar = mew / (1 - rho)  # z的期望值
    sigmaz = sigma / np.sqrt(1 - rho**2)  # z的标准差
    
    z_grid = zstar * np.ones(znum) + np.linspace(-Tauchen_q * sigmaz, Tauchen_q * sigmaz, znum)
    omega = z_grid[1] - z_grid[0]  # 网格点等距
    
    P = np.zeros((znum, znum))
    
    for i in range(znum):
        for j in range(znum):
            if j == 0:
                P[i, j] = norm.cdf((z_grid[j] + omega/2 - rho * z_grid[i] - mew) / sigma)
            elif j == znum - 1:
                P[i, j] = 1 - norm.cdf((z_grid[j] - omega/2 - rho * z_grid[i] - mew) / sigma)
            else:
                P[i, j] = norm.cdf((z_grid[j] + omega/2 - rho * z_grid[i] - mew) / sigma) - \
                         norm.cdf((z_grid[j] - omega/2 - rho * z_grid[i] - mew) / sigma)
    
    return z_grid, P

# 退休期间的值函数计算
def cocco_fun_valuefunc_retired(x, cash, nextV, gret, rf, ret_fac, gamma, beta, psi_1, psi_2, theta, gcash, survprob, weig):
    """
    退休期间的值函数计算
    
    参数:
    x - 决策变量 [消费, 风险资产投资比例]
    cash - 手中现金
    nextV - 下期值函数
    gret - 随机收益率
    rf - 无风险收益率
    ret_fac - 固定养老金
    gamma - 相对风险规避系数
    beta - 贴现因子
    psi_1, psi_2, theta - Epstein-Zin效用函数参数
    gcash - 现金网格
    survprob - 存活概率
    weig - 随机状态权重
    
    返回:
    v - 函数值
    """
    auxVV = 0
    sav = cash - x[0]
    
    u = -(x[0])**(1-gamma)
    
    # 下期的cash-on-hand
    cash_1 = (rf * (1 - x[1]) + gret * x[1]) * sav + ret_fac
    
    # 插值
    int_V = interp1d(gcash, nextV, kind='cubic', fill_value='extrapolate')(cash_1)
    
    # 计算期望值
    auxVV = auxVV + np.dot(weig, survprob * int_V)
    
    if auxVV == 0:
        v = -(u)
    else:
        v = -(u + beta * auxVV)
    
    return v

# 工作期间的值函数计算
def cocco_fun_valuefunc_work(x, cash, gyp, nextV, yh, gret, rf, gamma, beta, psi_1, psi_2, theta, gcash, survprob, n, weig):
    """
    工作期间的值函数计算
    
    参数:
    x - 决策变量 [消费, 风险资产投资比例]
    cash - 手中现金
    gyp - 劳动收入增长率
    nextV - 下期值函数
    yh - 暂时收入冲击
    gret - 随机收益率
    rf - 无风险收益率
    gamma - 相对风险规避系数
    beta - 贴现因子
    psi_1, psi_2, theta - Epstein-Zin效用函数参数
    gcash - 现金网格
    survprob - 存活概率
    n - 随机状态数量
    weig - 随机状态权重
    
    返回:
    v - 函数值
    """
    auxVV = 0
    
    u = -(x[0])**(1-gamma)
    
    # 创建索引网格
    i6_grid, i7_grid, i8_grid = np.meshgrid(np.arange(n), np.arange(n), np.arange(n), indexing='ij')
    
    # 计算所有可能的储蓄值
    sav_values = (cash - x[0]) / gyp[i6_grid, i8_grid]
    
    # 计算所有可能的现金值
    portfolio_return = rf * (1 - x[1]) + gret[i8_grid] * x[1]
    cash_1_values = portfolio_return * sav_values + yh[i7_grid, i8_grid]
    
    # 限制cash_1的范围
    cash_1_values = np.clip(cash_1_values, gcash[0], gcash[-1])
    
    # 创建插值函数
    interp_func = interp1d(gcash, nextV, kind='cubic', fill_value='extrapolate')
    
    # 对所有cash_1值应用插值
    int_V_values = interp_func(cash_1_values.flatten()).reshape(cash_1_values.shape)
    
    # 计算期望值
    auxVV = np.sum(weig * survprob * int_V_values * gyp[i6_grid, i8_grid])
    # 原循环
    # for i6 in range(n):  # z_t, 来自P_t=P_t-1+z_t
    # for i8 in range(n):  # epsilon_t, 即风险收益随机项
    #     for i7 in range(n):  # u_t, 来自logY_t = f_t + P_t + u_t
    #         # 下期的cash-on-hand
    #         sav = (cash - x[0]) / gyp[i6, i8]
    #         cash_1 = (rf * (1 - x[1]) + gret[i8] * x[1]) * sav + yh[i7, i8]
    #         cash_1 = max(min(cash_1, gcash[-1]), gcash[0])
            
    #         # 插值
    #         int_V = interp1d(gcash, nextV, kind='cubic', fill_value='extrapolate')(cash_1)
            
    #         # 计算期望值
    #         auxVV = auxVV + weig[i6, i7, i8] * survprob * (int_V * gyp[i6, i8])

    if auxVV == 0:
        v = -(u)
    else:
        v = -(u + beta * auxVV)
    
    return v

def solve_model():
    """
    求解模型的最优消费和投资决策
    """
    # 设置随机种子以获得可重复的结果
    np.random.seed(42)
    
    # 变量定义
    tb = 18    # 初始开始的年纪
    tr = 61    # 退休年龄
    td = 100   # 死亡年龄
    tn = td - tb + 1  # 总期数
    
    # 状态变量grid数量
    ncash = 51  # 手中现金
    n = 5       # 外生随机冲击的grid数量
    
    # 外生参数
    # 基础收入f的系数
    aa = (-2.170042 + 2.700381)
    b1 = 0.16818
    b2 = -0.0323371 / 10
    b3 = 0.0019704 / 100
    
    # 养老金相关
    ret_fac = 0.6827  # 退休后固定支付的工资（养老金）
    pension_pct = 0   # 工资扣除缴纳养老保险的比例
    
    # 随机冲击参数
    smay = np.sqrt(0.169993)  # 白噪声shock的标准差
    smav = np.sqrt(0.112572)  # 持续shock的标准差
    corr_z_epsilon = 0.0      # 工资收入白噪声与风险收益随机项的相关性
    corr_u_epsilon = 0.0      # 工资收入AR(1)随机项与风险收益随机项的相关性
    
    # 效用函数参数
    gamma = 3.84   # Epstein-Zin的相对风险规避系数
    beta = 0.95    # Epstein-Zin的discount factor
    psi = 0.15     # Epstein-Zin的跨期替代弹性
    
    # 资产收益率参数
    rf = 1.02    # 无风险总收入
    mu = 0.04    # 超额收益
    sigr = 0.27  # 风险资产收益率的标准差
    
    # 初始化数组
    survprob = np.zeros(tn - 1)
    beta2 = np.zeros(tn - 1)
    grid = np.zeros(n)
    weig = np.zeros(n)
    gret = np.zeros(n)
    ones_n_1 = np.ones(n)
    grid2 = np.zeros(n)
    yp = np.zeros((n, n))
    yh = np.zeros((n, n))
    nweig1 = np.zeros((n, n, n))
    f_y = np.zeros(tr - tb + 1)
    gy = np.zeros(tr - tb)
    gyp = np.zeros((n, n, tn - 1))
    gcash = np.zeros(ncash)
    lgcash = np.zeros(ncash)
    secd = np.zeros(ncash)
    C = np.zeros((ncash, tn))
    V = np.zeros((ncash, tn))
    A = np.ones((ncash, tn))
    
    # 正态分布的离散化近似
    gamma00 = 0  # AR自相关系数
    mew = 0      # AR的常数项
    sigma = 1    # 白噪声的标准差
    tauchenoptions = {'parallel': 0}
    grid, weig_matrix = discretizeAR1_Tauchen(mew, gamma00, sigma, n, 2, tauchenoptions)
    weig = np.diag(weig_matrix)
    
    # 条件生存概率 (从MATLAB代码中复制)
    survprob_data = [0.99975, 0.99974, 0.99973, 0.99972, 0.99971, 0.99969, 0.99968, 0.99966, 0.99963, 0.99961,
                     0.99958, 0.99954, 0.99951, 0.99947, 0.99943, 0.99938, 0.99934, 0.99928, 0.99922, 0.99916,
                     0.99908, 0.999, 0.99891, 0.99881, 0.9987, 0.99857, 0.99843, 0.99828, 0.99812, 0.99794,
                     0.99776, 0.99756, 0.99735, 0.99713, 0.9969, 0.99666, 0.9964, 0.99613, 0.99583, 0.99551,
                     0.99515, 0.99476, 0.99432, 0.99383, 0.9933, 0.9927, 0.99205, 0.99133, 0.99053, 0.98961,
                     0.98852, 0.98718, 0.98553, 0.98346, 0.98089, 0.97772, 0.97391, 0.96943, 0.96429, 0.95854,
                     0.95221, 0.94537, 0.93805, 0.93027, 0.92202, 0.91327, 0.90393, 0.89389, 0.88304, 0.87126,
                     0.85846, 0.84452, 0.82935, 0.81282, 0.79485, 0.77543, 0.75458, 0.7324, 0.70893, 0.68424,
                     0.68424, 0.68424]
    
    for i in range(min(len(survprob_data), len(survprob))):
        survprob[i] = survprob_data[i]
    
    # 风险资产收益率
    for i1 in range(n):
        gret[i1] = rf + mu + grid[i1] * sigr
    
    # 外生随机性的概率
    for i6 in range(n):
        for i7 in range(n):
            for i8 in range(n):
                nweig1[i6, i7, i8] = weig[i6] * weig[i7] * weig[i8]
    
    # Epstein-Zin效用函数参数
    theta = (1.0 - gamma) / (1.0 - 1.0 / psi)
    psi_1 = 1.0 - 1.0 / psi
    psi_2 = 1.0 / psi_1
    
    # cash-on-hand的grid
    maxcash = 100
    mincash = 0.25
    l_maxcash = np.log(maxcash)
    l_mincash = np.log(mincash)
    stepcash = (l_maxcash - l_mincash) / (ncash - 1)
    
    for i1 in range(ncash):
        lgcash[i1] = l_mincash + i1 * stepcash
    
    for i1 in range(ncash):
        gcash[i1] = np.exp(lgcash[i1])
    
    # 劳动收入（归一化）
    for i1 in range(n):
        grid2 = grid[i1] * corr_u_epsilon + grid * np.ones(n) * np.sqrt(1 - corr_u_epsilon**2)
        yh[:, i1] = np.exp(grid2 * smay)  # 白噪声随机项u_t的grid
    
    for i1 in range(n):
        grid2 = grid[i1] * corr_z_epsilon + grid * np.ones(n) * np.sqrt(1 - corr_z_epsilon**2)
        yp[:, i1] = grid2 * smav  # permanent labor p_t 的随机游走项白噪声z_t的grid
    
    # 随年龄变化的基础收入
    for i1 in range(tb, tr + 1):
        f_y[i1 - tb] = np.exp((aa + b1 * i1 + b2 * i1**2 + b3 * i1**3))
    
    # 工作年龄
    for i1 in range(tb, tr):
        gy[i1 - tb] = f_y[i1 - tb + 1] / f_y[i1 - tb] - 1.0  # 对基础工资收入进行normalization
        for i2 in range(n):
            gyp[:, i2, i1 - tb] = np.exp(gy[i1 - tb] * np.ones(n) + yp[:, i2])  # normalized工资收入+持久收入shock
    
    # 退休年龄
    for i1 in range(tr - tb, tn - 1):
        for i2 in range(n):
            gyp[:, i2, i1] = np.exp(0.0 * np.ones(n))
    
    # 终止期
    # 策略函数
    C[:, tn - 1] = gcash  # 最后一期的现金全部用来消费
    A[:, tn - 1] = 0.0    # 最后一期的投资全部为零
    delta = 0.2
    
    for i1 in range(ncash):
        V[i1, tn - 1] = -(C[i1, tn - 1])**(1 - gamma)
    
    # 值函数迭代
    # 退休期
    print("开始求解退休期值函数...")
    for i1 in tqdm(range(1,td - tr + 1), desc="退休期求解进度"):
        t = tn - i1 - 1
        # 确保survprob索引不会越界
        # t = min(t, len(survprob) - 1)
        
        for i3 in range(ncash):
            x0 = [0.2, 0.2]
            lb = [0, 0]
            ub = [0.999 * gcash[i3], 1]

            if i3 != 0:
                lb = [C[i3 - 1, t], 0]
            
            # 优化函数 - 使用next_t代替t+1，使用t代替t访问survprob
            result = minimize(
                lambda x: cocco_fun_valuefunc_retired(x, gcash[i3], V[:, t+1], gret, rf, ret_fac, gamma, beta, psi_1, psi_2, theta, gcash, survprob[t], weig),
                x0, bounds=list(zip(lb, ub)), method='L-BFGS-B'
            )

            policy = result.x
            C[i3, t] = policy[0]
            A[i3, t] = policy[1]
            V[i3, t] = -cocco_fun_valuefunc_retired(policy, gcash[i3], V[:, t+1], gret, rf, ret_fac, gamma, beta, psi_1, psi_2, theta, gcash, survprob[t], weig)
    
    # 工作期（退休前）
    print("开始求解工作期值函数...")
    for i1 in tqdm(range(1,tr - tb + 1), desc="工作期求解进度"):
        t = tr - tb - i1
        
        
        # 确保survprob索引不会越界
        
        for i3 in range(ncash):
            x0 = [0.1, 0.1]
            lb = [0, 0]
            ub = [0.999 * gcash[i3], 1]
            
            if i3 != 0:
                lb = [C[i3 - 1, t], 0]
            
            # 优化函数 - 使用next_t代替t+1，使用t代替t访问survprob
            # 确保初始值在边界范围内
            x0[0] = max(min(x0[0], ub[0]), lb[0])
            x0[1] = max(min(x0[1], ub[1]), lb[1])
            result = minimize(
                lambda x: cocco_fun_valuefunc_work(x, gcash[i3], gyp[:, :, t], V[:, t+1], yh, gret, rf, gamma, beta, psi_1, psi_2, theta, gcash, survprob[t], n, nweig1),
                x0, bounds=list(zip(lb, ub)), method='L-BFGS-B'
            )
            
            policy = result.x
            C[i3, t] = policy[0]
            A[i3, t] = policy[1]
            V[i3, t] = -cocco_fun_valuefunc_work(policy, gcash[i3], gyp[:, :, t], V[:, t+1], yh, gret, rf, gamma, beta, psi_1, psi_2, theta, gcash, survprob[t], n, nweig1)
    
    # 保存结果
    os.makedirs('result_cocco_python', exist_ok=True)
    pd.DataFrame(C).to_excel('result_cocco_python/C.xlsx')
    pd.DataFrame(A).to_excel('result_cocco_python/A.xlsx')
    pd.DataFrame(V).to_excel('result_cocco_python/V.xlsx')
    pd.DataFrame(gcash).to_excel('result_cocco_python/gcash.xlsx')

    print("模型求解完成，结果已保存到result_cocco_python文件夹")
    
    # 返回求解结果，供模拟使用
    model_results = {
        'C': C,
        'A': A,
        'V': V,
        'gcash': gcash,
        'params': {
            'tb': tb,
            'tr': tr,
            'td': td,
            'tn': tn,
            'n': n,
            'rf': rf,
            'mu': mu,
            'sigr': sigr,
            'smay': smay,
            'smav': smav,
            'ret_fac': ret_fac,
            'gy': gy
        }
    }
    
    return model_results

def simulate_model(model_results=None):
    """
    基于求解结果进行数值模拟
    
    参数：
        model_results: 模型求解结果。如果为None，会尝试从文件中读取
    """
    # 如果没有提供模型结果，尝试从文件中读取
    if model_results is None:
        try:
            # 检查文件是否存在
            if not os.path.exists('result_cocco_python/C.xlsx'):
                print("未找到模型求解结果文件，请先运行模型求解")
                return
                
            # 从文件中读取模型结果
            C = pd.read_excel('result_cocco_python/C.xlsx', index_col=0).values
            A = pd.read_excel('result_cocco_python/A.xlsx', index_col=0).values
            gcash = pd.read_excel('result_cocco_python/gcash.xlsx', index_col=0).values.flatten()
            
            # 设置默认参数
            tb = 18    # 初始开始的年纪
            tr = 61    # 退休年龄
            td = 100   # 死亡年龄
            tn = td - tb + 1  # 总期数
            n = 5       # 外生随机冲击的grid数量
            rf = 1.02    # 无风险总收入
            mu = 0.04    # 超额收益
            sigr = 0.27  # 风险资产收益率的标准差
            smay = np.sqrt(0.169993)  # 白噪声shock的标准差
            smav = np.sqrt(0.112572)  # 持续shock的标准差
            ret_fac = 0.6827  # 退休后固定支付的工资（养老金）
            
            # 计算gy
            aa = (-2.170042 + 2.700381)
            b1 = 0.16818
            b2 = -0.0323371 / 10
            b3 = 0.0019704 / 100
            f_y = np.zeros(tr - tb + 1)
            
            for i1 in range(tb, tr + 1):
                f_y[i1 - tb] = np.exp((aa + b1 * i1 + b2 * i1**2 + b3 * i1**3))
            
            gy = np.zeros(tr - tb)
            for i1 in range(tb, tr):
                gy[i1 - tb] = f_y[i1 - tb + 1] / f_y[i1 - tb] - 1.0
            
            model_results = {
                'C': C,
                'A': A,
                'gcash': gcash,
                'params': {
                    'tb': tb,
                    'tr': tr,
                    'td': td,
                    'tn': tn,
                    'n': n,
                    'rf': rf,
                    'mu': mu,
                    'sigr': sigr,
                    'smay': smay,
                    'smav': smav,
                    'ret_fac': ret_fac,
                    'gy': gy
                }
            }
            
        except Exception as e:
            print(f"读取模型求解结果出错: {e}")
            print("请先运行模型求解")
            return
    
    # 从模型结果中获取参数和数据
    C = model_results['C']
    A = model_results['A']
    gcash = model_results['gcash']
    
    # 获取参数
    params = model_results['params']
    tb = params['tb']
    tr = params['tr']
    td = params['td']
    tn = params['tn']
    n = params['n']
    rf = params['rf']
    mu = params['mu']
    sigr = params['sigr']
    smay = params['smay']
    smav = params['smav']
    ret_fac = params['ret_fac']
    gy = params['gy']
    
    print("开始数值模拟...")
    
    # 设置随机种子以获得可重复的结果
    np.random.seed(42)
    
    # 数值模拟参数
    nsim = 10000
    ones_nsim_1 = np.ones(nsim)
    meanY = np.zeros(tn)
    meanC = np.zeros(tn)
    meanW = np.zeros(tn)
    meanA = np.zeros(tn)
    meanS = np.zeros(tn)
    meanB = np.zeros(tn)
    meanWY = np.zeros(tn)
    meanalpha = np.zeros(tn)
    meanGPY = np.zeros(tn)
    cGPY = np.zeros(tn)
    meanYs = np.zeros(tn)
    meanCs = np.zeros(tn)
    meanWs = np.zeros(tn)
    simPY = np.zeros((tn, nsim))
    simGPY = np.zeros((tn, nsim))
    simY = np.zeros((tn, nsim))
    simC = np.zeros((tn, nsim))
    simW = np.zeros((tn, nsim))
    simA = np.zeros((tn, nsim))
    simS = np.zeros((tn, nsim))
    simS = np.zeros((tn, nsim))
    simB = np.zeros((tn, nsim))
    simW_Y = np.zeros((tn, nsim))
    simR = np.zeros((tn, nsim))
    
    eps_y = np.zeros((1, 1))
    simTY = np.zeros((1, 1))
    eps_r = np.zeros((1, 1))
    cash = np.zeros((tn, nsim))
    simC_pct = np.zeros((tn, nsim))
    
    # 1、模拟生成labor income
    for i1 in range(nsim // 2):  # 另外一半模拟完全对称
        # working period第一期
        eps_y[0, 0] = f_randn(1)  # N(0,1)
        simPY[0, i1] = eps_y[0, 0] * smav  # 初始的p
        simPY[0, nsim // 2 + i1] = -eps_y[0, 0] * smav
        simGPY[0, i1] = 1.0
        simGPY[0, nsim // 2 + i1] = 1.0
        simTY[0, 0] = f_randn(1)
        simY[0, i1] = np.exp(simTY[0, 0] * smay)
        simY[0, nsim // 2 + i1] = np.exp(-simTY[0, 0] * smay)
        
        # working period第2期~退休
        for i2 in range(1, tr - tb):
            w = i2 + tb - 1
            eps_y[0, 0] = f_randn(1)
            simPY[i2, i1] = eps_y[0, 0] * smav + simPY[i2 - 1, i1]
            simPY[i2, nsim // 2 + i1] = -eps_y[0, 0] * smav + simPY[i2 - 1, nsim // 2 + i1]
            simGPY[i2, i1] = np.exp(gy[i2 - 1]) * np.exp(simPY[i2, i1]) / np.exp(simPY[i2 - 1, i1])
            simGPY[i2, nsim // 2 + i1] = np.exp(gy[i2 - 1]) * np.exp(simPY[i2, nsim // 2 + i1]) / np.exp(simPY[i2 - 1, nsim // 2 + i1])
            simTY[0, 0] = f_randn(1)
            simY[i2, i1] = np.exp(simTY[0, 0] * smay)
            simY[i2, nsim // 2 + i1] = np.exp(-simTY[0, 0] * smay)
    
    # 退休期
    for t in range(tr - tb, tn):
        simY[t, :] = ret_fac
        simGPY[t, :] = 1.0
    
    # 2、模拟风险投资的收益率
    for t in range(tn):
        for i1 in range(nsim // 2):
            eps_r[0, 0] = f_randn(1)
            simR[t, i1] = mu + rf + sigr * eps_r[0, 0]
            simR[t, nsim // 2 + i1] = mu + rf - sigr * eps_r[0, 0]
    
    # 从第一期开始迭代，得到各控制变量的值
    simW[:, :] = 0
    
    for t in tqdm(range(tn), desc="模拟进度"):
        # 向量化计算所有模拟的财富收入比和现金持有量
        simW_Y[t, :] = simW[t, :] / simY[t, :]  # 上期财富-本期工资收入比
        cash[t, :] = simW[t, :] + simY[t, :]    # cash-on-hand
        
        # 创建插值函数
        interp_C = interp1d(gcash, C[:, t], kind='cubic', fill_value='extrapolate')
        interp_A = interp1d(gcash, A[:, t], kind='cubic', fill_value='extrapolate')
        
        # 向量化应用插值
        simC[t, :] = interp_C(cash[t, :])
        simA[t, :] = interp_A(cash[t, :])
        
        # 向量化确保约束条件满足
        simC[t, :] = np.maximum(np.minimum(simC[t, :], 0.9999 * cash[t, :]), 0)
        simA[t, :] = np.minimum(np.maximum(simA[t, :], 0), 1)
        
        # 向量化计算各种模拟变量
        simC_pct[t, :] = simC[t, :] / cash[t, :]
        sav = cash[t, :] - simC[t, :]  # 用于投资的金额
        simS[t, :] = simA[t, :] * sav   # 风险投资额
        simS[t, :] = np.minimum(simS[t, :], sav)
        simB[t, :] = sav - simS[t, :]   # 无风险投资额
        
        # 更新消费比例
        simC[t, :] = simC[t, :] / cash[t, :]
        
        # 计算下期财富
        if t < tn - 1:
            simW[t + 1, :] = (simB[t, :] * rf + simS[t, :] * simR[t, :]) / simGPY[t + 1, :]
    
    # 多次模拟path下变量平均值
    meanC = np.mean(simC, axis=1)
    meanC_pct = np.mean(simC_pct, axis=1)
    meanY = np.mean(simY, axis=1)
    meanW = np.mean(simW, axis=1)
    meanS = np.mean(simS, axis=1)
    meanB = np.mean(simB, axis=1)
    meanWY = np.mean(simW_Y, axis=1)
    meanalpha = np.mean(simA, axis=1)
    meanGPY = np.mean(simGPY, axis=1)
    
    # 保存模拟结果
    os.makedirs('result_cocco_python', exist_ok=True)
    sim_results = pd.DataFrame({
        'meanC': meanC,
        'meanC_pct': meanC_pct, 
        'meanY': meanY,
        'meanW': meanW,
        'meanS': meanS,
        'meanB': meanB,
        'meanWY': meanWY,
        'meanalpha': meanalpha,
        'meanGPY': meanGPY
    })
    sim_results.to_excel('result_cocco_python/simulation_results.xlsx')
    
    # 绘制结果
    plt.figure(figsize=(10, 6))
    plt.plot(meanC_pct, label='c')
    plt.plot(meanalpha, label=r'$\alpha$')
    plt.legend()
    plt.xlabel('Period')
    plt.ylabel('Value')
    plt.title('Consumption Ratio and Portfolio Allocation')
    plt.savefig('result_cocco_python/results_plot.png')
    plt.show()
    
    print("数值模拟完成，结果已保存到result_cocco_python文件夹")

def main(mode="both"):
    """
    主函数，可选择运行模式：
    - "solve": 只进行模型求解
    - "simulate": 只进行数值模拟（需要先有求解结果）
    - "both": 先求解再模拟（默认）
    """
    if mode.lower() == "solve" or mode.lower() == "both":
        print("开始求解模型...")
        model_results = solve_model()
        print("模型求解完成")
        
        if mode.lower() == "both":
            print("开始数值模拟...")
            simulate_model(model_results)
    
    elif mode.lower() == "simulate":
        print("开始数值模拟...")
        simulate_model()
    
    else:
        print("无效的运行模式，可选值: 'solve', 'simulate', 'both'")

if __name__ == "__main__":
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(description='生命周期模型求解与模拟')
    parser.add_argument('--mode', type=str, choices=['solve', 'simulate', 'both'], 
                        default='both', help='运行模式：求解(solve)，模拟(simulate)，或两者都运行(both)')
    
    # 解析命令行参数
    args = parser.parse_args()
    
    # 运行主函数
    main(args.mode) 