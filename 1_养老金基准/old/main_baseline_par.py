import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import minimize
from scipy.interpolate import interp1d
import pandas as pd
import os
from tqdm import tqdm
import argparse
import time
import functools
# 添加并行计算所需的库
from joblib import Parallel, delayed
import multiprocessing
from scipy.interpolate import RectBivariateSpline

# 创建一个字典来存储每个函数的统计信息
function_stats = {}

# 创建插值器缓存
interpolator_cache = {}

def get_cached_interpolator(grid_data, x_grid, y_grid, period=None, var_name=None, kx=1, ky=1):
    """
    获取或创建插值器并缓存
    
    参数:
    grid_data - 二维值网格
    x_grid - x轴网格点
    y_grid - y轴网格点
    period - 时期标识
    var_name - 变量名称
    kx, ky - 插值阶数
    
    返回:
    interpolator - 插值器对象
    """
    # 如果未提供period或var_name，则不使用缓存
    if period is None or var_name is None:
        return RectBivariateSpline(x_grid, y_grid, grid_data, kx=kx, ky=ky)
    
    cache_key = f"{var_name}_{period}"
    
    if cache_key not in interpolator_cache:
        # 创建并缓存插值器
        interpolator_cache[cache_key] = RectBivariateSpline(
            x_grid, y_grid, grid_data, kx=kx, ky=ky
        )
        
        # 如果缓存过大，可以清理不再需要的插值器
        if len(interpolator_cache) > 200:  # 限制缓存大小
            # 寻找可能不再需要的缓存项
            keys_to_remove = []
            for key in interpolator_cache:
                if key != cache_key:  # 不删除刚刚创建的
                    old_period = int(key.split('_')[1])
                    if old_period < period - 5:  # 只保留近期的
                        keys_to_remove.append(key)
            
            # 清理旧缓存
            for key in keys_to_remove:
                del interpolator_cache[key]
    
    return interpolator_cache[cache_key]

def monitor(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # 获取函数名称
        func_name = func.__name__
        
        # 如果函数不在统计字典中，则添加
        if func_name not in function_stats:
            function_stats[func_name] = {'calls': 0, 'total_time': 0, 'min_time': float('inf'), 'max_time': 0}
        
        # 增加调用次数
        function_stats[func_name]['calls'] += 1
        
        # 记录开始时间
        start_time = time.time()
        
        # 调用原始函数
        result = func(*args, **kwargs)
        
        # 计算执行时间
        execution_time = time.time() - start_time
        
        # 更新统计信息
        function_stats[func_name]['total_time'] += execution_time
        function_stats[func_name]['min_time'] = min(function_stats[func_name]['min_time'], execution_time)
        function_stats[func_name]['max_time'] = max(function_stats[func_name]['max_time'], execution_time)
        
        return result
    
    return wrapper

# 定义一个函数来打印统计信息
def print_function_stats():
    print("\nFunction Statistics:")
    print("-" * 80)
    print(f"{'Function Name':<30} {'Calls':<10} {'Total Time':<15} {'Avg Time':<15} {'Min Time':<15} {'Max Time':<15}")
    print("-" * 80)
    
    for func_name, stats in function_stats.items():
        avg_time = stats['total_time'] / stats['calls'] if stats['calls'] > 0 else 0
        print(f"{func_name:<30} {stats['calls']:<10} {stats['total_time']:<15.6f} {avg_time:<15.6f} {stats['min_time']:<15.6f} {stats['max_time']:<15.6f}")
        
# 辅助函数
def f_randn(size):
    """生成标准正态分布随机数"""
    return np.random.randn(size)


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
@monitor
def cocco_fun_valuefunc_retired(x, cash, fund, nextV, gret, rf, ret_fac, pension_pay, gamma, beta, psi_1, psi_2, theta, gcash, gfund, survprob, weig, t=None):
    """
    退休期间的值函数计算
    
    参数:
    x - 决策变量 [消费, 风险资产投资比例]
    cash - 手中现金
    fund - 养老基金账户余额
    nextV - 下期值函数
    nextF - 下期养老基金值函数
    gret - 随机收益率
    rf - 无风险收益率
    ret_fac - 固定基本养老金
    pension_pay - 个人养老金给付
    gamma - 相对风险规避系数
    beta - 贴现因子
    psi_1, psi_2, theta - Epstein-Zin效用函数参数
    gcash - 现金网格
    gfund - 养老基金网格
    survprob - 存活概率
    weig - 随机状态权重
    t - 当前时期（用于缓存）
    
    返回:
    v - 函数值
    """
    auxVV = 0
    sav = cash*(1 - x[0])
    
    # Epstein-Zin效用函数
    if gamma == 1:
        u = np.log(sav*x[0])  # 对数效用
    else:
        u = (cash*x[0]+1e-07)**psi_1  # CRRA效用
    
    # 下期的cash-on-hand
    cash_1 = (rf * (1 - x[1]) + gret * x[1]) * sav + pension_pay
    
    # 使用缓存的一维插值
    int_V = np.maximum(interp1d(gcash, nextV, kind='cubic', fill_value='extrapolate')(cash_1), 1e-20)
    
    # 计算期望值
    auxVV = auxVV + np.dot(weig, survprob * int_V**(1-gamma))
    
    if auxVV == 0:
        v = -((1-beta)*u)**(psi_2)  # 没有未来值函数时
    else:
        # Epstein-Zin递归效用
        v = -(((1-beta)*u + beta*((auxVV)**(1/theta)))**psi_2)
    
    return v

# 退休期间进行策略评估的函数

# 工作期间的值函数计算
@monitor
def cocco_fun_valuefunc_work(x, cash, fund, gyp, nextV, yh, gret, rf, rp, tau, gamma, beta, psi_1, psi_2, theta, gcash, gfund, survprob, n, weig, t=None):
    """
    工作期间的值函数计算
    
    参数:
    x - 决策变量 [消费, 风险资产投资比例, 个人养老金购买比例]
    cash - 手中现金
    fund - 养老基金账户余额
    gyp - 劳动收入增长率
    nextV - 下期值函数
    nextF - 下期养老基金值函数
    yh - 暂时收入冲击
    gret - 随机收益率
    rf - 无风险收益率
    rp - 个人养老基金收益率
    tau - 用于缴纳基本养老金的工资比例
    gamma - 相对风险规避系数
    beta - 贴现因子
    psi_1, psi_2, theta - Epstein-Zin效用函数参数
    gcash - 现金网格
    gfund - 养老基金网格
    survprob - 存活概率
    n - 随机状态数量
    weig - 随机状态权重
    t - 当前时期（用于缓存）
    
    返回:
    v - 函数值
    """
    auxVV = 0
        # Epstein-Zin效用函数
    if gamma == 1:
        u = np.log(cash*x[0])  # 对数效用
    else:
        u = (cash*x[0]+1e-07)**psi_1  # CRRA效用   

      
    # 创建索引网格
    # 'ij'表示使用矩阵索引顺序，而不是默认的笛卡尔坐标顺序
    # 这里创建三维网格，i6对应收益率状态，i7对应暂时收入冲击，i8对应持久收入冲击
    i6_grid, i7_grid, i8_grid = np.meshgrid(np.arange(n), np.arange(n), np.arange(n), indexing='ij')

    labor_income = yh[i7_grid, i8_grid]  # 下一期劳动收入  
    
    # 计算所有可能的储蓄值和个人养老金购买值
    pension_contrib = x[2] * cash * (1 - x[0])  # 个人养老金购买
    sav_values = cash*(1 - x[0]) *(1 - x[2])/ gyp[i6_grid, i8_grid]  # 剩余储蓄
    
    # 计算下期养老基金账户余额
    fund_1_values = (fund + pension_contrib) * rp / gyp[i6_grid, i8_grid] + tau * labor_income
    
    # 计算所有可能的现金值
    portfolio_return = rf * (1 - x[1]) + gret[i8_grid] * x[1]

    cash_1_values = portfolio_return * sav_values + (1-tau) * labor_income
    
    # 限制cash_1和fund_1的范围
    cash_1_values = np.clip(cash_1_values, gcash[0], gcash[-1])
    fund_1_values = np.clip(fund_1_values, gfund[0], gfund[-1])
    
    # 使用缓存的插值器而非重新创建
    if t is not None:
        interpolator_V = get_cached_interpolator(nextV, gcash, gfund, t+1, "V", kx=1, ky=1)
    else:
        # 兼容旧代码，如果未提供t则回退到原始方式
        interpolator_V = RectBivariateSpline(gcash, gfund, nextV, kx=1, ky=1)
    
    # cash_1_values和fund_1_values是3x3x3的矩阵，需要重塑后进行插值
    original_shape = cash_1_values.shape
    cash_1_flat = cash_1_values.flatten()
    fund_1_flat = fund_1_values.flatten()
    
    # 对扁平化后的数组进行插值
    int_V_flat = interpolator_V.ev(cash_1_flat, fund_1_flat)
    
    # 将结果重塑回原始形状
    int_V_values = int_V_flat.reshape(original_shape)
    
    # 计算期望值
    auxVV = np.sum(weig * survprob * (int_V_values * gyp[i6_grid, i8_grid])**(1-gamma))

    if auxVV == 0:
        v = -((1-beta)*u)**(psi_2)  # 没有未来值函数时
    else:
        # Epstein-Zin递归效用
        v = -(((1-beta)*u + beta*((auxVV)**(1/theta)))**psi_2)
    
    return v

# 添加并行优化任务函数 - 工作期
@monitor
def optimize_policy_work(i2, i3, x0, gcash, gfund, gyp, V_next, t, yh, gret, rf, rp, 
                        pension_pct, gamma, beta, psi_1, psi_2, theta, survprob, n, nweig1):
    """单个状态点的工作期策略优化函数，用于并行计算"""
    lb = [0, 0, 0]
    ub = [1, 1, 1]
    
    # 确保初始值在边界范围内
    x0[0] = max(min(x0[0], ub[0]), lb[0])
    x0[1] = max(min(x0[1], ub[1]), lb[1])
    x0[2] = max(min(x0[2], ub[2]), lb[2])
    
    # 优化函数
    result = minimize(
        lambda x: cocco_fun_valuefunc_work(
            x, gcash[i3], gfund[i2], 
            gyp, V_next, 
            yh, gret, rf, rp, pension_pct, 
            gamma, beta, psi_1, psi_2, theta, 
            gcash, gfund, survprob, n, nweig1, t  # 添加时期参数t
        ),
        x0, bounds=list(zip(lb, ub)),  method='SLSQP', 
        options={'ftol': 1e-7, 'maxiter': 200}
    )
    
    return (i2, i3, result.x)

# 添加并行优化任务函数 - 退休期
@monitor
def optimize_policy_retired(i2, i3, x0, gcash, gfund, V_next, t, gret, rf, ret_fac, pension_rate, 
                           gamma, beta, psi_1, psi_2, theta, survprob, weig):
    """单个状态点的退休期策略优化函数，用于并行计算"""
    lb = [0, 0]
    ub = [1, 1]
    
    # 确保初始值在边界范围内
    x0[0] = max(min(x0[0], ub[0]), lb[0])
    x0[1] = max(min(x0[1], ub[1]), lb[1])
    
    # 计算养老金给付
    pension_pay = pension_rate * gfund[i2]
    
    # 优化函数
    result = minimize(
        lambda x: cocco_fun_valuefunc_retired(
            x, gcash[i3], gfund[i2], 
            V_next, 
            gret, rf, ret_fac, pension_pay, 
            gamma, beta, psi_1, psi_2, theta, 
            gcash, gfund, survprob, weig, t  # 添加时期参数t
        ),
        x0, bounds=list(zip(lb, ub)), method='SLSQP', 
        options={'ftol': 1e-7, 'maxiter': 200}
    )
    
    return (i2, i3, result.x)

# 添加一个函数来检查值函数范围并在必要时进行缩放
def scale_value_function(V_matrix, min_range=0.01):
    """
    检查值函数的范围，如果范围小于阈值则进行缩放
    
    参数:
    V_matrix - 要检查的值函数矩阵
    min_range - 最小可接受的值函数范围
    
    返回:
    scaled_V - 缩放后的值函数
    scaling_factor - 使用的缩放因子
    scaled - 是否进行了缩放
    """
    v_min = np.min(V_matrix)
    v_max = np.max(V_matrix)
    v_range = v_max - v_min
    
    # 如果值函数范围小于阈值，则进行缩放
    if v_range < min_range:
        # 计算需要的缩放因子，使得缩放后的范围至少为min_range
        scaling_factor = min_range / v_range
        # 缩放值函数，保持相对关系不变
        scaled_V = v_min + (V_matrix - v_min) * scaling_factor
        return scaled_V, scaling_factor, True
    
    # 如果范围足够大，则不需要缩放
    return V_matrix, 1.0, False

def solve_model():
    """
    使用策略函数迭代方法求解模型的最优消费和投资决策
    """
    # 设置随机种子以获得可重复的结果
    np.random.seed(42)
    # 创建保存每期数据的目录
    os.makedirs('result_baseline_python_PFI', exist_ok=True)
    os.makedirs('result_baseline_python_PFI/slices', exist_ok=True)
    os.makedirs('result_baseline_python_PFI/slices/retire', exist_ok=True)
    os.makedirs('result_baseline_python_PFI/slices/work', exist_ok=True)
    # 变量定义
    tb = 18    # 初始开始的年纪
    tr = 61    # 退休年龄
    td = 100   # 死亡年龄
    tn = td - tb + 1  # 总期数
    # 状态变量grid数量
    ncash = 51  # 手中现金
    nfund = 51  # 养老基金余额
    n = 3       # 外生随机冲击的grid数量
    # 外生参数
    # 基础收入f的系数
    aa = (-2.170042 + 2.700381)
    b1 = 0.16818
    b2 = -0.0323371 / 10
    b3 = 0.0019704 / 100
    # 养老金相关
    ret_fac = 0.6827  # 退休后固定支付的工资（基本养老金）
    pension_pct = 0.08   # 工资扣除缴纳基本养老保险的比例（个人账户部分）
    rp = 1.04    # 个人养老基金的无风险收益率
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
    # PFI迭代参数
    max_iterations = 400  # 最大迭代次数
    tolerance = 1e-5      # 收敛容差    
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
    gfund = np.zeros(nfund)
    lgfund = np.zeros(nfund)
    secd = np.zeros(ncash)
    C = np.zeros((ncash, nfund, tn))
    V = np.zeros((ncash, nfund, tn))
    A = np.ones((ncash, nfund, tn))
    Q = np.zeros((ncash, nfund, tn))  # 个人养老金购买决策
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
        gcash[i1] = np.exp(l_mincash + i1 * stepcash)
            
    
    # fund的grid
    maxfund = 100
    minfund = 0.25
    l_maxfund = np.log(maxfund)
    l_minfund = np.log(minfund)
    stepfund = (l_maxfund - l_minfund) / (nfund - 1)
    
    for i1 in range(nfund):
        gfund[i1] = np.exp(l_minfund + i1 * stepfund)
        
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
    for i1 in range(ncash):
        for i2 in range(nfund):
            C[i1, i2, tn - 1] = 1 # 最后一期的现金全部用来消费
            A[i1, i2, tn - 1] = 0.0      # 最后一期的投资全部为零
            Q[i1, i2, tn - 1] = 0.0      # 最后一期不购买养老金
            V[i1, i2, tn - 1] = (C[i1, i2, tn - 1]*gcash[i1]) * (1 - beta)**psi_2
    
    # 计算个人养老金给付
    # 假设退休后按照剩余年份平均分配养老基金余额
    pension_period = td - tr + 1
    pension_rate = 1 / pension_period
    
    # 确定可用CPU核心数，留出一个核心给系统
    num_cores = max(1, multiprocessing.cpu_count())
    print(f"使用 {num_cores} 个CPU核心进行并行计算")
    
    # 策略函数迭代
    # 退休期
    print("开始求解退休期值函数（使用策略函数迭代）...")
    for i1 in tqdm(range(1, td - tr + 1), desc="退休期求解进度"):
        t = tn - i1 - 1
        
        # 初始化策略 - 可以使用合理的猜测
        for i2 in range(nfund):  # 遍历养老基金状态
            for i3 in range(ncash):  # 遍历现金状态
                # # 初始策略猜测
                # if i3 == 0:
                #     C[i3, i2, t] = 0.8
                #     A[i3, i2, t] = 0.5
                # else:
                #     # 使用上一个状态的策略作为初始猜测
                C[i3, i2, t] = C[i3, i2, t+1]
                A[i3, i2, t] = A[i3, i2, t+1]

                # 退休期不需要购买个人养老金
                Q[i3, i2, t] = 0.0        
        # 策略迭代循环
        for iter_count in range(max_iterations):
            # 策略评估 - 计算当前策略的值函数
            V_new = np.zeros((ncash, nfund))
            
            # 并行计算值函数 - 这部分可以不做并行，因为评估函数本身已经高度向量化
            for i2 in range(nfund):  # 遍历养老基金状态
                for i3 in range(ncash):  # 遍历现金状态
                    policy = np.array([C[i3, i2, t], A[i3, i2, t]])
                    
                    # 计算养老金给付
                    pension_pay = pension_rate * gfund[i2]
                    
                    V_new[i3, i2] = -cocco_fun_valuefunc_retired(
                        policy, gcash[i3], gfund[i2], 
                        V[:, i2, t+1],
                        gret, rf, ret_fac, pension_pay, 
                        gamma, beta, psi_1, psi_2, theta, 
                        gcash, gfund, survprob[t], weig, t
                    )
            
            # 检查并缩放值函数，如果范围太小
            V_new, scaling_factor, was_scaled = scale_value_function(V_new)
            if was_scaled:
                print(f"  退休期时期 {t}: 值函数范围过小，应用了缩放因子 {scaling_factor:.4f}")
            
            # 策略改进 - 基于当前值函数更新策略，使用并行计算
            C_new = np.zeros((ncash, nfund))
            A_new = np.zeros((ncash, nfund))
            
            # 修改并行策略：只在养老基金维度(i2)进行并行化
            # 每个并行任务处理一个养老基金水平的所有现金状态
            def optimize_fund_level(i2):
                """处理单个养老基金水平的所有现金状态"""
                results = []
                for i3 in range(ncash):
                    x0 = [C[i3, i2, t], A[i3, i2, t]]
                    # 计算养老金给付
                    pension_pay = pension_rate * gfund[i2]
                    
                    # 优化单个状态点
                    result = minimize(
                        lambda x: cocco_fun_valuefunc_retired(
                            x, gcash[i3], gfund[i2], 
                            V[:, i2, t+1],
                            gret, rf, ret_fac, pension_pay, 
                            gamma, beta, psi_1, psi_2, theta, 
                            gcash, gfund, survprob[t], weig, t
                        ),
                        x0, bounds=[(0,1), (0,1)], method='SLSQP', 
                        options={'ftol': 1e-7, 'maxiter': 200}
                    )
                    
                    results.append((i3, result.x))
                return i2, results
            
            # 并行执行养老基金维度的优化
            fund_results = Parallel(n_jobs=num_cores)(
                delayed(optimize_fund_level)(i2) for i2 in range(nfund)
            )
            
            # 处理并行计算结果
            for i2, cash_results in fund_results:
                for i3, policy in cash_results:
                    C_new[i3, i2] = policy[0]
                    A_new[i3, i2] = policy[1]
            
            # 计算策略变化量
            policy_change = (np.max(np.abs(C_new - C[:, :, t])) + 
                            np.max(np.abs(A_new - A[:, :, t])))
            
            # 自适应步长策略
            alpha = 0.1  # 接近收敛阶段使用较小步长
            
            # 更新策略
            C[:, :, t] = (1-alpha) * C[:, :, t] + alpha * C_new
            A[:, :, t] = (1-alpha) * A[:, :, t] + alpha * A_new
            V[:, :, t] = V_new
            
            # 检查收敛
            if policy_change < tolerance:
                print(f"  退休期时期 {t} 在 {iter_count+1} 次迭代后收敛")
                break
            
            # 如果达到最大迭代次数仍未收敛
            if iter_count == max_iterations - 1:
                print(f"  警告：退休期时期 {t} 在 {max_iterations} 次迭代后仍未收敛")
            
            # 保存当前时期的切片数据
            slice_dir = f'result_baseline_python_PFI/slices/retire/period_{t}'
            os.makedirs(slice_dir, exist_ok=True)
            pd.DataFrame(C[:, :, t]).to_excel(f'{slice_dir}/C_slice_t{t}.xlsx')
            pd.DataFrame(A[:, :, t]).to_excel(f'{slice_dir}/A_slice_t{t}.xlsx')
            pd.DataFrame(Q[:, :, t]).to_excel(f'{slice_dir}/Q_slice_t{t}.xlsx')
            pd.DataFrame(V[:, :, t]).to_excel(f'{slice_dir}/V_slice_t{t}.xlsx')
                        
            # 保存为numpy数组以便于后续加载（可选）
            np.save(f'{slice_dir}/C_t{t}.npy', C[:, :, t])
            np.save(f'{slice_dir}/A_t{t}.npy', A[:, :, t])
            np.save(f'{slice_dir}/Q_t{t}.npy', Q[:, :, t])
            np.save(f'{slice_dir}/V_t{t}.npy', V[:, :, t])

    # 工作期（退休前）
    print("开始求解工作期值函数（使用策略函数迭代）...")
    for i1 in tqdm(range(1, tr - tb + 1), desc="工作期求解进度"):
        t = tr - tb - i1
        
        # 初始化策略 - 使用更智能的初始猜测
        for i2 in range(nfund):  # 遍历养老基金状态
            for i3 in range(ncash):  # 遍历现金状态
                # 使用下一期的决策作为初始猜测
                C[i3, i2, t] = C[i3, i2, t+1]
                A[i3, i2, t] = A[i3, i2, t+1]
                Q[i3, i2, t] = 0.1
        
        # 策略迭代循环
        for iter_count in range(max_iterations):
            # 策略评估 - 计算当前策略的值函数
            V_new = np.zeros((ncash, nfund))
            
            for i2 in range(nfund):  # 遍历养老基金状态
                for i3 in range(ncash):  # 遍历现金状态
                    policy = np.array([C[i3, i2, t], A[i3, i2, t], Q[i3, i2, t]])
                    
                    V_new[i3, i2] = -cocco_fun_valuefunc_work(
                        policy, gcash[i3], gfund[i2], 
                        gyp[:, :, t], V[:, :, t+1],
                        yh, gret, rf, rp, pension_pct, 
                        gamma, beta, psi_1, psi_2, theta, 
                        gcash, gfund, survprob[t], n, nweig1, t  # 添加时期参数t
                    )
            
            # 检查并缩放值函数，如果范围太小
            V_new, scaling_factor, was_scaled = scale_value_function(V_new)
            if was_scaled:
                # print(f"  工作期时期 {t}: 值函数范围过小，应用了缩放因子 {scaling_factor:.4f}")
                
                # 对未来值函数进行状态依赖的缩放，使得值函数差异体现不同状态的重要性
                if t < 20:  # 早期工作期特别处理
                    # 保存原始值函数的拷贝，用于打印诊断信息
                    V_before = V_new.copy()
                    
                    # 应用状态依赖的缩放因子
                    for i3 in range(ncash):
                        cash_position = (gcash[i3] - gcash[0]) / (gcash[-1] - gcash[0])
                        for i2 in range(nfund):
                            fund_position = (gfund[i2] - gfund[0]) / (gfund[-1] - gfund[0])
                            # 状态依赖的缩放：现金和养老金越多，值函数差异越大
                            state_factor = 1.0 + 5.0 * (cash_position + fund_position)
                            V_new[i3, i2] *= state_factor
                    
                    # 打印诊断信息
                    v_range_before = np.max(V_before) - np.min(V_before)
                    v_range_after = np.max(V_new) - np.min(V_new)
                    # print(f"    状态依赖缩放: 范围从 {v_range_before:.6f} 增加到 {v_range_after:.6f}")
            
            # 策略改进 - 基于当前值函数更新策略，使用并行计算
            C_new = np.zeros((ncash, nfund))
            A_new = np.zeros((ncash, nfund))
            Q_new = np.zeros((ncash, nfund))
            
            # 修改并行策略：只在养老基金维度(i2)进行并行化
            # 每个并行任务处理一个养老基金水平的所有现金状态
            def optimize_fund_level_work(i2):
                """处理单个养老基金水平的所有现金状态 - 工作期"""
                results = []
                for i3 in range(ncash):
                    x0 = [C[i3, i2, t], A[i3, i2, t], Q[i3, i2, t]]
                    
                    # 优化单个状态点
                    result = minimize(
                        lambda x: cocco_fun_valuefunc_work(
                            x, gcash[i3], gfund[i2], 
                            gyp[:, :, t], V[:, :, t+1],
                            yh, gret, rf, rp, pension_pct, 
                            gamma, beta, psi_1, psi_2, theta, 
                            gcash, gfund, survprob[t], n, nweig1, t
                        ),
                        x0, bounds=[(0,1), (0,1), (0,1)], method='SLSQP', 
                        options={'ftol': 1e-7, 'maxiter': 200}
                    )
                    
                    results.append((i3, result.x))
                return i2, results
            
            # 并行执行养老基金维度的优化
            fund_results = Parallel(n_jobs=num_cores)(
                delayed(optimize_fund_level_work)(i2) for i2 in range(nfund)
            )
            
            # 处理并行计算结果
            for i2, cash_results in fund_results:
                for i3, policy in cash_results:
                    C_new[i3, i2] = policy[0]
                    A_new[i3, i2] = policy[1]
                    Q_new[i3, i2] = policy[2]
            
            # 计算策略变化量
            policy_change = (np.max(np.abs(C_new - C[:, :, t])) + 
                            np.max(np.abs(A_new - A[:, :, t])) + 
                           np.max(np.abs(Q_new - Q[:, :, t])))
            
            # 自适应步长策略
            alpha = 0.1  # 接近收敛阶段使用较小步长
            
            # 更新策略
            C[:, :, t] = (1-alpha) * C[:, :, t] + alpha * C_new
            A[:, :, t] = (1-alpha) * A[:, :, t] + alpha * A_new
            Q[:, :, t] = (1-alpha) * Q[:, :, t] + alpha * Q_new
            V[:, :, t] = V_new
            
            # 检查收敛
            if policy_change < tolerance:
                print(f"  工作期时期 {t} 在 {iter_count+1} 次迭代后收敛")
                break
            
            # 如果达到最大迭代次数仍未收敛
            if iter_count == max_iterations - 1:
                print(f"  警告：工作期时期 {t} 在 {max_iterations} 次迭代后仍未收敛")
            
            # 保存当前时期的切片数据
            slice_dir = f'result_baseline_python_PFI/slices/work/period_{t}'
            os.makedirs(slice_dir, exist_ok=True)
            pd.DataFrame(C[:, :, t]).to_excel(f'{slice_dir}/C_slice_t{t}.xlsx')
            pd.DataFrame(A[:, :, t]).to_excel(f'{slice_dir}/A_slice_t{t}.xlsx')
            pd.DataFrame(Q[:, :, t]).to_excel(f'{slice_dir}/Q_slice_t{t}.xlsx')
            pd.DataFrame(V[:, :, t]).to_excel(f'{slice_dir}/V_slice_t{t}.xlsx')
            
            
            # 保存为numpy数组以便于后续加载（可选）
            np.save(f'{slice_dir}/C_t{t}.npy', C[:, :, t])
            np.save(f'{slice_dir}/A_t{t}.npy', A[:, :, t])
            np.save(f'{slice_dir}/Q_t{t}.npy', Q[:, :, t])
            np.save(f'{slice_dir}/V_t{t}.npy', V[:, :, t])
    
    # 保存结果
    os.makedirs('result_baseline_python_PFI', exist_ok=True)
    
    # 保存完整的三维策略函数和值函数
    np.save('result_baseline_python_PFI/C.npy', C)
    np.save('result_baseline_python_PFI/A.npy', A)
    np.save('result_baseline_python_PFI/Q.npy', Q)
    np.save('result_baseline_python_PFI/V.npy', V)
    np.save('result_baseline_python_PFI/gcash.npy', gcash)
    np.save('result_baseline_python_PFI/gfund.npy', gfund)
    
    # 同时保存一个切片用于Excel查看（可选）
    pd.DataFrame(C[:, 0, :]).to_excel('result_baseline_python_PFI/C_slice.xlsx')
    pd.DataFrame(A[:, 0, :]).to_excel('result_baseline_python_PFI/A_slice.xlsx')
    pd.DataFrame(Q[:, 0, :]).to_excel('result_baseline_python_PFI/Q_slice.xlsx')
    pd.DataFrame(V[:, 0, :]).to_excel('result_baseline_python_PFI/V_slice.xlsx')
    pd.DataFrame(gcash).to_excel('result_baseline_python_PFI/gcash.xlsx')
    pd.DataFrame(gfund).to_excel('result_baseline_python_PFI/gfund.xlsx')

    print("模型求解完成，结果已保存到result_baseline_python_PFI文件夹")
    
    # 返回求解结果，供模拟使用
    model_results = {
        'C': C,
        'A': A,
        'Q': Q,
        'V': V,
        'gcash': gcash,
        'gfund': gfund,
        'params': {
            'tb': tb,
            'tr': tr,
            'td': td,
            'tn': tn,
            'n': n,
            'rf': rf,
            'rp': rp,
            'mu': mu,
            'sigr': sigr,
            'smay': smay,
            'smav': smav,
            'ret_fac': ret_fac,
            'pension_pct': pension_pct,
            'pension_rate': pension_rate,
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
            if not os.path.exists('result_baseline_python_PFI/C.npy'):
                print("未找到模型求解结果文件，请先运行模型求解")
                return
                
            # 从文件中读取模型结果
            C = np.load('result_baseline_python_PFI/C.npy')
            A = np.load('result_baseline_python_PFI/A.npy')
            Q = np.load('result_baseline_python_PFI/Q.npy')
            V = np.load('result_baseline_python_PFI/V.npy')
            gcash = np.load('result_baseline_python_PFI/gcash.npy')
            gfund = np.load('result_baseline_python_PFI/gfund.npy')
            
            # 设置默认参数
            tb = 18    # 初始开始的年纪
            tr = 61    # 退休年龄
            td = 100   # 死亡年龄
            tn = td - tb + 1  # 总期数
            n = 5       # 外生随机冲击的grid数量
            rf = 1.02    # 无风险总收入
            rp = 1.02    # 个人养老基金收益率
            mu = 0.04    # 超额收益
            sigr = 0.27  # 风险资产收益率的标准差
            smay = np.sqrt(0.169993)  # 白噪声shock的标准差
            smav = np.sqrt(0.112572)  # 持续shock的标准差
            ret_fac = 0.6827  # 退休后固定支付的工资（养老金）
            pension_pct = 0.08  # 工资扣除缴纳基本养老金的比例（个人账户部分）
            pension_period = td - tr + 1
            pension_rate = 1 / pension_period
            
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
                'Q': Q,
                'V': V,
                'gcash': gcash,
                'gfund': gfund,
                'params': {
                    'tb': tb,
                    'tr': tr,
                    'td': td,
                    'tn': tn,
                    'n': n,
                    'rf': rf,
                    'rp': rp,
                    'mu': mu,
                    'sigr': sigr,
                    'smay': smay,
                    'smav': smav,
                    'ret_fac': ret_fac,
                    'pension_pct': pension_pct,
                    'pension_rate': pension_rate,
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
    Q = model_results['Q']
    V = model_results['V']
    gcash = model_results['gcash']
    gfund = model_results['gfund']
    
    # 获取参数
    params = model_results['params']
    tb = params['tb']
    tr = params['tr']
    td = params['td']
    tn = params['tn']
    n = params['n']
    rf = params['rf']
    rp = params['rp']
    mu = params['mu']
    sigr = params['sigr']
    smay = params['smay']
    smav = params['smav']
    ret_fac = params['ret_fac']
    pension_pct = params['pension_pct']
    pension_rate = params['pension_rate']
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
    meanQ = np.zeros(tn)  # 平均个人养老金购买
    meanF = np.zeros(tn)  # 平均养老基金余额
    meanP = np.zeros(tn)  # 平均个人养老金给付
    cGPY = np.zeros(tn)
    meanYs = np.zeros(tn)
    meanCs = np.zeros(tn)
    meanWs = np.zeros(tn)
    simPY = np.zeros((tn, nsim))
    simGPY = np.zeros((tn, nsim))
    simY = np.zeros((tn, nsim))
    simC = np.zeros((tn, nsim))
    simW = np.zeros((tn, nsim))
    simS = np.zeros((tn, nsim))
    simA = np.zeros((tn, nsim))
    simB = np.zeros((tn, nsim))
    simS = np.zeros((tn, nsim))
    simW_Y = np.zeros((tn, nsim))
    simR = np.zeros((tn, nsim))
    simQ = np.zeros((tn, nsim))  # 个人养老金购买决策
    simQ_pct = np.zeros((tn, nsim))  # 个人养老金购买比例
    simF = np.zeros((tn, nsim))  # 养老基金账户余额
    simP = np.zeros((tn, nsim))  # 个人养老金给付
    eps_y = np.zeros((1, 1))
    simTY = np.zeros((1, 1))
    eps_r = np.zeros((1, 1))
    cash = np.zeros((tn, nsim))
    simC_pct = np.zeros((tn, nsim))
    

    time_start = time.time()
    # 1、模拟生成labor income
    for i1 in range(nsim // 2): # 另外一半模拟完全对称
        # working period第一期
        eps_y[0, 0] = f_randn(1) # N(0,1)
        simPY[0, i1] = eps_y[0, 0] * smav # 初始的p
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
    time_end = time.time()
    print(f"模拟风险投资的收益率时间: {time_end - time_start} 秒")
    
    # 从第一期开始迭代，得到各控制变量的值
    simW[:, :] = 0.2
    simF[:, :] = 0
    simF[0, :] = pension_pct * simY[0, :]
    # 添加对双线性插值函数的监控
    @monitor
    def get_bilinear_interp(value_grid, gcash, gfund, cash_value, fund_value):
        """
        对二维网格进行双线性插值
        
        参数:
        value_grid - 待插值的二维值函数网格
        gcash - 现金网格
        gfund - 养老基金网格
        cash_value - 要插值的现金值
        fund_value - 要插值的养老基金值
        
        返回:
        interpolated_value - 插值结果
        """
        # 限制输入值在网格范围内
        cash_value = max(min(cash_value, gcash[-1]), gcash[0])
        fund_value = max(min(fund_value, gfund[-1]), gfund[0])
        
        # 找到在gcash网格中的位置
        cash_idx = np.searchsorted(gcash, cash_value) - 1
        cash_idx = max(0, min(cash_idx, len(gcash) - 2))
        cash_weight = (cash_value - gcash[cash_idx]) / (gcash[cash_idx + 1] - gcash[cash_idx])
        
        # 找到在gfund网格中的位置
        fund_idx = np.searchsorted(gfund, fund_value) - 1
        fund_idx = max(0, min(fund_idx, len(gfund) - 2))
        fund_weight = (fund_value - gfund[fund_idx]) / (gfund[fund_idx + 1] - gfund[fund_idx])
        
        # 双线性插值
        v00 = value_grid[cash_idx, fund_idx]
        v01 = value_grid[cash_idx, fund_idx + 1]
        v10 = value_grid[cash_idx + 1, fund_idx]
        v11 = value_grid[cash_idx + 1, fund_idx + 1]
        
        interp_value = (1 - cash_weight) * (1 - fund_weight) * v00 + \
                      (1 - cash_weight) * fund_weight * v01 + \
                      cash_weight * (1 - fund_weight) * v10 + \
                      cash_weight * fund_weight * v11
        
        return interp_value

    # 开始模拟
    print("预先计算所有时期的插值器以提高模拟效率...")
    
    # 预先计算所有时期的插值器
    interpolator_cache.clear()  # 清空缓存，重新开始
    for t in range(tn):
        _ = get_cached_interpolator(C[:, :, t], gcash, gfund, t, "C")
        _ = get_cached_interpolator(A[:, :, t], gcash, gfund, t, "A")
        if t < tr - tb:  # 工作期才需要Q的插值器
            _ = get_cached_interpolator(Q[:, :, t], gcash, gfund, t, "Q")
    
    print(f"预计算完成，共缓存 {len(interpolator_cache)} 个插值器")
    
    for t in tqdm(range(tn), desc="模拟进度"):
        # 向量化计算所有模拟的财富收入比和现金持有量
        if t < tr - tb:  # 工作期
            simW_Y[t, :] = simW[t, :] / simY[t, :]  # 上期财富-本期工资收入比
            cash[t, :] = simW[t, :]  # 当期现金为上期的财富
        else:  # 退休期
            cash[t, :] = simW[t, :]
        
        # 根据当期的状态查找或插值获取策略
        for i in range(nsim):
            # 使用缓存的插值器获取当前状态下的策略
            interpolator_C = get_cached_interpolator(C[:, :, t], gcash, gfund, t, "C")
            interpolator_A = get_cached_interpolator(A[:, :, t], gcash, gfund, t, "A")
            
            # 使用RectBivariateSpline进行二维插值
            simC[t, i] = interpolator_C.ev(cash[t, i], simF[t, i])
            simA[t, i] = interpolator_A.ev(cash[t, i], simF[t, i])
            
            # 确保约束条件满足
            simC[t, i] = max(min(simC[t, i], 1), 0)
            simA[t, i] = max(min(simA[t, i], 1), 0)
            
            # 工作期，还需要获取养老金购买比例
            if t < tr - tb:
                interpolator_Q = get_cached_interpolator(Q[:, :, t], gcash, gfund, t, "Q")
                simQ[t, i] = interpolator_Q.ev(cash[t, i], simF[t, i])
                simQ[t, i] = max(min(simQ[t, i], 1), 0)
        
        # 记录消费占收入的比例
        if t < tr - tb:  # 工作期
            simC_pct[t, :] = simC[t, :]
            simQ_pct[t, :] = simQ[t, :] * (1-simC[t, :])
        else:  # 退休期
            simC_pct[t, :] = simC[t, :]
            simQ_pct[t, :] = 0

        if t == tn-1:
            simC_pct[t, :] = 1

        
        # 计算个人养老金购买金额和养老基金账户更新
        if t < tn - 1:
            if t < tr - tb:  # 工作期
                # 个人养老金购买
                pension_contrib = simQ[t, :] * cash[t, :] * (1-simC[t, :])
                
                # 计算养老基金账户更新
                simF[t+1, :] = (simF[t, :]+ pension_contrib)* rp + pension_pct * simY[t+1, :]
            else:  # 退休期
                # 计算个人养老金给付
                simP[t+1, :] = pension_rate * simF[tr-tb-1, :]        
                # 养老基金余额保持不变
                simF[t+1, :] = simF[t, :]
        
        # 剩余用于投资的资金
        sav = cash[t, :]*(1-simC[t, :])
        
        # 工作期扣除个人养老金购买金额
        if t < tr - tb:
            sav = (1-simQ[t, :]) * cash[t, :] * (1-simC[t, :])
        
        # 风险投资额
        simS[t, :] = simA[t, :] * sav
        
        # 无风险投资额
        simB[t, :] = sav - simS[t, :]
        
        # 计算下期财富
        if t < tn - 1:
            # 工作期
            if t < tr - tb - 1:
                simW[t + 1, :] = (simB[t, :] * rf + simS[t, :] * simR[t, :]) / simGPY[t + 1, :] + (1 - pension_pct) * simY[t + 1, :]
            # 退休期
            else:
                simW[t + 1, :] = (simB[t, :] * rf + simS[t, :] * simR[t, :]) + simP[t+1, :]
    
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
    meanQ = np.mean(simQ, axis=1)
    meanQ_pct = np.mean(simQ_pct, axis=1)
    meanF = np.mean(simF, axis=1)
    meanP = np.mean(simP, axis=1)
    # 保存状态变量路径，保存为excel文件
    os.makedirs('result_baseline_python_PFI/state_variables', exist_ok=True)
    pd.DataFrame(simW).to_excel('result_baseline_python_PFI/state_variables/simW.xlsx')
    pd.DataFrame(simF).to_excel('result_baseline_python_PFI/state_variables/simF.xlsx')
    

    # 保存模拟结果
    os.makedirs('result_baseline_python_PFI', exist_ok=True)
    sim_results = pd.DataFrame({
        'meanC': meanC,
        'meanC_pct': meanC_pct,
        'meanY': meanY,
        'meanW': meanW,
        'meanS': meanS,
        'meanB': meanB,
        'meanWY': meanWY,
        'meanalpha': meanalpha,
        'meanGPY': meanGPY,
        'meanQ': meanQ,
        'meanF': meanF,
        'meanP': meanP
    })
    sim_results.to_excel('result_baseline_python_PFI/simulation_results.xlsx')
    
    # 绘制结果
    # 添加中文支持
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    
    plt.figure(figsize=(15,6))
    
    # 消费比例、风险资产配置和个人养老金购买比例
    plt.subplot(1, 2, 1)
    plt.plot(meanC_pct, label='消费比例')
    plt.plot(meanalpha, label='风险资产配置')
    plt.plot(meanQ_pct, label='养老金购买比例')
    plt.axvline(x=tr-tb, color='k', linestyle='--', label='退休')
    plt.legend()
    plt.xlabel('期数')
    plt.ylabel('比例')
    plt.title('消费、风险资产和养老金购买比例')
    
    # 第二个图：财富和养老金综合展示
    plt.subplot(1, 2, 2)
    plt.plot(meanW, label='现金财富')
    plt.plot(meanF, label='养老基金余额')
    plt.plot(meanW + meanF, label='总财富')
    plt.plot(meanP, label='养老金给付')
    plt.axvline(x=tr-tb, color='k', linestyle='--', label='退休')
    plt.legend()
    plt.xlabel('期数')
    plt.ylabel('金额')
    plt.title('财富和养老金综合展示')
    

    
    plt.tight_layout()
    plt.savefig('result_baseline_python_PFI/results_plot_with_pension.png')
    plt.show()
    
    print("数值模拟完成，结果已保存到result_baseline_python_PFI文件夹")

def main(mode="both"):
    """
    主函数，可选择运行模式：
    "solve": 只进行模型求解
    "simulate": 只进行数值模拟（需要先有求解结果）
    "both": 先求解再模拟（默认）
    """
    # 在开始时重置统计信息字典
    global function_stats
    function_stats = {}
    
    start_time = time.time()
    
    if mode.lower() == "solve" or mode.lower() == "both":
        print("开始求解模型（策略函数迭代方法）...")
        model_results = solve_model()
        print("模型求解完成")
        
        # 打印求解阶段的性能统计信息
        print("\n=== 求解阶段性能统计信息 ===")
        print_function_stats()
        
        # 重置统计信息，以便单独记录模拟阶段的性能
        function_stats = {}
        
        if mode.lower() == "both":
            print("\n开始数值模拟...")
            simulate_model(model_results)
            
            # 打印模拟阶段的性能统计信息
            print("\n=== 模拟阶段性能统计信息 ===")
            print_function_stats()
    elif mode.lower() == "simulate":
        print("开始数值模拟...")
        simulate_model()
        
        # 打印模拟阶段的性能统计信息
        print("\n=== 模拟阶段性能统计信息 ===")
        print_function_stats()
    else:
        print("无效的运行模式，可选值: 'solve', 'simulate', 'both'")
    
    # 打印总运行时间
    total_time = time.time() - start_time
    print(f"\n总运行时间: {total_time:.2f} 秒 ({total_time/60:.2f} 分钟)")

if __name__ == "__main__":
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(description='生命周期模型求解与模拟（策略函数迭代方法）')
    parser.add_argument('--mode', type=str, choices=['solve', 'simulate', 'both'],
        default='both', help='运行模式：求解(solve)，模拟(simulate)，或两者都运行(both)')
    parser.add_argument('--profile', action='store_true', default=False, help='是否进行详细性能分析')
    
    # 解析命令行参数
    args = parser.parse_args()
    
    # 运行主函数
    if args.profile:
        import cProfile
        import pstats
        
        # 使用cProfile进行详细性能分析
        profiler = cProfile.Profile()
        profiler.enable()
        
        main(args.mode)
        
        profiler.disable()
        stats = pstats.Stats(profiler).sort_stats('cumtime')
        stats.print_stats(20)  # 打印前20个耗时最多的函数
        stats.dump_stats('profile_results.prof')  # 保存分析结果
        print("详细性能分析结果已保存到 profile_results.prof")
    else:
        main(args.mode)