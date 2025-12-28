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

# 添加JAX相关库
import jax
import jax.numpy as jnp
from jax import grad, jit, vmap, lax
# 修复导入错误：从jax.example_libraries导入optimizers，而不是从jax.experimental
from jax.example_libraries import optimizers
import jax.scipy as jscipy
from jax.scipy.optimize import minimize as jax_minimize
from jax import debug

# 设置JAX配置
jax.config.update("jax_enable_x64", True)  # 启用64位精度
print(f"JAX设备: {jax.devices()}")  # 打印可用设备

# 创建一个字典来存储每个函数的统计信息
function_stats = {}

# 创建插值器缓存
interpolator_cache = {}

# JAX版本的双线性插值函数
@jit
def bilinear_interpolation_jax(values, x_grid, y_grid, x_points, y_points):
    """
    使用JAX实现的双线性插值函数，增强数值稳定性以避免NaN
    
    参数:
    values - 二维网格上的值，形状为(len(x_grid), len(y_grid))
    x_grid - x方向的网格点，必须严格递增
    y_grid - y方向的网格点，必须严格递增
    x_points - 要插值的x坐标点
    y_points - 要插值的y坐标点
    
    返回:
    插值结果
    """
    # 确保输入是JAX数组
    values = jnp.asarray(values)
    x_grid = jnp.asarray(x_grid)
    y_grid = jnp.asarray(y_grid)
    x_points = jnp.asarray(x_points)
    y_points = jnp.asarray(y_points)
    
    # 确保插值点在网格范围内，避免越界
    x_points = jnp.clip(x_points, x_grid[0], x_grid[-1])
    y_points = jnp.clip(y_points, y_grid[0], y_grid[-1])
    
    # 找到x_grid中小于等于x_points的最大索引
    x_indices = jnp.searchsorted(x_grid, x_points, side='right') - 1
    # 确保索引在有效范围内
    x_indices = jnp.clip(x_indices, 0, len(x_grid) - 2)
    
    # 找到y_grid中小于等于y_points的最大索引
    y_indices = jnp.searchsorted(y_grid, y_points, side='right') - 1
    # 确保索引在有效范围内
    y_indices = jnp.clip(y_indices, 0, len(y_grid) - 2)
    
    # 获取相邻网格点
    x0 = x_grid[x_indices]
    x1 = x_grid[x_indices + 1]
    y0 = y_grid[y_indices]
    y1 = y_grid[y_indices + 1]
    
    # 计算插值权重，添加数值稳定性保护
    # 确保分母不为零
    x_denom = jnp.maximum(x1 - x0, 1e-10)
    y_denom = jnp.maximum(y1 - y0, 1e-10)
    
    x_weights = (x_points - x0) / x_denom
    y_weights = (y_points - y0) / y_denom
    
    # 确保权重在[0,1]范围内
    x_weights = jnp.clip(x_weights, 0.0, 1.0)
    y_weights = jnp.clip(y_weights, 0.0, 1.0)
    
    # 获取四个角点的值
    v00 = values[x_indices, y_indices]
    v01 = values[x_indices, y_indices + 1]
    v10 = values[x_indices + 1, y_indices]
    v11 = values[x_indices + 1, y_indices + 1]
    
    # 执行双线性插值
    v0 = v00 * (1 - x_weights) + v10 * x_weights  # 底边插值
    v1 = v01 * (1 - x_weights) + v11 * x_weights  # 顶边插值
    v = v0 * (1 - y_weights) + v1 * y_weights     # 最终插值
    
    # 检查结果是否为NaN，如果是则使用最近邻插值作为后备
    is_nan = jnp.isnan(v)
    
    # 计算最近邻插值作为后备
    nearest_x_idx = jnp.round(x_indices + x_weights).astype(jnp.int32)
    nearest_y_idx = jnp.round(y_indices + y_weights).astype(jnp.int32)
    nearest_x_idx = jnp.clip(nearest_x_idx, 0, len(x_grid) - 1)
    nearest_y_idx = jnp.clip(nearest_y_idx, 0, len(y_grid) - 1)
    nearest_v = values[nearest_x_idx, nearest_y_idx]
    
    # 如果双线性插值结果为NaN，则使用最近邻插值结果
    result = jnp.where(is_nan, nearest_v, v)
    
    return result

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

# 退休期间的值函数计算 - JAX版本


@jit
def cocco_fun_valuefunc_retired_jax(x, cash, fund, nextV, gret, rf, ret_fac, pension_pay, gamma, beta, psi_1, psi_2, theta, gcash, survprob, weig):
    """
    JAX版本的退休期间值函数计算
    
    参数:
    x - 决策变量 [消费, 风险资产投资比例]
    cash - 手中现金
    fund - 养老基金账户余额
    nextV - 下期值函数
    gret - 随机收益率
    rf - 无风险收益率
    ret_fac - 固定基本养老金
    pension_pay - 个人养老金给付
    gamma - 相对风险规避系数
    beta - 贴现因子
    psi_1, psi_2, theta - Epstein-Zin效用函数参数
    gcash - 现金网格
    survprob - 存活概率
    weig - 随机状态权重
    
    返回:
    当前状态下的值函数
    """
    # 提取决策变量
    c = x[0]  # 消费率
    alpha = x[1]  # 风险资产配置比例
    
    # 计算当期消费
    cons = jnp.maximum(c * cash, 1e-10)  # 确保消费为正
    
    # 计算当期效用 - 使用lax.cond替代if语句
    util = lax.cond(
        jnp.equal(gamma, 1.0),
        lambda _: jnp.log(cons),
        lambda _: cons**(1-gamma) / (1-gamma),
        operand=None
    )
    
    # 计算剩余资金
    sav = cash * (1-c)
    
    # 计算风险资产和无风险资产投资
    risky = alpha * sav
    safe = sav - risky
    
    # 计算下期财富
    next_cash = safe * rf + risky * gret + ret_fac + pension_pay
    
    # 计算下期期望效用
    next_util = jnp.zeros_like(gret)
    
    # 使用JAX的双线性插值计算下期期望效用
    for i in range(len(gret)):
        # 使用修复后的bilinear_interpolation_jax函数进行插值
        # 确保传递正确的参数格式
        next_util = next_util.at[i].set(
            bilinear_interpolation_jax(
                nextV,  # 值函数网格
                gcash,  # x网格
                jnp.array([fund]),  # y网格，确保是数组
                next_cash[i],  # x点
                fund  # y点
            )
        )
    
    # 计算期望效用
    exp_util = jnp.sum(weig * next_util)
    
    # 计算Epstein-Zin效用 - 使用lax.cond替代if语句
    value = lax.cond(
        jnp.equal(theta, 1.0),
        lambda _: util + beta * survprob * exp_util,
        lambda _: ((1-beta) * util**(1-psi_1) + beta * survprob * exp_util**(1-psi_2))**(1/(1-theta)),
        operand=None
    )
    
    # 返回负值函数（用于最小化）
    return -value

# 保留原始NumPy版本以便兼容
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
        u = np.log(cash*x[0])  # 对数效用
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

# 工作期间的值函数计算 - JAX版本


@jit
def cocco_fun_valuefunc_work_jax(x, cash, fund, gyp, nextV, yh, gret, rf, rp, tau, gamma, beta, psi_1, psi_2, theta, gcash, gfund, survprob, n, weig):
    """
    JAX版本的工作期间值函数计算
    
    参数:
    x - 决策变量 [消费率, 风险资产配置比例, 养老金购买比例]
    cash - 手中现金
    fund - 养老基金账户余额
    gyp - 收入过程
    nextV - 下期值函数
    yh - 收入冲击
    gret - 随机收益率
    rf - 无风险收益率
    rp - 个人养老基金收益率
    tau - 工资扣除缴纳基本养老金的比例
    gamma - 相对风险规避系数
    beta - 贴现因子
    psi_1, psi_2, theta - Epstein-Zin效用函数参数
    gcash - 现金网格
    gfund - 养老基金网格
    survprob - 存活概率
    n - 外生随机冲击的grid数量
    weig - 随机状态权重
    
    返回:
    当前状态下的值函数
    """
    # 提取决策变量
    c = x[0]  # 消费率
    alpha = x[1]  # 风险资产配置比例
    q = x[2]  # 养老金购买比例
    
    # 计算当期消费
    cons = jnp.maximum(c * cash, 1e-10)  # 确保消费为正
    
    # 计算当期效用 - 使用lax.cond替代if语句
    util = lax.cond(
        jnp.equal(gamma, 1.0),
        lambda _: jnp.log(cons),
        lambda _: cons**(1-gamma) / (1-gamma),
        operand=None
    )
    
    # 计算剩余资金
    sav = cash * (1-c)
    
    # 计算养老金购买
    pension_contrib = q * sav
    
    # 计算风险资产和无风险资产投资
    risky = alpha * (sav - pension_contrib)
    safe = (sav - pension_contrib) - risky
    
    # 初始化下期期望效用
    exp_util = 0.0
    
    # 计算下期期望效用
    for i1 in range(n):
        for i2 in range(n):
            # 计算下期财富
            next_cash = (safe * rf + risky * gret[i2]) / gyp[i1, i2] + (1 - tau) * jnp.exp(yh[i1])
            next_fund = (fund + pension_contrib) * rp + tau * jnp.exp(yh[i1])
            
            # 使用修复后的JAX双线性插值计算下期值函数
            next_value = bilinear_interpolation_jax(
                nextV,  # 值函数网格
                gcash,  # x网格
                gfund,  # y网格
                next_cash,  # x点
                next_fund  # y点
            )
            
            # 累加期望效用
            exp_util = exp_util + weig[i1] * weig[i2] * next_value
    
    # 计算Epstein-Zin效用 - 使用lax.cond替代if语句
    value = lax.cond(
        jnp.equal(theta, 1.0),
        lambda _: util + beta * survprob * exp_util,
        lambda _: ((1-beta) * util**(1-psi_1) + beta * survprob * exp_util**(1-psi_2))**(1/(1-theta)),
        operand=None
    )
    
    # 返回负值函数（用于最小化）
    return -value

# 保留原始NumPy版本以便兼容
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
def optimize_policy_work_jax(i2, i3, x0, gcash, gfund, gyp, V_next, yh, gret, rf, rp, 
                           pension_pct, gamma, beta, psi_1, psi_2, theta, 
                           survprob, n, nweig1):
    """
    使用JAX优化工作期策略
    """
    # 获取当前状态
    cash = gcash[i3]
    fund = gfund[i2]
    
    # 设置优化边界
    lb = [0.0, 0.0, 0.0]  # 下界：消费率、风险资产配置比例、养老金购买比例
    ub = [1.0, 1.0, 1.0]  # 上界：消费率、风险资产配置比例、养老金购买比例
    
    # 确保初始值在边界内
    x0 = np.clip(x0, lb, ub)
    
    # 转换为JAX数组
    x0_jax = jnp.array(x0)
    cash_jax = jnp.array(cash)
    fund_jax = jnp.array(fund)
    gyp_jax = jnp.array(gyp)
    V_next_jax = jnp.array(V_next)
    yh_jax = jnp.array(yh)
    gret_jax = jnp.array(gret)
    survprob_jax = jnp.array(survprob)
    nweig1_jax = jnp.array(nweig1)
    
    # 定义目标函数
    def objective(x):
        return cocco_fun_valuefunc_work_jax(
            x, cash_jax, fund_jax, gyp_jax, V_next_jax, yh_jax, gret_jax,
            rf, rp, pension_pct, gamma, beta, psi_1, psi_2, theta,
            jnp.array(gcash), jnp.array(gfund), survprob_jax, n, nweig1_jax
        )
    
    try:
        # 使用JAX的优化器，使用BFGS方法（JAX支持的方法）
        result = jax_minimize(
            objective, x0_jax, 
            method="BFGS",  # JAX支持的方法
            options={'gtol': 1e-7, 'maxiter': 200}
        )
        
        # 获取优化结果
        policy = np.array(result.x)
        
        # 确保结果在边界内
        policy = np.clip(policy, lb, ub)
        
    except Exception as e:
        print(f"JAX优化失败: {e}，回退到SciPy优化")
        # 回退到SciPy优化
        result = minimize(
            lambda x: float(objective(x)), 
            x0, 
            method='L-BFGS-B',
            bounds=list(zip(lb, ub)),
            options={'ftol': 1e-7, 'maxiter': 200}
        )
        policy = result.x
    
    return i2, i3, policy

# 保留原始NumPy版本以便兼容
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

def optimize_policy_retired_jax(i2, i3, x0, cash, fund, V_next, gret, rf, ret_fac, pension_rate,
                             gamma, beta, psi_1, psi_2, theta, gcash, survprob, nweig2):
    """
    使用JAX优化退休期策略
    """
    # 设置优化边界
    lb = [0.0, 0.0]  # 下界：消费率、风险资产配置比例
    ub = [1.0, 1.0]  # 上界：消费率、风险资产配置比例
    
    # 确保初始值在边界内
    x0 = np.clip(x0, lb, ub)
    
    # 转换为JAX数组
    x0_jax = jnp.array(x0)
    cash_jax = jnp.array(cash)
    fund_jax = jnp.array(fund)
    V_next_jax = jnp.array(V_next)
    gret_jax = jnp.array(gret)
    survprob_jax = jnp.array(survprob)
    nweig2_jax = jnp.array(nweig2)
    gcash_jax = jnp.array(gcash)
    # 定义目标函数
    def objective(x):
        return cocco_fun_valuefunc_retired_jax(
            x, cash_jax, fund_jax, V_next_jax, gret_jax, rf, ret_fac, pension_rate,
            gamma, beta, psi_1, psi_2, theta, gcash_jax, survprob_jax, nweig2_jax  # 确保传递nweig2_jax作为weig参数
        )
    
    try:
        # 使用JAX的优化器，使用BFGS方法（JAX支持的方法）
        result = jax_minimize(
            objective, x0_jax, 
            method="BFGS",  # JAX支持的方法
            options={'gtol': 1e-7, 'maxiter': 200}
        )
        
        # 获取优化结果
        policy = np.array(result.x)
        
        # 确保结果在边界内
        policy = np.clip(policy, lb, ub)
        
    except Exception as e:
        print(f"JAX优化失败: {e}，回退到SciPy优化")
        # 回退到SciPy优化
        result = minimize(
            lambda x: float(objective(x)), 
            x0, 
            method='L-BFGS-B',
            bounds=list(zip(lb, ub)),
            options={'ftol': 1e-7, 'maxiter': 200}
        )
        policy = result.x
    
    return i2, i3, policy

# 在solve_model函数中添加JAX版本的退休期策略迭代
def solve_retire_period_jax(t, C, A, V, gcash, gfund, gret, rf, ret_fac, pension_rate, 
                           gamma, beta, psi_1, psi_2, theta, 
                           survprob, weig, tolerance, max_iterations, num_cores):
    """
    使用JAX加速的退休期策略迭代
    """
    
    # 将NumPy数组转换为JAX数组
    gcash_jax = jnp.array(gcash)
    gfund_jax = jnp.array(gfund)
    V_next_jax = jnp.array(V[:, :, t+1])
    gret_jax = jnp.array(gret)
    survprob_jax = jnp.array(survprob[t])
    weig_jax = jnp.array(weig)
    
    # 初始化新的策略和值函数
    C_new = np.array(C[:, :, t])
    A_new = np.array(A[:, :, t])
    V_new = np.zeros_like(C_new)
    
    # 策略迭代循环
    for iter_count in range(max_iterations):
        # 策略评估 - 计算当前策略的值函数
        for i2 in range(len(gfund)):
            for i3 in range(len(gcash)):
                policy = np.array([C_new[i3, i2], A_new[i3, i2]])
                
                # 使用JAX版本的值函数计算
                V_new[i3, i2] = -cocco_fun_valuefunc_retired_jax(
                    policy, gcash_jax[i3], gfund_jax[i2],
                    V_next_jax, gret_jax, rf, ret_fac, pension_rate,
                    gamma, beta, psi_1, psi_2, theta,
                    gcash_jax, survprob_jax, weig_jax  # 确保传递weig_jax作为weig参数
                )
                raise
        
        # 检查并缩放值函数，如果范围太小
        V_new, scaling_factor, was_scaled = scale_value_function(V_new)
        if was_scaled:
            print(f"  退休期时期 {t}: 值函数范围过小，应用了缩放因子 {scaling_factor:.4f}")
        
        # 策略改进 - 基于当前值函数更新策略，使用并行计算
        C_opt = np.zeros_like(C_new)
        A_opt = np.zeros_like(A_new)
        
        # 准备并行计算的任务
        tasks = []
        for i2 in range(len(gfund)):
            for i3 in range(len(gcash)):
                x0 = [C_new[i3, i2], A_new[i3, i2]]
                tasks.append((i2, i3, x0))
        
        # 并行执行优化 - 使用JAX版本的优化函数
        results = Parallel(n_jobs=num_cores)(
            delayed(optimize_policy_retired_jax)(
                i2, i3, x0, gcash[i3], gfund[i2], 
                V[:, :, t+1], gret, rf, ret_fac, pension_rate, 
                gamma, beta, psi_1, psi_2, theta, gcash,
                survprob[t], weig  # 确保传递weig作为nweig2参数
            ) for i2, i3, x0 in tasks
        )
        
        # 处理并行计算结果
        for i2, i3, policy in results:
            C_opt[i3, i2] = policy[0]
            A_opt[i3, i2] = policy[1]
        
        # 计算策略变化量
        policy_change = (np.max(np.abs(C_opt - C_new)) + 
                        np.max(np.abs(A_opt - A_new)))
        
        # 自适应步长策略
        alpha = 0.1  # 接近收敛阶段使用较小步长
        
        # 更新策略
        C_new = (1-alpha) * C_new + alpha * C_opt
        A_new = (1-alpha) * A_new + alpha * A_opt
        
        # 检查收敛
        if policy_change < tolerance:
            print(f"  退休期时期 {t} 在 {iter_count+1} 次迭代后收敛")
            break
        
        # 如果达到最大迭代次数仍未收敛
        if iter_count == max_iterations - 1:
            print(f"  警告：退休期时期 {t} 在 {max_iterations} 次迭代后仍未收敛")
    
    # 更新策略和值函数
    C[:, :, t] = C_new
    A[:, :, t] = A_new
    V[:, :, t] = V_new
    
    return C, A, V

# 在solve_model函数中添加JAX版本的工作期策略迭代
def solve_work_period_jax(t, C, A, Q, V, gcash, gfund, gyp, yh, gret, rf, rp, 
                         pension_pct, gamma, beta, psi_1, psi_2, theta, 
                         survprob, n, nweig1, tolerance, max_iterations, num_cores):
    """使用JAX加速的工作期策略迭代"""
    
    # 将NumPy数组转换为JAX数组
    gcash_jax = jnp.array(gcash)
    gfund_jax = jnp.array(gfund)
    V_next_jax = jnp.array(V[:, :, t+1])
    gyp_jax = jnp.array(gyp[:, :, t])
    yh_jax = jnp.array(yh)
    gret_jax = jnp.array(gret)
    survprob_jax = jnp.array(survprob[t])
    nweig1_jax = jnp.array(nweig1)
    
    # 初始化新的策略和值函数
    C_new = np.array(C[:, :, t])
    A_new = np.array(A[:, :, t])
    Q_new = np.array(Q[:, :, t])
    V_new = np.zeros_like(C_new)
        
        # 策略迭代循环
    for iter_count in range(max_iterations):
            # 策略评估 - 计算当前策略的值函数
        for i2 in range(len(gfund)):
            for i3 in range(len(gcash)):
                policy = np.array([C_new[i3, i2], A_new[i3, i2], Q_new[i3, i2]])
                
                # 使用JAX版本的值函数计算
                V_new[i3, i2] = -cocco_fun_valuefunc_work_jax(
                    policy, gcash_jax[i3], gfund_jax[i2],
                    gyp_jax, V_next_jax, yh_jax, gret_jax,
                    rf, rp, pension_pct, gamma, beta, psi_1, psi_2, theta,
                    gcash_jax, gfund_jax, survprob_jax, n, nweig1_jax
                    )
            
            # 检查并缩放值函数，如果范围太小
            V_new, scaling_factor, was_scaled = scale_value_function(V_new)
            if was_scaled:
                print(f"  工作期时期 {t}: 值函数范围过小，应用了缩放因子 {scaling_factor:.4f}")
                
                # 对未来值函数进行状态依赖的缩放，使得值函数差异体现不同状态的重要性
                if t < 20:  # 早期工作期特别处理
                    # 保存原始值函数的拷贝，用于打印诊断信息
                    V_before = V_new.copy()
                    
                    # 应用状态依赖的缩放因子
                for i3 in range(len(gcash)):
                    cash_position = (gcash[i3] - gcash[0]) / (gcash[-1] - gcash[0])
                    for i2 in range(len(gfund)):
                            fund_position = (gfund[i2] - gfund[0]) / (gfund[-1] - gfund[0])
                            # 状态依赖的缩放：现金和养老金越多，值函数差异越大
                            state_factor = 1.0 + 5.0 * (cash_position + fund_position)
                            V_new[i3, i2] *= state_factor
                    
                    # 打印诊断信息
                    v_range_before = np.max(V_before) - np.min(V_before)
                    v_range_after = np.max(V_new) - np.min(V_new)
                print(f"    状态依赖缩放: 范围从 {v_range_before:.6f} 增加到 {v_range_after:.6f}")
            
            # 策略改进 - 基于当前值函数更新策略，使用并行计算
        C_opt = np.zeros_like(C_new)
        A_opt = np.zeros_like(A_new)
        Q_opt = np.zeros_like(Q_new)
            
            # 准备并行计算的任务
        tasks = []
        for i2 in range(len(gfund)):
            for i3 in range(len(gcash)):
                x0 = [C_new[i3, i2], A_new[i3, i2], Q_new[i3, i2]]
                tasks.append((i2, i3, x0))
            
        # 并行执行优化 - 使用JAX版本的优化函数
            results = Parallel(n_jobs=num_cores)(
            delayed(optimize_policy_work_jax)(
                    i2, i3, x0, gcash, gfund, 
                gyp[:, :, t], V[:, :, t+1], 
                    yh, gret, rf, rp, pension_pct, 
                    gamma, beta, psi_1, psi_2, theta, 
                    survprob[t], n, nweig1
                ) for i2, i3, x0 in tasks
            )
            
            # 处理并行计算结果
            for i2, i3, policy in results:
                C_opt[i3, i2] = policy[0]
                A_opt[i3, i2] = policy[1]
                Q_opt[i3, i2] = policy[2]
            
            # 计算策略变化量
        policy_change = (np.max(np.abs(C_opt - C_new)) + 
                        np.max(np.abs(A_opt - A_new)) + 
                       np.max(np.abs(Q_opt - Q_new)))
            
        # 自适应步长策略
        alpha = 0.1  # 接近收敛阶段使用较小步长
            
            # 更新策略
        C_new = (1-alpha) * C_new + alpha * C_opt
        A_new = (1-alpha) * A_new + alpha * A_opt
        Q_new = (1-alpha) * Q_new + alpha * Q_opt
            
        # 检查收敛
        if policy_change < tolerance:
            print(f"  工作期时期 {t} 在 {iter_count+1} 次迭代后收敛")
            break
        
        # 如果达到最大迭代次数仍未收敛
        if iter_count == max_iterations - 1:
            print(f"  警告：工作期时期 {t} 在 {max_iterations} 次迭代后仍未收敛")
            
    # 更新策略和值函数
    C[:, :, t] = C_new
    A[:, :, t] = A_new
    Q[:, :, t] = Q_new
    V[:, :, t] = V_new
    
    return C, A, Q, V
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
    使用JAX加速的生命周期模型求解函数
    """
    print("开始求解生命周期模型（使用JAX加速）...")
    start_time = time.time()
    
    # 设置参数
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
    
    # 效用函数参数
    gamma = 5.0
    beta = 0.96
    psi_1 = 0.5
    psi_2 = 0.5
    theta = 0.5
    
    # 生存概率
    survprob = np.zeros(tn - 1)
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
    
    # 设置状态空间
    # 现金网格
    gcash_min = 0.0
    gcash_max = 20.0
    gcash_n = 30
    gcash = np.linspace(gcash_min, gcash_max, gcash_n)
    
    # 养老基金网格
    gfund_min = 0.0
    gfund_max = 10.0
    gfund_n = 20
    gfund = np.linspace(gfund_min, gfund_max, gfund_n)
    
    # 初始化策略和值函数
    C = np.zeros((gcash_n, gfund_n, tn))  # 消费
    A = np.zeros((gcash_n, gfund_n, tn))  # 风险资产配置
    Q = np.zeros((gcash_n, gfund_n, tn))  # 个人养老金购买
    V = np.zeros((gcash_n, gfund_n, tn))  # 值函数
    
    # 设置最后一期的策略和值函数
    C[:, :, -1] = 1.0  # 最后一期全部消费
    A[:, :, -1] = 0.0  # 最后一期不投资
    Q[:, :, -1] = 0.0  # 最后一期不购买个人养老金
    
    # 计算最后一期的值函数
    for i3 in range(gcash_n):
        for i2 in range(gfund_n):
            cash = gcash[i3]
            fund = gfund[i2]
            # 最后一期的值函数
            V[i3, i2, -1] = cash**psi_1
    
    # 设置随机冲击
    # 收入冲击
    yh = np.zeros(n)
    nweig1 = np.zeros(n)
    
    # 使用Tauchen方法离散化AR(1)过程
    Tauchen_q = 3
    tauchenoptions = {'parallel': False}
    yh, nweig1 = discretizeAR1_Tauchen(0, 0, smay, n, Tauchen_q, tauchenoptions)
    
    # 风险资产收益率
    gret = np.zeros(n)
    nweig2 = np.zeros(n)
    gret, nweig2 = discretizeAR1_Tauchen(mu, 0, sigr, n, Tauchen_q, tauchenoptions)
    gret = gret + rf
    
    # 计算收入过程
    gyp = np.zeros((n, n, tr - tb))
    for i1 in range(tr - tb):
        for i2 in range(n):
            for i3 in range(n):
                gyp[i2, i3, i1] = np.exp(gy[i1] + yh[i2] + yh[i3])
    
    # 设置并行计算的核心数
    num_cores = multiprocessing.cpu_count() - 1
    print(f"使用 {num_cores} 个CPU核心进行并行计算")
    
    # 设置求解参数
    tolerance = 1e-4
    max_iterations = 20
    
    # 从退休期开始向前求解
    print("开始求解退休期...")
    for t in range(tn-2, tr-tb-1, -1):
        print(f"求解退休期时期 {t}...")
        C, A, V = solve_retire_period_jax(
            t, C, A, V, gcash, gfund, gret, rf, ret_fac, pension_rate,
            gamma, beta, psi_1, psi_2, theta, survprob, nweig2,
            tolerance, max_iterations, num_cores
        )
    
    # 求解工作期
    print("开始求解工作期...")
    for t in range(tr-tb-1, -1, -1):
        print(f"求解工作期时期 {t}...")
        C, A, Q, V = solve_work_period_jax(
            t, C, A, Q, V, gcash, gfund, gyp, yh, gret, rf, rp,
            pension_pct, gamma, beta, psi_1, psi_2, theta,
            survprob, n, nweig1, tolerance, max_iterations, num_cores
        )
    
    # 保存结果
    os.makedirs('result_baseline_jax_PFI', exist_ok=True)
    np.save('result_baseline_jax_PFI/C.npy', C)
    np.save('result_baseline_jax_PFI/A.npy', A)
    np.save('result_baseline_jax_PFI/Q.npy', Q)
    np.save('result_baseline_jax_PFI/V.npy', V)
    np.save('result_baseline_jax_PFI/gcash.npy', gcash)
    np.save('result_baseline_jax_PFI/gfund.npy', gfund)
    
    # 打印函数调用统计信息
    print_function_stats()
    
    end_time = time.time()
    print(f"模型求解完成，总耗时: {end_time - start_time:.2f} 秒")
    
    # 返回模型结果
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
            if not os.path.exists('result_baseline_jax_PFI/C.npy'):
                print("未找到模型求解结果文件，请先运行模型求解")
                return
                
            # 从文件中读取模型结果
            C = np.load('result_baseline_jax_PFI/C.npy')
            A = np.load('result_baseline_jax_PFI/A.npy')
            Q = np.load('result_baseline_jax_PFI/Q.npy')
            V = np.load('result_baseline_jax_PFI/V.npy')
            gcash = np.load('result_baseline_jax_PFI/gcash.npy')
            gfund = np.load('result_baseline_jax_PFI/gfund.npy')
            
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
    os.makedirs('result_baseline_jax_PFI/state_variables', exist_ok=True)
    pd.DataFrame(simW).to_excel('result_baseline_jax_PFI/state_variables/simW.xlsx')
    pd.DataFrame(simF).to_excel('result_baseline_jax_PFI/state_variables/simF.xlsx')
    

    # 保存模拟结果
    os.makedirs('result_baseline_jax_PFI', exist_ok=True)
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
    sim_results.to_excel('result_baseline_jax_PFI/simulation_results.xlsx')
    
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
    plt.savefig('result_baseline_jax_PFI/results_plot_with_pension.png')
    plt.show()
    
    print("数值模拟完成，结果已保存到result_baseline_jax_PFI文件夹")

def main(mode="both"):
    """
    主函数，根据模式执行模型求解和/或模拟
    
    参数：
        mode: 执行模式，可选值为"solve"（仅求解）、"simulate"（仅模拟）或"both"（求解和模拟）
    """
    # 解析命令行参数
    parser = argparse.ArgumentParser(description='生命周期模型求解与模拟（JAX加速版）')
    parser.add_argument('--mode', type=str, default=mode, choices=['solve', 'simulate', 'both'],
                        help='执行模式：solve（仅求解）、simulate（仅模拟）或both（求解和模拟）')
    args = parser.parse_args()
    
    # 根据模式执行相应操作
    if args.mode in ['solve', 'both']:
        # 求解模型
        model_results = solve_model()
        
        # 如果只求解不模拟，则直接返回
        if args.mode == 'solve':
            return model_results
    else:
        # 如果只模拟不求解，则model_results为None，将从文件中读取
        model_results = None
    
    # 模拟模型
    sim_results = simulate_model_jax(model_results)
    
    return sim_results


if __name__ == "__main__":
    main()



def simulate_model_jax(model_results=None):
    """
    基于求解结果进行数值模拟，使用JAX加速
    
    参数：
        model_results: 模型求解结果。如果为None，会尝试从文件中读取
    """
    # 如果没有提供模型结果，尝试从文件中读取
    if model_results is None:
        try:
            # 检查文件是否存在
            if not os.path.exists('result_baseline_jax_PFI/C.npy'):
                print("未找到模型求解结果文件，请先运行模型求解")
                return
                
            # 从文件中读取模型结果
            C = np.load('result_baseline_jax_PFI/C.npy')
            A = np.load('result_baseline_jax_PFI/A.npy')
            Q = np.load('result_baseline_jax_PFI/Q.npy')
            V = np.load('result_baseline_jax_PFI/V.npy')
            gcash = np.load('result_baseline_jax_PFI/gcash.npy')
            gfund = np.load('result_baseline_jax_PFI/gfund.npy')
            
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
    C = jnp.array(model_results['C'])
    A = jnp.array(model_results['A'])
    Q = jnp.array(model_results['Q'])
    V = jnp.array(model_results['V'])
    gcash = jnp.array(model_results['gcash'])
    gfund = jnp.array(model_results['gfund'])
    
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
    gy = jnp.array(params['gy'])
    
    print("开始数值模拟（使用JAX加速）...")
    
    # 设置随机种子以获得可重复的结果
    key = jax.random.PRNGKey(42)
    
    # 数值模拟参数
    nsim = 10000
    
    # 初始化模拟数组
    simPY = jnp.zeros((tn, nsim))
    simGPY = jnp.zeros((tn, nsim))
    simY = jnp.zeros((tn, nsim))
    simC = jnp.zeros((tn, nsim))
    simW = jnp.zeros((tn, nsim))
    simS = jnp.zeros((tn, nsim))
    simA = jnp.zeros((tn, nsim))
    simB = jnp.zeros((tn, nsim))
    simW_Y = jnp.zeros((tn, nsim))
    simR = jnp.zeros((tn, nsim))
    simQ = jnp.zeros((tn, nsim))  # 个人养老金购买决策
    simF = jnp.zeros((tn, nsim))  # 养老基金账户余额
    simP = jnp.zeros((tn, nsim))  # 个人养老金给付
    cash = jnp.zeros((tn, nsim))
    simC_pct = jnp.zeros((tn, nsim))
    
    # 1、模拟生成labor income
    # 使用JAX的随机数生成
    for i1 in range(nsim // 2): # 另外一半模拟完全对称
        # 生成随机数
        key, subkey = jax.random.split(key)
        eps_y = jax.random.normal(subkey)
        
        # working period第一期
        simPY = simPY.at[0, i1].set(eps_y * smay)
        simPY = simPY.at[0, nsim // 2 + i1].set(-eps_y * smay)
        simGPY = simGPY.at[0, i1].set(1.0)
        simGPY = simGPY.at[0, nsim // 2 + i1].set(1.0)
        
        key, subkey = jax.random.split(key)
        simTY = jax.random.normal(subkey)
        simY = simY.at[0, i1].set(jnp.exp(simTY * smay))
        simY = simY.at[0, nsim // 2 + i1].set(jnp.exp(-simTY * smay))
        
        # working period第2期~退休
        for i2 in range(1, tr - tb):
            w = i2 + tb - 1
            key, subkey = jax.random.split(key)
            eps_y = jax.random.normal(subkey)
            
            simPY = simPY.at[i2, i1].set(eps_y * smav + simPY[i2 - 1, i1])
            simPY = simPY.at[i2, nsim // 2 + i1].set(-eps_y * smav + simPY[i2 - 1, nsim // 2 + i1])
            
            simGPY = simGPY.at[i2, i1].set(jnp.exp(gy[i2 - 1]) * jnp.exp(simPY[i2, i1]) / jnp.exp(simPY[i2 - 1, i1]))
            simGPY = simGPY.at[i2, nsim // 2 + i1].set(jnp.exp(gy[i2 - 1]) * jnp.exp(simPY[i2, nsim // 2 + i1]) / jnp.exp(simPY[i2 - 1, nsim // 2 + i1]))
            
            key, subkey = jax.random.split(key)
            simTY = jax.random.normal(subkey)
            simY = simY.at[i2, i1].set(jnp.exp(simTY * smay))
            simY = simY.at[i2, nsim // 2 + i1].set(jnp.exp(-simTY * smay))
    
    # 退休期
    for t in range(tr - tb, tn):
        simY = simY.at[t, :].set(ret_fac)
        simGPY = simGPY.at[t, :].set(1.0)
    
    # 2、模拟风险投资的收益率
    for t in range(tn):
        for i1 in range(nsim // 2):
            key, subkey = jax.random.split(key)
            eps_r = jax.random.normal(subkey)
            simR = simR.at[t, i1].set(mu + rf + sigr * eps_r)
            simR = simR.at[t, nsim // 2 + i1].set(mu + rf - sigr * eps_r)
    
    # 从第一期开始迭代，得到各控制变量的值
    simW = simW.at[:, :].set(0.2)
    simF = simF.at[:, :].set(0)
    simF = simF.at[0, :].set(pension_pct * simY[0, :])
    
    # 预编译JAX函数以提高性能
    @jit
    def interpolate_policy(C, A, Q, gcash, gfund, cash_value, fund_value, t, is_working_period):
        """使用JAX的双线性插值获取策略"""
        # 确保输入在网格范围内
        cash_value = jnp.clip(cash_value, gcash[0], gcash[-1])
        fund_value = jnp.clip(fund_value, gfund[0], gfund[-1])
        
        # 找到在gcash中的位置
        cash_idx = jnp.searchsorted(gcash, cash_value, side='right') - 1
        cash_idx = jnp.clip(cash_idx, 0, len(gcash)-2)
        cash_weight = (cash_value - gcash[cash_idx]) / (gcash[cash_idx+1] - gcash[cash_idx])
        
        # 找到在gfund中的位置
        fund_idx = jnp.searchsorted(gfund, fund_value, side='right') - 1
        fund_idx = jnp.clip(fund_idx, 0, len(gfund)-2)
        fund_weight = (fund_value - gfund[fund_idx]) / (gfund[fund_idx+1] - gfund[fund_idx])
        
        # 双线性插值计算C和A
        v00_c = C[cash_idx, fund_idx, t]
        v01_c = C[cash_idx, fund_idx+1, t]
        v10_c = C[cash_idx+1, fund_idx, t]
        v11_c = C[cash_idx+1, fund_idx+1, t]
        
        v00_a = A[cash_idx, fund_idx, t]
        v01_a = A[cash_idx, fund_idx+1, t]
        v10_a = A[cash_idx+1, fund_idx, t]
        v11_a = A[cash_idx+1, fund_idx+1, t]
        
        c_interp = (1-cash_weight)*(1-fund_weight)*v00_c + cash_weight*(1-fund_weight)*v10_c + \
                  (1-cash_weight)*fund_weight*v01_c + cash_weight*fund_weight*v11_c
        
        a_interp = (1-cash_weight)*(1-fund_weight)*v00_a + cash_weight*(1-fund_weight)*v10_a + \
                  (1-cash_weight)*fund_weight*v01_a + cash_weight*fund_weight*v11_a
        
        # 如果是工作期，还需要计算Q
        q_interp = jnp.where(
            is_working_period,
            (1-cash_weight)*(1-fund_weight)*Q[cash_idx, fund_idx, t] + 
            cash_weight*(1-fund_weight)*Q[cash_idx+1, fund_idx, t] + 
            (1-cash_weight)*fund_weight*Q[cash_idx, fund_idx+1, t] + 
            cash_weight*fund_weight*Q[cash_idx+1, fund_idx+1, t],
            0.0
        )
        
        # 确保约束条件满足
        c_interp = jnp.clip(c_interp, 0.0, 1.0)
        a_interp = jnp.clip(a_interp, 0.0, 1.0)
        q_interp = jnp.clip(q_interp, 0.0, 1.0)
        
        return c_interp, a_interp, q_interp
    
    # 使用JAX的vmap函数向量化模拟过程
    @jit
    def simulate_single_period(t, simW, simF, C, A, Q, gcash, gfund, simY, simR, rf, rp, pension_pct, pension_rate, tr, tb):
        """模拟单个时期的状态转移"""
        # 判断是否为工作期
        is_working_period = t < tr - tb
        
        # 计算当期现金
        cash_t = simW[t, :]
        
        # 向量化插值计算
        c_interp, a_interp, q_interp = jax.vmap(
            lambda cash, fund: interpolate_policy(
                C, A, Q, gcash, gfund, cash, fund, t, is_working_period
            )
        )(cash_t, simF[t, :])
        
        # 记录策略
        simC_t = c_interp
        simA_t = a_interp
        simQ_t = q_interp
        
        # 计算个人养老金购买金额和养老基金账户更新
        if is_working_period:  # 工作期
            # 个人养老金购买
            pension_contrib = simQ_t * cash_t * (1-simC_t)
            
            # 计算养老基金账户更新
            simF_next = (simF[t, :] + pension_contrib) * rp + pension_pct * simY[t+1, :]
            
            # 个人养老金给付为0
            simP_next = jnp.zeros_like(simF_next)
        else:  # 退休期
            # 计算个人养老金给付
            simP_next = pension_rate * simF[tr-tb-1, :]
            
            # 养老基金余额保持不变
            simF_next = simF[t, :]
        
        # 剩余用于投资的资金
        sav = cash_t * (1-simC_t)
        
        # 工作期扣除个人养老金购买金额
        if is_working_period:
            sav = (1-simQ_t) * cash_t * (1-simC_t)
        
        # 风险投资额
        simS_t = simA_t * sav
        
        # 无风险投资额
        simB_t = sav - simS_t
        
        # 计算下期财富
        if is_working_period and t < tr - tb - 1:  # 工作期
            simW_next = (simB_t * rf + simS_t * simR[t, :]) / simGPY[t+1, :] + (1 - pension_pct) * simY[t+1, :]
        else:  # 退休期
            simW_next = (simB_t * rf + simS_t * simR[t, :]) + simP_next
        
        return simC_t, simA_t, simQ_t, simS_t, simB_t, simW_next, simF_next, simP_next
    
    # 开始模拟
    for t in tqdm(range(tn-1), desc="模拟进度"):
        # 向量化计算所有模拟的财富收入比和现金持有量
        if t < tr - tb:  # 工作期
            simW_Y = simW_Y.at[t, :].set(simW[t, :] / simY[t, :])  # 上期财富-本期工资收入比
            cash = cash.at[t, :].set(simW[t, :])  # 当期现金为上期的财富
        else:  # 退休期
            cash = cash.at[t, :].set(simW[t, :])
        
        # 使用JAX向量化计算单期模拟
        simC_t, simA_t, simQ_t, simS_t, simB_t, simW_next, simF_next, simP_next = simulate_single_period(
            t, simW, simF, C, A, Q, gcash, gfund, simY, simR, rf, rp, pension_pct, pension_rate, tr, tb
        )
        
        # 更新模拟结果
        simC = simC.at[t, :].set(simC_t)
        simA = simA.at[t, :].set(simA_t)
        simQ = simQ.at[t, :].set(simQ_t)
        simS = simS.at[t, :].set(simS_t)
        simB = simB.at[t, :].set(simB_t)
        simW = simW.at[t+1, :].set(simW_next)
        simF = simF.at[t+1, :].set(simF_next)
        simP = simP.at[t+1, :].set(simP_next)
        
        # 记录消费占收入的比例
        simC_pct = simC_pct.at[t, :].set(simC_t)
    
    # 最后一期全部消费
    simC = simC.at[tn-1, :].set(1.0)
    simC_pct = simC_pct.at[tn-1, :].set(1.0)
    
    # 将JAX数组转回NumPy数组进行后续处理
    simC_np = np.array(simC)
    simC_pct_np = np.array(simC_pct)
    simY_np = np.array(simY)
    simW_np = np.array(simW)
    simS_np = np.array(simS)
    simB_np = np.array(simB)
    simW_Y_np = np.array(simW_Y)
    simA_np = np.array(simA)
    simGPY_np = np.array(simGPY)
    simQ_np = np.array(simQ)
    simF_np = np.array(simF)
    simP_np = np.array(simP)
    
    # 多次模拟path下变量平均值
    meanC = np.mean(simC_np, axis=1)
    meanC_pct = np.mean(simC_pct_np, axis=1)
    meanY = np.mean(simY_np, axis=1)