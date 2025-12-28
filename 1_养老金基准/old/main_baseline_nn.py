import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import minimize
from scipy.interpolate import interp1d
import pandas as pd
import os
from tqdm import tqdm
import argparse
# 添加PyTorch相关库
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import torch.nn.functional as F
import random
# 设置环境变量以解决PyTorch在某些系统上的库加载问题
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

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
def cocco_fun_valuefunc_retired(x, cash, fund, nextV, nextF, gret, rf, ret_fac, pension_pay, gamma, beta, psi_1, psi_2, theta, gcash, gfund, survprob, weig):
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
    
    返回:
    v - 函数值
    """
    auxVV = 0
    sav = cash - x[0]
    
    u = -(x[0])**(1-gamma)
    
    # 下期的cash-on-hand
    cash_1 = (rf * (1 - x[1]) + gret * x[1]) * sav + ret_fac + pension_pay
    
    # 插值
    int_V = interp1d(gcash, nextV, kind='cubic', fill_value='extrapolate')(cash_1)
    
    # 计算期望值
    auxVV = auxVV + np.dot(weig, survprob * int_V)
    
    if auxVV == 0:
        v = -(u)
    else:
        v = -(u + beta * auxVV)
    
    return v


# 定义神经网络模型（退休期值函数网络 - 32层深度）
class ValueFunctionNetworkRetired(nn.Module):
    def __init__(self, input_dim=3, hidden_dim=32, depth=32):
        super(ValueFunctionNetworkRetired, self).__init__()
        self.input_layer = nn.Linear(input_dim, hidden_dim)
        
        # 创建多层隐藏层
        self.hidden_layers = nn.ModuleList()
        for _ in range(depth-2):  # 除去输入层和输出层
            self.hidden_layers.append(nn.Linear(hidden_dim, hidden_dim))
        
        self.output_layer = nn.Linear(hidden_dim, 1)
        
    def forward(self, x):
        x = F.relu(self.input_layer(x))
        
        for layer in self.hidden_layers:
            x = F.relu(layer(x))
        
        # 输出值函数（不需要特别的激活函数，因为值函数可以是任意实数）
        return self.output_layer(x)

# 计算特定状态下退休期最优值函数
def compute_optimal_value_retired(cash, fund, next_value_func, next_fund_func, 
                                 gret, rf, ret_fac, pension_pay, gamma, beta, 
                                 psi_1, psi_2, theta, gcash, gfund, survprob, weig):
    """
    计算退休期特定状态下的最优值函数
    
    通过数值优化找到最优的消费和资产配置决策，并计算对应的值函数
    """
    from scipy.optimize import minimize
    
    # 定义目标函数（最小化负的值函数）
    def objective(x):
        return cocco_fun_valuefunc_retired(
            x, cash, fund, next_value_func, next_fund_func, 
            gret, rf, ret_fac, pension_pay, gamma, beta, 
            psi_1, psi_2, theta, gcash, gfund, survprob, weig
        )
    
    # 设置初始猜测值和约束条件
    x0 = [0.5 * cash, 0.5]  # 初始猜测：消费一半现金，风险资产占比50%
    bounds = [(0, 0.999 * cash), (0, 1)]  # 约束条件
    
    # 使用优化算法找到最优策略
    result = minimize(objective, x0, bounds=bounds, method='SLSQP')
    
    if not result.success:
        print(f"优化未收敛，状态: cash={cash}, fund={fund}")
    
    # 计算并返回最优值函数（注意：目标函数最小化的是负的值函数）
    return -objective(result.x)

# 训练退休期统一值函数网络
def train_value_network_retired(network, cash_grid, fund_grid, V, F, 
                               gret, rf, ret_fac, pension_rate, gamma, beta, 
                               psi_1, psi_2, theta, gcash, gfund, survprob, weig, 
                               tn, tr, tb, epochs=100, batch_size=64, learning_rate=0.001):
    """
    训练单一的退休期值函数网络
    
    参数:
    network - 待训练的神经网络
    cash_grid - 现金状态网格
    fund_grid - 养老基金状态网格
    V - 值函数数组
    F - 养老基金数组
    其他参数同前
    
    返回:
    trained_network - 训练好的神经网络
    """
    # 准备训练数据 - 所有退休期的状态和对应值函数
    states = []
    target_values = []
    
    print("生成退休期训练数据...")
    # 遍历所有退休期时间点
    for t in range(tr-tb, tn):
        time_normalized = (t - (tr-tb)) / (tn - (tr-tb))
        
        for cash_idx, cash_val in enumerate(cash_grid):
            for fund_idx, fund_val in enumerate(fund_grid):
                states.append([cash_val, fund_val, time_normalized])
                
                # 如果是最后一期，直接使用终值效用
                if t == tn-1:
                    value = V[cash_idx, fund_idx, t]  # 已经在solve_model中设置
                else:
                    # 查看值函数是否已经计算
                    if V[cash_idx, fund_idx, t] != 0:
                        value = V[cash_idx, fund_idx, t]
                    else:
                        # 计算养老金给付
                        pension_pay = float(pension_rate * fund_val)
                        
                        # 计算最优值函数
                        value = compute_optimal_value_retired(
                            cash_val, fund_val, V[:, :, t+1], F[:, :, t+1], 
                            gret, rf, ret_fac, pension_pay, gamma, beta, psi_1, psi_2, 
                            theta, gcash, gfund, survprob[t], weig
                        )
                
                target_values.append(value)
    
    # 转换为PyTorch张量
    states_tensor = torch.FloatTensor(states)
    target_tensor = torch.FloatTensor(target_values).unsqueeze(1)
    
    # 创建数据集和数据加载器
    dataset = TensorDataset(states_tensor, target_tensor)
    dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)
    
    # 定义优化器和损失函数
    optimizer = optim.Adam(network.parameters(), lr=learning_rate)
    criterion = nn.MSELoss()
    
    # 训练网络
    print(f"开始训练退休期值函数网络...")
    for epoch in range(epochs):
        total_loss = 0
        for batch_states, batch_targets in dataloader:
            # 梯度清零
            optimizer.zero_grad()
            
            # 前向传播
            outputs = network(batch_states)
            
            # 计算损失
            loss = criterion(outputs, batch_targets)
            
            # 反向传播与优化
            loss.backward()
            optimizer.step()
            
            total_loss += loss.item()
        
        # 每10个epoch打印一次损失
        if (epoch + 1) % 10 == 0:
            print(f"Epoch [{epoch+1}/{epochs}], Loss: {total_loss/len(dataloader):.6f}")
    
    return network

# 使用神经网络计算最优策略（退休期）
def compute_optimal_policy_retired_with_nn(cash, fund, retired_value_net, gret, rf, ret_fac, 
                                          pension_pay, gamma, beta, time_normalized_next, weig, survprob):
    """
    使用训练好的神经网络计算退休期最优策略
    """
    from scipy.optimize import minimize
    
    # 定义目标函数（使用神经网络预测下期值函数）
    def objective(x):
        sav = cash - x[0]
        u = -(x[0])**(1-gamma)
        
        # 计算下期状态并使用神经网络预测值函数
        next_states = []
        for idx, g_ret in enumerate(gret):
            cash_next = (rf * (1 - x[1]) + g_ret * x[1]) * sav + ret_fac + pension_pay
            next_states.append([cash_next, fund, time_normalized_next])
        
        next_states_tensor = torch.FloatTensor(next_states)
        with torch.no_grad():
            next_values = retired_value_net(next_states_tensor).numpy().flatten()
        
        # 计算期望值
        auxVV = np.dot(weig, survprob * next_values)
        
        if auxVV == 0:
            v = -(u)
        else:
            v = -(u + beta * auxVV)
        
        return v
    
    # 设置初始猜测值和约束条件
    x0 = [0.5 * cash, 0.5]  # 初始猜测：消费一半现金，风险资产占比50%
    bounds = [(0, 0.999 * cash), (0, 1)]  # 约束条件
    
    # 使用优化算法找到最优策略
    result = minimize(objective, x0, bounds=bounds, method='SLSQP')
    
    if not result.success:
        print(f"策略优化未收敛，状态: cash={cash}, fund={fund}")
    
    return result.x, -objective(result.x)

# 添加工作期的值函数计算
def cocco_fun_valuefunc_work(x, cash, fund, gyp, nextV, nextF, yh, gret, rf, rp, pension_pct, gamma, beta, psi_1, psi_2, theta, gcash, gfund, survprob, n, nweig1):
    """
    工作期值函数计算
    
    参数:
    x - 决策变量 [消费, 风险资产投资比例, 个人养老金购买比例]
    cash - 手中现金
    fund - 养老基金账户余额
    gyp - 随机收入过程
    nextV - 下期值函数
    nextF - 下期养老基金函数
    yh - 劳动收入冲击
    gret - 随机收益率
    rf - 无风险收益率
    rp - 养老基金收益率
    pension_pct - 强制缴纳养老金比例
    gamma - 相对风险规避系数
    beta - 贴现因子
    psi_1, psi_2, theta - Epstein-Zin效用函数参数
    gcash - 现金网格
    gfund - 养老基金网格
    survprob - 存活概率
    n - 随机状态数量
    nweig1 - 随机状态权重
    
    返回:
    v - 函数值
    """
    auxVV = 0
    au = 0
    sav = cash - x[0] - x[2] * (1 - pension_pct)
    
    if sav <= 0:
        return 1e10  # 如果储蓄为负，返回一个大的惩罚值
    
    # 效用函数
    u = -(x[0])**(1-gamma)
    
    # 计算期望值函数
    for i1 in range(n):
        for i2 in range(n):
            for i3 in range(n):
                # 下期现金流
                pension_contribution = pension_pct + x[2] * (1 - pension_pct)
                cash_1 = (1 - pension_contribution) * gyp[i2, i3] * yh[i1, i2] + (rf * (1 - x[1]) + gret[i3] * x[1]) * sav
                
                # 下期养老基金
                fund_1 = fund * rp + pension_contribution * gyp[i2, i3] * yh[i1, i2]
                
                # 插值计算期望
                cash_1 = min(max(cash_1, gcash[0]), gcash[-1])
                fund_1 = min(max(fund_1, gfund[0]), gfund[-1])
                
                # 找到在现金网格中的位置
                cash_idx = f_ntoil_modified(cash_1, gcash, len(gcash))
                fund_idx = f_ntoil_modified(fund_1, gfund, len(gfund))
                
                # 双线性插值
                cash_weight = (cash_1 - gcash[cash_idx]) / (gcash[min(cash_idx+1, len(gcash)-1)] - gcash[cash_idx]) if cash_idx < len(gcash)-1 else 0
                fund_weight = (fund_1 - gfund[fund_idx]) / (gfund[min(fund_idx+1, len(gfund)-1)] - gfund[fund_idx]) if fund_idx < len(gfund)-1 else 0
                
                v00 = nextV[cash_idx, fund_idx]
                v01 = nextV[cash_idx, min(fund_idx+1, len(gfund)-1)] if fund_idx < len(gfund)-1 else v00
                v10 = nextV[min(cash_idx+1, len(gcash)-1), fund_idx] if cash_idx < len(gcash)-1 else v00
                v11 = nextV[min(cash_idx+1, len(gcash)-1), min(fund_idx+1, len(gfund)-1)] if cash_idx < len(gcash)-1 and fund_idx < len(gfund)-1 else v00
                
                int_V = (1-cash_weight)*(1-fund_weight)*v00 + (1-cash_weight)*fund_weight*v01 + cash_weight*(1-fund_weight)*v10 + cash_weight*fund_weight*v11
                
                auxVV = auxVV + nweig1[i1, i2, i3] * survprob * int_V
    
    if auxVV == 0:
        v = -(u)
    else:
        v = -(u + beta * auxVV)
    
    return v


# 定义工作期值函数网络模型（32层深度）
class ValueFunctionNetworkWorking(nn.Module):
    def __init__(self, input_dim=2, hidden_dim=32, depth=32):
        super(ValueFunctionNetworkWorking, self).__init__()
        self.input_layer = nn.Linear(input_dim, hidden_dim)
        
        # 创建多层隐藏层
        self.hidden_layers = nn.ModuleList()
        for _ in range(depth-2):  # 除去输入层和输出层
            self.hidden_layers.append(nn.Linear(hidden_dim, hidden_dim))
        
        self.output_layer = nn.Linear(hidden_dim, 1)
        
    def forward(self, x):
        x = F.relu(self.input_layer(x))
        
        for layer in self.hidden_layers:
            x = F.relu(layer(x))
        
        # 输出值函数（不需要特别的激活函数，因为值函数可以是任意实数）
        return self.output_layer(x)

# 计算特定状态下工作期最优值函数
def compute_optimal_value_working(cash, fund, gyp, nextV, nextF, yh, gret, rf, rp, 
                                 pension_pct, gamma, beta, psi_1, psi_2, theta, 
                                 gcash, gfund, survprob, n, nweig1):
    """
    计算工作期特定状态下的最优值函数
    
    通过数值优化找到最优的消费、资产配置和养老金购买决策，并计算对应的值函数
    """
    from scipy.optimize import minimize
    
    # 定义目标函数（最小化负的值函数）
    def objective(x):
        return cocco_fun_valuefunc_work(
            x, cash, fund, gyp, nextV, nextF, yh, gret, rf, rp, 
            pension_pct, gamma, beta, psi_1, psi_2, theta, 
            gcash, gfund, survprob, n, nweig1
        )
    
    # 设置初始猜测值和约束条件
    x0 = [0.5 * cash, 0.5, 0.1]  # 初始猜测：消费一半现金，风险资产占比50%，养老金购买比例10%
    bounds = [(0, 0.999 * cash), (0, 1), (0, 1)]  # 约束条件
    
    # 使用优化算法找到最优策略
    result = minimize(objective, x0, bounds=bounds, method='SLSQP')
    
    if not result.success:
        print(f"优化未收敛，状态: cash={cash}, fund={fund}")
    
    # 计算并返回最优值函数（注意：目标函数最小化的是负的值函数）
    return -objective(result.x), result.x

# 训练工作期值函数网络
def train_value_network_working(network, cash_grid, fund_grid, gyp, nextV, nextF, 
                              yh, gret, rf, rp, pension_pct, gamma, beta, psi_1, 
                              psi_2, theta, gcash, gfund, survprob, n, nweig1, 
                              epochs=100, batch_size=64, learning_rate=0.001):
    """
    训练工作期值函数网络
    
    参数:
    network - 待训练的神经网络
    cash_grid - 现金状态网格
    fund_grid - 养老基金状态网格
    gyp - 随机收入过程
    nextV - 下期值函数
    nextF - 下期养老基金函数
    其他参数与compute_optimal_value_working函数相同
    
    返回:
    trained_network - 训练好的神经网络
    optimal_policies - 每个状态下的最优策略 [C, A, Q]
    """
    # 准备训练数据
    states = []
    target_values = []
    optimal_policies = []
    
    # 生成训练数据
    print("生成训练数据...")
    for cash_idx, cash_val in enumerate(tqdm(cash_grid)):
        for fund_idx, fund_val in enumerate(fund_grid):
            states.append([cash_val, fund_val])
            
            # 计算最优值函数和策略
            value, policy = compute_optimal_value_working(
                cash_val, fund_val, gyp, nextV, nextF, yh, gret, rf, rp, 
                pension_pct, gamma, beta, psi_1, psi_2, theta, 
                gcash, gfund, survprob, n, nweig1
            )
            
            target_values.append(value)
            optimal_policies.append(policy)
    
    # 转换为PyTorch张量
    states_tensor = torch.FloatTensor(states)
    target_tensor = torch.FloatTensor(target_values).unsqueeze(1)
    
    # 创建数据集和数据加载器
    dataset = TensorDataset(states_tensor, target_tensor)
    dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)
    
    # 定义优化器和损失函数
    optimizer = optim.Adam(network.parameters(), lr=learning_rate)
    criterion = nn.MSELoss()
    
    # 训练网络
    print(f"开始训练值函数网络...")
    for epoch in range(epochs):
        total_loss = 0
        for batch_states, batch_targets in dataloader:
            # 梯度清零
            optimizer.zero_grad()
            
            # 前向传播
            outputs = network(batch_states)
            
            # 计算损失
            loss = criterion(outputs, batch_targets)
            
            # 反向传播与优化
            loss.backward()
            optimizer.step()
            
            total_loss += loss.item()
        
        # 每10个epoch打印一次损失
        if (epoch + 1) % 10 == 0:
            print(f"Epoch [{epoch+1}/{epochs}], Loss: {total_loss/len(dataloader):.6f}")
    
    # 返回训练好的网络和最优策略
    optimal_policies = np.array(optimal_policies)
    return network, optimal_policies

def solve_model():
    """
    使用神经网络拟合值函数方法求解模型的最优消费和投资决策
    """
    # 设置随机种子以获得可重复的结果
    np.random.seed(42)
    torch.manual_seed(42)
    random.seed(42)
    
    # 创建结果存储目录
    results_dir = 'result_baseline_python_nn'
    os.makedirs(results_dir, exist_ok=True)
    
    # 变量定义
    tb = 18    # 初始开始的年纪
    tr = 61    # 退休年龄
    td = 100   # 死亡年龄
    tn = td - tb + 1  # 总期数
    
    # 状态变量grid数量
    ncash = 21  # 手中现金
    nfund = 21  # 养老基金余额
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
    rp = 1.02    # 个人养老基金的无风险收益率
    
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
    gfund = np.zeros(nfund)
    lgfund = np.zeros(nfund)
    secd = np.zeros(ncash)
    C = np.zeros((ncash, nfund, tn))
    V = np.zeros((ncash, nfund, tn))
    A = np.ones((ncash, nfund, tn))
    Q = np.zeros((ncash, nfund, tn))  # 个人养老金购买决策
    F = np.zeros((ncash, nfund, tn))  # 养老基金余额状态变量
    
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
    
    # fund的grid
    maxfund = 100
    minfund = 0.25
    l_maxfund = np.log(maxfund)
    l_minfund = np.log(minfund)
    stepfund = (l_maxfund - l_minfund) / (nfund - 1)
    
    for i1 in range(nfund):
        lgfund[i1] = l_minfund + i1 * stepfund
    
    for i1 in range(nfund):
        gfund[i1] = np.exp(lgfund[i1])
    
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
            C[i1, i2, tn - 1] = gcash[i1]  # 最后一期的现金全部用来消费
            A[i1, i2, tn - 1] = 0.0      # 最后一期的投资全部为零
            Q[i1, i2, tn - 1] = 0.0      # 最后一期不购买养老金
            F[i1, i2, tn - 1] = gfund[i2] # 记录养老基金状态
            V[i1, i2, tn - 1] = -(C[i1, i2, tn - 1])**(1 - gamma)
    
    # 计算个人养老金给付
    # 假设退休后按照剩余年份平均分配养老基金余额
    pension_period = td - tr + 1
    pension_rate = 1 / pension_period
    
    # 退休期 - 使用单一神经网络值函数近似
    print("开始求解退休期值函数（使用单一32层神经网络拟合整个退休期值函数）...")
    
    # 初始化一个单一的神经网络模型 - 32层深度, 3个输入 (cash, fund, time)
    retired_value_net = ValueFunctionNetworkRetired(input_dim=3, hidden_dim=32, depth=32)
    
    # 从倒数第二期开始计算退休期的值函数，用于训练数据
    print("预计算退休期的值函数用于训练...")
    for t in tqdm(range(tn-2, tr-tb-1, -1), desc="预计算退休期值函数"):
        for i2 in range(nfund):  # 遍历养老基金状态
            for i3 in range(ncash):  # 遍历现金状态
                # 养老金给付
                pension_pay = pension_rate * gfund[i2]
                
                # 计算最优值函数和策略
                def objective(x):
                    return cocco_fun_valuefunc_retired(
                        x, gcash[i3], gfund[i2], 
                        V[:, :, t+1], F[:, :, t+1], 
                        gret, rf, ret_fac, pension_pay, 
                        gamma, beta, psi_1, psi_2, theta, 
                        gcash, gfund, survprob[t], weig
                    )
                
                x0 = [0.5 * gcash[i3], 0.5]  # 初始猜测
                bounds = [(0, 0.999 * gcash[i3]), (0, 1)]
                
                result = minimize(objective, x0, bounds=bounds, method='SLSQP')
                
                # 更新值函数和策略
                V[i3, i2, t] = -objective(result.x)
                C[i3, i2, t] = result.x[0]
                A[i3, i2, t] = result.x[1]
                F[i3, i2, t] = gfund[i2]
    
    # 训练统一的退休期值函数网络
    print("开始训练统一的退休期值函数网络...")
    retired_value_net = train_value_network_retired(
        retired_value_net, gcash, gfund, 
        V, F, gret, rf, ret_fac, pension_rate, 
        gamma, beta, psi_1, psi_2, theta, 
        gcash, gfund, survprob, weig, 
        tn, tr, tb, epochs=100, batch_size=64, learning_rate=0.001
    )
    
    # 保存训练好的退休期值函数网络
    torch.save(retired_value_net.state_dict(), f"{results_dir}/retired_value_net.pt")
    
    # 使用训练好的网络重新计算所有退休期的值函数和策略
    print("使用训练好的网络更新退休期值函数和策略...")
    for t in tqdm(range(tr-tb, tn), desc="更新退休期策略"):
        time_normalized = (t - (tr-tb)) / (tn - (tr-tb))
        time_normalized_next = ((t+1) - (tr-tb)) / (tn - (tr-tb)) if t < tn-1 else 1.0
        
        for i2 in range(nfund):  # 遍历养老基金状态
            for i3 in range(ncash):  # 遍历现金状态
                # 更新养老基金状态
                F[i3, i2, t] = gfund[i2]
                
                # 如果是最后一期，直接使用终值
                if t == tn-1:
                    C[i3, i2, t] = gcash[i3]
                    A[i3, i2, t] = 0.0
                    with torch.no_grad():
                        V[i3, i2, t] = -(gcash[i3])**(1-gamma)
                else:
                    # 养老金给付
                    pension_pay = pension_rate * gfund[i2]
                    
                    # 使用神经网络计算最优策略
                    policy, value = compute_optimal_policy_retired_with_nn(
                        gcash[i3], gfund[i2], retired_value_net, 
                        gret, rf, ret_fac, pension_pay, gamma, beta, 
                        time_normalized_next, weig, survprob[t]
                    )
                    
                    # 更新策略和值函数
                    C[i3, i2, t] = policy[0]
                    A[i3, i2, t] = policy[1]
                    
                    # 使用神经网络预测值函数
                    state = torch.FloatTensor([[gcash[i3], gfund[i2], time_normalized]])
                    with torch.no_grad():
                        V[i3, i2, t] = retired_value_net(state).item()

    # 工作期（退休前）- 使用神经网络值函数近似
    print("开始求解工作期值函数（使用32层神经网络拟合值函数）...")
    for i1 in tqdm(range(1, tr - tb + 1), desc="工作期求解进度"):
        t = tr - tb - i1
        
        # 初始化当前期的值函数神经网络模型 - 32层深度网络
        value_net = ValueFunctionNetworkWorking(input_dim=2, hidden_dim=32, depth=32)
        
        # 训练值函数神经网络并获取最优策略
        print(f"训练时期 {t} 的值函数网络...")
        value_net, optimal_policies = train_value_network_working(
            value_net, gcash, gfund, 
            gyp[:, :, t], V[:, :, t+1], F[:, :, t+1], 
            yh, gret, rf, rp, pension_pct, 
            gamma, beta, psi_1, psi_2, theta, 
            gcash, gfund, survprob[t], n, nweig1, 
            epochs=50, batch_size=32, learning_rate=0.001
        )
        
        # 使用训练好的网络计算所有状态的值函数
        policy_idx = 0
        for i2 in range(nfund):  # 遍历养老基金状态
            for i3 in range(ncash):  # 遍历现金状态
                # 构建输入状态
                state = torch.FloatTensor([[gcash[i3], gfund[i2]]])
                
                # 预测值函数
                with torch.no_grad():
                    V[i3, i2, t] = value_net(state).item()
                
                # 更新策略（使用优化过程中得到的最优策略）
                policy = optimal_policies[policy_idx]
                C[i3, i2, t] = policy[0]
                A[i3, i2, t] = policy[1]
                Q[i3, i2, t] = policy[2]
                
                policy_idx += 1
        
        # 保存训练好的网络模型
        torch.save(value_net.state_dict(), f"{results_dir}/value_net_t{t}.pt")
    
    # 保存结果
    os.makedirs(results_dir, exist_ok=True)
    
    # 保存完整的三维策略函数和值函数
    np.save(f'{results_dir}/C.npy', C)
    np.save(f'{results_dir}/A.npy', A)
    np.save(f'{results_dir}/Q.npy', Q)
    np.save(f'{results_dir}/V.npy', V)
    np.save(f'{results_dir}/F.npy', F)
    np.save(f'{results_dir}/gcash.npy', gcash)
    np.save(f'{results_dir}/gfund.npy', gfund)
    
    # 同时保存一个切片用于Excel查看（可选）
    pd.DataFrame(C[:, 0, :]).to_excel(f'{results_dir}/C_slice.xlsx')
    pd.DataFrame(A[:, 0, :]).to_excel(f'{results_dir}/A_slice.xlsx')
    pd.DataFrame(Q[:, 0, :]).to_excel(f'{results_dir}/Q_slice.xlsx')
    pd.DataFrame(V[:, 0, :]).to_excel(f'{results_dir}/V_slice.xlsx')
    pd.DataFrame(gcash).to_excel(f'{results_dir}/gcash.xlsx')
    pd.DataFrame(gfund).to_excel(f'{results_dir}/gfund.xlsx')

    print("模型求解完成，结果已保存到" + results_dir + "文件夹")
    
    # 返回求解结果，供模拟使用
    model_results = {
        'C': C,
        'A': A,
        'Q': Q,
        'V': V,
        'F': F,
        'gcash': gcash,
        'gfund': gfund,
        'retired_value_net': retired_value_net,  # 返回训练好的退休期值函数网络
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
            results_dir = 'result_baseline_python_nn'
            if not os.path.exists(f'{results_dir}/C.npy'):
                print("未找到模型求解结果文件，请先运行模型求解")
                return
                
            # 从文件中读取模型结果
            C = np.load(f'{results_dir}/C.npy')
            A = np.load(f'{results_dir}/A.npy')
            Q = np.load(f'{results_dir}/Q.npy')
            V = np.load(f'{results_dir}/V.npy')
            F = np.load(f'{results_dir}/F.npy')
            gcash = np.load(f'{results_dir}/gcash.npy')
            gfund = np.load(f'{results_dir}/gfund.npy')
            
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
                'F': F,
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
    
    # 使用模型结果进行模拟
    # 从模型结果中获取参数和数据
    C = model_results['C']
    A = model_results['A']
    Q = model_results['Q']
    V = model_results['V']
    F = model_results['F']
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
    simA = np.zeros((tn, nsim))
    simS = np.zeros((tn, nsim))
    simB = np.zeros((tn, nsim))
    simW_Y = np.zeros((tn, nsim))
    simR = np.zeros((tn, nsim))
    simQ = np.zeros((tn, nsim))  # 个人养老金购买决策
    simF = np.zeros((tn, nsim))  # 养老基金账户余额
    simP = np.zeros((tn, nsim))  # 个人养老金给付
    eps_y = np.zeros((1, 1))
    simTY = np.zeros((1, 1))
    eps_r = np.zeros((1, 1))
    cash = np.zeros((tn, nsim))
    simC_pct = np.zeros((tn, nsim))
    
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
    
    # 从第一期开始迭代，得到各控制变量的值
    simW[:, :] = 0
    simF[:, :] = 0
    
    # 创建二维插值函数
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
    for t in tqdm(range(tn), desc="模拟进度"):
        # 向量化计算所有模拟的财富收入比和现金持有量
        if t < tr - tb:  # 工作期
            simW_Y[t, :] = simW[t, :] / simY[t, :]  # 上期财富-本期工资收入比
            cash[t, :] = simW[t, :]  # 当期现金为上期的财富
        else:  # 退休期
            cash[t, :] = simW[t, :]
        
        # 根据当期的状态查找或插值获取策略
        for i in range(nsim):
            # 使用双线性插值获取当前状态下的策略
            simC[t, i] = get_bilinear_interp(C[:, :, t], gcash, gfund, cash[t, i], simF[t, i])
            simA[t, i] = get_bilinear_interp(A[:, :, t], gcash, gfund, cash[t, i], simF[t, i])
            
            # 确保约束条件满足
            simC[t, i] = max(min(simC[t, i], 0.9999 * cash[t, i]), 0)
            simA[t, i] = max(min(simA[t, i], 1), 0)
            
            # 工作期，还需要获取养老金购买比例
            if t < tr - tb:
                simQ[t, i] = get_bilinear_interp(Q[:, :, t], gcash, gfund, cash[t, i], simF[t, i])
                simQ[t, i] = max(min(simQ[t, i], 1), 0)
        
        # 记录消费占收入的比例
        if t < tr - tb:  # 工作期
            simC_pct[t, :] = simC[t, :] / simY[t, :]
        else:  # 退休期
            simC_pct[t, :] = simC[t, :] / (ret_fac + simP[t, :])
        
        # 计算个人养老金购买金额和养老基金账户更新
        if t < tr - tb:  # 工作期
            # 个人养老金购买
            pension_contrib = simQ[t, :] * (1 - pension_pct) * simY[t, :]
            
            # 计算养老基金账户更新
            if t == 0:
                simF[t, :] = pension_contrib + pension_pct * simY[t, :]
            else:
                simF[t, :] = simF[t-1, :] * rp + pension_contrib + pension_pct * simY[t, :]
        else:  # 退休期
            # 计算个人养老金给付
            simP[t, :] = pension_rate * simF[tr-tb-1, :]
            
            # 养老基金余额保持不变
            if t == tr - tb:
                simF[t, :] = simF[t-1, :]
            else:
                simF[t, :] = simF[t-1, :]
        
        # 剩余用于投资的资金
        sav = cash[t, :] - simC[t, :]
        
        # 工作期扣除个人养老金购买金额
        if t < tr - tb:
            sav = sav - simQ[t, :] * (1 - pension_pct) * simY[t, :]
        
        # 风险投资额
        simS[t, :] = simA[t, :] * sav
        simS[t, :] = np.minimum(simS[t, :], sav)
        
        # 无风险投资额
        simB[t, :] = sav - simS[t, :]
        
        # 计算下期财富
        if t < tn - 1:
            # 工作期
            if t < tr - tb - 1:
                simW[t + 1, :] = (simB[t, :] * rf + simS[t, :] * simR[t, :]) / simGPY[t + 1, :] + (1 - pension_pct) * simY[t + 1, :]
            # 退休期
            else:
                simW[t + 1, :] = (simB[t, :] * rf + simS[t, :] * simR[t, :]) + ret_fac + simP[t+1, :]
    
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
    meanF = np.mean(simF, axis=1)
    meanP = np.mean(simP, axis=1)
    
    # 保存模拟结果
    os.makedirs('result_baseline_python_nn', exist_ok=True)
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
    sim_results.to_excel('result_baseline_python_nn/simulation_results.xlsx')
    
    # 绘制结果
    plt.figure(figsize=(15, 10))
    
    # 消费比例、风险资产配置和个人养老金购买比例
    plt.subplot(2, 2, 1)
    plt.plot(meanC_pct, label='消费比例')
    plt.plot(meanalpha, label='风险资产配置')
    plt.plot(meanQ, label='养老金购买比例')
    plt.axvline(x=tr-tb, color='k', linestyle='--', label='退休')
    plt.legend()
    plt.xlabel('期数')
    plt.ylabel('比例')
    plt.title('消费、风险资产和养老金购买比例')
    
    # 第二个图：财富和养老金综合展示
    plt.subplot(2, 2, 2)
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
    plt.savefig('result_baseline_python_nn/results_plot_with_pension.png')
    plt.show()
    
    print("数值模拟完成，结果已保存到result_baseline_python_nn文件夹")

def main(mode="both"):
    """
    主函数，可选择运行模式：
    "solve": 只进行模型求解
    "simulate": 只进行数值模拟（需要先有求解结果）
    "both": 先求解再模拟（默认）
    """
    if mode.lower() == "solve" or mode.lower() == "both":
        print("开始求解模型（策略函数迭代方法）...")
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
    parser = argparse.ArgumentParser(description='生命周期模型求解与模拟（策略函数迭代方法）')
    parser.add_argument('--mode', type=str, choices=['solve', 'simulate', 'both'],
        default='both', help='运行模式：求解(solve)，模拟(simulate)，或两者都运行(both)')
    # 解析命令行参数
    args = parser.parse_args()
    # 运行主函数
    main(args.mode)