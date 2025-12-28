import scipy.io as sio
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import os
from matplotlib.font_manager import FontProperties
import matplotlib as mpl
import pandas as pd
from scipy.interpolate import interp2d, griddata, RegularGridInterpolator

# 设置中文字体支持
plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号
plt.rcParams['font.size'] = 12  # 增加默认字体大小
mpl.rcParams['font.family'] = 'SimHei'  # 确保全局字体设置为支持中文的字体

# 获取当前notebook所在目录的路径
current_dir = Path().absolute()

def simulate_model_cocco(mat_file_path):
# 首先打印文件路径，确认路径是否正确
    # print(f"尝试读取文件：{mat_file_path}")
    # print(f"文件是否存在：{mat_file_path.exists()}")

    # try:
    # 使用scipy.io替代h5py来读取.mat文件
    mat_contents = sio.loadmat(str(mat_file_path))

    # 获取model_results结构体中的数据
    model_results = mat_contents['model_results']

    # 读取各个数组
    # 对于MATLAB结构体，需要使用特殊的索引方式
    C = model_results['C'][0,0]  # 消费比例
    A = model_results['A'][0,0]  # 风险资产购买比例
    V = model_results['V'][0,0]  # 值函数
    tb = model_results['params'][0,0]['tb'][0,0][0,0]  # 时间起点
    td = model_results['params'][0,0]['td'][0,0][0,0]  # 时间终点
    tr = model_results['params'][0,0]['tr'][0,0][0,0]  # 退休时间点

    # 提取更多参数
    rf = model_results['params'][0,0]['rf'][0,0][0,0]  # 无风险收益率
    mu = model_results['params'][0,0]['mu'][0,0][0,0]  # 超额收益
    sigr = model_results['params'][0,0]['sigr'][0,0][0,0]  # 风险资产收益率的标准差
    smay = model_results['params'][0,0]['smay'][0,0][0,0]  # 白噪声shock的标准差
    smav = model_results['params'][0,0]['smav'][0,0][0,0]  # 持续shock的标准差
    ret_fac = model_results['params'][0,0]['ret_fac'][0,0][0,0]  # 退休后固定支付的工资（基本养老金）
    gy = model_results['params'][0,0]['gy'][0,0]  # 随年龄变化的工资增长率
    gcash = model_results['gcash'][0,0]  # 现金W_t的范围

    # 定义刻度


    # print("C策略函数图像已成功生成并保存到figs目录")

# 数值模拟函数

    """
    基于模型结果进行数值模拟，相当于MATLAB中的simulate_model函数
    """
    print('数值模拟开始...')
    
    # 设置随机种子以获得可重复的结果
    np.random.seed(42)
    
    # 获取模型参数
    ncash = len(gcash)
    tn = int(td - tb + 1)  # 总期数
    n = 5  # 外生随机冲击的grid数量，与MATLAB代码一致
    
    # 数值模拟参数
    nsim = 1000
    
    # 创建保存结果的数组
    meanY = np.zeros(tn)
    meanC = np.zeros(tn)
    meanW = np.zeros(tn)
    meanA = np.zeros(tn)
    meanS = np.zeros(tn)
    meanB = np.zeros(tn)

    meanWY = np.zeros(tn)
    meanalpha = np.zeros(tn)
    meanGPY = np.zeros(tn)
    
    simPY = np.zeros((tn, nsim))
    simGPY = np.zeros((tn, nsim))
    simY = np.zeros((tn, nsim))
    simC = np.zeros((tn, nsim))
    simW = np.zeros((tn, nsim))
    simS = np.zeros((tn, nsim))
    simB = np.zeros((tn, nsim))
    simA = np.zeros((tn, nsim))  # 风险资产投资比例 - 添加这个变量的初始化
    simW_Y = np.zeros((tn, nsim))
    simR = np.zeros((tn, nsim))
    cash = np.zeros((tn, nsim))
    simC_pct = np.zeros((tn, nsim))
    
    # 1、模拟生成labor income
    # print('生成劳动收入模拟数据...')
    for i1 in range(int(nsim/2)):  # 另外一半模拟完全对称
        # 工作期第一期
        eps_y = np.random.normal(0, 1)  # N(0,1)
        simPY[0, i1] = eps_y * smav  # 初始的p
        simPY[0, int(nsim/2) + i1] = -eps_y * smav
        simGPY[0, i1] = 1.0
        simGPY[0, int(nsim/2) + i1] = 1.0
        simTY = np.random.normal(0, 1)
        simY[0, i1] = np.exp(simTY * smay)
        simY[0, int(nsim/2) + i1] = np.exp(-simTY * smay)
        
        # 工作期第2期~退休
        for i2 in range(1, int(tr-tb)):
            w = i2 + tb
            eps_y = np.random.normal(0, 1)
            simPY[i2, i1] = eps_y * smav + simPY[i2-1, i1]
            simPY[i2, int(nsim/2) + i1] = -eps_y * smav + simPY[i2-1, int(nsim/2) + i1]
            simGPY[i2, i1] = np.exp(gy[i2-1]) * np.exp(simPY[i2, i1]) / np.exp(simPY[i2-1, i1])
            simGPY[i2, int(nsim/2) + i1] = np.exp(gy[i2-1]) * np.exp(simPY[i2, int(nsim/2) + i1]) / np.exp(simPY[i2-1, int(nsim/2) + i1])
            simTY = np.random.normal(0, 1)
            simY[i2, i1] = np.exp(simTY * smay)
            simY[i2, int(nsim/2) + i1] = np.exp(-simTY * smay)
    
    # 退休期
    for t in range(int(tr-tb), tn):
        simY[t, :] = ret_fac
        simGPY[t, :] = 1.0
    
    # 2、模拟风险投资的收益率
    # print('生成风险投资收益率模拟数据...')
    for t in range(tn):
        for i1 in range(int(nsim/2)):
            eps_r = np.random.normal(0, 1)
            simR[t, i1] = mu + rf + sigr * eps_r
            simR[t, int(nsim/2) + i1] = mu + rf - sigr * eps_r
    
    # 3、从第一期开始迭代，得到各控制变量的值
    # print('模拟控制变量...')
    simW[:, :] = 0  # 初始财富设置为0
    
    # 创建网格点
    
    for t in range(tn):
        # print(f'模拟进度: {t+1}/{tn}')
        
        # 向量化计算所有模拟的财富收入比和现金持有量
        if np.any(simY[t, :] == 0):
            # print(f'警告: 第{t+1}期存在零收入，已调整为极小值')
            simY[t, simY[t, :] == 0] = 1e-6
        
        simW_Y[t, :] = simW[t, :] / simY[t, :]  # 上期财富-本期工资收入比
        cash[t, :] = simW[t, :] + simY[t, :]  # cash-on-hand
        
        # 处理cash和fund超出范围的情况
        cash_t = np.clip(cash[t, :], gcash[0], gcash[-1])
        
        # 准备当前时期的策略函数值
        C_t = C[:,t].ravel()  # 展平为一维数组
        A_t = A[:,t].ravel()
        

        # 使用griddata进行插值，设置fill_value为0
        simC_t = griddata(gcash, C_t, cash_t, method='cubic')
        simA_t = griddata(gcash, A_t, cash_t, method='cubic')
        
        # 确保约束条件满足
        simC_t = np.clip(simC_t, 0, 1)  # 消费比例范围约束为0-1
        simA_t = np.clip(simA_t, 0, 1)  # 风险资产投资比例范围约束为0-1
        
        # 存储策略 
        simC[t, :] = simC_t * cash[t, :] # 具体消费额
        simC_pct[t, :] = simC_t  # 消费比例
        
        # 根据MATLAB代码的逻辑调整simA和simQ
        # 如果消费接近1，则风险资产和养老金购买都设为0
        # high_consumption_mask = simC_t > 0.99
        # simA_t[high_consumption_mask] = 0
        # simQ_t[high_consumption_mask] = 0
        
        simA[t, :] = simA_t
        
        # 计算各种模拟变量
        # 计算总储蓄
        sav = (1 - simC_pct[t, :]) * cash[t, :]  # 用于投资的金额 = (1-消费比例)*现金
        
        # 剩余投资
        remaining_sav = sav  # 除去养老金购买后剩余的储蓄
        
        # 风险资产投资
        simS[t, :] = simA[t, :] * remaining_sav  # 风险投资额
        simB[t, :] = remaining_sav - simS[t, :]  # 无风险投资额
        
        # 计算下期财富和养老基金余额
        if t < tn - 1:
            # 确保没有零除错误
            if np.any(simGPY[t+1, :] == 0):
                print(f'警告: 第{t+2}期存在零GPY，已调整为1')
                simGPY[t+1, simGPY[t+1, :] == 0] = 1
            
            # 更新财富
            simW[t+1, :] = (simB[t, :] * rf + simS[t, :] * simR[t, :]) / simGPY[t+1, :]
            
    
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
    
    # 计算平均实际消费值
    meanCash = meanW + meanY  # 平均现金持有量 = 平均财富 + 平均收入
    meanRealC = meanC * meanCash  # 平均实际消费值 = 消费比例 * 现金持有量
    
    
    # 创建表格保存模拟结果
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
        'meanRealC': meanRealC  # 添加实际消费值
    })

    return sim_results    

def simulate_model_baseline(mat_file_path):


    # try:
    # 使用scipy.io替代h5py来读取.mat文件
    mat_contents = sio.loadmat(str(mat_file_path))

    # 获取model_results结构体中的数据
    model_results = mat_contents['model_results']

    # 读取各个数组
    # 对于MATLAB结构体，需要使用特殊的索引方式
    C = model_results['C'][0,0]  # 消费比例
    A = model_results['A'][0,0]  # 风险资产购买比例
    Q = model_results['Q'][0,0]  # 养老金购买决策策略函数
    V = model_results['V'][0,0]  # 值函数
    tb = model_results['params'][0,0]['tb'][0,0][0,0]  # 时间起点
    td = model_results['params'][0,0]['td'][0,0][0,0]  # 时间终点
    tr = model_results['params'][0,0]['tr'][0,0][0,0]  # 退休时间点

    # 提取更多参数
    rf = model_results['params'][0,0]['rf'][0,0][0,0]  # 无风险收益率
    rp = model_results['params'][0,0]['rp'][0,0][0,0]  # 养老基金收益率
    mu = model_results['params'][0,0]['mu'][0,0][0,0]  # 超额收益
    sigr = model_results['params'][0,0]['sigr'][0,0][0,0]  # 风险资产收益率的标准差
    smay = model_results['params'][0,0]['smay'][0,0][0,0]  # 白噪声shock的标准差
    smav = model_results['params'][0,0]['smav'][0,0][0,0]  # 持续shock的标准差
    ret_fac = model_results['params'][0,0]['ret_fac'][0,0][0,0]  # 退休后固定支付的工资（基本养老金）
    pension_rate = model_results['params'][0,0]['pension_rate'][0,0][0,0]  # 养老金支付比例
    gy = model_results['params'][0,0]['gy'][0,0]  # 随年龄变化的工资增长率
    gcash = model_results['gcash'][0,0]  # 现金W_t的范围
    gfund = model_results['gfund'][0,0]   # 养老基金余额F_t的范围


    """
    基于模型结果进行数值模拟，相当于MATLAB中的simulate_model函数
    """
    print('数值模拟开始...')
    
    # 设置随机种子以获得可重复的结果
    np.random.seed(42)
    
    # 获取模型参数
    ncash = len(gcash)
    nfund = len(gfund)
    tn = int(td - tb + 1)  # 总期数
    n = 5  # 外生随机冲击的grid数量，与MATLAB代码一致
    
    # 数值模拟参数
    nsim = 1000
    
    # 创建保存结果的数组
    meanY = np.zeros(tn)
    meanC = np.zeros(tn)
    meanW = np.zeros(tn)
    meanA = np.zeros(tn)
    meanS = np.zeros(tn)
    meanB = np.zeros(tn)
    meanF = np.zeros(tn)  # 平均养老基金余额
    meanP = np.zeros(tn)  # 平均养老金给付
    meanQ = np.zeros(tn)  # 平均养老金购买比例
    meanWY = np.zeros(tn)
    meanalpha = np.zeros(tn)
    meanGPY = np.zeros(tn)
    
    simPY = np.zeros((tn, nsim))
    simGPY = np.zeros((tn, nsim))
    simY = np.zeros((tn, nsim))
    simC = np.zeros((tn, nsim))
    simW = np.zeros((tn, nsim))
    simS = np.zeros((tn, nsim))
    simB = np.zeros((tn, nsim))
    simF = np.zeros((tn, nsim))  # 模拟养老基金余额
    simP = np.zeros((tn, nsim))  # 模拟养老金给付
    simQ = np.zeros((tn, nsim))  # 模拟养老金购买比例
    simA = np.zeros((tn, nsim))  # 风险资产投资比例 - 添加这个变量的初始化
    simW_Y = np.zeros((tn, nsim))
    simR = np.zeros((tn, nsim))
    cash = np.zeros((tn, nsim))
    fund = np.zeros((tn, nsim))  # 每期养老基金账户余额
    simC_pct = np.zeros((tn, nsim))
    
    # 1、模拟生成labor income
    # print('生成劳动收入模拟数据...')
    for i1 in range(int(nsim/2)):  # 另外一半模拟完全对称
        # 工作期第一期
        eps_y = np.random.normal(0, 1)  # N(0,1)
        simPY[0, i1] = eps_y * smav  # 初始的p
        simPY[0, int(nsim/2) + i1] = -eps_y * smav
        simGPY[0, i1] = 1.0
        simGPY[0, int(nsim/2) + i1] = 1.0
        simTY = np.random.normal(0, 1)
        simY[0, i1] = np.exp(simTY * smay)
        simY[0, int(nsim/2) + i1] = np.exp(-simTY * smay)
        
        # 工作期第2期~退休
        for i2 in range(1, int(tr-tb)):
            w = i2 + tb
            eps_y = np.random.normal(0, 1)
            simPY[i2, i1] = eps_y * smav + simPY[i2-1, i1]
            simPY[i2, int(nsim/2) + i1] = -eps_y * smav + simPY[i2-1, int(nsim/2) + i1]
            simGPY[i2, i1] = np.exp(gy[i2-1]) * np.exp(simPY[i2, i1]) / np.exp(simPY[i2-1, i1])
            simGPY[i2, int(nsim/2) + i1] = np.exp(gy[i2-1]) * np.exp(simPY[i2, int(nsim/2) + i1]) / np.exp(simPY[i2-1, int(nsim/2) + i1])
            simTY = np.random.normal(0, 1)
            simY[i2, i1] = np.exp(simTY * smay)
            simY[i2, int(nsim/2) + i1] = np.exp(-simTY * smay)
    
    # 退休期
    for t in range(int(tr-tb), tn):
        simY[t, :] = ret_fac
        simGPY[t, :] = 1.0
    
    # 2、模拟风险投资的收益率
    # print('生成风险投资收益率模拟数据...')
    for t in range(tn):
        for i1 in range(int(nsim/2)):
            eps_r = np.random.normal(0, 1)
            simR[t, i1] = mu + rf + sigr * eps_r
            simR[t, int(nsim/2) + i1] = mu + rf - sigr * eps_r
    
    # 3、从第一期开始迭代，得到各控制变量的值
    # print('模拟控制变量...')
    simW[:, :] = 0  # 初始财富设置为0
    simF[:, :] = 0  # 初始养老基金账户余额设置为0
    
    # 创建网格点
    X, Y = np.meshgrid(gfund, gcash)
    grid_points = np.column_stack((X.ravel(), Y.ravel()))
    
    for t in range(tn):
        # print(f'模拟进度: {t+1}/{tn}')
        
        # 向量化计算所有模拟的财富收入比和现金持有量
        if np.any(simY[t, :] == 0):
            # print(f'警告: 第{t+1}期存在零收入，已调整为极小值')
            simY[t, simY[t, :] == 0] = 1e-6
        
        simW_Y[t, :] = simW[t, :] / simY[t, :]  # 上期财富-本期工资收入比
        cash[t, :] = simW[t, :] + simY[t, :]  # cash-on-hand
        fund[t, :] = simF[t, :]  # 养老基金账户余额
        
        # 处理cash和fund超出范围的情况
        cash_t = np.clip(cash[t, :], gcash[0], gcash[-1])
        fund_t = np.clip(fund[t, :], gfund[0], gfund[-1])
        
        # 准备当前时期的策略函数值
        C_t = C[:,:,t].ravel()  # 展平为一维数组
        A_t = A[:,:,t].ravel()
        Q_t = Q[:,:,t].ravel()
        
        # 确保没有None值
        # C_t = np.nan_to_num(C_t, nan=0.0)  # 将NaN替换为0
        # A_t = np.nan_to_num(A_t, nan=0.0)
        # Q_t = np.nan_to_num(Q_t, nan=0.0)
        
        # 准备查询点
        query_points = np.column_stack((fund_t, cash_t))
        
        # 使用griddata进行插值，设置fill_value为0
        simC_t = griddata(grid_points, C_t, query_points, method='cubic')
        simA_t = griddata(grid_points, A_t, query_points, method='cubic')
        simQ_t = griddata(grid_points, Q_t, query_points, method='cubic')
        
        # 确保约束条件满足
        simC_t = np.clip(simC_t, 0, 1)  # 消费比例范围约束为0-1
        simA_t = np.clip(simA_t, 0, 1)  # 风险资产投资比例范围约束为0-1
        simQ_t = np.clip(simQ_t, 0, 1)  # 养老金购买比例范围约束为0-1
        
        # 存储策略
        simC[t, :] = simC_t * cash[t, :]
        simC_pct[t, :] = simC_t  # 消费比例
        
        # 根据MATLAB代码的逻辑调整simA和simQ
        # 如果消费接近1，则风险资产和养老金购买都设为0
        # high_consumption_mask = simC_t > 0.99
        # simA_t[high_consumption_mask] = 0
        # simQ_t[high_consumption_mask] = 0
        
        simA[t, :] = simA_t
        simQ[t, :] = simQ_t
        
        # 计算各种模拟变量
        # 计算总储蓄
        sav = (1 - simC_pct[t, :]) * cash[t, :]  # 用于投资的金额 = (1-消费比例)*现金
        
        # 计算养老基金购买
        simP[t, :] = 0  # 默认养老金给付为0
        if t < (tr-tb):  # 工作期
            pension_purchase = simQ[t, :] * sav  # 养老金购买金额
        else:  # 退休期
            # 计算养老金给付
            simP[t, :] = fund[t, :] * pension_rate  # 每期领取养老基金的一定比例
            pension_purchase = np.zeros(nsim)  # 退休期不再购买养老金
        
        # 剩余投资
        remaining_sav = sav - pension_purchase  # 除去养老金购买后剩余的储蓄
        
        # 风险资产投资
        simS[t, :] = simA[t, :] * remaining_sav  # 风险投资额
        simB[t, :] = remaining_sav - simS[t, :]  # 无风险投资额
        
        # 计算下期财富和养老基金余额
        if t < tn - 1:
            # 确保没有零除错误
            if np.any(simGPY[t+1, :] == 0):
                print(f'警告: 第{t+2}期存在零GPY，已调整为1')
                simGPY[t+1, simGPY[t+1, :] == 0] = 1
            
            # 更新财富
            simW[t+1, :] = (simB[t, :] * rf + simS[t, :] * simR[t, :]) / simGPY[t+1, :] + simP[t, :]
            
            # 更新养老基金账户余额
            if t < (tr-tb):  # 工作期 - 养老基金增长
                simF[t+1, :] = (simF[t, :] + pension_purchase) * rp / simGPY[t+1, :]
            else:  # 退休期 - 养老基金减少，但至少保留原值（根据MATLAB代码的行为）
                simF[t+1, :] = simF[t, :]
    
    # 多次模拟path下变量平均值
    meanC = np.mean(simC, axis=1)
    meanC_pct = np.mean(simC_pct, axis=1)
    meanY = np.mean(simY, axis=1)
    meanW = np.mean(simW, axis=1)
    meanS = np.mean(simS, axis=1)
    meanB = np.mean(simB, axis=1)
    meanF = np.mean(simF, axis=1)  # 养老基金余额
    meanP = np.mean(simP, axis=1)  # 养老金给付
    meanQ = np.mean(simQ, axis=1)  # 养老金购买比例
    meanWY = np.mean(simW_Y, axis=1)
    meanalpha = np.mean(simA, axis=1)
    meanGPY = np.mean(simGPY, axis=1)
    
    # 计算平均实际消费值
    meanCash = meanW + meanY  # 平均现金持有量 = 平均财富 + 平均收入
    meanRealC = meanC * meanCash  # 平均实际消费值 = 消费比例 * 现金持有量
    
    # 保存模拟结果
    # sim_dir = current_dir / 'result_baseline_matlab_PFI' / 'simulation'
    # figs_dir = current_dir / 'figs'
    # if not os.path.exists(sim_dir):
    #     os.makedirs(sim_dir)
    
    # 创建表格保存模拟结果
    sim_results = pd.DataFrame({
        'meanC': meanC,
        'meanC_pct': meanC_pct,
        'meanY': meanY,
        'meanW': meanW,
        'meanS': meanS,
        'meanB': meanB,
        'meanF': meanF,
        'meanP': meanP,
        'meanQ': meanQ,
        'meanWY': meanWY,
        'meanalpha': meanalpha,
        'meanGPY': meanGPY,
        'meanRealC': meanRealC  # 添加实际消费值
    })
    # sim_results.to_csv(sim_dir / 'simulation_results.csv', index=False)

    return sim_results
    # 保存原始模拟数据
    # np.savez(sim_dir / 'raw_simulation.npz',
    #          simC=simC, simA=simA, simY=simY, simW=simW, simS=simS, simB=simB,
    #          simF=simF, simP=simP, simQ=simQ, simR=simR, simGPY=simGPY)
    
    
    


# 运行模拟函数
if __name__ == "__main__":

    mat_file_path = current_dir / 'result' / 'cocco.mat'
    cocco_sim_results = simulate_model_cocco(mat_file_path)

    mat_file_path = current_dir / 'result' / 'baseline.mat'
    baseline_sim_results = simulate_model_baseline(mat_file_path)

    mat_file_path = current_dir / 'result' / 'baseline_no_ret_fac.mat'
    baseline_no_ret_fac_sim_results = simulate_model_baseline(mat_file_path)  
    print(baseline_no_ret_fac_sim_results.meanY)
    # 获取时间参数 (使用baseline_sim_results的时间参数作为基准)
    tb = 18
    td = 100
    tr = 61
    age_range = np.arange(tb, td+1)

    # 创建保存图像的目录
    figs_dir = current_dir / 'figs'
    if not os.path.exists(figs_dir):
        os.makedirs(figs_dir)

    # 绘制比较图 - 2x2布局而非1x4布局
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))  # 增加图表尺寸
       
    # 第一个子图：比较三个模型的具体消费额度(meanC)
    axes[0,0].plot(age_range, cocco_sim_results['meanC'], linewidth=1.5)
    axes[0,0].plot(age_range, baseline_sim_results['meanC'], linewidth=1.5)
    axes[0,0].plot(age_range, baseline_no_ret_fac_sim_results['meanC'], linewidth=1.5)
    axes[0,0].axvline(x=tr, linestyle='--', color='gray', linewidth=1)
    axes[0,0].set_xlabel('年龄', fontsize=14)
    axes[0,0].set_ylabel('消费额度', fontsize=14)
    axes[0,0].set_title('(a) 具体消费额度', fontsize=16)
    axes[0,0].tick_params(axis='both', which='major', labelsize=12)
    axes[0,0].grid(True, alpha=0.3)
    axes[0,0].legend(['BP', 'EP+BP', 'EP', '退休年龄'], 
        loc='upper left', ncol=1, fontsize=10)   

    
    # 第二个子图：比较三个模型的风险资产比例(meanalpha)
    axes[0,1].plot(age_range, cocco_sim_results['meanalpha'], linewidth=1.5)
    axes[0,1].plot(age_range, baseline_sim_results['meanalpha'], linewidth=1.5)
    axes[0,1].plot(age_range, baseline_no_ret_fac_sim_results['meanalpha'], linewidth=1.5)
    axes[0,1].axvline(x=tr, linestyle='--', color='gray', linewidth=1)
    axes[0,1].set_xlabel('年龄', fontsize=14)
    axes[0,1].set_ylabel('比例', fontsize=14)
    axes[0,1].set_title('(b) 风险资产比例', fontsize=16)
    axes[0,1].tick_params(axis='both', which='major', labelsize=12)
    axes[0,1].grid(True, alpha=0.3)
    axes[0,1].legend(['BP', 'EP+BP', 'EP', '退休年龄'], 
        loc='lower left', ncol=1, fontsize=10)   

    # 第三个子图：比较三个模型的消费比例(meanC_pct)
    axes[1,0].plot(age_range, cocco_sim_results['meanC_pct'], linewidth=1.5)
    axes[1,0].plot(age_range, baseline_sim_results['meanC_pct'], linewidth=1.5)
    axes[1,0].plot(age_range, baseline_no_ret_fac_sim_results['meanC_pct'], linewidth=1.5)
    axes[1,0].axvline(x=tr, linestyle='--', color='gray', linewidth=1)
    axes[1,0].set_xlabel('年龄', fontsize=14)
    axes[1,0].set_ylabel('比例', fontsize=14)
    axes[1,0].set_title('(c) 消费比例', fontsize=16)
    axes[1,0].tick_params(axis='both', which='major', labelsize=12)
    axes[1,0].grid(True, alpha=0.3)
    axes[1,0].legend(['BP', 'EP+BP', 'EP', '退休年龄'], 
        loc='lower right', ncol=1, fontsize=10)   
    
    # 第四个子图：比较后两个模型的个人养老金决策比例(meanQ)
    l1 = axes[1,1].plot(age_range, baseline_sim_results['meanQ'], linewidth=1.5, label='EP+BP')
    l2 = axes[1,1].plot(age_range, baseline_no_ret_fac_sim_results['meanQ'], linewidth=1.5, label='EP')
    l3 = axes[1,1].axvline(x=tr, linestyle='--', color='gray', linewidth=1, label='退休年龄')
    axes[1,1].set_xlabel('年龄', fontsize=14)
    axes[1,1].set_ylabel('比例', fontsize=14)
    axes[1,1].set_title('(d) 养老金购买比例', fontsize=16)
    axes[1,1].tick_params(axis='both', which='major', labelsize=12)
    axes[1,1].grid(True, alpha=0.3)
    axes[1,1].legend([ 'EP+BP', 'EP', '退休年龄'], 
        loc='upper right', ncol=1, fontsize=10) 
    
    # 只在最后一个子图添加legend
    # 添加前三个子图的图例项
    # l0 = axes[1,1].plot([], [], linewidth=1.5, label='Cocco模型')[0]


    
    # 调整子图之间的间距和整体布局
    plt.tight_layout(rect=[0, 0.08, 1, 0.95])
    
    # 保存图片
    plt.savefig(figs_dir / 'p6_model_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()
    
    print(f'三个模型的比较图表已保存到 {figs_dir / "p6_model_comparison.png"}')