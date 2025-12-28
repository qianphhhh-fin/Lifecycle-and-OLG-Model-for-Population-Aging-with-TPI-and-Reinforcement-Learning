import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from scipy.stats import norm


# Variable Definitions
tb = 20  # initial starting age 
tr = 66  # retire age
td = 100  # death age 
tn = td - tb + 1  # total period

# Grid dimensions
na = 21      # asset grid num
ncash = 21   # cash grid num 
n = 3        # stochastic states
nc = 21      # consumption grid

# Parameters
maxcash = 200.0
mincash = 0.25

# Income coefficients
aa = -2.170042 + 2.700381
b1 = 0.16818
b2 = -0.0323371/10
b3 = 0.0019704/100

# Other parameters
ret_fac = 0.68212  # fixed pension
smay = 0.1    # std of transitory shock
smav = 0.1    # std of permanent shock
corr_v = 0.0  # correlation between wage and return shocks
corr_y = 0.0  # correlation between AR(1) wage and return shocks
rho = 10.0    # relative risk aversion
delta = 0.97  # discount factor
psi = 0.3     # elasticity of intertemporal substitution
r = 1.015     # risk-free return
mu = 0.04     # excess return
sigr = 0.2    # std of risky return

# Initialize arrays
survprob = np.zeros(tn-1)
delta2 = np.zeros(tn-1)
grid = np.zeros(n)
weig = np.zeros(n)
gret = np.zeros(n)
ones_n_1 = np.ones(n)
grid2 = np.zeros(n)
yp = np.zeros((n,n))
yh = np.zeros((n,n))
nweig1 = np.zeros((n,n,n))
f_y = np.zeros(tr-tb+1)
gy = np.zeros(tr-tb)
gyp = np.zeros((n,n,tn-1))
gcash = np.zeros(ncash)
lgcash = np.zeros(ncash)
ga = np.zeros(na)
riskret = np.zeros((na,n))
gc = np.zeros(nc)
auxV = np.zeros((na,nc))
vec_V = np.zeros(na*nc)
secd = np.zeros(ncash)
C = np.zeros((ncash,tn))
c = np.zeros((ncash,tn))
V = np.zeros((ncash,tn))
A = np.ones((ncash,tn))

# Simulation parameters
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
simPY = np.zeros((tn,nsim))
simGPY = np.zeros((tn,nsim))
simY = np.zeros((tn,nsim))
simC = np.zeros((tn,nsim))
simc = np.zeros((tn,nsim))
simW = np.zeros((tn,nsim))
simA = np.zeros((tn,nsim))
simS = np.zeros((tn,nsim))
simB = np.zeros((tn,nsim))
simW_Y = np.zeros((tn,nsim))
simR = np.zeros((tn,nsim))




# Part 2: 辅助函数定义
def f_randn(size):
    """生成标准正态随机数"""
    return np.random.normal(0, 1, size)

def discretize_AR1_Tauchen(mu, rho, sigma, n, m):
    """Tauchen方法离散化AR(1)过程"""
    step = 2 * m * sigma / np.sqrt(1 - rho**2) / (n-1)
    grid = np.zeros(n)
    for i in range(n):
        grid[i] = -m * sigma / np.sqrt(1 - rho**2) + i * step
    
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
    
    # 计算不变分布
    eigvals, eigvecs = np.linalg.eig(trans_mat.T)
    idx = np.argmin(np.abs(eigvals - 1))
    stationary_dist = eigvecs[:,idx].real
    stationary_dist = stationary_dist / np.sum(stationary_dist)
    
    return grid, stationary_dist

# Part 3: 初始化计算
# 生存概率设置
survprob = np.array([
0.99975, 0.99974, 0.99973, 0.99972, 0.99971, 0.99969, 0.99968, 0.99966,
0.99963, 0.99961, 0.99958, 0.99954, 0.99951, 0.99947, 0.99943, 0.99938,
0.99934, 0.99928, 0.99922, 0.99916, 0.99908, 0.999, 0.99891, 0.99881,
0.9987, 0.99857, 0.99843, 0.99828, 0.99812, 0.99794, 0.99776, 0.99756,
0.99735, 0.99713, 0.9969, 0.99666, 0.9964, 0.99613, 0.99583, 0.99551,
0.99515, 0.99476, 0.99432, 0.99383, 0.9933, 0.9927, 0.99205, 0.99133,
0.99053, 0.98961, 0.98852, 0.98718, 0.98553, 0.98346, 0.98089, 0.97772,
0.97391, 0.96943, 0.96429, 0.95854, 0.95221, 0.94537, 0.93805, 0.93027,
0.92202, 0.91327, 0.90393, 0.89389, 0.88304, 0.87126, 0.85846, 0.84452,
0.82935, 0.81282, 0.79485, 0.77543, 0.75458, 0.7324, 0.70893, 0.68424])  # 这里需要填入完整的生存概率数据

# 离散化正态分布
gamma00 = 0  # AR自相关系数
mew = 0     # AR常数项
sigma = 1    # 白噪声标准差
grid, weig = discretize_AR1_Tauchen(mew, gamma00, sigma, n, 2)
weig = np.diag(weig)

# 计算风险资产收益
for i in range(n):
    gret[i] = r + mu + grid[i] * sigr

# 计算外生随机性概率
for i in range(n):
    for j in range(n):
        for k in range(n):
            nweig1[i,j,k] = weig[i,i] * weig[j,j] * weig[k,k]

# Epstein-Zin参数
theta = (1.0-rho)/(1.0-1.0/psi)
psi_1 = 1.0-1.0/psi
psi_2 = 1.0/psi_1

# Part 4: 网格设置
# 风险资产持有比例网格
for i in range(na):
    ga[i] = (na-1-i)/(na-1.0)

# 投资组合收益率
for i in range(na):
    for j in range(n):
        riskret[i,j] = r*(1-ga[i]) + gret[j]*ga[i]

# cash-on-hand网格
l_maxcash = np.log(maxcash)
l_mincash = np.log(mincash)
stepcash = (l_maxcash-l_mincash)/(ncash-1)

for i in range(ncash):
    lgcash[i] = l_mincash + i*stepcash
    gcash[i] = np.exp(lgcash[i])

# Part 5: 劳动收入计算
# 工资收入shock
for i in range(n):
    grid2 = grid[i]*corr_y + grid*(1-corr_y**2)**0.5
    yh[:,i] = np.exp(grid2*smay)
    
    grid2 = grid[i]*corr_v + grid*(1-corr_v**2)**0.5
    yp[:,i] = grid2*smav

# 基础收入
for i in range(tb, tr+1):
    f_y[i-tb] = np.exp(aa + b1*i + b2*i**2 + b3*i**3)

# 工作期收入增长
for i in range(tb, tr):
    gy[i-tb] = f_y[i-tb+1]/f_y[i-tb] - 1.0
    for j in range(n):
        gyp[:,j,i-tb] = np.exp(gy[i-tb]*ones_n_1 + yp[:,j])

# 退休期收入
for i in range(tr-tb, tn-1):
    for j in range(n):
        gyp[:,j,i] = np.exp(np.zeros_like(ones_n_1))

# Part 6: 终末期计算
C[:,tn-1] = gcash
c[:,tn-1] = 1
A[:,tn-1] = 0.0

for i in range(ncash):
    V[i,tn-1] = C[i,tn-1]*((1.0-delta)**(psi/(psi-1.0)))

# Part 7: 值函数迭代
# 退休期
for t in range(tn-2, tr-tb-1, -1):
    print(t)
    
    for i3 in range(ncash):
        maxc = 1
        minc = 0
        stepc = (maxc-minc)/(nc-1)
        for i9 in range(nc):
            gc[i9] = minc + i9*stepc
            
        for i4 in range(nc):
            u = (1.0-delta)*((gc[i4]*gcash[i3])**psi_1)
            sav = gcash[i3]*(1-gc[i4])
            
            for i5 in range(na):
                auxVV = 0.0
                for i8 in range(n):
                    cash_1 = riskret[i5,i8]*sav + ret_fac
                    cash_1 = max(min(cash_1,gcash[ncash-1]),gcash[0])
                    
                    # 插值计算下期值函数
                    int_V = np.interp(cash_1, gcash, V[:,t+1])
                    auxVV += weig[i8,i8]*survprob[t]*(int_V**(1.0-rho))
                
                auxV[i5,i4] = (u + delta*(auxVV**(1.0/theta)))**psi_2
                
        vec_V = auxV.reshape(-1)
        V[i3,t], pt = np.max(vec_V), np.argmax(vec_V)
        aux2 = pt//na
        
        C[i3,t] = gc[aux2]*gcash[i3]
        c[i3,t] = gc[aux2] 
        A[i3,t] = ga[pt-na*aux2]

# 优化版本的工作期值函数迭代
for t in range(tr-tb-1, -1, -1):
    print(t)
    
    for i3 in range(ncash):
        maxc = 1
        minc = 0
        stepc = (maxc-minc)/(nc-1)
        for i9 in range(nc):
            gc[i9] = minc + i9*stepc
            
        for i4 in range(nc):
            u = (1.0-delta)*((gc[i4]*gcash[i3])**psi_1)
            sav = gcash[i3]*(1-gc[i4])
            
            for i5 in range(na):
                auxVV = 0.0
                for i6 in range(n):
                    for i8 in range(n):
                        for i7 in range(n):
                            cash_1 = riskret[i5,i8]*sav/gyp[i6,i8,t] + yh[i7,i8]
                            cash_1 = max(min(cash_1,gcash[ncash-1]),gcash[0])
                            int_V = np.interp(cash_1, gcash, V[:,t+1])
                            auxVV += nweig1[i6,i7,i8]*survprob[t]*((int_V*gyp[i6,i8,t])**(1.0-rho))
                
                auxV[i5,i4] = (u + delta*(auxVV**(1.0/theta)))**psi_2
                
        vec_V = auxV.reshape(-1)
        V[i3,t], pt = np.max(vec_V), np.argmax(vec_V)
        aux2 = pt//na
        
        C[i3,t] = gc[aux2]*gcash[i3]
        c[i3,t] = gc[aux2]
        A[i3,t] = ga[pt-na*aux2]

# Part 8: 模拟
nsim_half = nsim // 2
for i in range(nsim_half):
    # 第一期
    eps_y = f_randn(1)
    simPY[0,i] = eps_y*smav
    simPY[0,nsim_half+i] = -eps_y*smav
    simGPY[0,i] = 1.0
    simGPY[0,nsim_half+i] = 1.0
    simTY = f_randn(1)
    simY[0,i] = np.exp(simTY*smay)
    simY[0,nsim_half+i] = np.exp(-simTY*smay)
    
    # 工作期
    for t in range(1, tr-tb):
        eps_y = f_randn(1)
        simPY[t,i] = eps_y*smav + simPY[t-1,i]
        simPY[t,nsim_half+i] = -eps_y*smav + simPY[t-1,nsim_half+i]
        simGPY[t,i] = np.exp(gy[t-1])*np.exp(simPY[t,i])/np.exp(simPY[t-1,i])
        simGPY[t,nsim_half+i] = np.exp(gy[t-1])*np.exp(simPY[t,nsim_half+i])/np.exp(simPY[t-1,nsim_half+i])
        simTY = f_randn(1)
        simY[t,i] = np.exp(simTY*smay)
        simY[t,nsim_half+i] = np.exp(-simTY*smay)

# 退休期
simY[tr-tb:,:] = ret_fac
simGPY[tr-tb:,:] = 1.0

# 模拟风险资产收益
for t in range(tn):
    for i in range(nsim_half):
        eps_r = f_randn(1)
        simR[t,i] = mu + r + sigr*eps_r
        simR[t,nsim_half+i] = mu + r - sigr*eps_r

# 计算路径
simW.fill(0.0)
for t in range(tn):
    for i in range(nsim):
        simW_Y[t,i] = simW[t,i]/simY[t,i]
        cash = simW[t,i] + simY[t,i]
        
        # 插值计算最优控制
        ic = np.searchsorted(gcash, cash)
        ic = np.clip(ic, 1, ncash-1)
        ttc = (cash - gcash[ic-1])/(gcash[ic] - gcash[ic-1])
        ttc = np.clip(ttc, 0, 1)
        
        simc[t,i] = (1-ttc)*c[ic-1,t] + ttc*c[ic,t]
        simA[t,i] = (1-ttc)*A[ic-1,t] + ttc*A[ic,t]
        simc[t,i] = np.clip(simc[t,i], 0, 1)
        simC[t,i] = cash*simc[t,i]
        
        sav = cash*(1-simc[t,i])
        simS[t,i] = np.minimum(simA[t,i]*sav, sav)
        simB[t,i] = sav - simS[t,i]
        
        if t < tn-1:
            simW[t+1,i] = (simB[t,i]*r + simS[t,i]*simR[t,i])/simGPY[t+1,i]

# 计算均值
meanC = np.mean(simC, axis=1)
meanc = np.mean(simc, axis=1)
meanY = np.mean(simY, axis=1)
meanW = np.mean(simW, axis=1)
meanS = np.mean(simS, axis=1)
meanB = np.mean(simB, axis=1)
meanWY = np.mean(simW_Y, axis=1)
meanalpha = np.mean(simA, axis=1)
meanGPY = np.mean(simGPY, axis=1)

# 绘图
plt.figure(figsize=(10,6))
plt.plot(meanc, label='c')
plt.plot(meanalpha, label='α')
plt.legend()
plt.grid(True)
plt.show()