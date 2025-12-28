# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 14:54:23 2024

@author: qianp
"""
import numpy as np
import math
import scipy as sp
import scipy.stats
from scipy.interpolate import CubicSpline
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize
import pandas as pd
import os
cdir = os.path.dirname(os.path.abspath(__file__)) # 获取当前py文件所在目录

from gymnasium.envs.registration import register
register(
    id='cocco-v0',
    entry_point='cocco_env:CoccoEnv',
)

def tauchen(N, mu, rho, sigma, m=2):
    """
    Approximate an AR1 process by a finite markov chain using Tauchen's method.

    :param N: scalar, number of nodes for Z
    :param mu: scalar, unconditional mean of process
    :param rho: scalar
    :param sigma: scalar, std. dev. of epsilons
    :param m: max +- std. devs.
    :returns: Z, N*1 vector, nodes for Z. Zprob, N*N matrix, transition probabilities
    SJB: This is a port of Martin Floden's 1996 Matlab code to implement Tauchen 1986 Economic Letters method The following comments are Floden's. Finds a Markov chain whose sample paths approximate those of the AR(1) process z(t+1) = (1-rho)*mu + rho * z(t) + eps(t+1) where eps are normal with stddev sigma.
    """
    Z = np.zeros((N, 1))
    Zprob = np.zeros((N, N))
    a = (1 - rho) * mu

    Z[-1] = m * math.sqrt(sigma**2 / (1 - (rho**2)))
    Z[0] = -1 * Z[-1]
    zstep = (Z[-1] - Z[0]) / (N - 1)

    for i in range(1, N):
        Z[i] = Z[0] + zstep * (i)

    Z = Z + a / (1 - rho)

    for j in range(0, N):
        for k in range(0, N):
            if k == 0:
                Zprob[j, k] = sp.stats.norm.cdf(
                    (Z[0] - a - rho * Z[j] + zstep / 2) / sigma
                )
            elif k == (N - 1):
                Zprob[j, k] = 1 - sp.stats.norm.cdf(
                    (Z[-1] - a - rho * Z[j] - zstep / 2) / sigma
                )
            else:
                up = sp.stats.norm.cdf((Z[k] - a - rho * Z[j] + zstep / 2) / sigma)
                down = sp.stats.norm.cdf((Z[k] - a - rho * Z[j] - zstep / 2) / sigma)
                Zprob[j, k] = up - down

    return (Z, Zprob)

def f_spline(x, y, gam):
    # 计算导数 yp1
    yp1 = x[0] ** (-gam)
    
    # 使用CubicSpline进行插值，设置bc_type='natural'以匹配边界条件
    cs = CubicSpline(x.flatten(), y.flatten(), bc_type=((1, yp1), 'natural'))
    
    # 第二导数
    y2 = cs(x, 2)  # 求第二导数

    return y2

# 定义辅助函数，用于计算值函数
def cocco_fun_valuefunc(x, cash, nextV, gret, rf, income, delta, gcash, survprob, weig, riskaversion):
    auxVV = 0
    u = (x[0]*cash)**(1- riskaversion)/(1- riskaversion) if x[0]>0 else -9999

    # # 循环计算
    # for i6 in range(n):
    #     for i7 in range(n):
    #         # 下期的cash-on-hand
    #         sav = (cash*(1 - x[0]))
    #         cash_1 = (rf * (1 - x[1]) + gret[i6] * x[1]) * sav + income[i7]
    #         cash_1 = np.clip(cash_1, gcash[0], gcash[-1])  # 限制在gcash的范围内
    #         # 插值计算
    #         int_V = interp1d(gcash.squeeze(), nextV, kind='slinear')(cash_1) # , fill_value="extrapolate"
    #         auxVV += weig[i6, i7] * survprob * int_V

    
    # 向量化计算下期cash-on-hand
    sav = cash*(1 - x[0]) 
    # 创建网格点
    gret_grid, income_grid = np.meshgrid(gret, income)   
    # 向量化计算下期现金
    cash_1 = (rf * (1 - x[1]) + gret_grid.flatten() * x[1]) * sav + income_grid.flatten()
    cash_1 = np.clip(cash_1, gcash[0], gcash[-1])   
    # 向量化插值
    int_V = interp1d(gcash.squeeze(), nextV, kind='slinear')(cash_1) # , fill_value="extrapolate"
    # 重塑权重矩阵并计算期望值
    weig_flat = weig.flatten()
    auxVV = np.sum(weig_flat * survprob * int_V)

    v = -(u + delta * auxVV)
    return v

# def cocco_fun_valuefunc_work(x, cash, gyp, nextV, yh, gret, rf, rho, beta, gcash, survprob, n, weig):
#     auxVV = 0
#     u = (x[0]*cash)**(1-3.84)/(1-3.84)   
#     for i6 in range(n):
#         for i8 in range(n):
#             for i7 in range(n):
#                 # 下期的cash-on-hand
#                 sav = (cash*(1 - x[0])) / gyp[i6, i8]
#                 cash_1 = (rf * (1 - x[1]) + gret[i8] * x[1]) * sav + yh[i7, i8]
#                 cash_1 = np.clip(cash_1, gcash[0], gcash[-1])  # 限制在gcash的范围内
#                 # 插值计算
#                 int_V = interp1d(gcash.squeeze(), nextV, kind='cubic', fill_value="extrapolate")(cash_1)
#                 auxVV += weig[i6, i7, i8] * survprob * (int_V * gyp[i6, i8])   

#     v = -(u + beta * auxVV )
#     return v

tb     = 20 #  initial starting age
td     = 100 #  death age
tn     = td-tb+1 #  total period
ncash  = 51
n      = 5

maxcash     = 50.0
mincash     = 0.25
# 基础收入f的系数
aa          = -2.170042+2.700381
b1          = 0.16818
b2          = -0.0323371/10
b3          = 0.0019704/100

# ret_fac     = 0 #  固定养老金
smay        = 0.1 #  白噪声shock的标准差
delta       = 0.97 #  Epstein-Zin的discount factor, 对应Gomes(2020)中的beta
r           = 1.015 #  无风险总收入   
mu          = 0.04 #  超额收益
sigr        = 0.5 #  风险资产收益率的标准差
riskaversion = 3 # CRRA中的风险规避系数

# 离散化连续随机分布
gamma00 = 0 # AR自相关系数
mew = 0 # AR的常数项
sigma = 1 # 白噪声的标准差
(grid,weig) = tauchen(n, mew, gamma00, sigma, m=2)
weig = np.diag(weig)

# Conditional Survival Probabilities
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

# risky (total) return
gret = np.zeros((n,1))
for i1 in range(n):
    gret[i1,0] = r+mu+grid[i1,0]*sigr;


# 外生随机性的概率
nweig1 = np.zeros((n,n))
for i6 in range(n):
    for i7 in range(n):
            nweig1[i6,i7] = weig[i6]*weig[i7]



## Grids for the State Variables and for Portfolio Rule


# cash-on-hand的grid
l_maxcash = np.log(maxcash);
l_mincash = np.log(mincash);
stepcash = (l_maxcash-l_mincash)/(ncash-1);

lgcash = np.zeros((ncash,1))
for i1 in range(ncash):
   lgcash[i1,0] = l_mincash+i1*stepcash; 

gcash = np.zeros((ncash,1))
for i1 in range(ncash):
   gcash[i1,0] = np.exp(lgcash[i1,0]);

# Labor Income (normalization)

# 初始化变量
grid2 = np.zeros((n, 1))
ones_n_1 = np.ones((n, 1))
yh = np.zeros((n, 1))
f_y = np.zeros((tn, 1))

# 工资白噪声随机项u_t的grid
yh = np.exp(grid * smay) 

# 随年龄变化的基础收入 随年龄先递增后递减
tr = 65  # 设置退休年龄为65岁
for i1 in range(tb, td+1):
    if i1 < tr:
        # 工作期间的收入
        f_y[i1 - tb, 0] = np.exp(aa + b1 * i1 + b2 * i1**2 + b3 * i1**3)/100
    else:
        # 退休后的固定收入,设为退休前最后一年收入的60%
        f_y[i1 - tb, 0] = f_y[tr - tb - 1, 0] * 0.6


# In[1]: 最后一期
# policy function
C = np.zeros((ncash,tn))
A = np.ones((ncash,tn))
V = np.zeros((ncash,tn))

# 最后一期的现金全部用来消费
C[:, -1] = 1
# 最后一期的投资全部为零
A[:, -1] = 0.0

# 最后一期的值函数（Epstein-Zin效用函数）
for i1 in range(ncash):
    V[i1, -1] = (C[i1, -1]*gcash[i1,]) **(1-riskaversion)/(1-riskaversion)




# In[2]: 退休期
for i1 in range(1, tn):
    t = tn - i1
    print(t)
    
    for i3 in range(ncash):
        x0 = [0.5, 0.2]
        lb = [0, 0]
        ub = [1, 1]

        bounds = [(lb[i], ub[i]) for i in range(len(lb))]
        res = minimize(cocco_fun_valuefunc, x0, args=(gcash[i3], V[:, t], gret, r, f_y[t] + yh, delta,\
                                                       gcash, survprob[t-1], nweig1, riskaversion),
                        bounds=bounds, method='Nelder-Mead', options={'disp': False})

        C[i3, t-1] = res.x[0]
        if res.x[0]==1:
            A[i3, t-1] = 0
        else:
            A[i3, t-1] = res.x[1]
        # print(np.mean(A[:,-1]))
        V[i3, t-1] = -cocco_fun_valuefunc(res.x, gcash[i3], V[:, t], gret, r, f_y[t] + yh, delta,\
                                           gcash, survprob[t-1], nweig1, riskaversion)
np.savez(cdir + '\\datas\\cocco.npz', A=A, V=V, C=C, gcash=gcash) 
# 将C,A,V保存为excel

# 转换为DataFrame
df_C = pd.DataFrame(C)
df_A = pd.DataFrame(A) 
df_V = pd.DataFrame(V)
# 保存为excel文件
df_C.to_excel(cdir + '\\datas\\C.xlsx', index=False)
df_A.to_excel(cdir + '\\datas\\A.xlsx', index=False)
df_V.to_excel(cdir + '\\datas\\V.xlsx', index=False)


