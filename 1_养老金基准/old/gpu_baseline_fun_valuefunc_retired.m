function v = gpu_baseline_fun_valuefunc_retired(x, cash, fund, nextV, gret, rf, ret_fac, pension_pay, gamma, beta, psi_1, psi_2, theta, gcash, gfund, survprob, weig, t)
% 退休期间的值函数计算（GPU版本）
%
% 参数:
% x - 决策变量 [消费, 风险资产投资比例]
% cash - 手中现金
% fund - 养老基金账户余额
% nextV - 下期值函数
% gret - 随机收益率
% rf - 无风险收益率
% ret_fac - 固定基本养老金
% pension_pay - 个人养老金给付
% gamma - 相对风险规避系数
% beta - 贴现因子
% psi_1, psi_2, theta - Epstein-Zin效用函数参数
% gcash - 现金网格
% gfund - 养老基金网格
% survprob - 存活概率
% weig - 随机状态权重
% t - 当前时期（用于缓存，可选）
%
% 返回:
% v - 函数值

% 确保输入是GPU数组
if ~isgpuarray(x)
    x = gpuArray(x);
end
if ~isgpuarray(cash)
    cash = gpuArray(cash);
end
if ~isgpuarray(fund)
    fund = gpuArray(fund);
end
if ~isgpuarray(nextV)
    nextV = gpuArray(nextV);
end
if ~isgpuarray(gret)
    gret = gpuArray(gret);
end
if ~isgpuarray(gcash)
    gcash = gpuArray(gcash);
end
if ~isgpuarray(gfund)
    gfund = gpuArray(gfund);
end
if ~isgpuarray(weig)
    weig = gpuArray(weig);
end

% 计算储蓄
sav = cash * (1 - x(1));

% Epstein-Zin效用函数
u = (cash * x(1) + 0.00001).^psi_1;  % CRRA效用

% 下期的cash-on-hand
cash_1 = (rf * (1 - x(2)) + gret .* x(2)) .* sav + pension_pay;

% 使用MATLAB的一维插值（GPU版本）
try
    int_V = interp1(gcash, nextV, cash_1, 'spline');
catch
    % 如果GPU插值失败，尝试在CPU上进行插值
    int_V = interp1(gather(gcash), gather(nextV), gather(cash_1), 'spline');
    int_V = gpuArray(int_V);
end
int_V = max(int_V, 1e-20); % 确保插值结果为正

% 计算期望值
auxVV = weig' * (survprob * int_V.^(1-gamma));

% Epstein-Zin递归效用
v = -(((1-beta)*u + beta*((auxVV)^(1/theta)))^psi_2);

% 如果是用于网格搜索，确保返回标量
if isscalar(v)
    v = gather(v);
end
end