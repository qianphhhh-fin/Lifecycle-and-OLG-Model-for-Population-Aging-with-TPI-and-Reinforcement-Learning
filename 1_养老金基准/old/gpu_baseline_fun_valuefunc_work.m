function v = gpu_baseline_fun_valuefunc_work(x, cash, fund, gyp, nextV, yh, gret, rf, rp, tau, gamma, beta, psi_1, psi_2, theta, gcash, gfund, survprob, n, weig, t)
% 工作期间的值函数计算（GPU版本）
%
% 参数:
% x - 决策变量 [消费, 风险资产投资比例, 个人养老金购买比例]
% cash - 手中现金
% fund - 养老基金账户余额
% gyp - 劳动收入增长率
% nextV - 下期值函数
% yh - 暂时收入冲击
% gret - 随机收益率
% rf - 无风险收益率
% rp - 个人养老基金收益率
% tau - 用于缴纳基本养老金的工资比例
% gamma - 相对风险规避系数
% beta - 贴现因子
% psi_1, psi_2, theta - Epstein-Zin效用函数参数
% gcash - 现金网格
% gfund - 养老基金网格
% survprob - 存活概率
% n - 随机状态数量
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
if ~isgpuarray(gyp)
    gyp = gpuArray(gyp);
end
if ~isgpuarray(nextV)
    nextV = gpuArray(nextV);
end
if ~isgpuarray(yh)
    yh = gpuArray(yh);
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

% Epstein-Zin效用函数
u = (cash * x(1) + 1e-07)^psi_1;  % CRRA效用

% 计算个人养老金购买
pension_contrib = x(3) * cash * (1 - x(1));  % 个人养老金购买

% 获取维度
n1 = min(n, size(gyp, 1));  % z_t维度
n2 = min(n, size(gyp, 2));  % eta_t维度
n3 = min(n, size(gyp, 3));  % epsilon_t维度

% 使用meshgrid创建网格
[I6, I7, I8] = meshgrid(1:n1, 1:n2, 1:n3);
I6 = permute(I6, [2, 1, 3]); % 调整维度顺序为(n1,n2,n3)
I7 = permute(I7, [2, 1, 3]);
I8 = permute(I8, [2, 1, 3]);

% 提取相应的gyp值
gyp_values = gyp(sub2ind(size(gyp), I6, I7, I8));

% 计算储蓄值 - 矩阵运算
sav_matrix = cash * (1 - x(1)) * (1 - x(3)) ./ gyp_values;

% 计算投资组合收益率
portfolio_return = rf * (1 - x(2)) + gret(I8) * x(2);

% 提取yh值
yh_values = yh(sub2ind(size(yh), I7, I8));

% 计算下期现金值
cash_1_matrix = portfolio_return .* sav_matrix + (1-tau) * yh_values;

% 限制cash_1的范围
cash_1_matrix = max(min(cash_1_matrix, gcash(end)), gcash(1));

% 计算下期养老基金账户余额
fund_1_matrix = (fund + pension_contrib) * rp ./ gyp_values + tau * yh_values;

% 限制fund_1的范围
fund_1_matrix = max(min(fund_1_matrix, gfund(end)), gfund(1));

% 使用meshgrid为插值准备网格
[X, Y] = meshgrid(gcash, gfund);

% 准备值函数网格
V_grid = permute(nextV, [2, 1]);

% 使用interp2进行插值
% 将3D矩阵重塑为2D以便批量插值
cash_1_flat = reshape(cash_1_matrix, [], 1);
fund_1_flat = reshape(fund_1_matrix, [], 1);

% 批量插值
try
    int_V_flat = interp2(X, Y, V_grid, cash_1_flat, fund_1_flat, 'spline');
catch
    % 如果GPU插值失败，尝试在CPU上进行插值
    int_V_flat = interp2(gather(X), gather(Y), gather(V_grid), gather(cash_1_flat), gather(fund_1_flat), 'spline');
    int_V_flat = gpuArray(int_V_flat);
end

% 重塑回3D矩阵
int_V_matrix = reshape(int_V_flat, size(cash_1_matrix));
int_V_matrix = max(int_V_matrix, 1e-20); % 确保插值结果为正

% 提取权重
weig_values = weig(sub2ind(size(weig), I6, I7, I8));

% 计算期望值
expected_value = weig_values .* survprob .* ((int_V_matrix .* gyp_values).^(1-gamma));
auxVV = sum(expected_value(:));

% Epstein-Zin递归效用
v = -(((1-beta)*u + beta*((auxVV)^(1/theta)))^psi_2);

% 如果是用于网格搜索，确保返回标量
if isscalar(v)
    v = gather(v);
end
end

% =================== 循环版本 =====================

% function v = baseline_fun_valuefunc_work(x, cash, fund, gyp, nextV, yh, gret, rf, rp, tau, gamma, beta, psi_1, psi_2, theta, gcash, gfund, survprob, n, weig, t)
% % 工作期间的值函数计算
% %
% % 参数:
% % x - 决策变量 [消费, 风险资产投资比例, 个人养老金购买比例]
% % cash - 手中现金
% % fund - 养老基金账户余额
% % gyp - 劳动收入增长率
% % nextV - 下期值函数
% % yh - 暂时收入冲击
% % gret - 随机收益率
% % rf - 无风险收益率
% % rp - 个人养老基金收益率
% % tau - 用于缴纳基本养老金的工资比例
% % gamma - 相对风险规避系数
% % beta - 贴现因子
% % psi_1, psi_2, theta - Epstein-Zin效用函数参数
% % gcash - 现金网格
% % gfund - 养老基金网格
% % survprob - 存活概率
% % n - 随机状态数量
% % weig - 随机状态权重
% % t - 当前时期（用于缓存，可选）
% %
% % 返回:
% % v - 函数值
% 
% auxVV = 0;
% 
% % Epstein-Zin效用函数
% u = (cash * x(1) + 1e-07)^psi_1;  % CRRA效用
% 
% % 计算个人养老金购买
% pension_contrib = x(3) * cash * (1 - x(1));  % 个人养老金购买
% 
% % 使用三层循环
% for i6 = 1:min(n, size(gyp, 1))  % z_t, 来自P_t=P_t-1+z_t
%     for i7 = 1:min(n, size(gyp, 2))  % eta_t, 即暂时收入冲击
%         for i8 = 1:min(n, size(gyp, 3))  % epsilon_t, 即风险收益随机项
%             % 计算储蓄值
%             sav = cash * (1 - x(1)) * (1 - x(3)) / gyp(i6, i7, i8);
% 
%             % 计算下期现金值
%             portfolio_return = rf * (1 - x(2)) + gret(i8) * x(2);
%             cash_1 = portfolio_return * sav + (1-tau) * yh(i7, i8);
% 
%             % 限制cash_1的范围
%             cash_1 = max(min(cash_1, gcash(end)), gcash(1));
% 
%             % 计算下期养老基金账户余额
%             fund_1 = (fund + pension_contrib) * rp / gyp(i6, i7, i8) + tau * yh(i7, i8);
% 
%             % 限制fund_1的范围
%             fund_1 = max(min(fund_1, gfund(end)), gfund(1));
% 
%             % 使用MATLAB的插值函数进行二维插值
%             int_V = interp2(gcash, gfund, permute(nextV, [2, 1]), cash_1, fund_1, 'spline');
%             int_V = max(int_V, 1e-20); % 确保插值结果为正
% 
%             % 累加期望值
%             auxVV = auxVV + weig(i6, i7, i8) * survprob * ((int_V * gyp(i6, i7, i8))^(1-gamma));
%         end
%     end
% end
% 
% % Epstein-Zin递归效用
% v = -(((1-beta)*u + beta*((auxVV)^(1/theta)))^psi_2);
% 
% 
% end

% ============= 测试函数 ===================
% % 测试函数，用于测试baseline_fun_valuefunc_work的功能
% fprintf('正在执行baseline_fun_valuefunc_work.m的测试...\n');
%
% % 设置测试参数
% x = [0.3, 0.4, 0.1];  % 消费比例30%，风险资产投资比例40%，个人养老金购买比例10%
% cash = 10;  % 手中现金
% fund = 5;   % 养老基金账户余额
%
% % 创建简单的测试网格
% n = 3;  % 随机状态数量
% gyp = ones(n, n, n) * 1.02;  % 劳动收入增长率
% yh = ones(n, n) * 2;  % 暂时收入冲击
% gret = [1.01, 1.05, 1.10];  % 随机收益率
%
% % 创建简单的下期值函数网格
% gcash = linspace(1, 20, 10);
% gfund = linspace(1, 10, 10);
% [X, Y] = meshgrid(gcash, gfund);
% nextV = X + Y;  % 简单的值函数，随现金和基金余额线性增长
%
% % 其他参数
% rf = 1.02;  % 无风险收益率
% rp = 1.02;  % 个人养老基金收益率
% tau = 0.08;  % 用于缴纳基本养老金的工资比例
% gamma = 2;  % 相对风险规避系数
% beta = 0.96;  % 贴现因子
% psi_1 = 1 - 1/gamma;  % Epstein-Zin效用函数参数
% psi_2 = 1/(1 - 1/gamma);  % Epstein-Zin效用函数参数
% theta = (1 - gamma)/(1 - 1/psi_1);  % Epstein-Zin效用函数参数
% survprob = 0.98;  % 存活概率
%
% % 创建权重
% weig = ones(n, n, n) / (n^3);  % 均匀权重
%
% % 调用函数
% tic;
% v = baseline_fun_valuefunc_work(x, cash, fund, gyp, nextV, yh, gret, rf, rp, tau, gamma, beta, psi_1, psi_2, theta, gcash, gfund, survprob, n, weig, 1);
% elapsed_time = toc;
%
% % 输出结果
% fprintf('测试结果:\n');
% fprintf('值函数结果: %f\n', v);
% fprintf('计算时间: %f 秒\n', elapsed_time);
%
% % 测试不同的消费比例
% fprintf('\n测试不同的消费比例对值函数的影响:\n');
% consumption_ratios = [0.1, 0.2, 0.3, 0.4, 0.5];
% values = zeros(size(consumption_ratios));
%
% for i = 1:length(consumption_ratios)
%     x_test = x;
%     x_test(1) = consumption_ratios(i);
%     values(i) = baseline_fun_valuefunc_work(x_test, cash, fund, gyp, nextV, yh, gret, rf, rp, tau, gamma, beta, psi_1, psi_2, theta, gcash, gfund, survprob, n, weig, 1);
%     fprintf('消费比例 %.1f: 值函数 = %f\n', consumption_ratios(i), values(i));
% end
%
% % 绘制结果
% figure;
% plot(consumption_ratios, values, 'o-', 'LineWidth', 2);
% title('消费比例对值函数的影响（越小越大）');
% xlabel('消费比例');
% ylabel('值函数');
% grid on;
%
% fprintf('\n测试完成!\n');
