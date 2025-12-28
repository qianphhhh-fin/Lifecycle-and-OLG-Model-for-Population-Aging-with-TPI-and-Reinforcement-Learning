clear
% 基于main_baseline_PFI.m保存的结果进行数值模拟
% 返回模拟结果

% 开始计时
simu_start_time = tic;
fprintf('开始数值模拟...\n');

% 检查结果文件是否存在
if ~exist('result_baseline_matlab_PFI/model_results.mat', 'file')
    error('未找到模型求解结果文件，请先运行main_baseline_PFI.m');
end

% 加载模型结果
fprintf('加载模型求解结果...\n');
load('result_baseline_matlab_PFI/model_results.mat');

% 设置随机种子以获得可重复的结果
rng(42);

% 数值模拟参数
nsim = 10000;
ones_nsim_1 = ones(nsim, 1);
meanY = zeros(tn, 1);
meanC = zeros(tn, 1);
meanW = zeros(tn, 1);
meanA = zeros(tn, 1);
meanS = zeros(tn, 1);
meanB = zeros(tn, 1);
meanWY = zeros(tn, 1);
meanalpha = zeros(tn, 1);
meanGPY = zeros(tn, 1);
meanQ = zeros(tn, 1);  % 平均个人养老金购买
meanF = zeros(tn, 1);  % 平均养老基金余额
meanP = zeros(tn, 1);  % 平均个人养老金给付
cGPY = zeros(tn, 1);
meanYs = zeros(tn, 1);
meanCs = zeros(tn, 1);
meanWs = zeros(tn, 1);
simPY = zeros(tn, nsim);
simGPY = zeros(tn, nsim);
simY = zeros(tn, nsim);
simC = zeros(tn, nsim);
simW = zeros(tn, nsim);
simS = zeros(tn, nsim);
simA = zeros(tn, nsim);
simB = zeros(tn, nsim);
simW_Y = zeros(tn, nsim);
simR = zeros(tn, nsim);
simQ = zeros(tn, nsim);  % 个人养老金购买决策
simQ_pct = zeros(tn, nsim);  % 个人养老金购买比例
simF = zeros(tn, nsim);  % 养老基金账户余额
simP = zeros(tn, nsim);  % 个人养老金给付
eps_y = zeros(1, 1);
simTY = zeros(1, 1);
eps_r = zeros(1, 1);
cash = zeros(tn, nsim);
simC_pct = zeros(tn, nsim);

% 1、模拟生成labor income
fprintf('模拟生成劳动收入...\n');
labor_start_time = tic;

for i1 = 1:(nsim/2)  % 另外一半模拟完全对称
    % working period第一期
    eps_y(1, 1) = randn(1);  % N(0,1)
    simPY(1, i1) = eps_y(1, 1) * smav;  % 初始的p
    simPY(1, nsim/2 + i1) = -eps_y(1, 1) * smav;
    simGPY(1, i1) = 1.0;
    simGPY(1, nsim/2 + i1) = 1.0;
    simTY(1, 1) = randn(1);
    simY(1, i1) = exp(simTY(1, 1) * smay);
    simY(1, nsim/2 + i1) = exp(-simTY(1, 1) * smay);
    
    % working period第2期~退休
    for i2 = 2:(tr-tb)
        w = i2 + tb - 1;
        eps_y(1, 1) = randn(1);
        simPY(i2, i1) = eps_y(1, 1) * smav + simPY(i2-1, i1);
        simPY(i2, nsim/2 + i1) = -eps_y(1, 1) * smav + simPY(i2-1, nsim/2 + i1);
        simGPY(i2, i1) = exp(gy(i2-1)) * exp(simPY(i2, i1)) / exp(simPY(i2-1, i1));
        simGPY(i2, nsim/2 + i1) = exp(gy(i2-1)) * exp(simPY(i2, nsim/2 + i1)) / exp(simPY(i2-1, nsim/2 + i1));
        simTY(1, 1) = randn(1);
        simY(i2, i1) = exp(simTY(1, 1) * smay);
        simY(i2, nsim/2 + i1) = exp(-simTY(1, 1) * smay);
    end
end

% 退休期
for t = (tr-tb+1):tn
    simY(t, :) = ret_fac;
    simGPY(t, :) = 1.0;
end

% 2、模拟风险投资的收益率
fprintf('模拟风险投资的收益率...\n');
for t = 1:tn
    for i1 = 1:(nsim/2)
        eps_r(1, 1) = randn(1);
        simR(t, i1) = mu + rf + sigr * eps_r(1, 1);
        simR(t, nsim/2 + i1) = mu + rf - sigr * eps_r(1, 1);
    end
end

labor_time = toc(labor_start_time);
fprintf('模拟劳动收入和风险收益率完成，耗时 %.2f 秒\n', labor_time);

% 从第一期开始迭代，得到各控制变量的值
simW(:, :) = 0.2;
simF(:, :) = 0;
simF(1, :) = pension_pct * simY(1, :);

% 创建网格矩阵用于插值
[gW_mat, gF_mat] = meshgrid(gcash, gfund);

% 开始模拟
fprintf('开始模拟策略...\n');
simu_policy_start_time = tic;

% 第一期的初始现金
cash(1, :) = simY(1, :);

% 确定可用CPU核心数
num_cores = 14; %feature('numcores');
fprintf('使用最多 %d 个CPU核心进行并行计算\n', num_cores);

% 设置并行池
if isempty(gcp('nocreate'))
    parpool('local', num_cores);
end

% 从第一期开始迭代
for t = 1:tn
    fprintf('模拟时期 %d/%d\n', t, tn);
    
    % 当前时期的策略函数
    C_t = C(:, :, t);
    A_t = A(:, :, t);
    Q_t = Q(:, :, t);
    
    % 临时存储变量 - 为parfor准备
    temp_C_pct = zeros(1, nsim);
    temp_C = zeros(1, nsim);
    temp_A = zeros(1, nsim);
    temp_Q = zeros(1, nsim);
    temp_Q_pct = zeros(1, nsim);
    temp_S = zeros(1, nsim);
    temp_W = zeros(1, nsim);
    temp_B = zeros(1, nsim);
    temp_P = zeros(1, nsim);
    temp_F_next = zeros(1, nsim);
    temp_cash_next = zeros(1, nsim);
    
    % 对每个模拟样本并行处理
    parfor i = 1:nsim
        % 当前状态
        cash_i = cash(t, i);
        fund_i = simF(t, i);
        
        % 限制在网格范围内
        cash_i = min(max(cash_i, gcash(1)), gcash(end));
        fund_i = min(max(fund_i, gfund(1)), gfund(end));
        
        % 使用interp2进行双线性插值（与PFI文件匹配）
        temp_C_pct(i) = interp2(gW_mat, gF_mat, C_t', cash_i, fund_i, 'spline');
        temp_A(i) = interp2(gW_mat, gF_mat, A_t', cash_i, fund_i, 'spline');
        temp_Q(i) = interp2(gW_mat, gF_mat, Q_t', cash_i, fund_i, 'spline');
        
        % 计算消费和储蓄
        temp_C(i) = temp_C_pct(i) * cash_i;
        temp_S(i) = cash_i - temp_C(i);
        
        % 计算个人养老金购买
        temp_Q_pct(i) = temp_Q(i);
        Q_amt = temp_Q_pct(i) * cash_i;
        
        % 计算风险资产投资（确保与PFI逻辑一致）
        % 总投资 = 现金 - 消费 - 养老金购买
        total_investment = cash_i - temp_C(i) - Q_amt;
        temp_W(i) = temp_A(i) * total_investment;
        temp_B(i) = total_investment - temp_W(i);
        
        % 计算养老金给付（退休期）
        if t >= (tr - tb + 1)
            temp_P(i) = pension_rate * fund_i;
        else
            temp_P(i) = 0;
        end
        
        % 计算下一期的状态（如果不是最后一期）
        if t < tn
            % 下一期的养老基金余额
            if t < (tr - tb)
                % 工作期
                temp_F_next(i) = (fund_i + Q_amt) * rp / simGPY(t, i) + pension_pct * simY(t+1, i);
            else
                % 退休期
                temp_F_next(i) = (fund_i - temp_P(i)) * rp;
            end
            
            % 下一期的现金
            if t < (tr - tb)
                % 工作期
                temp_cash_next(i) = (rf * temp_B(i) + simR(t, i) * temp_W(i)) / simGPY(t, i) + (1 - pension_pct) * simY(t+1, i);
            else
                % 退休期
                temp_cash_next(i) = rf * temp_B(i) + simR(t, i) * temp_W(i) + temp_P(i);
            end
        end
    end
    
    % 将临时变量赋值回全局数组
    simC_pct(t, :) = temp_C_pct;
    simC(t, :) = temp_C;
    simA(t, :) = temp_A;
    simQ(t, :) = temp_Q;
    simQ_pct(t, :) = temp_Q_pct;
    simS(t, :) = temp_S;
    simW(t, :) = temp_W;
    simB(t, :) = temp_B;
    simP(t, :) = temp_P;
    
    % 更新下一期状态
    if t < tn
        simF(t+1, :) = temp_F_next;
        cash(t+1, :) = temp_cash_next;
    end
    
    % 计算平均值
    meanY(t) = mean(simY(t, :));
    meanC(t) = mean(simC(t, :));
    meanW(t) = mean(simW(t, :));
    meanA(t) = mean(simA(t, :));
    meanS(t) = mean(simS(t, :));
    meanB(t) = mean(simB(t, :));
    meanWY(t) = mean(simW(t, :) ./ simY(t, :));
    meanalpha(t) = mean(simW(t, :) ./ (simW(t, :) + simB(t, :) + 1e-10)); % 避免除以零
    meanGPY(t) = mean(simGPY(t, :));
    meanQ(t) = mean(simQ(t, :));
    meanF(t) = mean(simF(t, :));
    meanP(t) = mean(simP(t, :));
end

simu_policy_time = toc(simu_policy_start_time);
fprintf('模拟策略完成，耗时 %.2f 秒\n', simu_policy_time);

% 计算标准化的平均值
for t = 1:tn
    cGPY(t) = prod(meanGPY(1:t));
    meanYs(t) = meanY(t) / cGPY(t);
    meanCs(t) = meanC(t) / cGPY(t);
    meanWs(t) = meanW(t) / cGPY(t);
end

% 保存模拟结果 - 修改为使用PFI目录
fprintf('保存模拟结果...\n');
if ~exist('result_baseline_matlab_PFI/simulation', 'dir')
    mkdir('result_baseline_matlab_PFI/simulation');
end

% 保存主要结果
save('result_baseline_matlab_PFI/simulation/sim_results.mat', 'simY', 'simC', 'simW', 'simA', 'simS', 'simB', ...
    'simW_Y', 'simR', 'simQ', 'simF', 'simP', 'meanY', 'meanC', 'meanW', 'meanA', 'meanS', 'meanB', ...
    'meanWY', 'meanalpha', 'meanGPY', 'meanQ', 'meanF', 'meanP', 'cGPY', 'meanYs', 'meanCs', 'meanWs');

% 保存CSV文件用于Excel查看
writematrix(meanY, 'result_baseline_matlab_PFI/simulation/meanY.csv');
writematrix(meanC, 'result_baseline_matlab_PFI/simulation/meanC.csv');
writematrix(meanW, 'result_baseline_matlab_PFI/simulation/meanW.csv');
writematrix(meanA, 'result_baseline_matlab_PFI/simulation/meanA.csv');
writematrix(meanS, 'result_baseline_matlab_PFI/simulation/meanS.csv');
writematrix(meanB, 'result_baseline_matlab_PFI/simulation/meanB.csv');
writematrix(meanWY, 'result_baseline_matlab_PFI/simulation/meanWY.csv');
writematrix(meanalpha, 'result_baseline_matlab_PFI/simulation/meanalpha.csv');
writematrix(meanQ, 'result_baseline_matlab_PFI/simulation/meanQ.csv');
writematrix(meanF, 'result_baseline_matlab_PFI/simulation/meanF.csv');
writematrix(meanP, 'result_baseline_matlab_PFI/simulation/meanP.csv');
writematrix(meanYs, 'result_baseline_matlab_PFI/simulation/meanYs.csv');
writematrix(meanCs, 'result_baseline_matlab_PFI/simulation/meanCs.csv');
writematrix(meanWs, 'result_baseline_matlab_PFI/simulation/meanWs.csv');

% 绘制结果图表
fprintf('绘制结果图表...\n');
if ~exist('result_baseline_matlab_PFI/simulation/figures', 'dir')
    mkdir('result_baseline_matlab_PFI/simulation/figures');
end

% 年龄轴
age = (tb:td)';

% 图1：状态变量 - 收入、消费、财富和养老基金余额
figure('Position', [100, 100, 900, 600]);
subplot(1, 1, 1);
plot(age, meanYs, 'b-', 'LineWidth', 2);
hold on;
plot(age, meanCs, 'r-', 'LineWidth', 2);
plot(age, meanWs, 'g-', 'LineWidth', 2);
plot(age, meanF, 'm-', 'LineWidth', 2);
xline(tr, 'k--', 'LineWidth', 1.5);
hold off;
title('状态变量随年龄变化');
xlabel('年龄');
ylabel('标准化值');
legend('收入', '消费', '财富', '养老基金余额', '退休年龄', 'Location', 'best');
grid on;
saveas(gcf, 'result_baseline_matlab_PFI/simulation/figures/state_variables.png');

% 图2：策略变量 - 消费比例、风险资产配置比例、养老金购买比例
figure('Position', [100, 100, 900, 600]);
subplot(1, 1, 1);

% 创建左侧y轴
yyaxis left;
plot(age, mean(simC_pct, 2), 'b-', 'LineWidth', 2);
ylabel('消费比例');
ylim([0, 1]);

% 创建右侧y轴
yyaxis right;
plot(age, meanalpha, 'r-', 'LineWidth', 2);
hold on;
plot(age, meanQ, 'g-', 'LineWidth', 2);
xline(tr, 'k--', 'LineWidth', 1.5);
hold off;
ylabel('比例');
ylim([0, 1]);

title('策略变量随年龄变化');
xlabel('年龄');
legend('消费比例', '风险资产配置比例', '养老金购买比例', '退休年龄', 'Location', 'best');
grid on;
saveas(gcf, 'result_baseline_matlab_PFI/simulation/figures/policy_variables.png');

% 图3：养老金给付
figure('Position', [100, 100, 800, 600]);
plot(age, meanP, 'b-', 'LineWidth', 2);
hold on;
xline(tr, 'k--', 'LineWidth', 1.5);
hold off;
title('养老金给付随年龄变化');
xlabel('年龄');
ylabel('养老金给付');
grid on;
saveas(gcf, 'result_baseline_matlab_PFI/simulation/figures/pension_payment.png');

% 总耗时
total_time = toc(simu_start_time);
fprintf('数值模拟完成，总耗时 %.2f 秒 (%.2f 分钟)\n', total_time, total_time/60);

% 返回模拟结果
model_simu = struct();
model_simu.simY = simY;
model_simu.simC = simC;
model_simu.simW = simW;
model_simu.simA = simA;
model_simu.simS = simS;
model_simu.simB = simB;
model_simu.simW_Y = simW_Y;
model_simu.simR = simR;
model_simu.simQ = simQ;
model_simu.simF = simF;
model_simu.simP = simP;
model_simu.meanY = meanY;
model_simu.meanC = meanC;
model_simu.meanW = meanW;
model_simu.meanA = meanA;
model_simu.meanS = meanS;
model_simu.meanB = meanB;
model_simu.meanWY = meanWY;
model_simu.meanalpha = meanalpha;
model_simu.meanGPY = meanGPY;
model_simu.meanQ = meanQ;
model_simu.meanF = meanF;
model_simu.meanP = meanP;
model_simu.cGPY = cGPY;
model_simu.meanYs = meanYs;
model_simu.meanCs = meanCs;
model_simu.meanWs = meanWs;
