clear
% 三期生命周期模型数值模拟
% 基于main_3stages_PFI.m保存的结果进行数值模拟
% 所有决策变量都是比例形式：
% C - 消费比例（占当前财富）
% Q - 养老金购买比例（占当前财富）
% A - 风险资产投资比例（占总投资）

% 开始计时
simu_start_time = tic;
fprintf('开始三期生命周期模型数值模拟...\n');

% 检查结果文件是否存在
if ~exist('result_3stages_PFI/model_results.mat', 'file')
    error('未找到模型求解结果文件，请先运行main_3stages_PFI.m');
end

% 加载模型结果
fprintf('加载模型求解结果...\n');
load('result_3stages_PFI/model_results.mat');

% 设置随机种子以获得可重复的结果
rng(42);

% 数值模拟参数
nsim = 10000;  % 模拟次数
nperiods = 3;  % 时期数

% 初始化存储变量
simY = zeros(nperiods, nsim);          % 劳动收入
simC = zeros(nperiods, nsim);          % 消费
simC_pct = zeros(nperiods, nsim);      % 消费比例
simW = zeros(nperiods, nsim);          % 现金财富
simA = zeros(nperiods, nsim);          % 风险资产配置比例
simA_amt = zeros(nperiods, nsim);      % 风险资产投资金额
simB = zeros(nperiods, nsim);          % 无风险资产投资金额
simQ = zeros(nperiods, nsim);          % 养老金购买比例
simQ_amt = zeros(nperiods, nsim);      % 养老金购买金额
simF = zeros(nperiods, nsim);          % 养老金余额
simR = zeros(nperiods, nsim);          % 风险资产收益率
simS = zeros(nperiods, nsim);          % 总储蓄金额

% 均值存储变量
meanY = zeros(nperiods, 1);          % 平均劳动收入
meanC = zeros(nperiods, 1);          % 平均消费
meanC_pct = zeros(nperiods, 1);      % 平均消费比例
meanW = zeros(nperiods, 1);          % 平均现金财富
meanA = zeros(nperiods, 1);          % 平均风险资产配置比例
meanA_amt = zeros(nperiods, 1);      % 平均风险资产投资金额
meanB = zeros(nperiods, 1);          % 平均无风险资产投资金额
meanQ = zeros(nperiods, 1);          % 平均养老金购买比例
meanQ_amt = zeros(nperiods, 1);      % 平均养老金购买金额
meanF = zeros(nperiods, 1);          % 平均养老金余额
meanS = zeros(nperiods, 1);          % 平均总储蓄金额

% 1. 生成随机劳动收入
fprintf('生成随机劳动收入...\n');
labor_start_time = tic;

% 第1期和第2期的劳动收入（工作期）
for i = 1:nsim
    simY(1, i) = exp(normrnd(0, sigma_eps));
    simY(2, i) = exp(normrnd(0, sigma_eps));
end
% 第3期为退休期，收入为养老金
simY(3, :) = 0;

% 生成风险资产收益率
for t = 1:2  % 只需要前两期的投资收益率
    simR(t, :) = Rs_mean + normrnd(0, sigma_s, 1, nsim);
end

labor_time = toc(labor_start_time);
fprintf('随机劳动收入和收益率生成完成，耗时 %.2f 秒\n', labor_time);

% 2. 初始化状态变量
% 第1期初始现金为第1期劳动收入
simW(1, :) = simY(1, :);
% 第1期初始养老金余额为0
simF(1, :) = 0;

% 创建网格矩阵用于插值
[gW_mat, gF_mat] = meshgrid(gW, gF);

% 3. 开始模拟决策
fprintf('开始模拟决策...\n');
simu_policy_start_time = tic;

for t = 1:nperiods
    fprintf('模拟第 %d 期决策...\n', t);
    
    % 当前时期的策略函数
    C_t = C(:, :, t);
    A_t = A(:, :, t);
    Q_t = Q(:, :, t);
    
    % 对每个模拟样本计算决策
    for i = 1:nsim
        % 当前状态
        cash_i = simW(t, i);
        fund_i = simF(t, i);
        
        % 确保状态在网格范围内
        cash_i = min(max(cash_i, gW(1)), gW(end));
        fund_i = min(max(fund_i, gF(1)), gF(end));
        
        % 使用interp2函数获取策略（消费比例、风险资产配置比例、养老金购买比例）
        simC_pct(t, i) = interp2(gW_mat, gF_mat, C_t', cash_i, fund_i, 'linear');
        simA(t, i) = interp2(gW_mat, gF_mat, A_t', cash_i, fund_i, 'linear');
        simQ(t, i) = interp2(gW_mat, gF_mat, Q_t', cash_i, fund_i, 'linear');
        
        % 计算绝对值
        simC(t, i) = simC_pct(t, i) * cash_i;  % 消费金额
        simQ_amt(t, i) = simQ(t, i) * cash_i;  % 养老金购买金额
        
        % 计算储蓄和投资
        simS(t, i) = cash_i - simC(t, i) - simQ_amt(t, i);  % 总储蓄
        simA_amt(t, i) = simA(t, i) * cash_i;  % 风险资产投资金额
        simB(t, i) = cash_i - simC(t, i) - simQ_amt(t, i) - simA_amt(t, i);  % 无风险资产投资金额
        
        % 计算下一期状态（如果不是最后一期）
        if t < nperiods
            % 下一期养老金余额
            simF(t+1, i) = (simF(t, i) + simQ_amt(t, i)) * Rp;
            
            % 下一期现金财富
            if t == 1  % 第二期还有劳动收入
                simW(t+1, i) = simA_amt(t, i) * simR(t, i) + simB(t, i) * Rf + simY(t+1, i);
            else  % 第三期是退休期，收入只有养老金
                simW(t+1, i) = simA_amt(t, i) * simR(t, i) + simB(t, i) * Rf;
            end
        end
    end
    
    % 计算当期均值
    meanY(t) = mean(simY(t, :));
    meanC(t) = mean(simC(t, :));
    meanC_pct(t) = mean(simC_pct(t, :));
    meanW(t) = mean(simW(t, :));
    meanA(t) = mean(simA(t, :));
    meanA_amt(t) = mean(simA_amt(t, :));
    meanB(t) = mean(simB(t, :));
    meanQ(t) = mean(simQ(t, :));
    meanQ_amt(t) = mean(simQ_amt(t, :));
    meanF(t) = mean(simF(t, :));
    meanS(t) = mean(simS(t, :));
end

simu_policy_time = toc(simu_policy_start_time);
fprintf('决策模拟完成，耗时 %.2f 秒\n', simu_policy_time);

% 4. 保存模拟结果
fprintf('保存模拟结果...\n');
if ~exist('result_3stages_PFI/simulation', 'dir')
    mkdir('result_3stages_PFI/simulation');
end

% 保存MAT文件
save('result_3stages_PFI/simulation/sim_results.mat', 'simY', 'simC', 'simC_pct', ...
    'simW', 'simA', 'simA_amt', 'simB', 'simQ', 'simQ_amt', 'simF', 'simR', 'simS', ...
    'meanY', 'meanC', 'meanC_pct', 'meanW', 'meanA', 'meanA_amt', 'meanB', ...
    'meanQ', 'meanQ_amt', 'meanF', 'meanS');

% 保存CSV文件用于Excel查看
writematrix(meanY, 'result_3stages_PFI/simulation/meanY.csv');
writematrix(meanC, 'result_3stages_PFI/simulation/meanC.csv');
writematrix(meanC_pct, 'result_3stages_PFI/simulation/meanC_pct.csv');
writematrix(meanW, 'result_3stages_PFI/simulation/meanW.csv');
writematrix(meanA, 'result_3stages_PFI/simulation/meanA.csv');
writematrix(meanA_amt, 'result_3stages_PFI/simulation/meanA_amt.csv');
writematrix(meanB, 'result_3stages_PFI/simulation/meanB.csv');
writematrix(meanQ, 'result_3stages_PFI/simulation/meanQ.csv');
writematrix(meanQ_amt, 'result_3stages_PFI/simulation/meanQ_amt.csv');
writematrix(meanF, 'result_3stages_PFI/simulation/meanF.csv');
writematrix(meanS, 'result_3stages_PFI/simulation/meanS.csv');

% 5. 计算理论解析解
fprintf('计算理论解析解...\n');

% 风险溢价参数
mu = Rs_mean - Rf;  % 风险溢价

% 第3期（退休期）的解析解
% 消费比例：C_3* = (W_3 + F_3) / W_3
C3_theory = 1;  
% 风险资产配置比例：α_3* = 0
A3_theory = 0;  
% 养老金购买比例：Q_3* = 0（退休期无需购买养老金）
Q3_theory = 0;  

% 第2期（中年期）的解析解
% 均值估计值
expected_Y3 = 0;  % 第3期无劳动收入
expected_F3 = (meanF(2) + meanQ_amt(2)) * Rp;  % 预期第3期养老金余额

% 消费比例近似值（基于生命周期假说和解析解的简化版本）
denominator2 = 1 + (beta * Rf^(1-gamma))^(1/gamma);
C2_theory = 1 / denominator2;

% 风险资产配置比例：α_2* = (μ/(γ*σ_s^2)) * (W_2 - C_2* - Q_2*)/(W_2 + (F_2 + Q_2*)*R_p/R_f)
A2_theory = (mu / (gamma * sigma_s^2)) * (meanW(2) - meanC(2) - meanQ_amt(2)) / ...
    (meanW(2) + (meanF(2) + meanQ_amt(2)) * Rp / Rf);
A2_theory = min(max(A2_theory, 0), 1);  % 确保在[0,1]范围内

% 养老金购买比例估计（根据模型的简化版本）
zeta2 = ((gamma * sigma_s^2 * Rf) / (mu * (meanW(2) * Rf + meanF(2) * Rp)))^(1/(1+gamma));
Q2_numerator = ((beta * Rp^(1-gamma))^(1/gamma)) * meanW(2) - zeta2 * meanF(2);
Q2_denominator = ((beta * Rp^(1-gamma))^(1/gamma)) + zeta2;
Q2_theory = Q2_numerator / Q2_denominator / meanW(2);  % 转为比例
Q2_theory = min(max(Q2_theory, 0), 1);  % 确保在[0,1]范围内

% 第1期（青年期）的解析解
% 预期第2期劳动收入（假设为均值）
expected_Y2 = exp(0.5 * sigma_eps^2);  % E[exp(N(0,sigma_eps^2))]

% 消费比例：C_1* = Y_1/(1 + (β*R_f^(1-γ))^(1/γ) + (β*R_p^(1-γ))^(1/γ))
denominator1 = 1 + (beta * Rf^(1-gamma))^(1/gamma) + (beta * Rp^(1-gamma))^(1/gamma);
C1_theory = 1 / denominator1;

% 养老金购买比例：Q_1* = ((β*R_p^(1-γ))^(1/γ)*Y_1)/(1 + (β*R_p^(1-γ))^(-1/γ) + (β*R_f^(1-γ))^(-1/γ))
Q1_theory = ((beta * Rp^(1-gamma))^(1/gamma)) / denominator1;

% 风险资产配置比例：α_1* = μ/(γ*σ_s^2) * (Y_1 - C_1* - Q_1*)/(Y_1 + E[Y_2]/R_f + Q_1*R_p/R_f^2)
A1_numerator = meanW(1) - meanC(1) - meanQ_amt(1);
A1_denominator = meanW(1) + expected_Y2/Rf + meanQ_amt(1)*Rp/(Rf^2);
A1_theory = (mu / (gamma * sigma_s^2)) * (A1_numerator / A1_denominator);
A1_theory = min(max(A1_theory, 0), 1);  % 确保在[0,1]范围内

% 组合理论解
C_theory = [C1_theory; C2_theory; C3_theory];
A_theory = [A1_theory; A2_theory; A3_theory];
Q_theory = [Q1_theory; Q2_theory; Q3_theory];

% 6. 绘制结果图表
fprintf('绘制结果图表...\n');
if ~exist('result_3stages_PFI/simulation/figures', 'dir')
    mkdir('result_3stages_PFI/simulation/figures');
end

% 消费比例图
figure('Position', [100, 100, 800, 600]);
plot(1:nperiods, meanC_pct, 'b-o', 'LineWidth', 2, 'DisplayName', '数值模拟');
hold on;
plot(1:nperiods, C_theory, 'r--*', 'LineWidth', 2, 'DisplayName', '理论解析解');
hold off;
title('消费比例(C_t/W_t)');
xlabel('时期');
ylabel('比例');
legend('Location', 'best');
grid on;
ylim([0, 1]);
saveas(gcf, 'result_3stages_PFI/simulation/figures/consumption_ratio.png');

% 风险资产配置比例图
figure('Position', [100, 100, 800, 600]);
plot(1:nperiods, meanA, 'b-o', 'LineWidth', 2, 'DisplayName', '数值模拟');
hold on;
plot(1:nperiods, A_theory, 'r--*', 'LineWidth', 2, 'DisplayName', '理论解析解');
hold off;
title('风险资产配置比例(\alpha_t)');
xlabel('时期');
ylabel('比例');
legend('Location', 'best');
grid on;
ylim([0, 1]);
saveas(gcf, 'result_3stages_PFI/simulation/figures/risky_asset_ratio.png');

% 养老金购买比例图
figure('Position', [100, 100, 800, 600]);
plot(1:nperiods, meanQ, 'b-o', 'LineWidth', 2, 'DisplayName', '数值模拟');
hold on;
plot(1:nperiods, Q_theory, 'r--*', 'LineWidth', 2, 'DisplayName', '理论解析解');
hold off;
title('养老金购买比例(Q_t/W_t)');
xlabel('时期');
ylabel('比例');
legend('Location', 'best');
grid on;
ylim([0, 1]);
saveas(gcf, 'result_3stages_PFI/simulation/figures/pension_ratio.png');

% 策略变量对比图
figure('Position', [100, 100, 900, 600]);
subplot(1, 3, 1);
plot(1:nperiods, meanC_pct, 'b-o', 'LineWidth', 2);
hold on;
plot(1:nperiods, C_theory, 'r--*', 'LineWidth', 2);
hold off;
title('消费比例(C_t/W_t)');
xlabel('时期');
ylabel('比例');
legend('数值模拟', '理论解析解', 'Location', 'best');
grid on;
ylim([0, 1]);

subplot(1, 3, 2);
plot(1:nperiods, meanA, 'b-o', 'LineWidth', 2);
hold on;
plot(1:nperiods, A_theory, 'r--*', 'LineWidth', 2);
hold off;
title('风险资产配置比例(\alpha_t)');
xlabel('时期');
legend('数值模拟', '理论解析解', 'Location', 'best');
grid on;
ylim([0, 1]);

subplot(1, 3, 3);
plot(1:nperiods, meanQ, 'b-o', 'LineWidth', 2);
hold on;
plot(1:nperiods, Q_theory, 'r--*', 'LineWidth', 2);
hold off;
title('养老金购买比例(Q_t/W_t)');
xlabel('时期');
legend('数值模拟', '理论解析解', 'Location', 'best');
grid on;
ylim([0, 1]);

saveas(gcf, 'result_3stages_PFI/simulation/figures/policy_comparison.png');

% 状态变量图
figure('Position', [100, 100, 800, 600]);
plot(1:nperiods, meanW, 'g-o', 'LineWidth', 2, 'DisplayName', '现金财富');
hold on;
plot(1:nperiods, meanF, 'm-o', 'LineWidth', 2, 'DisplayName', '养老金余额');
plot(1:nperiods, meanC, 'r-o', 'LineWidth', 2, 'DisplayName', '消费');
plot(1:nperiods, meanY, 'b-o', 'LineWidth', 2, 'DisplayName', '收入');
hold off;
title('状态变量');
xlabel('时期');
ylabel('金额');
legend('Location', 'best');
grid on;
saveas(gcf, 'result_3stages_PFI/simulation/figures/state_variables.png');

% 总耗时
total_time = toc(simu_start_time);
fprintf('三期生命周期模型数值模拟完成，总耗时 %.2f 秒 (%.2f 分钟)\n', total_time, total_time/60);

% 创建模拟结果结构体
model_simu = struct();
model_simu.simY = simY;
model_simu.simC = simC;
model_simu.simC_pct = simC_pct;
model_simu.simW = simW;
model_simu.simA = simA;
model_simu.simA_amt = simA_amt;
model_simu.simB = simB;
model_simu.simQ = simQ;
model_simu.simQ_amt = simQ_amt;
model_simu.simF = simF;
model_simu.simR = simR;
model_simu.simS = simS;
model_simu.meanY = meanY;
model_simu.meanC = meanC;
model_simu.meanC_pct = meanC_pct;
model_simu.meanW = meanW;
model_simu.meanA = meanA;
model_simu.meanA_amt = meanA_amt;
model_simu.meanB = meanB;
model_simu.meanQ = meanQ;
model_simu.meanQ_amt = meanQ_amt;
model_simu.meanF = meanF;
model_simu.meanS = meanS;
model_simu.theory = struct('C', C_theory, 'A', A_theory, 'Q', Q_theory);

% 自定义violin图函数 (如果MATLAB没有violin包)
function violin(data)
    % 简单的替代violin plot的函数
    % 如果没有violin plot包，可以使用这个简单的boxplot替代
    try
        % 尝试使用violin plot (如果已安装)
        violinplot(data);
    catch
        % 否则使用boxplot
        boxplot(cell2mat(data));
    end
end 