% Huggett (1996) Replication
% Multi-Period OLG Model with Earnings and Lifetime Uncertainty
% Yanran Guo
% 8/7/2018
% 修改为使用HuggettUtils类

%--------------------------------------------------------------------
%{
In this exercise, there are several points need to be noticed

(1). Grid for k
     Number of grids can be set arbitrarily
     kMin is set to be 0 due to borrowing constraint
     kMax is set by using function kgrid_olgm

(2). Labor endowment process
     This process is calibrated by using function EarningProcess_olgm
     - Use Tauchen (1986) method to produce finte state Markov approximation 
       of the AR(1) process
     - For the initial labor endowment, approximate a normal distribution
       on the grids provided by Tauchen method

(3). Age efficiency profile 
     When calibrating labor endowment process, function 'EarningProcess_olgm'
     generates pure labor endowment states and is net of the age efficiency
     profile. 
     So we need to generate a rough linear approximation of Huggett (1996), 
     by physical age levels

%}


%% Environment Setup
clc;                     % Clear screen
clear;                   % Clear memory
close all;               % Close open windows
% addpath(genpath(pwd));   % Include subfolders of the current folder


%% Fixed Model Parameters
cS = HuggettUtils.ParameterValues_Fixed();

% Precalibrate labor endowment process
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = HuggettUtils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV);

% Age efficiency profile
% A rough linear approximation of Huggett (1996), by physical age levels
ageEffV = zeros(100, 1);
ageEffV(20:72) = [linspace(0.3, 1.5, 36-20+1), 1.5 .* ones(1, 47-37+1), ...
                  linspace(1.5, 0.2, 65-48+1), linspace(0.18, 0, 72-66+1)];
paramS.ageEffV = ageEffV(cS.age1 : cS.ageLast);


%% Solve HH Problem with Calibration 
% Step 1. Precompute the aggregate labor supply, L. Since age efficiency is
%         fixed and there is no labor supply decision, L can be precalibrated 
%         without solving household problem
% Step 2. Guess total capital KGuess and lump-sum transfer of accidental bequests TGuess
% Step 3. Given KGuess and L, compute prices R, w, and b
% Step 4. Given prices R, w, b, and given TGuess, solve HH problem and 
%         get policy function for saving
% Step 5. Given solution to HH problem, compute the aggregate saving K and 
%         accidental bequests T
% Step 6. devV = [KGuess - K; TGuess - T]
%         If devV is above the tolerance level, update guess for K and T,
%         and go back to Step 2.


% 步骤1：预计算总劳动供给
eIdxM = HuggettUtils.LaborEndowSimulation_olgm(cS, paramS);
[~, L] = HuggettUtils.LaborSupply_Huggett(eIdxM, cS, paramS);

% 初始猜测K和T的值
KGuess = 56.2229;
TGuess = 0.9039;

% 设置迭代参数
maxIter = 100;          % 最大迭代次数
tolLevel = 1e-6;        % 收敛容差
dampK = 0.3;            % K的阻尼系数
dampT = 0.3;            % T的阻尼系数
iter = 0;               % 迭代计数器

fprintf('迭代开始寻找均衡 (K*, T*)\n');
fprintf('迭代   K         T         KDev       TDev       偏差范数\n');
fprintf('------------------------------------------------------\n');

for iter = 1:maxIter
    % 步骤3：给定KGuess和L，计算价格R, w和b
    [~, R, w, b] = HuggettUtils.HHPrices_Huggett(KGuess, L, cS);
    bV = [zeros(1, cS.aR), ones(1, cS.aD - cS.aR) .* b];
    
    % 步骤4：给定价格和TGuess，解决家户问题
    [cPolM, kPolM, valueM] = HuggettUtils.HHSolution_VFI_Huggett(R, w, TGuess, bV, paramS, cS);
    
    % 步骤5：计算模型中的总资本和意外遗产
    [kHistM, ~] = HuggettUtils.HHSimulation_olgm(kPolM, cPolM, eIdxM, cS);
    
    % 计算模型中的总资本
    KModel = mean(kHistM,1) * cS.ageMassV'; 
    KModel = max(0.01, KModel);
    
    % 计算意外遗产
    kprimeHistM = [kHistM(:,2:end), zeros(size(kHistM,1),1)];
    ageDeathMassV = cS.ageMassV .* cS.d;
    TModel = (mean(kprimeHistM * R, 1) * ageDeathMassV') / (1 - cS.popGrowth);
    
    % 步骤6：计算偏差
    KDev = KGuess - KModel;
    TDev = TGuess - TModel;
    devNorm = sqrt(KDev^2 + TDev^2);
    
    % 输出当前迭代结果
    fprintf('%3d    %8.4f  %8.4f  %10.6f  %10.6f  %10.6f\n', iter, KGuess, TGuess, KDev, TDev, devNorm);
    
    % 检查收敛
    if devNorm < tolLevel
        fprintf('均衡已收敛! K* = %8.4f, T* = %8.4f\n', KGuess, TGuess);
        break;
    end
    
    % 更新猜测值（使用阻尼系数更新）
    KGuess = KGuess - dampK * KDev;
    TGuess = TGuess - dampT * TDev;
end

if iter == maxIter && devNorm >= tolLevel
    fprintf('警告：达到最大迭代次数，但未收敛。\n');
end

% 使用最终的均衡值计算最终结果
K = KGuess;
T = TGuess;

% 获取家户面临的价格
[~, R, w, b] = HuggettUtils.HHPrices_Huggett(K, L, cS);
bV = [zeros(1, cS.aR), ones(1, cS.aD - cS.aR) .* b];

% 再次计算最终的政策函数
[cPolM, kPolM, valueM] = HuggettUtils.HHSolution_VFI_Huggett(R, w, T, bV, paramS, cS);


%% Figure
% Figure 1. 
% Policy function -- optimal consumption
figure;
plot(cS.kGridV, cPolM(:, 1, 23), 'k', 'LineWidth', 1);
hold on
plot(cS.kGridV, cPolM(:, 10, 23), 'k--', 'LineWidth', 1);
plot(cS.kGridV, cPolM(:, 18, 23), 'r', 'LineWidth', 1);
hold off
xlabel('wealth k');
ylabel('c(k,e,23)');
title('optimal consumption at model age 23')
legend('e=1', 'e=10', 'e=18', 'Location', 'northwest');

% Figure 2.
% Policy function -- optimal capital level
figure;
plot(cS.kGridV, kPolM(:, 1, 23), 'k', 'LineWidth', 1);
hold on
plot(cS.kGridV, kPolM(:, 10, 23), 'k--', 'LineWidth', 1);
plot(cS.kGridV, kPolM(:, 18, 23), 'r', 'LineWidth', 1);
hold off
xlabel('wealth k');
ylabel('kprime(k,e,23)');
title('optimal capital level at model age 23')
legend('e=1', 'e=10', 'e=18', 'Location', 'northwest');


%% 打印最终的均衡结果
fprintf('\n最终均衡结果\n');
fprintf('总资本K = %8.4f\n', K);
fprintf('意外遗产T = %8.4f\n', T);
fprintf('利率r = %8.4f\n', R - 1);
fprintf('工资w = %8.4f\n', w);
fprintf('社会保障福利b = %8.4f\n', b); 