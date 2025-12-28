%% Main script for the Lifecycle Model - GOMES REPLICATION
%  This script sets up the model, solves it using value function iteration,
%  simulates lifecycle paths, and plots the results.

clear all;
close all;
clc;

fprintf('Starting the lifecycle model solution and simulation (Gomes Replication)...\n\n');

% -------------------------------------------------------------------------
% 1. Set Model Parameters
% -------------------------------------------------------------------------
fprintf('1. Setting model parameters...\n');
tic;
% =================================================================
% == 功能: 设置模型所有参数和网格 (V2.1 - 使用Gomes的离散化)
% == 核心: 所有随机过程统一使用 Gomes(2020) 的Gauss-Hermite节点
% =================================================================

% --- 1. 生命周期与人口结构 ---
cS.tb_age = 20;         % 工作起始年龄
cS.tr_age = 66;         % 退休年龄
cS.td_age = 100;        % 最高寿命
cS.tb = 1;              % 模型起始期
cS.tr = cS.tr_age - cS.tb_age + 1; % 模型退休期
cS.tn = cS.td_age - cS.tb_age + 1; % 模型总期数

% 条件生存概率
survprob_data = [0.99845, 0.99839, 0.99833, 0.9983, 0.99827, 0.99826, 0.99824, 0.9982, 0.99813, 0.99804, 0.99795, 0.99785, 0.99776, 0.99766, 0.99755, 0.99743, 0.9973, 0.99718, 0.99707, 0.99696, 0.99685, 0.99672, 0.99656, 0.99635, 0.9961, 0.99579, 0.99543, 0.99504, 0.99463, 0.9942, 0.9937, 0.99311, 0.99245, 0.99172, 0.99091, 0.99005, 0.98911, 0.98803, 0.9868, 0.98545, 0.98409, 0.9827, 0.98123, 0.97961, 0.97786, 0.97603, 0.97414, 0.97207, 0.9697, 0.96699, 0.96393, 0.96055, 0.9569, 0.9531, 0.94921, 0.94508, 0.94057, 0.9357, 0.93031, 0.92424, 0.91717, 0.90922, 0.90089, 0.89282, 0.88503, 0.87622, 0.86576, 0.8544, 0.8423, 0.82942, 0.8154, 0.80002, 0.78404, 0.76842, 0.75382, 0.73996, 0.72464, 0.71057, 0.6961, 0.6809];
cS.pi_pathV = ones(cS.tn, 1);
len_surv = length(survprob_data);
cS.pi_pathV(2:min(cS.tn, len_surv+1)) = survprob_data(1:min(cS.tn-1, len_surv));

% --- 2. 偏好 (Epstein-Zin) ---
cS.beta = 0.97;
cS.gamma = 10.0;
cS.psi = 0.5;
cS.theta = (1.0 - cS.gamma) / (1.0 - 1.0/cS.psi);
cS.psi_1 = 1.0 - 1.0/cS.psi;
cS.psi_2 = 1.0 / cS.psi_1;

% --- 3. 收入过程 ---
aa = -2.170042 + 2.700381; b1 = 0.16818; b2 = -0.0323371/10; b3 = 0.0019704/100;
g_t = zeros(cS.tr, 1);
for t = 1:cS.tr
    age = t + cS.tb_age - 1;
    g_t(t) = aa + b1*age + b2*age^2 + b3*age^3;
end
cS.G_pathV = ones(cS.tn, 1);
cS.G_pathV(1:cS.tr-1) = exp(g_t(2:end) - g_t(1:end-1));
cS.sigma_u = 0.1;
cS.sigma_z = 0.1;

% --- 4. 养老金与税收 ---
cS.tau_y = 0.06;
cS.tau_q = 0.03;
cS.Q_max = 0.1;
cS.Y_bar_frac = 0.68212;

% --- 5. 资产与回报率 ---
cS.R_f = 1.015;
cS.mu = 0.04;
cS.sigma_eps = 0.2;
cS.R_p = cS.R_f;

% --- 6. 数值求解网格 ---
cS.nW = 51;
cS.nF = 31;
cS.w_min = 0.25; cS.w_max = 200;
cS.f_min = 0.0; cS.f_max = 100;
cS.wGridV = exp(linspace(log(cS.w_min), log(cS.w_max), cS.nW))'-cS.Y_bar_frac;

cS.fGridV = exp(linspace(log(cS.f_min+1e-6), log(cS.f_max), cS.nF-1))';
cS.fGridV = [0; cS.fGridV];
cS.nC = 21;
cS.nQ = 21;
cS.nAlpha = 51;
cS.alphaGridV = linspace(0, 1, cS.nAlpha)';
cS.qGridV = linspace(0, cS.Q_max, cS.nQ)';
cS.savingsFracGridV = linspace(0.001, 0.999, cS.nC)';

% --- [!!! 核心修改: 使用 Gomes(2020) 的 Gauss-Hermite 节点 !!!] ---
cS.nShocks = 5;

% 直接从 gomes2020.m 复制节点和权重
shock_nodes_std = [-2.85697001387280; -1.35562617997427; 0.0; 1.35562617997426; 2.85697001387280];
shock_probs_std = [0.01125741132772; 0.22207592200561; 0.53333333333333; 0.22207592200561; 0.01125741132772];

% 确保概率和为1
cS.shockProbs = shock_probs_std / sum(shock_probs_std);

% 根据各自的标准差缩放标准正态节点
cS.zNodes = shock_nodes_std * cS.sigma_z;
cS.uNodes = shock_nodes_std * cS.sigma_u;
cS.epsNodes = shock_nodes_std * cS.sigma_eps;

% 基于离散化的冲击节点计算风险回报率
cS.R_shock_V = cS.R_f + cS.mu + cS.epsNodes;
% --- [!!! 修正结束 !!!] ---

fprintf('   Parameters set. Elapsed time: %.2f seconds.\n\n', toc);


% -------------------------------------------------------------------------
% 2. Solve the Model using Value Function Iteration (No PPS)
% -------------------------------------------------------------------------
fprintf('2. Solving the household problem via VFI (No PPS)...\n');
tic;
% [polS, valS] = utils.value_function_iteration_matrix(cS); % <-- Comment out original
[polS, valS] = utils.value_function_iteration_matrix_nopps(cS); % <-- Use new function
fprintf('   VFI solved. Elapsed time: %.2f seconds.\n\n', toc);

% Save the results to avoid re-computation
output_dir = 'debug';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
save(fullfile(output_dir, 'vfi_results_gomes_replication.mat'), 'polS', 'valS', 'cS');
fprintf('   VFI results saved to %s.\n\n', fullfile(output_dir, 'vfi_results_gomes_replication.mat'));


% -------------------------------------------------------------------------
% 3. Simulate Lifecycle Paths (No PPS)
% -------------------------------------------------------------------------
fprintf('3. Simulating lifecycle paths for a panel of households (No PPS)...\n');
tic;
% simS = utils.simulation(polS, cS); % <-- Comment out original
simS = utils.simulation_nopps(polS, cS); % <-- Use new function
fprintf('   Simulation complete. Elapsed time: %.2f seconds.\n\n', toc);


% -------------------------------------------------------------------------
% 4. Plot Simulation Results
% -------------------------------------------------------------------------
fprintf('4. Generating and saving result plots...\n');
tic;
utils.plot_results(simS, cS);
fprintf('   Plots generated and saved to the /fig directory. Elapsed time: %.2f seconds.\n\n', toc);

fprintf('Lifecycle model execution (Gomes Replication) finished successfully.\n');