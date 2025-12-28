%% ========================================================================
%  == 调试脚本: 对比 Gomes(2020) 与复刻代码在最后两期的VFI解
%  == 目标: 隔离 t=T 和 t=T-1 的计算，验证复刻代码的贝尔曼方程
%  ==       是否与 Gomes 的原始实现一致。
%  ========================================================================

clear all;
close all;
clc;

fprintf('开始进行 t=T 和 t=T-1 的VFI对比...\n\n');

%% ========================================================================
%  == Section 1: 提取并运行 Gomes(2020) 的原始逻辑
%  ========================================================================
fprintf('--- 1. 运行 Gomes(2020) 原始代码逻辑 ---\n');

% --- 1.1 设定与Gomes完全一致的参数和网格 ---
% 经济参数
gomes_rho   = 10.0;     % 对应 cS.gamma
gomes_delta = 0.97;     % 对应 cS.beta
gomes_psi   = 0.5;
gomes_r     = 1.015;    % 对应 cS.R_f
gomes_mu    = 0.04;
gomes_sigr  = 0.2;      % 对应 cS.sigma_eps
gomes_ret_fac = 0.68212;% 对应 cS.Y_bar_frac

% 模型维度
gomes_tn    = 81;       % 总期数 (100-20+1)
gomes_ncash = 51;       % 状态变量网格点数
gomes_na    = 51;       % alpha网格点数
gomes_nc    = 21;       % 消费选择网格点数
gomes_n     = 5;        % 冲击节点数

% 状态变量网格 (Cash-on-Hand)
maxcash = 200.0; mincash = 0.25;
lgcash_gomes = linspace(log(mincash), log(maxcash), gomes_ncash)';
gcash_gomes  = exp(lgcash_gomes);

% 决策变量网格 (alpha)
ga_gomes = linspace(1, 0, gomes_na)'; % Gomes的alpha网格是反向的

% 冲击离散化 (Gauss-Hermite)
grid_gomes = [-2.85697001387280; -1.35562617997427; 0.0; 1.35562617997426; 2.85697001387280];
weig_gomes = [0.01125741132772; 0.22207592200561; 0.53333333333333; 0.22207592200561; 0.01125741132772];
gret_gomes = gomes_r + gomes_mu + grid_gomes * gomes_sigr;
riskret_gomes = (1 - ga_gomes) * gomes_r + ga_gomes * gret_gomes';

% 生存概率 (仅需 t=T-1 时刻的)
survprob_gomes_T_minus_1 = 0.6809; % survprob(80,1)

% --- 1.2 计算 t=T 时的价值函数 ---
C_gomes_T = gcash_gomes;
psi_1_gomes = 1.0 - 1.0/gomes_psi;
V_gomes_T = C_gomes_T * ((1.0 - gomes_delta)^(1.0/psi_1_gomes));
secd_gomes = f_spline(gcash_gomes, V_gomes_T, gomes_ncash, 1.0); % 价值函数的二阶导

% --- 1.3 计算 t=T-1 时的策略 (精确复刻 life_cycle.m:317-360) ---
t = gomes_tn - 1;
V_gomes_T_minus_1 = zeros(gomes_ncash, 1);
A_gomes_T_minus_1 = zeros(gomes_ncash, 1);
C_gomes_T_minus_1 = zeros(gomes_ncash, 1);

theta_gomes = (1.0 - gomes_rho) / (1.0 - 1.0/gomes_psi);
psi_2_gomes = 1.0 / psi_1_gomes;

for i3 = 1:gomes_ncash
    % 消费网格搜索范围
    if i3 == 1
        maxc = C_gomes_T(i3); minc = maxc/2;
    else
        minc = C_gomes_T_minus_1(i3-1);
        mpc = 0.1;
        if i3 > 9
           mpc = max((C_gomes_T_minus_1(i3-1)-C_gomes_T_minus_1(i3-9))/(gcash_gomes(i3-1) - gcash_gomes(i3-9)), 0.1);
        end
        maxc = minc + mpc*(gcash_gomes(i3) - gcash_gomes(i3-1));
    end
    gc_gomes = linspace(minc, maxc, gomes_nc)';
    
    auxV = zeros(gomes_na, gomes_nc);
    for i4 = 1:gomes_nc % Loop over consumption choices
        u = (1.0 - gomes_delta) * (gc_gomes(i4)^psi_1_gomes);
        sav = gcash_gomes(i3) - gc_gomes(i4);
        if sav < 0, auxV(:, i4) = -inf; continue; end
        
        for i5 = 1:gomes_na % Loop over alpha choices
            auxVV = 0.0;
            for i8 = 1:gomes_n % Loop over return shocks
                cash_1 = riskret_gomes(i5, i8) * sav + gomes_ret_fac;
                cash_1 = max(min(cash_1, gcash_gomes(end)), gcash_gomes(1));
                int_V  = f_sc_splint(gcash_gomes, V_gomes_T, secd_gomes, gomes_ncash, cash_1);
                auxVV = auxVV + weig_gomes(i8) * survprob_gomes_T_minus_1 * (int_V^(1.0 - gomes_rho));
            end
            auxV(i5, i4) = (u + gomes_delta * (auxVV^(1.0/theta_gomes)))^psi_2_gomes;
        end
    end
    
    [V_gomes_T_minus_1(i3), pt] = max(auxV(:));
    [pt_alpha, pt_c] = ind2sub(size(auxV), pt);
    C_gomes_T_minus_1(i3) = gc_gomes(pt_c);
    A_gomes_T_minus_1(i3) = ga_gomes(pt_alpha);
end
fprintf('   Gomes 原始逻辑计算完成。\n\n');


%% ========================================================================
%  == Section 2: 运行复刻代码的逻辑
%  ========================================================================
fprintf('--- 2. 运行复刻代码逻辑 ---\n');

% --- 2.1 设定参数和网格 (与Gomes对齐) ---
cS.beta = gomes_delta;
cS.gamma = gomes_rho;
cS.psi = gomes_psi;
cS.R_f = gomes_r;
cS.mu = gomes_mu;
cS.Y_bar_frac = gomes_ret_fac;
cS.pi_pathV = ones(gomes_tn,1); cS.pi_pathV(t+1) = survprob_gomes_T_minus_1; % t=80, t+1=81
cS.psi_1 = 1 - 1/cS.psi;
cS.psi_2 = 1/cS.psi_1;

% 网格
cS.nW = gomes_ncash;
cS.nAlpha = gomes_na;
cS.nC = gomes_nc;
cS.nShocks = gomes_n;

cS.wGridV = gcash_gomes; % [关键] 使用完全相同的状态网格
cS.w_min = min(cS.wGridV); cS.w_max = max(cS.wGridV);
cS.alphaGridV = linspace(0, 1, cS.nAlpha)'; % 复刻代码的alpha是正向的
cS.savingsFracGridV = linspace(0.001, 0.999, cS.nC)';
cS.shockProbs = weig_gomes;
cS.R_shock_V = gomes_r + gomes_mu + grid_gomes * gomes_sigr;

% --- 2.2 计算 t=T 时的价值函数 ---
final_c_mycode = cS.wGridV + cS.Y_bar_frac;
V_mycode_T = final_c_mycode * ((1 - cS.beta)^cS.psi_2);
V_interp_mycode = griddedInterpolant(cS.wGridV, V_mycode_T, 'spline', 'spline');

% --- 2.3 计算 t=T-1 时的策略 (精确复刻 value_function_iteration_matrix_nopps 退休期逻辑) ---
V_mycode_T_minus_1 = zeros(cS.nW, 1);
A_mycode_T_minus_1 = zeros(cS.nW, 1);
C_mycode_T_minus_1 = zeros(cS.nW, 1);

for iW = 1:cS.nW
    w_t = cS.wGridV(iW);
    coh = w_t + cS.Y_bar_frac;
    
    max_v = -inf;
    best_c = coh;
    best_alpha = 0;

    for iAlpha = 1:cS.nAlpha
        alpha_t = cS.alphaGridV(iAlpha);
        for iSav = 1:cS.nC
            sav = coh * cS.savingsFracGridV(iSav);
            c_t = coh - sav;
            if c_t <= 1e-6, continue; end

            util_c_term = (1 - cS.beta) * c_t^cS.psi_1;
            
            EV_inner = 0;
            for iEps = 1:cS.nShocks
                R_port = (1-alpha_t)*cS.R_f + alpha_t*cS.R_shock_V(iEps);
                w_next = max(cS.w_min, min(cS.w_max, sav * R_port)); % Gomes的cash_1包含了下一期的确定性收入
                v_prime = V_interp_mycode(w_next);
                EV_inner = EV_inner + cS.shockProbs(iEps) * cS.pi_pathV(t+1) * v_prime^(1 - cS.gamma);
            end
            
            certainty_equiv_term = EV_inner^((1-1/cS.psi)/(1-cS.gamma));
            current_v = (util_c_term + cS.beta * certainty_equiv_term)^cS.psi_2;
            
            if real(current_v) > max_v
                max_v = real(current_v);
                best_c = c_t;
                best_alpha = alpha_t;
            end
        end
    end
    V_mycode_T_minus_1(iW) = max_v;
    C_mycode_T_minus_1(iW) = best_c;
    A_mycode_T_minus_1(iW) = best_alpha;
end
fprintf('   复刻代码逻辑计算完成。\n\n');


%% ========================================================================
%  == Section 3: 结果对比
%  ========================================================================
fprintf('--- 3. 对比结果 ---\n');

% --- 3.1 图形对比 ---
figure('Name', 'Debug: Alpha(t=T-1) Policy Function Comparison', 'Position', [100, 100, 800, 600]);
plot(gcash_gomes, A_gomes_T_minus_1, 'b-o', 'LineWidth', 2, 'DisplayName', 'Gomes Original');
hold on;
plot(cS.wGridV, A_mycode_T_minus_1, 'r--x', 'LineWidth', 2, 'DisplayName', 'My Code (Corrected)');
hold off;
title('Policy Function for \alpha at t=T-1 (Age 99)');
xlabel('Cash-on-Hand / Wealth Grid');
ylabel('Risky Share (\alpha)');
legend('location', 'best');
grid on;
xlim([min(gcash_gomes), max(gcash_gomes)]);
set(gca, 'XScale', 'log'); % 使用对数坐标轴以便观察低财富区域

% --- 3.2 数值对比 ---
comparison_table = table(gcash_gomes, A_gomes_T_minus_1, A_mycode_T_minus_1);
comparison_table.Difference = comparison_table.A_gomes_T_minus_1 - comparison_table.A_mycode_T_minus_1;

fprintf('策略函数 alpha(t=T-1) 的数值对比:\n');
disp(comparison_table);

fprintf('\n\n对比分析:\n');
fprintf('1. 请观察生成的图表。如果两条线基本重合，说明我们的逻辑已经与Gomes的实现对齐。\n');
fprintf('2. 请检查命令行输出的表格。"Difference"列如果接近于0，则验证了我们的一致性。\n');
fprintf('3. 如果仍有显著差异，请检查第2部分的 `w_next` 定义。我已根据Gomes的 `cash_1` 定义进行了调整。\n');