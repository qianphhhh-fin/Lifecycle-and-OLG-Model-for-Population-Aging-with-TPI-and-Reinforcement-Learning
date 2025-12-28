% --- START OF FILE test_ss0_market_clearing.m ---
%
% 描述:
%   本脚本专门用于测试初始稳态求解器 (solve_steady_state_with_fund) 的
%   市场出清特性。它会在求解器的每一次迭代中，详细计算并打印出
%   商品市场 (Y-C-I-G) 的残差，以确保在收敛过程中，宏观会计恒等式
%   始终得到满足。
%
% 使用方法:
%   1. 确保此文件与 main_olg_v10_utils.m 在同一目录下。
%   2. 直接在MATLAB命令窗口运行此脚本: >> test_ss0_market_clearing
%
% ---

clear; close all;
addpath(pwd);
fprintf('=== 初始稳态求解器市场出清检验脚本 ===\n');

%% 1. 初始化所有必要的参数 (与 main_olg_v10.m 完全一致)
fprintf('\n--- 1. 初始化参数 ---\n');
cS = main_olg_v10_utils.ParameterValues_HuggettStyle();
cS.nSim = 5000;
paramS = struct();

% --- V10 固定的政策参数 ---
cS.rho_prime_payg = 0.4; 
cS.theta_payg = 0.15; 
B_p_Y_ratio_target = 0.05; % 目标：初始稳态养老金基金占GDP的5%

% --- 家庭偏好参数 ---
cS.sigma = 3.0;
cS.beta = 0.97;
cS.cFloor = 0.05;
cS.nSim = 1000;
cS.phi_bequest = 3; 
cS.sigma_bequest = cS.sigma;

% --- 生产技术参数 ---
cS.A = 2;
cS.alpha = 0.36;
cS.ddk = 0.3; % 1 - (1 - 0.08)^5;

% --- 政府财政参数 ---
cS.tau_k = 0.20;
cS.tau_c = 0.02;
cS.tau_l = 0.10;
cS.gov_exp_frac_Y = 0.15;
cS.gov_debt_frac_Y = 0.60;

% --- 网格参数 ---
ngrid = 50;
cS.nk = ngrid;
cS.nkpps = ngrid;
cS.nkprime = ngrid;
cS.npps = ngrid;
    
% --- PPS制度参数 (在初始稳态中会被禁用) ---
cS.pps_active = true;
cS.pps_tax_rate_withdrawal = 0.03;
cS.pps_return_rate_premium = 0.01;
cS.pps_withdrawal_rate = 0.15;
cS.pps_in_K = true;
cS.pps_bequeathable = true;
cS.pps_contrib_limit = 9999;
cS.pps_max_contrib_frac = 0.1;

% 重新生成网格
cS = main_olg_v10_utils.generateGrids(cS);
fprintf('参数初始化完成。\n');

%% 2. 计算外生路径 (仅为获取初始人口分布)
fprintf('\n--- 2. 计算外生路径 ---\n');
popS = main_olg_v10_utils.initPopulation(cS);
popS = main_olg_v10_utils.populationDynamics(popS, cS);
Z_path_norm = popS.ageDist; % 使用 ageDist 获取归一化后的人口分布
Z_0 = Z_path_norm(:, 1); % 初始稳态人口分布

[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v10_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));
paramS.ageEffV_new = cS.ageEffV_new;
eIdxM = main_olg_v10_utils.LaborEndowSimulation_olgm_AgeGroup(cS, paramS);
fprintf('人口和劳动冲击路径计算完成。\n');

%% 3. [核心] 带有市场出清检验的稳态求解过程
fprintf('\n--- 3. 开始求解初始稳态 (带迭代过程中的市场出清检验) ---\n');

% --- 求解器参数 ---
K_guess = 1.8816; % 初始猜测总资本
max_iter = 100;
tol = 1e-4;
damp = 0.5;
eq_found = false;

% --- 打印表头 ---
fprintf('\n%s\n', repmat('=', 1, 110));
fprintf('%5s | %12s | %12s | %12s | %12s | %12s | %12s | %15s\n', ...
    'Iter', 'K_guess', 'K_model', 'K_pvt', 'B_p', 'K_pps', 'K_error', 'Goods Mkt Resid');
fprintf('%s\n', repmat('-', 1, 110));

% --- 迭代循环 ---
for iter = 1:max_iter
    % --- 步骤 1: 基于 K_guess 计算当期价格和宏观量 ---
    paramS_ss = paramS;
    paramS_ss.ageMassV = Z_0;
    M_ss_prices = main_olg_v10_utils.get_prices_at_t(K_guess, Z_0, cS, paramS_ss, eIdxM);

    % --- 步骤 2: 计算当期政策变量并构建家庭面临的环境 M_ss ---
    M_ss = M_ss_prices;
    mass_workers_ss = sum(Z_0(1:cS.aR_new));
    avg_wage_ss = M_ss.w_t * M_ss.L_t / mass_workers_ss;
    M_ss.b_t = cS.rho_prime_payg * avg_wage_ss;
    M_ss.current_t = -1; % 标记为稳态求解

    % --- 步骤 3: 计算目标基金规模 ---
    B_p_target = B_p_Y_ratio_target * M_ss.Y_t;

    % --- 步骤 4: 求解家庭问题 (VFI) 并聚合得到下一期资本和当期消费 ---
    % 标记这是初始稳态求解，以禁用PPS
    paramS_ss.is_initial_steady_state = true;
    
    % [关键] simulate_private_capital_forward 同时返回下一期资本和当期消费
    [K_pvt_model, K_pps_model, C_model, ~] = main_olg_v10_utils.simulate_private_capital_forward(M_ss, Z_0, cS, paramS_ss, eIdxM);
    
    % --- 步骤 5: 计算模型的总资本需求 ---
    K_total_model = K_pvt_model + B_p_target + K_pps_model;
    
    % --- 步骤 6: 计算资本市场误差 ---
    K_error = K_guess - K_total_model;
    
    % --- 步骤 7: [核心检验] 计算商品市场出清残差 ---
    % a. 计算当期总产出 Y
    Y_t = M_ss.Y_t;
    
    % b. 获取当期总消费 C
    C_t = C_model;
    
    % c. 计算当期政府支出 G
    G_t = cS.gov_exp_frac_Y * Y_t;
    
    % d. 计算当期总投资 I
    %    在稳态下，投资等于资本折旧
    %    I = δ * K
    Investment_t = cS.ddk * K_guess; % 投资必须弥补当前资本存量(K_guess)的折旧
    
    % e. 计算商品市场残差
    goods_market_residual = Y_t - C_t - Investment_t - G_t;
    
    % --- 步骤 8: 打印所有检验结果 ---
    fprintf('%5d | %12.4f | %12.4f | %12.4f | %12.4f | %12.4f | %12.3e | %15.4e\n', ...
        iter, K_guess, K_total_model, K_pvt_model, B_p_target, K_pps_model, K_error, goods_market_residual);

    % --- 步骤 9: 检查收敛 ---
    if abs(K_error) < tol
        fprintf('\n✅ 资本市场均衡收敛！\n');
        eq_found = true;
        
        % 最终检查商品市场
        if abs(goods_market_residual) < 1e-5 % 使用一个较小的容忍度
            fprintf('✅ 商品市场在均衡点近似出清 (残差: %.4e)。\n', goods_market_residual);
        else
            fprintf('⚠️ 警告: 资本市场收敛，但商品市场残差较大 (残差: %.4e)。模型可能存在会计不一致问题。\n', goods_market_residual);
        end
        break;
    end
    
    % --- 步骤 10: 更新猜测 ---
    K_guess = (1 - damp) * K_guess + damp * K_total_model;
end

if ~eq_found
    fprintf('\n❌ 稳态均衡未在最大迭代次数内收敛。\n');
end

fprintf('\n--- 检验脚本执行完毕 ---\n');

% --- END OF FILE test_ss0_market_clearing.m ---