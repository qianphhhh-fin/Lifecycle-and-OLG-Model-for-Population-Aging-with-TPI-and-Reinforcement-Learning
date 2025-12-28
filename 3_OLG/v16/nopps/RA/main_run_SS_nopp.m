% =========================================================================
% == SCRIPT: main_run_SS_nopps.m
% == 版本: [v2.2 - 最终检验版]
% ==
% == 核心修正:
% ==   - 采用最终的、理论一致的检验逻辑。
% =========================================================================
clear; close all; clear classes;
addpath(pwd);
fprintf('=== OLG模型 (无PPS情景) 稳态求解与数据准备脚本 ===\n\n');

%% --- 1. 全局设置与初始化 ---
fprintf('--- 1. 全局设置与初始化 ---\n');

ngrid = 100;       % 常规资产网格密度
ngrid_pps = 1;    % 无PPS模式下，此参数无效但为保持结构完整性而设置

% 加载模型物理参数
fprintf('   加载模型物理参数...\n');
cS = model_setup_utils_bgp.ParameterValues();

% --- [v1.7 核心重构] 在顶层统一生成全局共享的参数 ---
fprintf('   设定并生成全局共享的资产网格...\n');
cS.nk = ngrid; 
cS.nkpps = ngrid_pps; 
cS.nkprime = ngrid; 
cS.npps = ngrid_pps;
cS = model_setup_utils_bgp.generateGrids(cS);
fprintf('   ✅ 全局网格(kGridV等)已生成并存入cS。\n');
% --- [重构结束] ---

% --- 设定全局政策情景开关 ---
cS.endogenous_theta_mode = true;  % true: DB模式(内生theta); false: DC模式(外生theta)
cS.pps_active = false;            % **[关键]** 全局强制关闭PPS制度

% 设定初始年份和模拟起点
cS.ss0_year = 2023;
cS.start_year = 2023;

% --- 如果是DB模式，在此处设定全局的政策参数 ---
if cS.endogenous_theta_mode
    cS.payg_replacement_rate = 0.2; % 目标替代率
    cS.theta_max = 0.99;             % 缴费率上限
end


%% --- 2. [核心] 生成完整的外生路径 ---
fprintf('\n--- 2. 生成完整的外生路径 (内生决定模拟长度) ---\n');
[Z_path,Z_path_raw, A_path, cS] = model_setup_utils_bgp.generate_exo_paths(cS, false); % false: 不画人口图
cS = model_setup_utils_bgp.calcaulte_theta_payg_path(cS, false);
theta_path = cS.theta_path;

ss0 = struct(); % 初始化以防求解失败
ssF = struct(); % 初始化以防求解失败

%% --- 3. 求解初始稳态 (ss0, t=0) ---
% 此部分保持注释

%% --- 4. 求解终期稳态 (ssF, t=T) ---
fprintf('\n\n--- 4. 求解终期BGP稳态 (ssF, t=T) ---\n');

cSF = cS; % 创建一个独立的副本
cSF.pps_active = false; 
cSF.s_pathV = cS.s_pathV(:,end); 
% cSF.n_ss= -0.005;
fprintf('   特征: 长期稳态人口(n=%.3f), 长期增长(g=%.3f), **PPS激活=false**\n', cSF.n_ss, cSF.g_A_ss);

if cSF.endogenous_theta_mode
    fprintf('   [养老金模式] DB模式激活, 目标替代率=%.2f%%, 上限=%.2f%%\n', cSF.payg_replacement_rate*100, cSF.theta_max*100);
    theta_for_ssF = [];
else
    fprintf('   [养老金模式] DC模式激活, 外生缴费率=%.4f\n', theta_path(end));
    theta_for_ssF = theta_path(end);
end

paramSF = struct();
[paramSF.leGridV, paramSF.TrProbM_by_age, paramSF.leProb1V, cSF.nw_expanded] = model_setup_utils_bgp.EarningProcess_AgeDependent(cSF);
params_for_ssF = struct('Z', Z_path(:,end), 'A', 1.0, 'theta', theta_for_ssF, 'g_A_ss', cSF.g_A_ss, 'n_ss', cSF.n_ss);

tic;
fprintf('   ⚙️  启动统一稳态求解器 (求解 ssF)...\n');
[ssF_result, distF, polF, valF] = main_steady_state_utils_bgp.solve_steady_state_unified(cSF, paramSF, params_for_ssF, true, [], 'lsqnonlin');
toc;

if isempty(ssF_result), warning('终期稳态(ssF)求解失败！程序将终止。'); return;
else, ssF = ssF_result; fprintf('✅ 终期稳态(ssF)求解成功！ r=%.4f, w=%.4f\n', ssF.r_mkt, ssF.w_hat); end

%% --- 5. 保存转轨所需数据 ---
if ~isempty(fieldnames(ssF))
    fprintf('\n\n--- 5. 保存转轨分析所需数据 ---\n');
    if ~exist('SS', 'dir'), mkdir('SS'); end
    save_filename = 'SS/data_for_transition_nopps.mat';
    
    fprintf('   将关键外生路径和模型维度参数打包进cS结构体...\n');
    cS.A_path = A_path;
    cS.Z_path = Z_path;
    cS.Z_path_raw = Z_path_raw;
    cS.theta_path = theta_path;
    cS.nw_expanded = cSF.nw_expanded; 
    
    data_for_transition = struct(...
        'ss0', ss0, 'ssF', ssF, 'dist0', [], 'distF', distF, 'pol0', [], ...
        'polF', valF, 'val0', [], 'valF', valF, 'cS', cS, 'paramS0', [], ...
        'paramSF', paramSF);
    
    save(save_filename, 'data_for_transition', '-v7.3');
    fprintf('✅ 所有转轨所需数据已成功保存至:\n   %s\n', save_filename);
else
    fprintf('\n\n数据未保存，因为终期稳态求解失败。\n');
end

%% --- 6. 结果分析与一致性检验 (最终版) ---
if ~isempty(fieldnames(ssF))
    fprintf('\n\n============================================================\n');
    fprintf('===            终期稳态 (ssF) 结果分析与检验           ===\n');
    fprintf('============================================================\n');
    
    % --- [核心] 资本存量稳态一致性检验 ---
    fprintf('--- [核心] 资本存量稳态一致性检验 ---\n');
    
    % fprintf('   检验0: 微观预算是否加总平衡？\n');
    % fprintf('   聚合后微观预算总误差             : %15.4e\n', ssF.micro_budget_balance_error);
    % if abs(ssF.micro_budget_balance_error) > 1e-9 % 使用一个合理的数值容忍度
    %      fprintf('   -> 警告: 微观预算约束在聚合层面不为零，核算口径不一致。\n');
    % else
    %      fprintf('   -> 通过: 微观预算在聚合层面是平衡的。\n');
    % end

    fprintf('\n   检验1: 稳态存量是否自洽 (期初总量 == 期末总量)？\n');
    fprintf('   K_p_hat_begin (求解器均衡点)       : %15.8f\n', ssF.K_private_begin_hat);
    fprintf('   K_p_hat_end (模型聚合期末选择)     : %15.8f\n', ssF.K_private_hat);
    ratio = ssF.K_private_hat / ssF.K_private_begin_hat;
    fprintf('   -> 两者之比 (End / Begin)           : %15.8f\n', ratio);

    consistency_error = ratio - 1;
    fprintf('   -> 一致性误差 ((比率)-1)          : %14.2f%%\n', consistency_error * 100);

    if abs(consistency_error) < 1e-4 % 容忍度设为 0.01%
        fprintf('   \n   检验结果: ✅ 通过。资本存量在BGP稳态下完全自洽。\n');
    else
        fprintf('   \n   检验结果: ⚠️  警告。资本存量不自洽。\n');
    end

    fprintf('============================================================\n');
else
    fprintf('\n\n结果分析无法执行，因为终期稳态(ssF)求解失败。\n');
end

fprintf('\n--- 无PPS稳态数据准备脚本执行完毕 ---\n');