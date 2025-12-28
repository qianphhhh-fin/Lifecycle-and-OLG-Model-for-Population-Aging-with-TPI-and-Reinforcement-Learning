% =========================================================================
% == SCRIPT: main_run_SS.m
% == 版本: [v2.4 - 转轨数据准备最终版]
% ==
% == 核心修改:
% ==   - 在保存转轨数据时，明确地将终期稳态的 pps_active 状态
% ==     赋给全局参数 cS，确保转轨路径求解器能正确识别PPS模式。
% =========================================================================
clear; close all; clear classes;
addpath(pwd);
fprintf('=== OLG模型 (无PPS情景) 稳态求解与数据准备脚本 ===\n\n');

%% --- 1. 全局设置与初始化 ---
fprintf('--- 1. 全局设置与初始化 ---\n');

ngrid = 40;
ngrid_pps = 1; 

fprintf('   加载模型物理参数...\n');
cS = utils.ParameterValues();

fprintf('   设定并生成全局共享的资产网格...\n');
cS.nk = ngrid; 
cS.nkpps = ngrid_pps; 
cS.nkprime = ngrid; 
cS.npps = ngrid_pps;
cS = utils.generateGrids(cS);
fprintf('   ✅ 全局网格(kGridV等)已生成并存入cS。\n');

cS.endogenous_theta_mode = true;
cS.pps_active = false;

cS.ss0_year = 2023;
cS.start_year = 2023;

if cS.endogenous_theta_mode
    cS.payg_replacement_rate = 0.4;
    cS.theta_max = 0.99;            
end


%% --- 2. [核心] 生成完整的外生路径 ---
fprintf('\n--- 2. 生成完整的外生路径 (内生决定模拟长度) ---\n');
[Z_path,Z_path_raw, A_path, cS] = utils.generate_exo_paths(cS, false);
cS = utils.calcaulte_theta_payg_path(cS, false);
theta_path = cS.theta_path;

ss0 = struct();
ssF = struct();

% --- [新增] 定义报告文件名 ---
report_file_initial = 'SS/初期稳态NIPS_noPPS.txt';
report_file_final = 'SS/终期稳态NIPS_PPS.txt';


%% --- 3. 求解初始稳态 (ss0, t=0) ---
fprintf('\n\n--- 3. 求解初始伪稳态 (ss0, t=0) ---\n');
fprintf('   特征: %d年人口结构, 无长期增长(g=0), 无PPS制度\n', cS.ss0_year);

cS0 = cS;
cS0.pps_active = false;
cS0.g_A_ss = 0.05;
cS0.s_pathV = cS.s_pathV(:,1); 
Z0 = utils.get_calibration_inputs(cS.ss0_year, cS0);

cS0.nk = ngrid; cS0.nkpps = 1; cS0.nkprime = ngrid; cS0.npps = 1;
cS0 = utils.generateGrids(cS0);
paramS0 = struct();
[paramS0.leGridV, paramS0.TrProbM_by_age, paramS0.leProb1V, cS0.nw_expanded] = ...
    utils.EarningProcess_AgeDependent(cS0);
paramS0.leLogGridV = log(paramS0.leGridV(1:cS0.nw));

params_for_ss0 = struct(...
    'Z', Z0, 'A', 1.0, 'theta', theta_path(1), 'g_A_ss', cS0.g_A_ss);

tic;
fprintf('   ⚙️  启动统一稳态求解器 (求解 ss0)...\n');
[ss0, dist0, ~, ~] = ...
    SS.solve_steady_state(cS0, paramS0, params_for_ss0, true, 'lsqnonlin');
toc;

if isempty(ss0)
    error('初始稳态(ss0)求解失败！请检查代码逻辑或参数。');
else
    fprintf('✅ 初始稳态(ss0)求解成功！\n');
    SS.display_national_accounts(ss0, cS0, paramS0, Z0, report_file_initial);
    fprintf('   利率和工资: r=%.4f, w=%.4f\n', ss0.r_mkt, ss0.w_hat);
end


%% --- 4. 求解终期稳态 (ssF, t=T) ---
fprintf('\n\n--- 4. 求解终期BGP稳态 (ssF, t=T) ---\n');

cSF = cS;
% *************************************************************************
% [修改演示] 在这里决定是否为终期稳态激活PPS
cSF.pps_active = true; % 改为 false 以匹配你当前的 "nopps" 设定
if cSF.pps_active
    cSF.nkpps = ngrid; % 如果激活，则使用完整的PPS网格
    cSF = utils.generateGrids(cSF, 'k_max', 40, 'kpps_max', 20);
else
    cSF.nkpps = 1; % 未激活，则设为1
    cSF = utils.generateGrids(cSF);
end
% *************************************************************************
cSF.s_pathV = cS.s_pathV(:,end); 
fprintf('   特征: 长期稳态人口(n=%.3f), 长期增长(g=%.3f), **PPS激活=%s**\n', cSF.n_ss, cSF.g_A_ss, string(cSF.pps_active));

if cSF.endogenous_theta_mode
    fprintf('   [养老金模式] DB模式激活, 目标替代率=%.2f%%, 上限=%.2f%%\n', cSF.payg_replacement_rate*100, cSF.theta_max*100);
else
    fprintf('   [养老金模式] DC模式激活, 外生缴费率=%.4f\n', theta_path(end));
end

paramSF = struct();
[paramSF.leGridV, paramSF.TrProbM_by_age, paramSF.leProb1V, cSF.nw_expanded] = utils.EarningProcess_AgeDependent(cSF);
params_for_ssF = struct('Z', Z_path(:,end), 'A', 1.0, 'g_A_ss', cSF.g_A_ss, 'n_ss', cSF.n_ss);

tic;
fprintf('   ⚙️  启动统一稳态求解器 (求解 ssF)...\n');
[ssF_result, distF, polF, valF] = SS.solve_steady_state(cSF, paramSF, params_for_ssF, true, 'lsqnonlin');
toc;

if isempty(ssF_result)
    warning('终期稳态(ssF)求解失败！程序将终止。'); 
    return;
else
    ssF = ssF_result; 
    SS.display_national_accounts(ssF, cSF, paramSF,  Z_path(:,end), report_file_final);
    fprintf('✅ 终期稳态(ssF)求解成功！ r=%.4f, w=%.4f\n', ssF.r_mkt, ssF.w_hat); 
end


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

    % ======================== [核心修正] ========================
    % 将最终稳态的PPS激活状态，赋给全局cS参数集，
    % 这样转轨路径求解器就知道该使用哪种模式。
    cS.pps_active = cSF.pps_active;
    fprintf('   已将终期PPS状态 (pps_active = %s) 保存至cS，用于转轨路径分析。\n', string(cS.pps_active));
    % ==========================================================
    
    data_for_transition = struct(...
        'ss0', ss0, 'ssF', ssF, ...
        'dist0', dist0, ...
        'distF', distF, 'pol0', [], 'polF', polF, 'val0', [], 'valF', valF, ...
        'cS', cS, 'paramS0', paramS0, 'paramSF', paramSF);
    
    save(save_filename, 'data_for_transition', '-v7.3');
    fprintf('✅ 所有转轨所需数据已成功保存至:\n   %s\n', save_filename);
else
    fprintf('\n\n数据未保存，因为终期稳态求解失败。\n');
end



%% --- 7. [新增] 调用封装好的比较函数进行最终分析 ---
if ~isempty(fieldnames(ss0)) && ~isempty(fieldnames(ssF))
    utils.compare_steady_states_from_files(report_file_initial, report_file_final);
else
    fprintf('\n\n比较分析无法执行，因为至少一个稳态求解失败。\n');
end

fprintf('\n--- 稳态数据准备与比较脚本执行完毕 ---\n');