% =========================================================================
% == SCRIPT: main_run_SS.m
% == 版本: [v3.0 - 真实数据校准最终版]
% ==
% == 核心修改:
% ==   - 统一外生路径生成，确保 ss0 和 ssF 基于同一套宏观背景。
% ==   - 移除过时的内生theta(替代率)逻辑，与新的PAYG求解机制完全对齐。
% ==   - 确保 ss0 的人口分布与转轨路径的 t=1 时刻完全一致。
% =========================================================================
clear; close all; clear classes;
addpath(pwd);
fprintf('=== OLG模型稳态求解与转轨数据准备脚本 (v3.0) ===\n\n');

%% --- 1. 全局设置与初始化 ---
fprintf('--- 1. 全局设置与初始化 ---\n');

ngrid = 40;
ngrid_pps = 1; % 为包含PPS的终期稳态预设网格

fprintf('   加载模型物理参数...\n');
cS = utils.ParameterValues();
cS.nk = ngrid; cS.nkpps = 1; cS.nkprime = ngrid;
cS = utils.generateGrids(cS);


% 设定模拟的起始年份
cS.ss0_year = 2023;
cS.start_year = 2023;

% --- 报告文件名 ---
report_file_initial = 'SS/初期稳态NIPS.txt';
report_file_final = 'SS/终期稳态NIPS.txt';

%% --- 2. [核心] 生成完整的外生路径 ---
fprintf('\n--- 2. 生成完整的外生路径 (内生决定模拟长度) ---\n');
[Z_path, Z_path_raw, A_path, cS] = utils.generate_exo_paths(cS, false); % 可视化检查
cS = utils.calcaulte_theta_payg_path(cS, false); % 可视化检查
theta_path = cS.theta_path;

[~, ~, ~, cS.nw_expanded] = utils.EarningProcess_AgeDependent(cS);


% 将核心路径存入cS，以便传递
cS.A_path = A_path;
cS.Z_path = Z_path; % 归一化分布
cS.Z_path_raw = Z_path_raw; % 绝对人口

%% --- 3. 求解初始稳态 (ss0, t=0) ---
fprintf('\n\n--- 3. 求解初始伪稳态 (ss0, t=0) ---\n');
fprintf('   特征: %d年人口结构, 无PPS制度\n', cS.ss0_year);

cS0 = cS;
cS0.pps_active = false; % 初始稳态不含PPS
cS0.g_A_ss = ((A_path(2)/A_path(1))^(1/cS.time_Step))-1; % 使用路径初始期的增长率
cS0.n_ss = ((sum(Z_path_raw(:,2))/sum(Z_path_raw(:,1)))^(1/cS.time_Step))-1;
cS0.s_pathV = cS.s_pathV(:,1); % 使用第一期的生存率
cS0.theta_path = theta_path(1);


% [关键] 使用与转轨完全一致的 t=1 人口分布 (绝对量)
Z0_abs = Z_path_raw(:,1);
Z0_norm = Z0_abs / sum(Z0_abs); % 求解器需要归一化分布

cS0.nk = ngrid; cS0.nkpps = 1; cS0.nkprime = ngrid;
cS0 = utils.generateGrids(cS0);
paramS0 = struct();
[paramS0.leGridV, paramS0.TrProbM_by_age, paramS0.leProb1V, ~] = utils.EarningProcess_AgeDependent(cS0);

% 准备传递给求解器的外部参数
params_for_ss0 = struct(...
    'Z', Z0_norm, ...
    'A', 1.0, ...
    'g_A_ss', cS0.g_A_ss, ...
    'n_ss', cS0.n_ss);

tic;
fprintf('   ⚙️  启动统一稳态求解器 (求解 ss0)...\n');
[ss0, dist0, ~, ~] = SS.solve_steady_state(cS0, paramS0, params_for_ss0, true, 'lsqnonlin');
toc;

if isempty(ss0)
    error('初始稳态(ss0)求解失败！请检查代码逻辑或参数。');
else
    fprintf('✅ 初始稳态(ss0)求解成功！\n');
    % 使用绝对人口分布来报告有意义的宏观总量
    dist0_abs = dist0 .* sum(Z0_abs); 
    SS.display_national_accounts(ss0, cS0, paramS0, Z0_norm, report_file_initial);
end


%% --- 4. 求解终期稳态 (ssF, t=T) ---
fprintf('\n\n--- 4. 求解终期BGP稳态 (ssF, t=T) ---\n');

cSF = cS;
cSF.pps_active = false; % 在这里决定终期稳态是否激活PPS
cSF.theta_path = theta_path;


if cSF.pps_active
    cSF.nk = ngrid; cSF.nkpps = ngrid_pps; cSF.nkprime = ngrid;
    cSF = utils.generateGrids(cSF, 'k_max', 40, 'kpps_max', 20);
    fprintf('   特征: 长期稳态人口(n=%.3f), 长期增长(g=%.3f), **PPS已激活**\n', cSF.n_ss, cSF.g_A_ss);
else
    cSF.nk = ngrid; cSF.nkpps = 1; cSF.nkprime = ngrid;
    cSF = utils.generateGrids(cSF);
    fprintf('   特征: 长期稳态人口(n=%.3f), 长期增长(g=%.3f), **PPS未激活**\n', cSF.n_ss, cSF.g_A_ss);
end

cSF.s_pathV = cS.s_pathV(:,end); % 使用最后一期的生存率
ZF_norm = Z_path(:,end); % 使用最后一期的归一化人口分布

paramSF = struct();
[paramSF.leGridV, paramSF.TrProbM_by_age, paramSF.leProb1V, ~] = utils.EarningProcess_AgeDependent(cSF);

params_for_ssF = struct(...
    'Z', ZF_norm, ...
    'A', 1.0, ...
    'g_A_ss', cSF.g_A_ss, ...
    'n_ss', cSF.n_ss);

tic;
fprintf('   ⚙️  启动统一稳态求解器 (求解 ssF)...\n');
[ssF, distF, polF, valF] = SS.solve_steady_state(cSF, paramSF, params_for_ssF, true, 'lsqnonlin');
toc;

if isempty(ssF)
    warning('终期稳态(ssF)求解失败！程序将终止。');
    return;
else
    fprintf('✅ 终期稳态(ssF)求解成功！\n');
    % 报告时，可以乘以一个任意的总人口（如1）来查看人均等价的宏观量
    distF_scaled = distF * sum(Z_path_raw(:,end));
    SS.display_national_accounts(ssF, cSF, paramSF, ZF_norm, report_file_final);
end


%% --- 5. 保存转轨所需数据 ---
if ~isempty(fieldnames(ssF))
    fprintf('\n\n--- 5. 保存转轨分析所需数据 ---\n');
    if ~exist('SS', 'dir'), mkdir('SS'); end
    save_filename = 'SS/data_for_transition.mat';
    
    % 将最终稳态的PPS激活状态，赋给全局cS参数集
    cS.pps_active = cSF.pps_active;
    fprintf('   已将终期PPS状态 (pps_active = %s) 保存至cS，用于转轨路径分析。\n', string(cS.pps_active));
    
    % 确保dist0是绝对分布，与Z_path_raw(:,1)的规模一致
    dist0_for_trans = dist0 * sum(Z_path_raw(:,1));
    
    data_for_transition = struct(...
        'ss0', ss0, 'ssF', ssF, ...
        'dist0', dist0_for_trans, ... % 使用与Z_path_raw对齐的绝对分布
        'distF', distF, 'polF', polF, 'valF', valF, ... % 稳态求解器本身输出的是归一化分布和策略
        'cS', cS, 'paramS0', paramS0, 'paramSF', paramSF);
    
    save(save_filename, 'data_for_transition', '-v7.3');
    fprintf('✅ 所有转轨所需数据已成功保存至:\n   %s\n', save_filename);
else
    fprintf('\n\n数据未保存，因为终期稳态求解失败。\n');
end

%% --- 6. [新增] 调用封装好的比较函数进行最终分析 ---
if ~isempty(fieldnames(ss0)) && ~isempty(fieldnames(ssF))
    utils.compare_steady_states_from_files(report_file_initial, report_file_final);
else
    fprintf('\n\n比较分析无法执行，因为至少一个稳态求解失败。\n');
end

fprintf('\n--- 稳态数据准备与比较脚本执行完毕 ---\n');