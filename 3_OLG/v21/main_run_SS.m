% =========================================================================
% == SCRIPT: main_run_SS.m
% == 版本: [v4.0 - 新PAYG制度与内生基金版]
% ==
% == 核心修改:
% ==   - 新增 K_payg_to_Y_ratio 校准参数。
% ==   - solve_steady_state 调用接口更新，传入 is_bgp_ss 标志。
% ==   - display_national_accounts 调用接口更新，传入完整的 dist 和 polS。
% ==   - 统一 is_final_ss 和 is_bgp_ss 为一个标志 is_bgp_ss。
% =========================================================================
clear; close all; clear classes;
addpath(pwd);
fprintf('=== OLG异质性模型稳态求解脚本 (v4.0) ===\n\n');

%% --- 1. 全局设置与初始化 ---
fprintf('--- 1. 全局设置与初始化 ---\n');
fprintf('   加载模型物理参数 (含异质性设定)...\n');
cS = utils.ParameterValues();

cS.nk = 50;
cS.g_A_ss = 0.015;

% --- [新设定] PAYG 基金校准目标 ---
cS.K_payg_to_Y_ratio_ss0 = 0.05; % 初始稳态基金/GDP目标
cS.K_payg_to_Y_ratio_ssF = 0.0;  % 终期稳态基金/GDP目标 (强制为0)

% --- 人口路径设定 ---
cS.ss0_year = 2023;
cS.start_year = 2023;
cS.pop_data_last_year = 2040;

% --- 报告文件名 ---
report_file_initial = 'SS/初期稳态NIPS_het.txt';
report_file_final = 'SS/终期稳态NIPS_het.txt';

%% --- 2. 生成完整的外生路径 ---
fprintf('\n--- 2. 生成完整的外生路径 (内生决定模拟长度) ---\n');
[Z_path, Z_path_raw, cS] = population.generate_Z_path(cS, false);
T_sim = size(Z_path_raw, 2);
g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
cS.A_path = (1 + g_A_period).^(0:(T_sim-1));

% --- 核心异质性设定 ---
cS.num_hh_types = 4;
nH = cS.num_hh_types;
[~, ~, ~, cS.nw_expanded] = utils.EarningProcess_AgeDependent(cS);
cS = utils.calcaulte_theta_payg_path(cS, false);
cS.Z_path = Z_path;
cS.Z_path_raw = Z_path_raw;

%% --- 3. 求解初始伪稳态 (ss0, t=1) ---
fprintf('\n\n--- 3. 求解初始伪稳态 (ss0, t=1) ---\n');
fprintf('   特征: %d年人口结构, K_payg/Y=%.2f, 非BGP均衡\n', cS.ss0_year, cS.K_payg_to_Y_ratio_ss0);

cS0 = cS;
cS0.pps_active = false; % 为简化，此处保持无PPS
cS0.nkpps = 1; cS0.n_pps_rate_grid = 1;
cS0 = utils.generateGrids(cS0);

cS0.n_ss = ((sum(Z_path_raw(:,2))/sum(Z_path_raw(:,1)))^(1/cS.time_Step))-1;
cS0.s_pathV = cS.s_pathV(:,1);
cS0.theta_path_h = cS.theta_path_h(:,1);

Z0_total_abs = Z_path_raw(:,1);
Z0_norm_vec = Z0_total_abs / sum(Z0_total_abs);
Z0_norm = Z0_norm_vec * cS.type_weights';

paramS0 = struct();
[paramS0.leGridV, paramS0.TrProbM_by_age, paramS0.leProb1V, ~] = utils.EarningProcess_AgeDependent(cS0);

params_for_ss0 = struct(...
    'Z', Z0_norm, ...
    'A', 1.0, ...
    'g_A_ss', cS0.g_A_ss, ...
    'n_ss', cS0.n_ss);

tic;
fprintf('   ⚙️  启动统一稳态求解器 (求解 ss0, 非BGP模式)...\n');
is_bgp_ss0 = false;
[ss0, dist0_h, polS0_h, valS0_h] = SS.solve_steady_state(cS0, paramS0, params_for_ss0, true, is_bgp_ss0);
toc;

if isempty(ss0)
    error('初始稳态(ss0)求解失败！请检查代码逻辑或参数。');
else
    fprintf('✅ 初始稳态(ss0)求解成功！\n');
    SS.display_national_accounts(ss0, cS0, paramS0, Z0_norm, dist0_h, polS0_h, report_file_initial, true, is_bgp_ss0);
end

%% --- 4. 求解终期稳态 (ssF, t=T) ---
fprintf('\n\n--- 4. 求解终期BGP稳态 (ssF, t=T) ---\n');
fprintf('   特征: 长期稳定人口结构, K_payg=0, 真实BGP均衡\n');

cSF = cS;
cSF.pps_active =false;
cSF.nkpps = 1; cSF.n_pps_rate_grid = 1;
cSF = utils.generateGrids(cSF);

cSF.theta_path_h = cS.theta_path_h(:,end);
cSF.s_pathV = cS.s_pathV(:,end);

ZF_total_abs = Z_path_raw(:,end);
ZF_norm_vec = ZF_total_abs / sum(ZF_total_abs);
ZF_norm = ZF_norm_vec * cS.type_weights';

paramSF = struct();
[paramSF.leGridV, paramSF.TrProbM_by_age, paramSF.leProb1V, ~] = utils.EarningProcess_AgeDependent(cSF);

params_for_ssF = struct(...
    'Z', ZF_norm, ...
    'A', 1.0, ...
    'g_A_ss', cSF.g_A_ss, ...
    'n_ss', cSF.n_ss);

tic;
fprintf('   ⚙️  启动统一稳态求解器 (求解 ssF, BGP模式)...\n');
is_bgp_ssF = true;
[ssF, distF_h, polSF_h, valSF_h] = SS.solve_steady_state(cSF, paramSF, params_for_ssF, true, is_bgp_ssF);
toc;

if isempty(ssF)
    warning('终期稳态(ssF)求解失败！程序将终止。');
    return;
else
    fprintf('✅ 终期稳态(ssF)求解成功！\n');
    SS.display_national_accounts(ssF, cSF, paramSF, ZF_norm, distF_h, polSF_h, report_file_final, true, is_bgp_ssF);
end

%% --- 5. 保存转轨所需数据 ---
if ~isempty(fieldnames(ssF))
    fprintf('\n\n--- 5. 保存转轨分析所需数据 ---\n');
    if ~exist('SS', 'dir'), mkdir('SS'); end
        if cS.pps_active
        save_filename = 'SS/data_for_het_transition_pps.mat';
    else
        save_filename = 'SS/data_for_het_transition_nopps.mat';
        end
    
    
    cS.pps_active = cSF.pps_active;
    fprintf('   已将终期PPS状态 (pps_active = %s) 保存至cS，用于转轨路径分析。\n', string(cS.pps_active));
    
    dist0_abs_h = cell(nH, 1);
    for h = 1:nH
        pop_h_by_age_t1 = cS.Z_path_raw(:, 1) * cS.type_weights(h);
        dist_cond_h = zeros(size(dist0_h{h}));
        for ia = 1:cS.aD_new
            pop_mass_norm_hia = sum(dist0_h{h}(:,:,:,ia), 'all');
            if pop_mass_norm_hia > 1e-12
                dist_cond_h(:,:,:,ia) = dist0_h{h}(:,:,:,ia) / pop_mass_norm_hia;
            end
        end
        dist0_abs_h{h} = dist_cond_h .* reshape(pop_h_by_age_t1, [1, 1, 1, cS.aD_new]);
    end
    
    cS_for_trans = cSF;
    cS_for_trans.theta_path_h = cS.theta_path_h;
    cS_for_trans.s_pathV = cS.s_pathV;
    cS_for_trans.Z_path_raw = cS.Z_path_raw;
    cS_for_trans.type_weights = cS.type_weights;
    cS_for_trans.ageEffV_new_h = cS.ageEffV_new_h;

    data_for_transition = struct();
    data_for_transition.ss0 = ss0;
    data_for_transition.ssF = ssF;
    data_for_transition.dist0_h = dist0_abs_h;
    data_for_transition.distF_h = distF_h;
    data_for_transition.polF_h = polSF_h;   % [修改] 保存终期策略
    data_for_transition.valF_h = valSF_h;   % [修改] 保存终期价值函数
    data_for_transition.cS = cS_for_trans;
    data_for_transition.paramSF = paramSF;
    

    save(save_filename, 'data_for_transition', "-v7.3");
    fprintf('✅ 所有转轨所需数据已成功保存至:\n   %s\n', save_filename);
else
    fprintf('\n\n数据未保存，因为终期稳态求解失败。\n');
end
%% --- 6. 比较初始与终期稳态的国民核算结果 ---
fprintf('\n--- 6. 比较初始稳态(ss0)与终期稳态(ssF)的国民核算 ---\n');
utils.compare_steady_states_from_files(report_file_initial, report_file_final);

fprintf('\n--- 稳态数据准备与比较脚本执行完毕 ---\n');