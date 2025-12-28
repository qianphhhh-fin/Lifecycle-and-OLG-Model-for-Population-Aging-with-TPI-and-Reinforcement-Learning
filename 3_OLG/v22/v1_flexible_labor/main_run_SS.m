% =========================================================================
% == SCRIPT: main_run_SS.m
% == 版本: [v3.1 - 异质性框架完全对齐版]
% ==
% == 核心修改:
% ==   - [!!!] 修正了 ss0 和 ssF 的参数构建逻辑，使其与异质性框架完全对齐。
% ==   - 确保将 Z_h (aD x nH 矩阵) 和 cS 中包含 _h 后缀的异质性参数
% ==     (如 theta_path_h, ageEffV_new_h) 正确传递给求解器。
% ==   - 统一外生路径生成，确保 ss0 和 ssF 基于同一套宏观背景。
% ==   - 确保 ss0 的人口分布与转轨路径的 t=1 时刻完全一致。
% =========================================================================
clear; close all; clear classes;
addpath(pwd);
fprintf('=== OLG异质性模型稳态求解与转轨数据准备脚本 (v3.1) ===\n\n');





%% --- 1. 全局设置与初始化 ---
fprintf('--- 1. 全局设置与初始化 ---\n');

fprintf('   加载模型物理参数 (含异质性设定)...\n');
cS = utils.ParameterValues();

cS.nk = 50;
% cS.g_A_ss = 0.015; % 长期技术增长率

% --- 人口路径设定 ---
cS.ss0_year = 2023;
cS.start_year = 2023;
cS.pop_data_last_year = 2070; % 人口数据在2050年后进入稳态

cS.g_A_ss = 0.015;                  % 长期稳态下技术进步的年化增长率 (TFP growth rate)。
cS.transition_period_years = 100;    % 外生技术进步路径从数据点过渡到长期稳态所需的年数。

% --- 报告文件名 ---
report_file_initial = 'SS/初期稳态NIPS_het.txt';
report_file_final = 'SS/终期稳态NIPS_het.txt';

%% --- 2. [核心] 生成完整的外生路径 ---
fprintf('\n--- 2. 生成完整的外生路径 (内生决定模拟长度) ---\n');
[Z_path,Z_path_raw, cS] = population.generate_Z_path(cS, false);
T_sim = size(Z_path_raw, 2);

[cS.A_path,g_path_annual_output,~] = utils.generate_tfp_path(cS,'baseline',false); % 'optimistic'  'pessimistic' 'baseline' 'constant'

% --- [!!! 核心异质性设定 !!!] ---
cS.num_hh_types = 4;
nH = cS.num_hh_types;
% [核心修改] cS 现在已包含 cS.theta_path_h 和 cS.ageEffV_new_h
[~, ~, ~, cS.nw_expanded] = utils.EarningProcess_AgeDependent(cS);
cS = utils.theta_payg_path(cS, 'baseline',false); % 可视化检查

cS = utils.retire_age_path(cS, 'mild', Z_path_raw,false); % 'mild' 'moderate' ''aggressive'

% 将核心路径存入cS，以便传递
cS.Z_path = Z_path; % 归一化分布
cS.Z_path_raw = Z_path_raw; % 绝对人口

cS = utils.labor_setting(cS); % false 表示不显示绘图
cS = utils.generate_PPS_path(cS, 'conservative',false); % 'conservative' 'moderate' 'ambitious'


%% --- 3. 求解初始伪稳态 (ss0, t=0) ---
fprintf('\n\n--- 3. 求解初始伪稳态 (ss0, t=0) ---\n');
fprintf('   特征: %d年人口结构, 无PPS制度\n', cS.ss0_year);

cS0 = cS;
cS0.pps_active = false;
cS0.nkpps = 1; cS0.n_pps_rate_grid = 1;
cS0 = utils.generateGrids(cS0);


cS0.aR_new = cS.aR_new_path(1); 
fprintf('      [参数设定] 初始稳态退休年龄索引(aR_new)设置为: %d\n', cS0.aR_new);

% 使用转轨第一期的参数
cS0.n_ss = ((sum(Z_path_raw(:,2))/sum(Z_path_raw(:,1)))^(1/cS.time_Step))-1;
cS0.s_pathV = cS.s_pathV(:,1);
cS0.theta_path_h = cS.theta_path_h(:,1); % 使用第一期的缴费率向量
cS0.g_A_ss = g_path_annual_output(1); 

% [核心修改] 构建 t=1 的异质性人口分布矩阵
Z0_total_abs = Z_path_raw(:,1);
Z0_norm = Z0_total_abs / sum(Z0_total_abs); % 求解器需要归一化分布
Z0_norm = Z0_norm * cS.type_weights'; % [aD x 1] * [1 x nH] -> [aD x nH]

paramS0 = struct();
[paramS0.leGridV, paramS0.TrProbM_by_age, paramS0.leProb1V, ~] = utils.EarningProcess_AgeDependent(cS0);

params_for_ss0 = struct(...
    'Z', Z0_norm, ... % [!!!] 传递矩阵形式的相对人口分布
    'A', 1.0, ...
    'g_A_ss', 0, ... %  
    'n_ss', 0); % 

tic;
fprintf('   ⚙️  启动统一稳态求解器 (求解 ss0, 异质性模式)...\n');
[ss0, dist0_h, ~, ~] = SS.solve_steady_state(cS0, paramS0, params_for_ss0, true, false);
toc;

if isempty(ss0)
    error('初始稳态(ss0)求解失败！请检查代码逻辑或参数。');
else
    fprintf('✅ 初始稳态(ss0)求解成功！\n');
    SS.display_national_accounts(ss0, cS0, paramS0, Z0_norm, report_file_initial, true, false);
end

%% --- 4. 求解终期稳态 (ssF, t=T) ---
fprintf('\n\n--- 4. 求解终期BGP稳态 (ssF, t=T) ---\n');

cSF = cS;
cSF.pps_active = false; % 在这里决定终期稳态是否激活PPS
cSF.nkpps = 1; cSF.n_pps_rate_grid = 1;
cSF.pps_fixed = cSF.pps_fixed_path(end);

% === 核心修改: 使用路径的【终点】作为终期稳态的退休年龄 ===
cSF.aR_new = cS.aR_new_path(end);
fprintf('      [参数设定] 终期稳态退休年龄索引(aR_new)设置为: %d\n', cSF.aR_new);
cSF.g_A_ss = g_path_annual_output(end); 
cSF = utils.generateGrids(cSF);

% 使用转轨最后一期的参数
cSF.theta_path_h = cS.theta_path_h(:,end);
cSF.s_pathV = cS.s_pathV(:,end);

% [核心修改] 构建 t=T 的异质性人口分布矩阵
ZF_total_abs = Z_path_raw(:,end);
ZF_norm = ZF_total_abs / sum(ZF_total_abs); % 求解器需要归一化分布
ZF_norm = ZF_norm * cS.type_weights'; % [aD x 1] * [1 x nH] -> [aD x nH]


paramSF = struct();
[paramSF.leGridV, paramSF.TrProbM_by_age, paramSF.leProb1V, ~] = utils.EarningProcess_AgeDependent(cSF);

params_for_ssF = struct(...
    'Z', ZF_norm, ... % [!!!] 传递矩阵形式的绝对人口分布
    'A', 1.0, ...
    'g_A_ss', cSF.g_A_ss , ... % 
    'n_ss', cSF.n_ss);

tic;
fprintf('   ⚙️  启动统一稳态求解器 (求解 ssF, 异质性模式)...\n');
[ssF, distF_h, polF_h, valF_h] = SS.solve_steady_state(cSF, paramSF, params_for_ssF, true, true);
toc;

if isempty(ssF)
    warning('终期稳态(ssF)求解失败！程序将终止。');
    return;
else
    fprintf('✅ 终期稳态(ssF)求解成功！\n');
    SS.display_national_accounts(ssF, cSF, paramSF, ZF_norm, report_file_final, true, true);
end

%% --- 5. 保存转轨所需数据 ---
if ~isempty(fieldnames(ssF))
    fprintf('\n\n--- 5. 保存转轨分析所需数据 ---\n');
    if ~exist('SS', 'dir'), mkdir('SS'); end

if cSF.pps_active
save_filename = ['SS/data_for_het_transition_pps.mat'];
else
save_filename = ['SS/data_for_het_transition_nopps.mat'];
end
    
    cS.pps_active = cSF.pps_active;
    fprintf('   已将终期PPS状态 (pps_active = %s) 保存至cS，用于转轨路径分析。\n', string(cS.pps_active));
    
    % =========================================================================
    % ==            [!!! 核心会计修正 !!!]
    % == 正确构造转轨第一期所需的【绝对人口分布】 dist0_abs_h
    % =========================================================================
    nH = cS.nTypes;
    dist0_abs_h = cell(nH, 1);

    for h = 1:nH
        % 1. 获取类型h在 t=1 时的【分年龄绝对人口】向量
        pop_h_by_age_t1 = cS.Z_path_raw(:, 1) * cS.type_weights(h); % [aD x 1]
        
        % 2. 从ss0返回的分布 dist0_h 中，计算出【条件分布】
        %    dist0_h{h} 是基于 Z0_norm(:,h) 计算的，Z0_norm(:,h) 是归一化人口乘以权重
        %    所以 dist0_h{h}(:,:,:,ia) 的加总之和等于 Z0_norm(ia, h)
        dist_cond_h = zeros(size(dist0_h{h}));
        for ia = 1:cS.aD_new
            pop_mass_norm_hia = sum(dist0_h{h}(:,:,:,ia), 'all');
            if pop_mass_norm_hia > 1e-12
                dist_cond_h(:,:,:,ia) = dist0_h{h}(:,:,:,ia) / pop_mass_norm_hia;
            end
        end
        
        % 3. 将条件分布乘以【分年龄绝对人口】，得到最终的绝对分布
        dist0_abs_h{h} = dist_cond_h .* reshape(pop_h_by_age_t1, [1, 1, 1, cS.aD_new]);
    end
    
    % [!!! 确保转轨参数是完整的路径 !!!]
    cS_for_trans = cSF; 
    cS_for_trans.theta_path_h = cS.theta_path_h;
    cS_for_trans.s_pathV = cS.s_pathV;
    cS_for_trans.Z_path_raw = cS.Z_path_raw;
    cS_for_trans.type_weights = cS.type_weights;
    cS_for_trans.ageEffV_new_h = cS.ageEffV_new_h; % 确保效率剖面也传递
    
    % === 核心修改: 将退休年龄路径也存入最终的 cS 结构体 ===
    cS_for_trans.aR_new_path = cS.aR_new_path;

    % [!!! 核心结构修正: 创建一个1x1的struct，其字段为元胞数组 !!!]
    data_for_transition = struct(); % 初始化一个1x1的struct
    data_for_transition.ss0 = ss0;
    data_for_transition.ssF = ssF;
    data_for_transition.dist0_h = dist0_abs_h; % 这是一个 cell array, 包含【绝对】人口分布
    data_for_transition.distF_h = distF_h;     % 这是一个 cell array, 包含【归一化】人口分布
    data_for_transition.polF_h = polF_h;       % 这是一个 cell array
    data_for_transition.valF_h = valF_h;       % 这是一个 cell array
    data_for_transition.cS = cS_for_trans;
    data_for_transition.paramSF = paramSF;
    save(save_filename, 'data_for_transition', '-v7.3');

    fprintf('✅ 所有转轨所需数据已成功保存至:\n   %s\n', save_filename);
else
    fprintf('\n\n数据未保存，因为终期稳态求解失败。\n');
end
%% --- 6. 比较初始与终期稳态的国民核算结果 ---
fprintf('\n--- 6. 比较初始稳态(ss0)与终期稳态(ssF)的国民核算 ---\n');
utils.compare_steady_states_from_files(report_file_initial, report_file_final);

fprintf('\n--- 稳态数据准备与比较脚本执行完毕 ---\n');