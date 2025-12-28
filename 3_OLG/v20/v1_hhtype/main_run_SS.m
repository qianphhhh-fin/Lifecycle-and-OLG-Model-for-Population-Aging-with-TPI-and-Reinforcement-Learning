% =========================================================================
% == SCRIPT: main_run_SS.m
% == 版本: [v4.0 - 异质性类型求解版]
% ==
% == 核心修改:
% ==   - 脚本完全适配包含4类异质性家庭的模型框架。
% ==   - 调用新的、支持异质性类型的 `ParameterValues` 和 `age_efficiency`。
% ==   - 正确处理返回的多类型策略 `polS` (cell数组) 和5维分布 `dist`。
% ==   - [!!!] 暂时注释掉了对 `display_national_accounts` 的调用，
% ==     因为它尚未被修改以处理新的 `polS` 和 `dist` 结构。
% =========================================================================
clear; close all; clear classes;
addpath(pwd);
fprintf('=== OLG模型稳态求解脚本 (v4.0 - 异质性类型版) ===\n\n');

%% --- 1. 全局设置与初始化 ---
fprintf('--- 1. 全局设置与初始化 ---\n');

fprintf('   加载模型物理参数 (支持异质性类型)...\n');
cS = utils.ParameterValues();

cS.nk = 50;
cS.g_A_ss = 0.02; 



cS.pop_data_last_year = 2030;
% --- 报告文件名 ---
report_file_initial = 'SS/初期稳态NIPS.txt';
report_file_final = 'SS/终期稳态NIPS.txt';

%% --- 2. [核心] 生成完整的外生路径 ---
fprintf('\n--- 2. 生成完整的外生路径 (内生决定模拟长度) ---\n');
[Z_path, Z_path_raw, cS] = population.generate_Z_path(cS, false); 
T_sim = size(Z_path,2);
g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
cS.A_path = (1 + g_A_period).^(0:(T_sim-1));

[~, ~, ~, cS.nw_expanded] = utils.EarningProcess_AgeDependent(cS);

% [修改] 调用修改后的PAYG费率函数
cS = utils.calcaulte_theta_payg_path(cS, false);

% 将核心路径存入cS，以便传递
cS.Z_path = Z_path; % 归一化分布
cS.Z_path_raw = Z_path_raw; % 绝对人口

%% --- 3. 求解初始稳态 (ss0, t=1) ---
% fprintf('\n\n--- 3. 求解初始伪稳态 (ss0, t=1) ---\n');
% fprintf('   特征: 2023年人口结构, 资本市场无需满足BGP积累, %d类家庭\n', cS.nTypes); % 更新描述
% 
% cS0 = cS;
% cS0.pps_active = false;
% cS0.nkpps = 1; cS0.n_pps_rate_grid = 1;
% cS0 = utils.generateGrids(cS0);
% 
% cS0.n_ss = ((sum(Z_path_raw(:,2))/sum(Z_path_raw(:,1)))^(1/cS.time_Step))-1;
% cS0.s_pathV = cS.s_pathV(:,1);
% cS0.theta_path_urban = cS.theta_path_urban(1);
% cS0.theta_path_resident = cS.theta_path_resident(1);
% 
% Z0_norm = Z_path(:,1);
% 
% paramS0 = struct();
% [paramS0.leGridV, paramS0.TrProbM_by_age, paramS0.leProb1V, ~] = utils.EarningProcess_AgeDependent(cS0);
% 
% params_for_ss0 = struct(...
%     'Z', Z0_norm, ...
%     'A', 1.0, ...
%     'g_A_ss', cS0.g_A_ss, ...
%     'n_ss', cS0.n_ss);
% 
% tic;
% fprintf('   ⚙️  启动【初始时期专用】求解器 (求解 ss0)...\n');
% % [!!! 核心修改: 调用新函数 !!!]
% [ss0, dist0, ~, ~] = SS.solve_steady_state_ss0(cS0, paramS0, params_for_ss0);
% toc;
% 
% if isempty(ss0)
%     error('初始稳态(ss0)求解失败！请检查代码逻辑或参数。');
% else
%     fprintf('✅ 初始稳态(ss0)求解成功！\n');
%     % 注意: display_national_accounts 中对投资缺口的解释现在对ss0有了新的意义
%     SS.display_national_accounts_ss0(ss0, cS0, paramS0, Z0_norm, report_file_initial, true);
% end
% % 在 main_run_SS.m 中
% % 调用函数
% utils.check_k_grid_limit(dist0, cS0.kGridV, 'ss0', 'k_private', true, [0.25, 0.5, 0.75, 0.9]);


%% --- 4. 求解终期稳态 (ssF, t=T) ---
fprintf('\n\n--- 4. 求解终期BGP稳态 (ssF, t=T) ---\n');

cSF = cS;
% cSF.n_ss=0;    
cSF.pps_active = false; % 在这里决定终期稳态是否激活PPS
cSF.nkpps = 1; cSF.n_pps_rate_grid = 1;
cSF = utils.generateGrids(cSF);

if cSF.pps_active
    fprintf('   特征: 长期稳态人口(n=%.3f), 长期增长(g=%.3f), %d类家庭, **PPS已激活**\n', cSF.n_ss, cSF.g_A_ss, cS.nTypes);
else
    fprintf('   特征: 长期稳态人口(n=%.3f), 长期增长(g=%.3f), %d类家庭, **PPS未激活**\n', cSF.n_ss, cSF.g_A_ss, cS.nTypes);
end

cSF.s_pathV = cS.s_pathV(:,end); % 使用最后一期的生存率
ZF_norm = Z_path(:,end); % 使用最后一期的归一化人口分布

% 选择正确的缴费率路径 (稳态使用最后一个值)
cSF.theta_path_urban = cS.theta_path_urban(end);
cSF.theta_path_resident = cS.theta_path_resident(end);

paramSF = struct();
[paramSF.leGridV, paramSF.TrProbM_by_age, paramSF.leProb1V, ~] = utils.EarningProcess_AgeDependent(cSF);

params_for_ssF = struct(...
    'Z', ZF_norm, ...
    'A', 1.0, ...
    'g_A_ss', cSF.g_A_ss, ...
    'n_ss', cSF.n_ss);

tic;
fprintf('   ⚙️  启动统一稳态求解器 (求解 ssF)...\n');
[ssF, distF, polF, valF] = SS.solve_steady_state(cSF, paramSF, params_for_ssF, true, true);
toc;
if isempty(ssF)
    warning('终期稳态(ssF)求解失败！程序将终止。');
    return;
else
    fprintf('✅ 终期稳态(ssF)求解成功！\n');
    % [!!!] 报告函数尚未更新，暂时注释
    SS.display_national_accounts(ssF, cSF, paramSF, ZF_norm, report_file_final, true, true);
end

%% --- 5. 保存转轨所需数据 ---
if ~isempty(ssF) && ~isempty(ss0)
    fprintf('\n\n--- 5. 保存转轨分析所需数据 ---\n');
    if ~exist('SS', 'dir'), mkdir('SS'); end
    if cSF.pps_active
        save_filename = 'SS/data_for_transition_pps_het_types.mat';
    else
        save_filename = 'SS/data_for_transition_nopps_het_types.mat';
    end

    cS.pps_active = cSF.pps_active;
    fprintf('   已将终期PPS状态 (pps_active = %s) 保存至cS，用于转轨路径分析。\n', string(cS.pps_active));
    
    % [核心会计] 构造转轨第一期需要的绝对人口分布 dist0 (5维)
    dist_cond_ss0_by_type = cell(cS.nTypes, 1);
    dist0_abs = zeros(cS0.nk, cS0.nkpps, cS0.nw_expanded, cS0.aD_new, cS.nTypes);

    for i_type = 1:cS.nTypes
        Z0_type_norm = dist0(:,:,:,:,i_type) / sum(dist0(:,:,:,:,i_type), 'all');
        dist_cond_ss0_by_type{i_type} = Z0_type_norm;
        
        Z0_abs_type = Z_path_raw(:,1) * cS.type_weights(i_type);
        dist0_abs(:,:,:,:,i_type) = dist_cond_ss0_by_type{i_type} .* reshape(Z0_abs_type, [1,1,1,cS.aD_new]);
    end

    % 准备用于转轨的 cS 结构体
    cS_for_trans = cSF;
    % 传递完整的路径
    cS_for_trans.theta_path_urban = cS.theta_path_urban;
    cS_for_trans.theta_path_resident = cS.theta_path_resident;
    cS_for_trans.s_pathV = cS.s_pathV; 

    data_for_transition = struct(...
        'ss0', ss0, 'ssF', ssF, ...
        'dist0', dist0_abs, ... % 使用与Z_path_raw对齐的5维绝对分布
        'distF', distF, 'polF', polF, 'valF', valF, ... 
        'cS', cS_for_trans, 'paramS0', paramS0, 'paramSF', paramSF);
    
    save(save_filename, 'data_for_transition', '-v7.3');
    fprintf('✅ 所有转轨所需数据已成功保存至:\n   %s\n', save_filename);
else
    fprintf('\n\n数据未保存，因为稳态求解失败。\n');
end


%% --- 6. [新增] 比较初始与终期稳态的国民核算结果 ---
fprintf('\n--- 6. 比较初始稳态(ss0)与终期稳态(ssF)的国民核算 ---\n');

% [核心] 调用已更新的、支持异质性类型的报告函数
% 格式: (ss, cS, paramS, Z_ss_norm, report_filename, print_flag, is_bgp_ss)
report0 = SS.display_national_accounts(ss0, cS0, paramS0, Z0_norm, report_file_initial, true, false);
reportF = SS.display_national_accounts(ssF, cSF, paramSF, ZF_norm, report_file_final, true, true);

% --- 定义要比较的变量和格式 ---
vars_to_compare = {
    '--- 宏观总量与比率 ---', '', '', '';
    '国内生产总值 (Y)',           'Y',          '%.4f',  'abs';
    '私人资本存量 (K_p)',       'K_p',        '%.4f',  'abs';
    '  其中: 常规资本 (K_k)',    'K_k',        '%.4f', 'abs';
    '  其中: PPS 资本 (K_pps)',  'K_pps',      '%.4f', 'abs';
    '有效劳动 (L)',                 'L',          '%.4f',  'abs';
    '私人消费 (C)',                 'C',          '%.4f',  'abs';
    '私人投资 (I_p, BGP)',        'I_p_bgp',    '%.4f',  'abs';
    '私人资本/产出比 (K_p/Y)',  'Kp_Y_ratio', '%.3f',  'diff';
    '--- 价格 ---', '', '', '';
    '年化真实利率 (r, %%)',       'r_annual',   '%.4f',  'diff';
    '真实工资率 (w_hat)',           'w',          '%.4f',  'abs';
    '--- 一般政府与PAYG (占Y的比例) ---', '', '', '';
    '一般税收收入 (Tax_General/Y)',  'Tax_General', '%.4f%%', 'ratio_y';
    '转移支付 (TR/Y)',             'TR',          '%.4f%%', 'ratio_y';
    'PAYG缴费 (Contrib/Y)',         'PAYG_contrib','%.4f%%', 'ratio_y';
    'PAYG福利 (Benefit/Y)',         'PAYG_benefit','%.4f%%', 'ratio_y';
    '--- PPS体系流量 (占Y的比例) ---', '', '', '';
    'PPS缴费 (Contrib/Y)',          'PPS_contrib',  '%.4f%%', 'ratio_y';
    'PPS提取 (Withdrawal/Y)',       'PPS_withdrawal','%.4f%%', 'ratio_y';
    'PPS提取税 (Tax/Y)',            'PPS_tax',      '%.4f%%', 'ratio_y';
    '--- 检验项 (占Y的比例) ---', '', '', '';
    '投资缺口 (%%)',              'Invest_Gap_pct','%.4f', 'diff';
    };

% --- 打印表头 ---
fprintf('\n\n========================================================================================================\n');
fprintf('###                          稳 态 核 心 变 量 对 比 (ss0 vs. ssF)                          ###\n');
fprintf('========================================================================================================\n');
fprintf('%-32s | %20s | %20s | %20s\n', '变量', '初始稳态 (ss0)', '终期稳态 (ssF)', '变化');
fprintf('%s\n', repmat('-', 1, 95));

% --- 循环打印表格内容 (增加稳健性检查) ---
for i = 1:size(vars_to_compare, 1)
    label = vars_to_compare{i, 1};
    var_name = vars_to_compare{i, 2};

    if isempty(var_name)
        fprintf('%s\n', label);
        continue;
    end

    % 稳健性检查：如果变量不存在，则设为NaN
    if isfield(report0, var_name), val0 = report0.(var_name); else, val0 = NaN; end
    if isfield(reportF, var_name), valF_val = reportF.(var_name); else, valF_val = NaN; end % Renamed to avoid conflict
    
    format_spec = vars_to_compare{i, 3};
    change_type = vars_to_compare{i, 4};

    % 根据类型计算变化值
    change_val_str = 'N/A';
    if ~isnan(val0) && ~isnan(valF_val)
        switch change_type
            case 'abs' % 绝对变化百分比
                if abs(val0) > 1e-9, change_val = (valF_val / val0 - 1) * 100; else, change_val = NaN; end
                if ~isnan(change_val), change_val_str = sprintf('%.2f %%', change_val); end
            case 'diff' % 百分点差异
                change_val = valF_val - val0;
                change_val_str = sprintf('%.4f pp', change_val);
            case 'ratio_y' % 占Y比重的百分点差异
                if isfield(report0, 'Y') && report0.Y ~= 0
                    val0_ratio = (val0 / report0.Y) * 100;
                else
                    val0_ratio = NaN;
                end
                if isfield(reportF, 'Y') && reportF.Y ~= 0
                    valF_ratio = (valF_val / reportF.Y) * 100;
                else
                    valF_ratio = NaN;
                end

                if ~isnan(val0_ratio) && ~isnan(valF_ratio)
                    change_val = valF_ratio - val0_ratio;
                    change_val_str = sprintf('%.4f pp', change_val);
                end
                val0 = val0_ratio;
                valF_val = valF_ratio;
        end
    end
    
    % 如果值为NaN，则打印N/A，否则按格式打印
    if isnan(val0)
        val0_str = 'N/A';
    else
        val0_str = sprintf(format_spec, val0);
    end
    
    if isnan(valF_val)
        valF_str = 'N/A';
    else
        valF_str = sprintf(format_spec, valF_val);
    end

    fprintf('%-32s | %20s | %20s | %20s\n', label, val0_str, valF_str, change_val_str);
end
fprintf('%s\n', repmat('-', 1, 95));

% 调用网格检查函数
utils.check_k_grid_limit(distF, cSF.kGridV, 'ssF', 'k_private');
if cSF.pps_active
    utils.check_k_grid_limit(dist0, cS0.kppsGridV, 'ss0', 'k_pps');
    utils.check_k_grid_limit(distF, cSF.kppsGridV, 'ssF', 'k_pps');
end

fprintf('\n--- 稳态数据准备与比较脚本执行完毕 ---\n');