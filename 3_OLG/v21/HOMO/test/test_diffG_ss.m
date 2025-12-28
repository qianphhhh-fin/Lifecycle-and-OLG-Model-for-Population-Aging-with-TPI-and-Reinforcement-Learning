% =========================================================================
% == SCRIPT: test_diffG_ss.m
% == 版本: [v1.0 - 技术冲击环境生成]
% ==
% == 目的:
% ==   1. 构建一个受控的技术冲击 (g_A shock) 实验环境。
% ==   2. 求解冲击前的初始稳态 (ss0)，其技术增长率为 g_A_ss_0。
% ==   3. 求解冲击后的终期稳态 (ssF)，其技术增长率为 g_A_ss_F。
% ==   4. 生成一个从 g_A_ss_0 平滑过渡到 g_A_ss_F 的技术路径 (A_path)，
% ==      以及一个保持恒定增长的人口路径 (Z_path)。
% ==   5. 保存所有TPI所需的数据，为 test_diffG_trans.m 做准备。
% =========================================================================
clear; close all; clear classes;
addpath(pwd);
fprintf('=== 受控技术动态(g_A Shock)稳态求解脚本 (v1.0) ===\n\n');

%% --- 1. 实验设置与参数初始化 ---
fprintf('--- 1. 实验设置与参数初始化 ---\n');
output_filename = 'SS/data_for_diff_G_transition.mat';

cS = utils.ParameterValues();
cS.pps_active = false;
cS.nk = 100; cS.nkpps = 1; cS.n_pps_rate_grid = 1;
cS = utils.generateGrids(cS);
[~,~,~, cS.nw_expanded] = utils.EarningProcess_AgeDependent(cS);
% cS.theta_path = 0.2;
% --- [核心设定] 技术增长率冲击 ---
g_A_ss0_annual = 0.01; % 初始稳态的技术年增长率
g_A_ssF_annual = 0.015; % 终期稳态的技术年增长率 (永久性提高)
n_ss_annual_constant = 0.01; % 人口年增长率保持不变

fprintf('   技术冲击设定: g_A 从 %.4f 永久性变为 %.4f\n', g_A_ss0_annual, g_A_ssF_annual);
fprintf('   人口设定: n_ss 保持恒定在 %.4f\n', n_ss_annual_constant);


%% --- 2. [核心] 构建外生路径 (恒定人口增长 + 技术路径转轨) ---
fprintf('\n--- 2. 构建外生路径 ---\n');
% --- 2a. 构建恒定增长的人口路径 ---
T_sim_guess = 200; % 先给一个较长的模拟期猜测
cS.T_sim = T_sim_guess;
s_pathV_ss = cS.s_pathV(:, end); % 使用稳态生存率
Z_ss_dist_const = population.compute_theoretical_ss_dist(s_pathV_ss, n_ss_annual_constant, cS.time_Step, cS.aD_new);
n_period_const = (1 + n_ss_annual_constant)^cS.time_Step - 1;

total_pop_path = (1 + n_period_const).^(0:T_sim_guess-1);
cS.Z_path_raw = Z_ss_dist_const * total_pop_path;
cS.Z_path = repmat(Z_ss_dist_const, 1, T_sim_guess);
cS.s_pathV = repmat(s_pathV_ss, 1, T_sim_guess);
cS.theta_path = repmat(cS.theta_path(end), 1, T_sim_guess);
fprintf('   ✅ 已构建基于 n_ss=%.4f 的恒定人口结构和增长路径 (T=%d)。\n', n_ss_annual_constant, T_sim_guess);

% --- 2b. 构建技术转轨路径 A_path ---
g_A_period_0 = (1 + g_A_ss0_annual)^cS.time_Step - 1;
g_A_period_F = (1 + g_A_ssF_annual)^cS.time_Step - 1;

transition_start_period = 10; % 技术冲击在第3期开始显现
transition_duration_periods = 12; % 技术增长率用10个时期(50年)完成过渡

g_A_period_path = ones(1, T_sim_guess) * g_A_period_0;
transition_end_period = transition_start_period + transition_duration_periods -1;
trans_path_g = utils.smooth_transition(g_A_period_0, g_A_period_F, transition_duration_periods, 1);
g_A_period_path(transition_start_period : transition_end_period) = trans_path_g;
g_A_period_path(transition_end_period+1 : end) = g_A_period_F;

cS.A_path = cumprod([1, 1 + g_A_period_path(1:end-1)]);
fprintf('   ✅ 已构建从 g_A=%.4f 过渡到 g_A=%.4f 的技术水平路径 A_path。\n', g_A_ss0_annual, g_A_ssF_annual);


%% --- 3. 求解初始稳态 (ss0) ---
fprintf('\n--- 3. 求解初始稳态 (ss0) @ g_A = %.4f ---\n', g_A_ss0_annual);
cS0 = cS;
cS0.g_A_ss = g_A_ss0_annual;
cS0.n_ss = n_ss_annual_constant;
cS0.s_pathV = s_pathV_ss;
Z0_norm = Z_ss_dist_const;

paramS0 = struct();
[paramS0.leGridV, paramS0.TrProbM_by_age, paramS0.leProb1V, cS0.nw_expanded] = utils.EarningProcess_AgeDependent(cS0);

params_for_ss0 = struct(...
    'Z', Z0_norm, 'A', 1.0, 'g_A_ss', cS0.g_A_ss, 'n_ss', cS0.n_ss);

fprintf('   ⚙️  启动稳态求解器 (求解 ss0)...\n');
[ss0, dist0_norm, ~, ~] = SS.solve_steady_state(cS0, paramS0, params_for_ss0, true, 'lsqnonlin');
if isempty(ss0), error('初始稳态(ss0)求解失败！'); else fprintf('✅ 初始稳态(ss0)求解成功！\n'); end

% [核心会计] 构造转轨第一期需要的绝对人口分布 dist0
dist_cond_ss0 = dist0_norm ./ reshape(Z0_norm, [1,1,1,cS.aD_new]);
dist0_abs = dist_cond_ss0 .* reshape(cS.Z_path_raw(:,1), [1,1,1,cS.aD_new]);


%% --- 4. 求解终期稳态 (ssF) ---
fprintf('\n--- 4. 求解终期稳态 (ssF) @ g_A = %.4f ---\n', g_A_ssF_annual);
cSF = cS;
cSF.g_A_ss = g_A_ssF_annual;
cSF.n_ss = n_ss_annual_constant;
cSF.s_pathV = s_pathV_ss;
ZF_norm = Z_ss_dist_const;

paramSF = struct();
[paramSF.leGridV, paramSF.TrProbM_by_age, paramSF.leProb1V, cSF.nw_expanded] = utils.EarningProcess_AgeDependent(cSF);

params_for_ssF = struct(...
    'Z', ZF_norm, 'A', 1.0, 'g_A_ss', cSF.g_A_ss, 'n_ss', cSF.n_ss);

fprintf('   ⚙️  启动稳态求解器 (求解 ssF)...\n');
[ssF, ~, polF, valF] = SS.solve_steady_state(cSF, paramSF, params_for_ssF, true, 'lsqnonlin');
if isempty(ssF), error('终期稳态(ssF)求解失败！'); else fprintf('✅ 终期稳态(ssF)求解成功！\n'); end


%% --- 5. 保存测试所需数据 ---
fprintf('\n--- 5. 保存转轨测试所需的数据 ---\n');
if ~exist('SS', 'dir'), mkdir('SS'); end

data_for_diff_G_transition = struct(...
    'ss0', ss0, 'ssF', ssF, ...
    'dist0', dist0_abs, ...
    'polF', polF, 'valF', valF, ...
    'cS', cS, ... % 保存包含完整路径的 cS
    'paramS0', paramS0, 'paramSF', paramSF);

save(output_filename, 'data_for_diff_G_transition', '-v7.3');
fprintf('✅ 所有技术冲击实验数据已成功保存至:\n   %s\n', output_filename);

fprintf('\n--- 受控技术动态(g_A Shock)稳态求解脚本执行完毕 ---\n');


%% --- 6. [新增] 比较初始与终期稳态的国民核算结果 ---
fprintf('\n--- 6. 比较初始稳态(ss0)与终期稳态(ssF)的国民核算 ---\n');
report0 = SS.display_national_accounts(ss0, cS0, paramS0, Z0_norm, '稳态报告_ss0.txt');
reportF = SS.display_national_accounts(ssF, cSF, paramSF, ZF_norm, '稳态报告_ssF.txt');

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
    if isfield(reportF, var_name), valF = reportF.(var_name); else, valF = NaN; end
    
    format_spec = vars_to_compare{i, 3};
    change_type = vars_to_compare{i, 4};

    % 根据类型计算变化值
    change_val_str = 'N/A';
    if ~isnan(val0) && ~isnan(valF)
        switch change_type
            case 'abs' % 绝对变化百分比
                if abs(val0) > 1e-9, change_val = (valF / val0 - 1) * 100; else, change_val = NaN; end
                if ~isnan(change_val), change_val_str = sprintf('%.2f %%', change_val); end
            case 'diff' % 百分点差异
                change_val = valF - val0;
                change_val_str = sprintf('%.4f pp', change_val);
            case 'ratio_y' % 占Y比重的百分点差异
                if isfield(report0, 'Y') && report0.Y ~= 0
                    val0_ratio = (val0 / report0.Y) * 100;
                else
                    val0_ratio = NaN;
                end
                if isfield(reportF, 'Y') && reportF.Y ~= 0
                    valF_ratio = (valF / reportF.Y) * 100;
                else
                    valF_ratio = NaN;
                end

                if ~isnan(val0_ratio) && ~isnan(valF_ratio)
                    change_val = valF_ratio - val0_ratio;
                    change_val_str = sprintf('%.4f pp', change_val);
                end
                val0 = val0_ratio;
                valF = valF_ratio;
        end
    end
    
    % 如果值为NaN，则打印N/A，否则按格式打印
    if isnan(val0)
        val0_str = 'N/A';
    else
        val0_str = sprintf(format_spec, val0);
    end
    
    if isnan(valF)
        valF_str = 'N/A';
    else
        valF_str = sprintf(format_spec, valF);
    end

    fprintf('%-32s | %20s | %20s | %20s\n', label, val0_str, valF_str, change_val_str);
end
fprintf('%s\n', repmat('-', 1, 95));

fprintf('\n--- test_diffG_ss.m 脚本执行完毕 ---\n');