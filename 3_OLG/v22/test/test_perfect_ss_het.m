% =========================================================================
% == SCRIPT: test_perfect_ss_het.m
% == 版本: [v1.0 - 异质性黄金标准生成]
% ==
% == 目的:
% ==   1. 创建两个内部完全自洽的异质性稳态 (ss0_perfect, ssF_perfect)。
% ==   2. 模拟一个从低缴费率(ss0)到高缴费率(ssF)的政策冲击。
% ==   3. 确保两个稳态之间的人口结构、生存率等都是理论上一致的。
% ==   4. 保存所有数据，为 test_perfect_trans_het.m 提供干净的输入。
% =========================================================================
clear; close all; clear classes;
addpath(pwd); 
fprintf('=== 异质性模型理想化环境稳态求解脚本 (v1.0) ===\n\n');

%% --- 1. 全局设置 ---
output_data_filename = 'SS/data_for_perfect_het_transition.mat';
report_file_initial = 'SS/SS0_perfect_het_report.txt';
report_file_final = 'SS/SSF_perfect_het_report.txt';

fprintf('   加载模型物理参数...\n');
cS = utils.ParameterValues();

% 禁用PPS以简化问题
cS.pps_active = false;
cS.nkpps = 1;
cS.nk = 50;
% --- [!!! 核心异质性设定 !!!] ---
cS.num_hh_types = 4;
nH = cS.num_hh_types;

%% --- 2. 构建理想化的宏观与人口环境 ---
fprintf('\n--- 2. 构建理想化的宏观与人口环境 ---\n');
% 这些参数在 ss0 和 ssF 之间保持不变
cS_common = cS; 
cS_common.g_A_ss = 0.015;
cS_common.n_ss = 0.01; % 设定一个统一的、不变的人口年增长率
fprintf('   [核心设定] 理想环境的人口年增长率 n_ss 被强制设为: %.4f\n', cS_common.n_ss);

% 使用最后一年（稳态）的生存率数据作为统一生存率
s_pathV_ss = cS_common.s_pathV(:, end);
cS_common.s_pathV = s_pathV_ss;

% 根据生存率和增长率，计算一个完全自洽的理论稳态人口分布
Z_theory = population.compute_theoretical_ss_dist(cS_common.s_pathV, cS_common.n_ss, cS_common.time_Step, cS_common.aD_new);
% 将其扩展到异质性维度
Z_theory_h = Z_theory * cS_common.type_weights'; % [aD x 1] * [1 x nH] -> [aD x nH]
fprintf('   ✅ 已根据 n_ss 和 s_pathV 计算出通用的理论人口分布 Z_theory_h。\n');

% 通用的劳动过程参数
paramS_common = struct();
[paramS_common.leGridV, paramS_common.TrProbM_by_age, paramS_common.leProb1V, cS_common.nw_expanded] = utils.EarningProcess_AgeDependent(cS_common);

%% --- 3. 求解初始完美稳态 (ss0_perfect) ---
fprintf('\n\n--- 3. 求解初始完美稳态 (ss0_perfect) ---\n');
cS0 = cS_common;
% [政策设定] 初始稳态使用较低的养老金缴费率
theta_low = 0.12;
cS0.theta_path_h = repmat(theta_low, cS0.nTypes, 1);
fprintf('   [政策设定] ss0 的所有类型theta统一设为: %.2f\n', theta_low);
cS0 = utils.generateGrids(cS0);

params_for_ss0 = struct(...
    'Z', Z_theory_h, ...
    'A', 1.0, ...
    'g_A_ss', cS0.g_A_ss, ...
    'n_ss', cS0.n_ss);

tic;
fprintf('   ⚙️  启动稳态求解器 (求解 ss0)...\n');
[ss0_perfect, dist0_perfect_h, pol0_perfect_h, ~] = SS.solve_steady_state(cS0, paramS_common, params_for_ss0, true, true);
toc;
if isempty(ss0_perfect), error('初始稳态(ss0)求解失败！'); end
SS.display_national_accounts(ss0_perfect, cS0, paramS_common, Z_theory_h, report_file_initial, true, true);


% --- 这是 test_vfi_preservation.m 的新第4和第5部分 ---

%% --- 4. 求解终期完美稳态 (ssF_perfect) 作为基准 ---
fprintf('\n\n--- 4. 求解终期完美稳态 (ssF_perfect) 作为基准 ---\n');
cSF = cS_common;
theta_high = 0.12; 
cSF.theta_path_h = repmat(theta_high, cSF.nTypes, 1);
cSF = utils.generateGrids(cSF);

params_for_ssF = struct(...
    'Z', Z_theory_h, ...
    'A', 1.0, ...
    'g_A_ss', cSF.g_A_ss, ...
    'n_ss', cSF.n_ss);

tic;
fprintf('   ⚙️  启动稳态求解器 (求解 ssF)...\n');
[ssF_perfect, ~, polF_perfect_h, valF_perfect_h] = SS.solve_steady_state(cSF, paramS_common, params_for_ssF, true, true);
toc;
if isempty(ssF_perfect), error('终期稳态(ssF)求解失败！'); end

% --- 这是 test_perfect_ss_het.m 的第 5 部分 ---

% --- 这是 test_perfect_ss_het.m 的第 5 部分 (已修正维度错误) ---

%% --- 5. 保存所有必要数据到 .mat 文件 ---
fprintf('\n\n--- 5. 保存结果用于完美转轨测试 ---\n');

% 构建一个在转轨期间保持不变的 cS 结构体
cS_for_trans = cSF; % 以终期为基础
T_trans = 160; % 设定一个转轨时长
cS_for_trans.T_sim = T_trans;

% [!!! 关键修正: 构造一个完全稳态的人口和技术路径 !!!]
% 1. 生存率路径: 保持稳态生存率不变
cS_for_trans.s_pathV = repmat(cS_common.s_pathV, 1, T_trans); 

% 2. 养老金缴费率路径: 保持稳态缴费率不变 (因为ss0和ssF的theta相同)
cS_for_trans.theta_path_h = repmat(cSF.theta_path_h, 1, T_trans);

% 3. 绝对人口路径: 必须与稳态的定义完全一致
%    Z_theory 是 [aD x 1] 的归一化总人口年龄结构
%    我们首先计算 t=1 时的绝对人口分布
Z_theory_total_abs_t1 = Z_theory * sum(Z_theory_h, 'all');

%    然后计算每一期的增长因子路径
g_n_period = (1 + cS_common.n_ss)^cS_common.time_Step - 1;
pop_growth_path = (1 + g_n_period) .^ (0:(T_trans-1)); % size [1 x T]

% [!!! 维度修正 !!!]
%    使用矩阵乘法将 [aD x 1] 的分布乘以 [1 x T] 的增长路径，得到 [aD x T] 的绝对人口路径
Z_path_raw_perfect = Z_theory_total_abs_t1 * pop_growth_path; 
cS_for_trans.Z_path_raw = Z_path_raw_perfect;
fprintf('   [核心修正] 已构造一个与BGP完全一致的绝对人口路径。\n');


% 4. 技术路径: 必须与BGP一致
g_A_period = (1 + cS_common.g_A_ss)^cS_common.time_Step - 1;
A_path_perfect = (1 + g_A_period) .^ (0:(T_trans-1)); 
cS_for_trans.A_path = A_path_perfect;

data_for_perfect_transition = struct();
data_for_perfect_transition.ss0=ss0_perfect;
data_for_perfect_transition.ssF=ssF_perfect;
data_for_perfect_transition.dist0_h=dist0_perfect_h; 
data_for_perfect_transition.valF_h=valF_perfect_h;
data_for_perfect_transition.polF_h=polF_perfect_h;
data_for_perfect_transition.cS=cS_for_trans;
data_for_perfect_transition.paramSF=paramS_common;

save(output_data_filename, 'data_for_perfect_transition', '-v7.3');
fprintf('✅ 所有完美转轨所需的数据已成功保存至: %s\n', output_data_filename);

fprintf('\n--- 理想化环境稳态求解脚本执行完毕 ---\n');