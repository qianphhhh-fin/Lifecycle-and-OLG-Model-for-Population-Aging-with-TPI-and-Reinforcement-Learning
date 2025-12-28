% =========================================================================
% == SCRIPT: test_falseSS_ss.m
% == 版本: [v1.0 - 伪稳态实验构建]
% ==
% == 目的:
% ==   - 创建一个受控实验，其转轨路径的起点是一个【伪稳态】，终点是
% ==     一个【真稳态】，用于模拟更真实的转轨场景。
% ==   - 伪稳态 (ss0): 使用一个【人为构造的非均衡人口分布】进行求解。
% ==   - 真稳态 (ssF): 使用理论上的稳态人口分布进行求解。
% ==   - 保存所有结果，为 test_falseSS_trans.m 提供输入。
% =========================================================================
clear; close all; clc;
addpath(pwd);
fprintf('=== 伪稳态 -> 真稳态 转轨实验数据准备脚本 (v1.0) ===\n\n');

%% --- 1. 全局设置与初始化 ---
fprintf('--- 1. 全局设置与初始化 ---\n');

fprintf('   加载模型物理参数...\n');
cS = utils.ParameterValues();
% 设定模拟的起始年份
cS.ss0_year = 2023;
cS.start_year = 2023;
cS.pop_data_last_year = 2023;


% 为了测试，使用较小的网格
cS.nk = 40;
cS.g_A_ss = 0.01; % 设定一个统一的技术进步率
cS.time_Step = 5;

% 设定报告文件名
report_file_initial = 'SS/伪稳态报告_ss0.txt';
report_file_final = 'SS/伪稳态实验_真稳态报告_ssF.txt';
save_filename = 'SS/data_for_falseSS_transition.mat';

% 生成一个时不变的生存率路径和一个基础的劳动过程
% 注意：我们将为ss0和ssF使用相同的生存率和劳动过程以隔离人口分布的影响
cS.pop_data_last_year = 2023; % 假设死亡率不再变化
[~, ~, cS] = population.generate_Z_path(cS, false); 
cS.s_pathV_ss = cS.s_pathV(:,1); % 取第一期的生存率作为基准
[paramS.leGridV, paramS.TrProbM_by_age, paramS.leProb1V, cS.nw_expanded] = utils.EarningProcess_AgeDependent(cS);


%% --- 2. 求解伪稳态 (ss0) ---
fprintf('\n\n--- 2. 求解伪稳态 (ss0) ---\n');

cS0 = cS;
cS0.pps_active = false;
cS0.nkpps = 1;
cS0 = utils.generateGrids(cS0);
cS0.theta_path = cS.theta_path(1);

% [!!! 核心: 人为构造一个非稳态初始人口分布 !!!]
% 例如，我们构造一个非常年轻化的人口结构
fprintf('   特征: 使用人为构造的、非均衡的年轻化人口分布。\n');
Z0_fake = ones(cS0.aD_new, 1);
Z0_fake(1:5) = linspace(0.15, 0.08, 5); % 年轻人占比很高
Z0_fake = Z0_fake / sum(Z0_fake); % 归一化

% [!!! 核心: 计算这个非稳态人口对应的“伪”增长率 !!!]
% 我们需要模拟一期来得到下一期的人口，从而计算增长率
Z1_fake_survivors = zeros(cS0.aD_new, 1);
Z1_fake_survivors(2:end) = Z0_fake(1:end-1) .* cS0.s_pathV_ss(1:end-1);
% 假设一个粗出生率来关闭模型，例如等于总死亡率
total_death_rate = 1 - sum(Z1_fake_survivors);
Z1_fake_newborns = total_death_rate;
Z1_fake = [Z1_fake_newborns; Z1_fake_survivors];
Z1_fake = Z1_fake / sum(Z1_fake);
% 计算这一期的伪人口增长率
n_ss0_fake_period = sum(Z1_fake) / sum(Z0_fake) - 1; 
cS0.n_ss = (1 + n_ss0_fake_period)^(1/cS0.time_Step) - 1;
fprintf('   计算出的伪稳态年化增长率 n_ss = %.4f%%\n', cS0.n_ss * 100);

params_for_ss0 = struct(...
    'Z', Z0_fake, ...
    'A', 1.0, ...
    'g_A_ss', cS0.g_A_ss, ...
    'n_ss', cS0.n_ss);

tic;
fprintf('   ⚙️  启动统一稳态求解器 (求解伪稳态 ss0)...\n');
[ss0, dist0, ~, ~] = SS.solve_steady_state(cS0, paramS, params_for_ss0, true, 'lsqnonlin');
toc;

if isempty(ss0)
    error('伪稳态(ss0)求解失败！');
else
    fprintf('✅ 伪稳态(ss0)求解成功！\n');
    SS.display_national_accounts(ss0, cS0, paramS, Z0_fake, report_file_initial);
end


%% --- 3. 求解终期真稳态 (ssF) ---
fprintf('\n\n--- 3. 求解终期BGP真稳态 (ssF) ---\n');

cSF = cS;
cSF.pps_active = false; % 假设终点也无PPS，以隔离冲击
cSF.nkpps = 1;
cSF = utils.generateGrids(cSF);
cSF.theta_path = cS.theta_path(end); % 使用不同的养老金费率

% [!!! 核心: 使用理论稳态人口分布 !!!]
% 假设一个长期的、较低的人口年增长率
n_ssF_annual = -0.005; % 例如，长期人口年均萎缩0.5%
cSF.n_ss = n_ssF_annual;
[ZF_theory, ~] = population.compute_theoretical_ss_dist(cSF.s_pathV_ss, cSF.n_ss, cSF.time_Step, cSF.aD_new);
fprintf('   特征: 使用理论计算的BGP稳态人口分布 (长期 n = %.4f%%)。\n', cSF.n_ss * 100);

params_for_ssF = struct(...
    'Z', ZF_theory, ...
    'A', 1.0, ...
    'g_A_ss', cSF.g_A_ss, ...
    'n_ss', cSF.n_ss);

tic;
fprintf('   ⚙️  启动统一稳态求解器 (求解真稳态 ssF)...\n');
[ssF, distF, polF, valF] = SS.solve_steady_state(cSF, paramS, params_for_ssF, true, 'lsqnonlin');
toc;

if isempty(ssF)
    error('真稳态(ssF)求解失败！');
else
    fprintf('✅ 真稳态(ssF)求解成功！\n');
    SS.display_national_accounts(ssF, cSF, paramS, ZF_theory, report_file_final);
end


%% --- 4. 保存转轨所需数据 ---
if ~isempty(fieldnames(ssF))
    fprintf('\n\n--- 4. 保存转轨分析所需数据 ---\n');
    if ~exist('SS', 'dir'), mkdir('SS'); end
    
    % [!!! 核心: 构造初始绝对人口分布 !!!]
    % 假设转轨的起点总人口为1000 (或任何标准化数值)
    total_pop_t1 = 1000; 
    dist0_abs = dist0 * total_pop_t1;
    
    % 在cS中设置一个模拟的人口路径，简单地从Z0_fake线性过渡到ZF_theory
    T_sim_fake = 40; % 假设一个40期的转轨时长
    cS.T_sim = T_sim_fake;
    Z_path_fake = zeros(cS.aD_new, T_sim_fake);
    for ia = 1:cS.aD_new
        Z_path_fake(ia,:) = linspace(Z0_fake(ia), ZF_theory(ia), T_sim_fake);
    end
    cS.Z_path = Z_path_fake;
    cS.Z_path_raw = Z_path_fake * total_pop_t1; % 假设总人口保持不变，仅结构变化
    cS.A_path = (1 + (1+cS.g_A_ss)^cS.time_Step - 1).^(0:(T_sim_fake-1));

    
    data_for_transition = struct(...
        'ss0', ss0, 'ssF', ssF, ...
        'dist0', dist0_abs, ...
        'distF', distF, 'polF', polF, 'valF', valF, ...
        'cS', cSF, ... % 保存与终点一致的参数
        'paramS0', paramS, 'paramSF', paramS); % 假设劳动过程不变
    
    save(save_filename, 'data_for_false_SS_transition', '-v7.3');
    fprintf('✅ 所有伪稳态转轨所需数据已成功保存至:\n   %s\n', save_filename);
else
    fprintf('\n\n数据未保存，因为终期稳态求解失败。\n');
end