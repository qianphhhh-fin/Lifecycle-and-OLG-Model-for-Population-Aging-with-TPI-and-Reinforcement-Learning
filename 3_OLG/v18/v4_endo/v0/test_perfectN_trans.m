% =========================================================================
% == SCRIPT: test_perfectN_trans.m
% == 版本: [v1.8 - 遗赠路径索引修正最终版]
% ==
% == 目的:
% ==   - 修正了 TPI 核心引擎的调用，使其与 v4.1 版本的新接口匹配。
% ==   - 确保了在最完美的 BGP 环境下，所有计算环节的时间索引都是自洽的。
% =========================================================================
clear; close all; clear classes;
addpath(pwd);
fprintf('=== TPI求解器零冲击诊断脚本 (v1.8 - 遗赠路径索引修正最终版) ===\n\n');

%% --- 1. 加载“黄金标准”稳态数据 ---
fprintf('--- 1. 加载黄金标准稳态数据 ---\n');
data_filename = 'SS/data_for_perfect_transition.mat';
if ~exist(data_filename, 'file'), error('黄金标准稳态数据文件 "%s" 不存在。', data_filename); end
load(data_filename, 'data_for_perfect_transition');
fprintf('   ✅ 已加载数据: %s\n', data_filename);

ssF = data_for_perfect_transition.ssF;
distF = data_for_perfect_transition.distF;
valF = data_for_perfect_transition.valF;
polF = data_for_perfect_transition.polF;
cS = data_for_perfect_transition.cS;
paramSF = data_for_perfect_transition.paramSF;

cS.T_sim = 20;
T = cS.T_sim;
fprintf('   [快速诊断模式] T已被临时设为: %d\n', T);

%% --- 2. 手动构建一个完美的、自洽的、平稳增长的人口环境 ---
fprintf('\n--- 2. 手动构建完美的、平稳增长的BGP环境 ---\n');
ss0 = ssF;
dist0 = distF;
fprintf('   ✅ 初始稳态(ss0)和分布(dist0)已设为与终期稳态(ssF)相同。\n');

Z_ss_abs = sum(distF, [1,2,3]); 
Z_ss_abs = Z_ss_abs(:);
g_n_period = (1 + cS.n_ss)^cS.time_Step - 1;
pop_growth_factors = (1 + g_n_period).^(0:T-1);
cS.Z_path_raw = Z_ss_abs * pop_growth_factors;
cS.Z_path = cS.Z_path_raw ./ sum(cS.Z_path_raw, 1);
fprintf('   ✅ 总人口路径 (Z_path_raw) 已被构建为以 n_ss=%.4f 恒定增长。\n', cS.n_ss);

s_pathV_ss = cS.s_pathV; 
if size(s_pathV_ss, 2) > 1, s_pathV_ss = s_pathV_ss(:, end); end
cS.s_pathV = repmat(s_pathV_ss, 1, T);

g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
A_path_flat = (1 + g_A_period).^(0:T-1);
cS.A_path = A_path_flat;

if isfield(cS, 'theta_path'), cS.theta_path = repmat(cS.theta_path(end), 1, T);
else, cS.theta_path = zeros(1, T); end

%% --- 3. 构造初始猜测：平坦的人均/价格路径 ---
fprintf('\n--- 3. 构造平坦的人均/价格初始猜测路径 ---\n');
r_path_guess = ones(1, T) * ss0.r_mkt;
w_hat_path_guess = ones(1, T) * ss0.w_hat;
b_hat_path_guess = ones(1, T) * ss0.b_hat;

total_pop0 = sum(cS.Z_path_raw(:,1));
tr0_per_hh_hat = ss0.TR_distributed_agg / total_pop0;
TR_per_hh_path_guess = ones(1, T) * tr0_per_hh_hat;

mass_newborns0 = cS.Z_path_raw(1,1);
beq0_per_newborn_hat = ss0.Bequest_generated_agg / mass_newborns0; % [修正] 使用产生的遗赠来初始化
Bequest_per_newborn_path_guess = ones(1, T) * beq0_per_newborn_hat;
fprintf('   ✅ 所有价格和人均路径已设为平坦的稳态值。\n');

%% --- 4. 运行TPI循环 ---
fprintf('\n--- 4. 启动TPI循环进行诊断 ---\n');
max_iter = 10;
tolerance = 1e-7;
damping.r = 0.2; damping.w = 0.2; damping.beq = 0.2; damping.tr = 0.2; damping.b_hat = 0.2;
fprintf('   阻尼因子: r=%.2f, w=%.2f, beq=%.2f, tr=%.2f, b_hat=%.2f\n', ...
    damping.r, damping.w, damping.beq, damping.tr, damping.b_hat);

r_path_iter = r_path_guess;
w_hat_path_iter = w_hat_path_guess;
Bequest_per_newborn_path_iter = Bequest_per_newborn_path_guess;
TR_per_hh_path_iter = TR_per_hh_path_guess;
b_hat_path_iter = b_hat_path_guess;

for iter = 1:max_iter
    [target_paths, errors, ~, ~, ~] = TPI.calculate_paths_and_errors(...
        r_path_iter, w_hat_path_iter, Bequest_per_newborn_path_iter, ...
        TR_per_hh_path_iter, b_hat_path_iter, ...
        ss0, ssF, cS, paramSF, valF, polF, dist0);

    fprintf('Iter [%2d/%2d] | Error: %.3e (r:%.2e, w:%.2e, bq:%.2e, tr:%.2e, b_hat:%.2e)\n', ...
        iter, max_iter, errors.total, errors.r, errors.w, errors.beq, errors.tr, errors.b_hat);

    if errors.total < tolerance, fprintf('\n✅ 诊断成功！\n'); break; end
    if iter == max_iter, fprintf('\n❌ 诊断失败！\n'); end

    % 阻尼更新
    r_path_iter(2:T) = (1 - damping.r) * r_path_iter(2:T) + damping.r * target_paths.r_path(2:T);
    w_hat_path_iter(2:T) = (1 - damping.w) * w_hat_path_iter(2:T) + damping.w * target_paths.w_hat_path(2:T);
    Bequest_per_newborn_path_iter(2:T) = (1 - damping.beq) * Bequest_per_newborn_path_iter(2:T) + damping.beq * target_paths.Bequest_per_newborn_path(2:T);
    TR_per_hh_path_iter(2:T) = (1 - damping.tr) * TR_per_hh_path_iter(2:T) + damping.tr * target_paths.TR_per_hh_path(2:T);
    b_hat_path_iter(2:T) = (1 - damping.b_hat) * b_hat_path_iter(2:T) + damping.b_hat * target_paths.b_hat_path(2:T);
end

%% --- 5. 结果分析 ---
if errors.total < tolerance
    fprintf('\n--- 5. 最终路径检验 ---\n');
    max_dev_r = max(abs(target_paths.r_path - ssF.r_mkt));
    max_dev_w = max(abs(target_paths.w_hat_path - ssF.w_hat));
    fprintf('   利率路径与稳态的最大偏差: %.3e\n', max_dev_r);
    fprintf('   工资路径与稳态的最大偏差: %.3e\n', max_dev_w);
    if max_dev_r < 1e-5 && max_dev_w < 1e-5
        fprintf('   ✅ 结论：求解出的路径确实是平坦的稳态路径。TPI核心逻辑通过检验。\n');
    else
        fprintf('   ⚠️ 警告：路径收敛了，但并非平坦路径。可能存在微小的不一致性。\n');
    end
end

fprintf('\n--- 黄金标准诊断脚本执行完毕 ---\n');