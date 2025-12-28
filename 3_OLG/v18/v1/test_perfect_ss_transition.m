% =========================================================================
% == SCRIPT: test_perfect_ss_transition.m
% == 版本: [v1.5 - 终极自洽环境构建版]
% ==
% == 目的:
% ==   - 彻底根除因外生数据不一致导致的人口动态问题。
% ==   - 不再信任加载数据中的人口路径，而是基于稳态值手动构建一个
% ==     完美的、平稳的、内部完全自洽的人口环境。
% =========================================================================
clear; close all; clear classes;
addpath(pwd);
fprintf('=== TPI求解器零冲击诊断脚本 (v1.5 - 终极自洽环境构建版) ===\n\n');

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

cS.T_sim = 10;
T = cS.T_sim;
fprintf('   [快速诊断模式] T已被临时设为: %d\n', T);

%% --- 2. [核心修正] 手动构建一个完美的、自洽的、平稳的人口环境 ---
fprintf('\n--- 2. 手动构建完美的零冲击环境 ---\n');
ss0 = ssF;
dist0 = distF;
fprintf('   ✅ 初始稳态(ss0)和分布(dist0)已设为与终期稳态(ssF)相同。\n');

% a. 提取稳态的、自洽的人口分布（绝对数量）和生存率
Z_ss_abs = sum(distF, [1,2,3]); Z_ss_abs = Z_ss_abs(:);
s_pathV_ss = cS.s_pathV; % 从黄金标准文件加载的cS.s_pathV现在是单一的稳态向量
if size(s_pathV_ss, 2) > 1
    warning('黄金标准文件中的cS.s_pathV不是一个向量，正在使用其最后一列作为稳态值。');
    s_pathV_ss = s_pathV_ss(:, end);
end

% b. 强制将整个模拟期间的人口分布和生存率路径设为恒定的稳态值
cS.Z_path_raw = repmat(Z_ss_abs, 1, T);
cS.Z_path = cS.Z_path_raw ./ sum(cS.Z_path_raw, 1);
cS.s_pathV = repmat(s_pathV_ss, 1, T);
fprintf('   ✅ 人口路径 (Z_path_raw) 和生存率路径 (s_pathV) 已被强制设为恒定的、自洽的稳态值。\n');

% c. 构建其他平稳路径
g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
A_path_flat = (1 + g_A_period).^(0:T-1);
cS.A_path = A_path_flat / A_path_flat(1);
fprintf('   ✅ 技术路径 (A_path) 已被设为以 g_A_ss 恒定增长。\n');

if isfield(cS, 'theta_path'), cS.theta_path = repmat(cS.theta_path(end), 1, T);
else, cS.theta_path = zeros(1, T); end
fprintf('   ✅ 政策路径 (theta_path) 已被设为常数。\n');

%% --- 2b. [新增] 前置稳态聚合一致性检验 ---
fprintf('\n--- 2b. 前置稳态聚合一致性检验 ---\n');

% 使用稳态策略(polF)和稳态分布(distF)手动聚合下一期的总资本
K_p_eop_from_ss_objects = 0;
g_n_ss_period = (1 + cS.n_ss)^cS.time_Step - 1;

for ia = 1:cS.aD_new
    k_prime_slice_ss = polF(ia).k_prime;
    if cS.pps_active
        k_prime_slice_ss = k_prime_slice_ss + polF(ia).kpps_prime;
    end
    mass_ia_ss_abs = sum(distF(:,:,:,ia), 'all');
    K_p_eop_from_ss_objects = K_p_eop_from_ss_objects + sum(k_prime_slice_ss .* distF(:,:,:,ia), 'all');
end

% 根据BGP关系，计算稳态资本存量 K_hat
K_hat_recomputed_from_ss = K_p_eop_from_ss_objects / (1 + g_n_ss_period);

fprintf('   从 ssF 结构体中读取的稳态资本 (K_private_hat): %.6f\n', ssF.K_private_hat);
fprintf('   根据 polF 和 distF 手动重新计算的稳态资本: %.6f\n', K_hat_recomputed_from_ss);

recomputation_error = abs(K_hat_recomputed_from_ss - ssF.K_private_hat);
fprintf('   => 两者之间的绝对误差: %.4e\n', recomputation_error);

if recomputation_error > 1e-6
    error('稳态对象不一致！ssF.K_private_hat 与从distF和polF聚合出的值不匹配。请重新运行 test_perfectN_SS.m 生成黄金标准文件。');
else
    fprintf('   ✅ 稳态对象一致性检验通过。\n');
end

%% --- 3. 构造初始猜测：一条完美的平坦路径 ---
fprintf('\n--- 3. 构造完美的平坦路径作为初始猜测 ---\n');
r_path_guess = ones(1, T) * ss0.r_mkt;
w_hat_path_guess = ones(1, T) * ss0.w_hat;
Bequest_hat_path_guess = ones(1, T) * ss0.Bequest_distributed_agg;
TR_hat_path_guess = ones(1, T) * ss0.TR_distributed_agg;
fprintf('   ✅ 初始猜测路径已设为稳态值。\n');

%% --- 4. 运行TPI循环 ---
% (TPI循环代码保持不变，但现在运行在一个真正自洽的环境中)
fprintf('\n--- 4. 启动TPI循环进行诊断 ---\n');
max_iter = 10;
tolerance = 1e-7;
damping.r = 0.1; damping.w = 0.1; damping.beq = 0.1; damping.tr = 0.1;
fprintf('   阻尼因子: r=%.2f, w=%.2f, beq=%.2f, tr=%.2f\n', damping.r, damping.w, damping.beq, damping.tr);

r_path_iter = r_path_guess;
w_hat_path_iter = w_hat_path_guess;
Bequest_hat_path_iter = Bequest_hat_path_guess;
TR_hat_path_iter = TR_hat_path_guess;

for iter = 1:max_iter
    [target_paths, errors, Pol_path, Dist_path, aggr_supply] = TPI.calculate_paths_and_errors(...
        r_path_iter, w_hat_path_iter, Bequest_hat_path_iter, TR_hat_path_iter, ...
        ss0, ssF, cS, paramSF, valF, polF, dist0);

    fprintf('Iter [%2d/%2d] | Error: %.3e (r:%.2e, w:%.2e, bq:%.2e, tr:%.2e)\n', ...
        iter, max_iter, errors.total, errors.r, errors.w, errors.beq, errors.tr);

    fprintf('   [DEBUG] Iter %d:\n', iter);
    max_c_dev = 0; max_k_dev = 0;
    for t = 1:(T-1)
        c_dev = max(abs(Pol_path{t}(1).c(:) - Pol_path{t+1}(1).c(:)));
        k_dev = max(abs(Pol_path{t}(1).k_prime(:) - Pol_path{t+1}(1).k_prime(:)));
        if c_dev > max_c_dev, max_c_dev = c_dev; end
        if k_dev > max_k_dev, max_k_dev = k_dev; end
    end
    fprintf('      -> Policy Func Consistency: Max C dev(t,t+1)=%.3e | Max K'' dev(t,t+1)=%.3e\n', max_c_dev, max_k_dev);
    pop_mass_sim = squeeze(sum(Dist_path, [1,2,3,4]));
    pop_mass_exo = sum(cS.Z_path_raw, 1)';
    mass_dev = max(abs(pop_mass_sim - pop_mass_exo));
    fprintf('      -> Distribution Mass Conservation: Max Abs Pop dev=%.3e\n', mass_dev);
    fprintf('         - Simulated Pop Path: '); fprintf('%.4f ', pop_mass_sim); fprintf('\n');
    fprintf('         - Exogenous Pop Path: '); fprintf('%.4f ', pop_mass_exo); fprintf('\n');
    fprintf('      -> Agg. Private K_hat Path (Supply): '); fprintf('%.4f ', aggr_supply.K_p_hat_path); fprintf('\n');
    fprintf('      -> Agg. Labor Path (Supply)       : '); fprintf('%.4f ', aggr_supply.L_path); fprintf('\n');
    fprintf('      -> Target Interest Rate Path (r)    : '); fprintf('%.4f ', target_paths.r_path); fprintf('\n');
    fprintf('      -> Target Wage_hat Path (w_hat)   : '); fprintf('%.4f ', target_paths.w_hat_path); fprintf('\n');
    fprintf('   ------------------------------------------------------------\n');

    if errors.total < tolerance, fprintf('\n✅ 诊断成功！\n'); break; end
    if iter == max_iter, fprintf('\n❌ 诊断失败！\n'); end

    r_path_iter(2:T) = (1 - damping.r) * r_path_iter(2:T) + damping.r * target_paths.r_path(2:T);
    w_hat_path_iter(2:T) = (1 - damping.w) * w_hat_path_iter(2:T) + damping.w * target_paths.w_hat_path(2:T);
    Bequest_hat_path_iter(2:T) = (1 - damping.beq) * Bequest_hat_path_iter(2:T) + damping.beq * target_paths.Bequest_hat_path(2:T);
    TR_hat_path_iter(2:T) = (1 - damping.tr) * TR_hat_path_iter(2:T) + damping.tr * target_paths.TR_hat_path(2:T);
end


%% --- 5. 结果分析 ---
% (代码保持不变)
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
else
    fprintf('\n--- 5. 失败分析 ---\n');
    fprintf('   既然在最理想的环境下都无法收敛，问题几乎100%%在于TPI的某个计算环节。\n');
    % figure('Name', '黄金标准诊断：失败的路径', 'Position', [200, 200, 1000, 400]);
    % time_axis = 1:T;
    % subplot(1,2,1); plot(time_axis, r_path_iter, 'r-o', 'DisplayName', '迭代路径 (r)'); hold on; plot(time_axis, target_paths.r_path, 'b--x', 'DisplayName', '目标路径 (r)'); yline(ssF.r_mkt, 'k:', 'DisplayName', '理论稳态 (r)'); title('利率路径'); legend; grid on;
    % subplot(1,2,2); plot(time_axis, w_hat_path_iter, 'r-o', 'DisplayName', '迭代路径 (w_hat)'); hold on; plot(time_axis, target_paths.w_hat_path, 'b--x', 'DisplayName', '目标路径 (w_hat)'); yline(ssF.w_hat, 'k:', 'DisplayName', '理论稳态 (w_hat)'); title('工资路径'); legend; grid on;
    % sgtitle('TPI在黄金标准环境下的发散行为');
end

fprintf('\n--- 黄金标准诊断脚本执行完毕 ---\n');