% =========================================================================
% == SCRIPT: test_perfectN_trans.m
% == 版本: [v2.4 - BGP猜测最终修正版]
% ==
% == 目的:
% ==   - 修正了对 TR (转移支付) 初始猜测路径的构造错误，这是导致
% ==     Pol_Err 和剩余宏观误差的根源。
% ==   - 明确 TR 总量在稳态中未经BGP调整，因此其人均猜测路径应直接
% ==     由 ss0.TR_distributed_agg 计算得出，无需任何还原操作。
% ==   - 对遗赠路径的构造也做了简化和澄清，确保所有猜测路径的会计
% ==     口径都与 TPI 引擎的预期完全匹配。
% =========================================================================
clear; close all; clear classes;
addpath(pwd);
fprintf('=== TPI求解器零冲击诊断脚本 (v2.4 - BGP猜测最终修正版) ===\n\n');

%% --- 1. 加载“黄金标准”稳态数据 ---
fprintf('--- 1. 加载黄金标准稳态数据 ---\n');
data_filename = 'SS/data_for_perfect_transition.mat';
if ~exist(data_filename, 'file'), error('黄金标准稳待数据文件 "%s" 不存在。', data_filename); end
load(data_filename, 'data_for_perfect_transition');
fprintf('   ✅ 已加载数据: %s\n', data_filename);

ssF = data_for_perfect_transition.ssF;
distF = data_for_perfect_transition.distF;
valF = data_for_perfect_transition.valF;
polF = data_for_perfect_transition.polF;
cS = data_for_perfect_transition.cS;
paramSF = data_for_perfect_transition.paramSF;

cS.T_sim = 50;
T = cS.T_sim;
fprintf('   [快速诊断模式] T已被临时设为: %d\n', T);

%% --- 2. 手动构建一个完美的、自洽的、平稳增长的人口环境 ---
fprintf('\n--- 2. 手动构建完美的、平稳增长的BGP环境 ---\n');
ss0 = ssF;
dist0 = distF;
fprintf('   ✅ 初始稳态(ss0)和分布(dist0)已设为与终期稳态(ssF)相同。\n');

Z_ss_abs_dist = sum(distF, [1,2,3]);
Z_ss_norm = Z_ss_abs_dist(:) / sum(Z_ss_abs_dist(:));
total_pop_ss0 = sum(dist0, 'all');

g_n_period = (1 + cS.n_ss)^cS.time_Step - 1;
pop_growth_factors = (1 + g_n_period).^(0:T-1);
cS.Z_path_raw = (total_pop_ss0 * Z_ss_norm) * pop_growth_factors;
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

%% --- 3. 构造初始猜测：平坦的价格路径和BGP一致的流量路径 ---
fprintf('\n--- 3. 构造BGP一致的猜测路径 ---\n');
r_path_guess = ones(1, T) * ss0.r_mkt;
w_hat_path_guess = ones(1, T) * ss0.w_hat;
b_hat_path_guess = ones(1, T) * ss0.b_hat;

% [!!! 核心修正: TR路径 !!!]
% ss0.TR_distributed_agg 是未经BGP调整的当期总量。
% 其人均猜测路径应该是平坦的。
tr0_hat_pc = ss0.TR_distributed_agg / total_pop_ss0;
TR_hat_pc_path_guess = ones(1, T) * tr0_hat_pc;

% [!!! 核心修正: Bequest路径 !!!]
% ss0.Bequest_generated_agg 是经过 `/(1+n)` 调整的量。
% 我们先将其还原为未经调整的“期末”原始总量。
Bequest_hat_total_raw_ss0 = ss0.Bequest_generated_agg * (1 + g_n_period);
% 基于这个原始总量，构造BGP路径
Bequest_hat_total_path_guess = Bequest_hat_total_raw_ss0 .* pop_growth_factors;
fprintf('   ✅ 所有价格和流量路径已构造完毕。\n');


%% --- 4. 运行TPI循环 ---
fprintf('\n--- 4. 启动TPI循环进行诊断 ---\n');
max_iter = 10;
tolerance = 1e-7;
damping.r = 0.2; damping.w = 0.2; damping.beq = 0.2; damping.tr = 0.2; damping.b_hat = 0.2;

r_path_iter = r_path_guess;
w_hat_path_iter = w_hat_path_guess;
Bequest_hat_total_path_iter = Bequest_hat_total_path_guess;
TR_hat_pc_path_iter = TR_hat_pc_path_guess;
b_hat_path_iter = b_hat_path_guess;

for iter = 1:max_iter
    [target_paths, errors, Pol_path, ~, ~] = TPI.calculate_paths_and_errors(...
        r_path_iter, w_hat_path_iter, Bequest_hat_total_path_iter, ...
        TR_hat_pc_path_iter, b_hat_path_iter, ...
        ss0, ssF, cS, paramSF, valF, polF, dist0);

    policy_error = 0;
    % 我们只关心第一次迭代的策略误差，之后会因宏观误差而自然偏离
    if iter == 1
        pol_t1 = Pol_path{1};
        for ia = 1:cS.aD_new
            err_c = max(abs(pol_t1(ia).c - polF(ia).c), [], 'all');
            err_k_prime = max(abs(pol_t1(ia).k_prime - polF(ia).k_prime), [], 'all');
            policy_error = max([policy_error, err_c, err_k_prime]);
            if isfield(polF, 'kpps_prime')
                 err_kpps_prime = max(abs(pol_t1(ia).kpps_prime - polF(ia).kpps_prime), [], 'all');
                 policy_error = max(policy_error, err_kpps_prime);
            end
        end
    end

    fprintf('Iter [%2d/%2d] | Error: %.2e (r:%.1e, w:%.1e, bq:%.1e, tr:%.1e, b:%.1e) | Pol_Err(t=1): %.2e\n', ...
        iter, max_iter, errors.total, errors.r, errors.w, errors.beq, errors.tr, errors.b_hat, policy_error);

    if errors.total < tolerance
        fprintf('\n✅ 诊断成功！宏观路径收敛。\n');
        break;
    end
    if iter == max_iter, fprintf('\n❌ 诊断失败！未在最大迭代次数内收敛。\n'); end

    r_path_iter(2:T) = (1 - damping.r) * r_path_iter(2:T) + damping.r * target_paths.r_path(2:T);
    w_hat_path_iter(2:T) = (1 - damping.w) * w_hat_path_iter(2:T) + damping.w * target_paths.w_hat_path(2:T);
    Bequest_hat_total_path_iter(2:T) = (1 - damping.beq) * Bequest_hat_total_path_iter(2:T) + damping.beq * target_paths.Bequest_hat_total_path(2:T);
    TR_hat_pc_path_iter(2:T) = (1 - damping.tr) * TR_hat_pc_path_iter(2:T) + damping.tr * target_paths.TR_hat_pc_path(2:T);
    b_hat_path_iter(2:T) = (1 - damping.b_hat) * b_hat_path_iter(2:T) + damping.b_hat * target_paths.b_hat_path(2:T);
end

%% --- 5. 结果分析 ---
if errors.total < tolerance
    fprintf('\n--- 5. 最终路径检验 ---\n');
    max_dev_r = max(abs(target_paths.r_path - ssF.r_mkt));
    max_dev_w = max(abs(target_paths.w_hat_path - ssF.w_hat));
    fprintf('   利率路径与稳态的最大偏差: %.3e\n', max_dev_r);
    fprintf('   工资路径与稳态的最大偏差: %.3e\n', max_dev_w);
    
    Bequest_hat_total_path_target_bgp = Bequest_hat_total_raw_ss0 .* pop_growth_factors;
    max_dev_beq_total = max(abs(target_paths.Bequest_hat_total_path - Bequest_hat_total_path_target_bgp));
    fprintf('   总量遗赠路径与BGP路径的最大偏差: %.3e\n', max_dev_beq_total);

    if max_dev_r < 1e-5 && max_dev_w < 1e-5 && max_dev_beq_total < 1e-5
        fprintf('   ✅ 结论：求解出的路径确实是平坦的稳态路径。TPI核心逻辑通过检验。\n');
    else
        fprintf('   ⚠️ 警告：路径收敛了，但并非平坦路径。可能存在微小的不一致性。\n');
    end
end

fprintf('\n--- 黄金标准诊断脚本执行完毕 ---\n');