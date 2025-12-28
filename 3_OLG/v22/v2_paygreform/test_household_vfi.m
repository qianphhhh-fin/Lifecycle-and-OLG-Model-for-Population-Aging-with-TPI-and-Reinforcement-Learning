% =========================================================================
% == 脚本: test_vfi_by_age_pps.m
% == 版本: [v2.1 - 结构体数组赋值修正版]
% ==
% == 核心修改:
% ==   - [!!! 关键BUG修复 !!!] 修正了结构体数组赋值错误。
% ==   - 使用 cell 数组 (polS_cell) 在循环中暂存策略结构体。
% ==   - 循环结束后，将 cell 数组一次性转换为标准的结构体数组 polS。
% ==     这是在MATLAB中处理此类问题的标准、稳健做法。
% =========================================================================

clear; clc; close all;
addpath(pwd); % 确保所有类定义都在路径上

fprintf('=== 启动 household 完整生命周期模拟测试脚本 (v2.1) ===\n\n');

%% --- 1. 参数与环境配置 (不变) ---
fprintf('--- 1. 正在配置测试参数和宏观环境...\n');

cS = utils.ParameterValues();
cS.nk = 40;
cS.nkpps = 20;
cS.nw = 1;
cS.pps_active = true;
cS.pps_fixed = 0.10;
cS.n_pps_rate_grid = 11;
cS.pps_withdrawal_rate = 0.05;
cS.pps_tax_rate_withdrawal = 0.03;
cS.g_A_ss = 0.015;
cS.ageEffV_new = cS.ageEffV_new_h(:,1);
cS.s_pathV = cS.s_pathV(:,1);

fprintf('   正在生成资产网格和收入过程...\n');
cS = utils.generateGrids(cS, 'k_max', 80, 'kpps_max', 40);
paramS = struct();
[paramS.leGridV, paramS.TrProbM_by_age, paramS.leProb1V, cS.nw_expanded] = utils.EarningProcess_AgeDependent(cS);

M_age = struct();
M_age.r_mkt_t = 0.04;
M_age.w_t = 1.5;
M_age.tr_per_hh = 0.1;
M_age.theta_t = 0.20;
M_age.g_A_t_plus_1 = (1 + cS.g_A_ss)^cS.time_Step - 1;

fprintf('   ✅ 参数配置完成。\n\n');

%% --- 2. [核心修改] 求解所有年龄的策略函数 (修正赋值逻辑) ---
fprintf('--- 2. 正在执行完整的VFI反向迭代以求解所有年龄的策略...\n');
tic;

valS = -Inf(cS.nk, cS.nkpps, cS.nw_expanded, cS.aD_new);
% [!!! 核心修正 1/3 !!!] 使用cell数组进行预分配
polS_cell = cell(cS.aD_new, 1);

b_age_val_retire = 0.8;

for a_idx = cS.aD_new : -1 : 1
    vPrime_kkppse_next = [];
    if a_idx < cS.aD_new
        vPrime_kkppse_next = valS(:,:,:,a_idx+1);
    end
    
    b_age_val = 0;
    if a_idx > cS.aR_new
        b_age_val = b_age_val_retire;
    end
    
    [val_age, pol_age] = household.VFI_by_age_PPS(a_idx, vPrime_kkppse_next, M_age, b_age_val, paramS, cS);
    
    valS(:,:,:,a_idx) = val_age;
    % [!!! 核心修正 2/3 !!!] 将结构体存入cell单元中
    polS_cell{a_idx} = pol_age;
end

% [!!! 核心修正 3/3 !!!] 循环结束后，将cell数组转换为结构体数组
polS = [polS_cell{:}]';

vfi_time = toc;
fprintf('   ✅ VFI求解完成，耗时: %.2f 秒。\n\n', vfi_time);


%% --- 3. 模拟一个代理人的完整生命周期路径 (不变) ---
fprintf('--- 3. 正在模拟单个代理人的生命周期路径...\n');

k_path = zeros(cS.aD_new, 1);
kpps_path = zeros(cS.aD_new, 1);
c_path = zeros(cS.aD_new, 1);
pps_rate_path = zeros(cS.aD_new, 1);

ie_path = round(cS.nw_expanded / 2); 
fprintf('   (为简化，假设代理人一生都处于劳动效率状态 ie = %d)\n', ie_path);

for a_idx = 1:(cS.aD_new - 1)
    k_now = k_path(a_idx);
    kpps_now = kpps_path(a_idx);
    
    pol_c_slice    = polS(a_idx).c(:,:,ie_path);
    pol_k_p_slice  = polS(a_idx).k_prime(:,:,ie_path);
    pol_kpps_p_slice = polS(a_idx).kpps_prime(:,:,ie_path);
    pol_rate_slice = polS(a_idx).pps_contrib_rate(:,:,ie_path);

    interp_c = griddedInterpolant({cS.kGridV, cS.kppsGridV}, pol_c_slice');
    interp_k_p = griddedInterpolant({cS.kGridV, cS.kppsGridV}, pol_k_p_slice');
    interp_kpps_p = griddedInterpolant({cS.kGridV, cS.kppsGridV}, pol_kpps_p_slice');
    interp_rate = griddedInterpolant({cS.kGridV, cS.kppsGridV}, pol_rate_slice');

    c_path(a_idx) = interp_c(k_now, kpps_now);
    pps_rate_path(a_idx) = interp_rate(k_now, kpps_now);
    
    k_path(a_idx + 1) = interp_k_p(k_now, kpps_now);
    kpps_path(a_idx + 1) = interp_kpps_p(k_now, kpps_now);
end

k_final = k_path(cS.aD_new);
kpps_final = kpps_path(cS.aD_new);
pol_c_final_slice = polS(cS.aD_new).c(:,:,ie_path);
interp_c_final = griddedInterpolant({cS.kGridV, cS.kppsGridV}, pol_c_final_slice');
c_path(cS.aD_new) = interp_c_final(k_final, kpps_final);

fprintf('   ✅ 前向模拟完成。\n\n');

%% --- 4. 可视化完整生命周期路径 (不变) ---
fprintf('--- 4. 正在生成生命周期路径可视化图表...\n');

age_axis = cS.age1_orig + (0:(cS.aD_new-1)) * cS.time_Step;
retirement_age_phys = cS.age1_orig + (cS.aR_new) * cS.time_Step;

figure('Name', '单一代理人完整生命周期路径模拟', 'Position', [100, 100, 1800, 500]);

subplot(1, 3, 1);
plot(age_axis, k_path, 'b-o', 'LineWidth', 2, 'DisplayName', '常规资产 (k)');
hold on;
plot(age_axis, kpps_path, 'r-s', 'LineWidth', 2, 'DisplayName', 'PPS资产 (kpps)');
plot(age_axis, k_path + kpps_path, 'k--', 'LineWidth', 1.5, 'DisplayName', '总资产');
xline(retirement_age_phys, 'm:', 'LineWidth', 2, 'DisplayName', '退休年龄');
title('生命周期资产路径');
xlabel('年龄 (岁)'); ylabel('资产水平');
legend('show', 'Location', 'best');
grid on;

subplot(1, 3, 2);
plot(age_axis, c_path, 'g-d', 'LineWidth', 2, 'DisplayName', '消费 (c)');
hold on;
xline(retirement_age_phys, 'm:', 'LineWidth', 2, 'DisplayName', '退休年龄');
title('生命周期消费路径');
xlabel('年龄 (岁)'); ylabel('消费水平');
legend('show', 'Location', 'best');
grid on;

subplot(1, 3, 3);
bar(age_axis, pps_rate_path * 100, 'FaceColor', [0.8 0.8 0.2]);
hold on;
xline(retirement_age_phys, 'm:', 'LineWidth', 2, 'DisplayName', '退休年龄');
title('生命周期PPS缴费率路径');
xlabel('年龄 (岁)'); ylabel('PPS缴费率 (%)');
ylim([-1, cS.pps_fixed * 100 + 1]);
legend('show', 'Location', 'best');
grid on;

fprintf('   ✅ 可视化完成。\n\n');
fprintf('=== 测试脚本执行完毕 ===\n');