% =========================================================================
% == SCRIPT: Interpolator Creation Debugger
% == 目的：精确检验从VFI求解器输出的策略矩阵kPolM，到创建的
% ==       插值器kPolInterp，信息是否被正确传递，特别是在k=0点。
% =========================================================================
clear; close all;
addpath(pwd);
fprintf('=== 插值器创建专项调试脚本 ===\n\n');

%% 1. 初始化一个与主程序完全一致的环境
fprintf('--- 1. 初始化一个与主程序完全一致的测试环境 ---\n');
cS = main_olg_v14_utils.ParameterValues_HuggettStyle();
paramS = struct();
% --- 使用之前确定最稳健的参数 ---
cS.beta = 0.985; cS.sigma = 3.0; cS.phi_bequest = 2.0; cS.sigma_bequest = cS.sigma;
cS.alpha = 0.4; cS.cFloor = 0.05; cS.nSim = 5;
cS.time_Step = 5; cS.ddk = 1 - (1 - 0.05)^cS.time_Step;
cS.tau_k = 0.25; cS.tau_l = 0.20; cS.tau_c = 0.10;
ngrid = 20; cS.nk = ngrid; cS.nkpps = 1; cS.nkprime = ngrid; cS.npps = 1;
cS.A = 1.0; cS = main_olg_v14_utils.generateGrids(cS);
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v14_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));

% --- 创建一个固定的、合理的宏观价格 M_vfi ---
M_vfi = struct();
M_vfi.w_t = 1.5;
M_vfi.r_mkt_period = 0.04 * cS.time_Step;
M_vfi.r_net_period = M_vfi.r_mkt_period * (1 - cS.tau_k);
M_vfi.b_t = 0.5;
cS.theta_t = 0.10; cS.pps_active = false;
fprintf('✅ 测试环境初始化完成。\n\n');


%% 2. 运行VFI求解器，获取“原始”策略矩阵 kPolM
fprintf('--- 2. 运行 HHSolution_VFI_Huggett 获取策略矩阵 kPolM ---\n');
tr_per_capita_vfi = 0; % 假设无遗赠

% [核心调用]
[~, kPolM, ~, ~] = main_olg_v14_utils.HHSolution_VFI_Huggett(M_vfi, tr_per_capita_vfi, paramS, cS);
fprintf('   VFI求解完毕。\n\n');


%% 3. [关键检验] 检查原始kPolM在k=0处的值
fprintf('--- 3. 检查【原始策略矩阵 kPolM】在 k=0 处的值 ---\n');
a_check = 5; % 我们关心第5个年龄组
ie_check = round(cS.nw / 2); % 一个中间的效率状态

% k=0 对应的是资产网格的第一个点，即索引 1
% k_pps=0 对应的是pps网格的第一个点，即索引 1
k_prime_at_zero_asset_raw = kPolM(1, 1, ie_check, a_check);

fprintf('   在 a=%d, e=%d, k=0 时，kPolM(1,1,ie,a) 的原始值是: %.6f\n', ...
    a_check, ie_check, k_prime_at_zero_asset_raw);

if k_prime_at_zero_asset_raw > 0
    fprintf('   ✅ 证据1: VFI求解器确实为 k=0 计算出了一个【正】的储蓄决策。\n');
else
    fprintf('   ❌ 证据1: VFI求解器为 k=0 计算出的储蓄决策是【零或负数】！问题源于VFI！\n');
end


%% 4. [关键检验] 手动创建插值器并立即测试
fprintf('\n--- 4. 手动创建【单个插值器】并立即测试其在 k=0 处的值 ---\n');

% 提取与 a_check, ie_check 对应的策略向量
kPol_vector = squeeze(kPolM(:, 1, ie_check, a_check));

fprintf('   用于创建插值器的数据向量 (前5个值): \n');
disp(kPol_vector(1:5)');

% 手动创建插值器
kPolInterp_test = griddedInterpolant(cS.kGridV, kPol_vector, 'linear', 'none');

% 立即测试 k=0 点
k_prime_from_interpolator = kPolInterp_test(0);

fprintf('\n   新创建的插值器在 k=0 点返回的值是: %.6f\n', k_prime_from_interpolator);

if abs(k_prime_from_interpolator - k_prime_at_zero_asset_raw) < 1e-9
    fprintf('   ✅ 证据2: 插值器正确地复现了 kPolM 在 k=0 处的值。\n');
else
    fprintf('   ❌ 证据2: 插值器【未能】正确复现 kPolM 在 k=0 处的值！\n');
    fprintf('      这意味着 griddedInterpolant 的行为存在问题或误用。\n');
end

%% 5. 最终结论
fprintf('\n--- 5. 最终诊断 ---\n');
if k_prime_at_zero_asset_raw > 0 && abs(k_prime_from_interpolator - k_prime_at_zero_asset_raw) < 1e-9
    fprintf('✅✅✅ 最终结论：VFI求解和插值器创建【均无问题】！\n');
    fprintf('   kPolM(1,...) 本身是正数，插值器也能正确返回这个正数。\n');
    fprintf('   如果 HHSimulation_olgm 仍然得到 k''(0)=0，\n');
    fprintf('   那唯一的解释是 HHSimulation_olgm 在调用插值器时传入的【上下文参数】(ie, a_idx)有误，\n');
    fprintf('   或者它调用了另一个被错误初始化的插值器实例。\n');
    fprintf('   请仔细检查 HHSimulation_olgm 中 kPolInterp{ie, a_idx} 的索引是否正确。\n');
else
    fprintf('❌❌❌ 最终结论：我们在证据1或证据2中发现了根本性错误。\n');
    fprintf('   请根据上面的诊断信息，定位到具体的出错环节。\n');
end
fprintf('--- 调试结束 ---\n');