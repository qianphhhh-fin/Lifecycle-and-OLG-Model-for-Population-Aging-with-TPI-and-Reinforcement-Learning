% =========================================================================
% == VFI求解器性能对比脚本
% == 比较 main_steady_state_utils_bgp_fmincon.m 和 main_steady_state_utils_bgp.m
% == 中初始稳态VFI求解器的输出一致性和求解速度
% =========================================================================

clear; clc; close all;
fprintf('========================================================================\n');
fprintf('== VFI求解器性能对比测试 (完整生命周期迭代)\n');
fprintf('== 原版BGP vs FMINCON版本\n');
fprintf('========================================================================\n\n');

%% 第一步：设置相同的参数和输入（复制主脚本的完整初始化流程）

fprintf('📋 步骤1: 初始化模型参数 (按照主脚本流程)...\n');

% === 1.1 基础设定参数（复制主脚本） ===
MODEL_START_YEAR = 2023;         % 数据的起始年份（稳态年份）
TIME_STEP = 5;                   % 模型每期代表的年数
T_SIM_MAX_PERIODS = 40;          % 模拟期数 (40期 * 5年/期 = 200年)

% 设置基础参数
cS = model_setup_utils_bgp.ParameterValues();

% === 1.2 时间相关参数设置（与主脚本完全一致） ===
cS.time_Step = TIME_STEP;
cS.ss0_year = MODEL_START_YEAR;
cS.start_year = MODEL_START_YEAR+1;
cS.end_year = cS.start_year + (T_SIM_MAX_PERIODS - 1) * cS.time_Step;
cS.T_sim = T_SIM_MAX_PERIODS;

% === 1.3 网格设定（复制主脚本） ===
% 为快速测试适当减少规模，但保持主脚本的比例关系
ngrid = 25;  % 原主脚本40，减少到25以加快测试
cS.nk = ngrid; 
cS.nkpps = 15;  % 原主脚本20，减少到15
cS.nkprime = ngrid; 
cS.npps = 15;  % 原主脚本20，减少到15

% 确保非PPS模式
cS.pps_active = false;
cS.g_A_ss = 0.0;  % 初始设置为无技术进步，避免复杂性

% === 1.4 生成网格（必须在EarningProcess之前） ===
cS = model_setup_utils_bgp.generateGrids(cS);

% === 1.5 设置收入过程（完全按照主脚本顺序） ===
paramS = struct();
[paramS.leGridV, paramS.TrProbM_by_age, paramS.leProb1V, cS.nw_expanded] = ...
    model_setup_utils_bgp.EarningProcess_AgeDependent(cS);
paramS.leLogGridV = log(paramS.leGridV(1:cS.nw));

% === 1.6 设置关键的缺失字段（复制主脚本逻辑） ===
% 计算theta路径（需要在VFI调用之前设置）
cS = model_setup_utils_bgp.calcaulte_theta_payg_path(cS, false);
% 设置当前期的theta值（VFI求解器内部需要）
if ~isfield(cS, 'theta_t') || isempty(cS.theta_t)
    cS.theta_t = cS.theta_path(1);
end

% === 1.7 创建相同的市场价格环境（基于主脚本） ===
M_test = struct();
M_test.r_mkt_t = 0.04;                    % 4%市场利率
M_test.w_t = 1.0;                         % 标准化工资
M_test.b_t = 0.3;                         % 养老金
M_test.K_p = 3.5;
M_test.K_g = 1.0;
M_test.L_t = 0.3;
M_test.Y_t = M_test.K_p^cS.alpha * M_test.L_t^(1-cS.alpha);

fprintf('   ✅ 参数设置完成 (按照主脚本流程)\n');
fprintf('      生命周期: %d期, 退休期: %d期\n', cS.aD_new, cS.aR_new);
fprintf('      资本网格: %d点, 收入状态: %d+%d\n', cS.nk, cS.nw, cS.nw_expanded - cS.nw);
fprintf('      技术增长率: %.3f\n', cS.g_A_ss);
fprintf('      模拟时期: %d期 (%d-%d年)\n', cS.T_sim, cS.start_year, cS.end_year);

%% 第二步：测试原版BGP求解器

fprintf('\n📊 步骤2: 测试原版BGP求解器...\n');

% 记录开始时间
tic_bgp = tic;

try
    % 调用原版BGP的VFI求解器
    [cPolM_bgp, kPolM_bgp, cPpsPolM_bgp, valM_bgp] = ...
        main_steady_state_utils_bgp.HHSolution_VFI(M_test, paramS, cS);
    
    time_bgp = toc(tic_bgp);
    success_bgp = true;
    
    fprintf('   ✅ 原版BGP求解成功\n');
    fprintf('      求解时间: %.3f秒\n', time_bgp);
    fprintf('      输出矩阵维度: cPol %s, kPol %s, Val %s\n', ...
        mat2str(size(cPolM_bgp)), mat2str(size(kPolM_bgp)), mat2str(size(valM_bgp)));
    
catch ME
    time_bgp = toc(tic_bgp);
    success_bgp = false;
    fprintf('   ❌ 原版BGP求解失败: %s\n', ME.message);
    
    % 创建空矩阵以便后续比较
    cPolM_bgp = zeros(cS.nk, cS.nw_expanded, cS.aD_new);
    kPolM_bgp = zeros(cS.nk, cS.nw_expanded, cS.aD_new);
    cPpsPolM_bgp = zeros(cS.nk, cS.nw_expanded, cS.aD_new);
    valM_bgp = -Inf(cS.nk, cS.nw_expanded, cS.aD_new);
end

%% 第三步：测试FMINCON版本求解器

fprintf('\n🔧 步骤3: 测试FMINCON版本求解器...\n');

% 记录开始时间  
tic_fmincon = tic;

try
    % 调用FMINCON版本的VFI求解器（非PPS模式）
    [cPolM_fmincon, kPolM_fmincon, TaxPolM_fmincon, ShockPolM_fmincon, valM_fmincon] = ...
        main_steady_state_utils_bgp_fmincon.HHSolution_VFI_fmincon(M_test, paramS, cS);
    
    % 创建兼容性变量（非PPS模式下cPpsPolM为零）
    cPpsPolM_fmincon = zeros(size(cPolM_fmincon));
    
    time_fmincon = toc(tic_fmincon);
    success_fmincon = true;
    
    fprintf('   ✅ FMINCON版本求解成功\n');
    fprintf('      求解时间: %.3f秒\n', time_fmincon);
    fprintf('      输出矩阵维度: cPol %s, kPol %s, Val %s\n', ...
        mat2str(size(cPolM_fmincon)), mat2str(size(kPolM_fmincon)), mat2str(size(valM_fmincon)));
    fprintf('      额外输出: TaxPol %s, ShockPol %s\n', ...
        mat2str(size(TaxPolM_fmincon)), mat2str(size(ShockPolM_fmincon)));
    
catch ME
    time_fmincon = toc(tic_fmincon);
    success_fmincon = false;
    fprintf('   ❌ FMINCON版本求解失败: %s\n', ME.message);
    
    % 创建空矩阵以便后续比较
    cPolM_fmincon = zeros(cS.nk, cS.nw_expanded, cS.aD_new);
    kPolM_fmincon = zeros(cS.nk, cS.nw_expanded, cS.aD_new);
    cPpsPolM_fmincon = zeros(cS.nk, cS.nw_expanded, cS.aD_new);
    valM_fmincon = -Inf(cS.nk, cS.nw_expanded, cS.aD_new);
    TaxPolM_fmincon = zeros(cS.nk, cS.nw_expanded, cS.aD_new);
    ShockPolM_fmincon = zeros(cS.nk, cS.nw_expanded, cS.aD_new);
end

%% 第四步：性能比较

fprintf('\n⚡ 步骤4: 性能比较结果\n');
fprintf('========================================================================\n');

if success_bgp && success_fmincon
    speedup_ratio = time_bgp / time_fmincon;
    if speedup_ratio > 1
        fprintf('🚀 FMINCON版本更快: %.2fx 加速\n', speedup_ratio);
    else
        fprintf('🐌 原版BGP更快: %.2fx 加速\n', 1/speedup_ratio);
    end
elseif success_bgp
    fprintf('⚠️  只有原版BGP成功求解\n');
elseif success_fmincon
    fprintf('⚠️  只有FMINCON版本成功求解\n');
else
    fprintf('❌ 两个版本都失败了\n');
end

fprintf('\n时间详细对比:\n');
fprintf('   原版BGP:      %.3f 秒\n', time_bgp);
fprintf('   FMINCON版本:  %.3f 秒\n', time_fmincon);
fprintf('   时间差:       %.3f 秒\n', abs(time_bgp - time_fmincon));

%% 第五步：输出一致性比较

if success_bgp && success_fmincon
    fprintf('\n🔍 步骤5: 输出一致性分析\n');
    fprintf('========================================================================\n');
    
    % 比较消费策略矩阵
    cPol_diff = abs(cPolM_bgp - cPolM_fmincon);
    cPol_max_diff = max(cPol_diff(:));
    cPol_mean_diff = mean(cPol_diff(:));
    cPol_rel_diff = max(cPol_diff(:)) / max(abs(cPolM_bgp(:)));
    
    fprintf('📊 消费策略 (cPolM) 比较:\n');
    fprintf('   最大绝对差异:  %.6e\n', cPol_max_diff);
    fprintf('   平均绝对差异:  %.6e\n', cPol_mean_diff);
    fprintf('   最大相对差异:  %.6e (%.4f%%)\n', cPol_rel_diff, cPol_rel_diff*100);
    
    % 比较储蓄策略矩阵
    kPol_diff = abs(kPolM_bgp - kPolM_fmincon);
    kPol_max_diff = max(kPol_diff(:));
    kPol_mean_diff = mean(kPol_diff(:));
    kPol_rel_diff = max(kPol_diff(:)) / max(abs(kPolM_bgp(:)));
    
    fprintf('\n🏦 储蓄策略 (kPolM) 比较:\n');
    fprintf('   最大绝对差异:  %.6e\n', kPol_max_diff);
    fprintf('   平均绝对差异:  %.6e\n', kPol_mean_diff);
    fprintf('   最大相对差异:  %.6e (%.4f%%)\n', kPol_rel_diff, kPol_rel_diff*100);
    
    % 比较价值函数矩阵（排除-Inf值）
    valid_mask = isfinite(valM_bgp) & isfinite(valM_fmincon);
    if any(valid_mask(:))
        val_diff = abs(valM_bgp - valM_fmincon);
        val_max_diff = max(val_diff(valid_mask));
        val_mean_diff = mean(val_diff(valid_mask));
        val_rel_diff = max(val_diff(valid_mask)) / max(abs(valM_bgp(valid_mask)));
        
        fprintf('\n💎 价值函数 (valM) 比较:\n');
        fprintf('   最大绝对差异:  %.6e\n', val_max_diff);
        fprintf('   平均绝对差异:  %.6e\n', val_mean_diff);
        fprintf('   最大相对差异:  %.6e (%.4f%%)\n', val_rel_diff, val_rel_diff*100);
        fprintf('   有效比较点数:  %d / %d\n', sum(valid_mask(:)), numel(valM_bgp));
    else
        fprintf('\n💎 价值函数 (valM) 比较:\n');
        fprintf('   ⚠️  没有有效的比较点\n');
    end
    
    % 一致性判断
    tolerance_loose = 1e-3;   % 宽松容差
    tolerance_strict = 1e-6;  % 严格容差
    
    fprintf('\n🎯 一致性判断:\n');
    
    is_consistent_loose = (cPol_max_diff < tolerance_loose) && ...
                         (kPol_max_diff < tolerance_loose) && ...
                         (exist('val_max_diff', 'var') && val_max_diff < tolerance_loose);
    
    is_consistent_strict = (cPol_max_diff < tolerance_strict) && ...
                          (kPol_max_diff < tolerance_strict) && ...
                          (exist('val_max_diff', 'var') && val_max_diff < tolerance_strict);
    
    if is_consistent_strict
        fprintf('   ✅ 高度一致 (误差 < %.0e)\n', tolerance_strict);
    elseif is_consistent_loose
        fprintf('   ✅ 基本一致 (误差 < %.0e)\n', tolerance_loose);
    else
        fprintf('   ⚠️  存在显著差异 (误差 > %.0e)\n', tolerance_loose);
    end
    
    % 生成详细的差异分布图
    if max(cPol_max_diff, kPol_max_diff) > 1e-8
        fprintf('\n📈 生成差异分布图...\n');
        
        figure('Name', 'VFI求解器比较', 'Position', [100, 100, 1200, 800]);
        
        subplot(2,3,1);
        imagesc(squeeze(cPolM_bgp(:,:,end-5))); colorbar;
        title('原版BGP - 消费策略 (某年龄组)');
        
        subplot(2,3,2);
        imagesc(squeeze(cPolM_fmincon(:,:,end-5))); colorbar;
        title('FMINCON版本 - 消费策略 (某年龄组)');
        
        subplot(2,3,3);
        imagesc(squeeze(cPol_diff(:,:,end-5))); colorbar;
        title('消费策略差异');
        
        subplot(2,3,4);
        imagesc(squeeze(kPolM_bgp(:,:,end-5))); colorbar;
        title('原版BGP - 储蓄策略 (某年龄组)');
        
        subplot(2,3,5);
        imagesc(squeeze(kPolM_fmincon(:,:,end-5))); colorbar;
        title('FMINCON版本 - 储蓄策略 (某年龄组)');
        
        subplot(2,3,6);
        imagesc(squeeze(kPol_diff(:,:,end-5))); colorbar;
        title('储蓄策略差异');
        
        sgtitle('VFI求解器输出比较', 'FontSize', 16, 'FontWeight', 'bold');
    end
    
else
    fprintf('\n❌ 步骤5: 无法进行一致性比较（至少一个求解器失败）\n');
end

%% 第六步：FMINCON版本特有功能测试

if success_fmincon
    fprintf('\n🔧 步骤6: FMINCON版本特有功能分析\n');
    fprintf('========================================================================\n');
    
    % 分析额外的税收和冲击支出矩阵
    tax_total = sum(TaxPolM_fmincon(:));
    shock_total = sum(ShockPolM_fmincon(:));
    tax_max = max(TaxPolM_fmincon(:));
    shock_max = max(ShockPolM_fmincon(:));
    
    fprintf('💰 税收策略矩阵 (TaxPolM) 分析:\n');
    fprintf('   总税收:        %.6f\n', tax_total);
    fprintf('   最大税收:      %.6f\n', tax_max);
    fprintf('   非零元素数:    %d / %d\n', sum(TaxPolM_fmincon(:) > 0), numel(TaxPolM_fmincon));
    
    fprintf('\n💥 冲击支出矩阵 (ShockPolM) 分析:\n');
    fprintf('   总冲击支出:    %.6f\n', shock_total);
    fprintf('   最大冲击支出:  %.6f\n', shock_max);
    fprintf('   非零元素数:    %d / %d\n', sum(ShockPolM_fmincon(:) > 0), numel(ShockPolM_fmincon));
    
    fprintf('\n🎯 FMINCON版本优势:\n');
    fprintf('   ✅ 事前计算：VFI阶段直接计算存储所有会计流量\n');
    fprintf('   ✅ 直接聚合：聚合阶段无需反解计算\n');
    fprintf('   ✅ 连续优化：使用fmincon替代离散网格搜索\n');
    if tax_total > 0
        fprintf('   ✅ 完整税收信息：平均税负 %.2f%%\n', tax_max/max(cPolM_fmincon(:))*100);
    end
    if shock_total > 0
        fprintf('   ✅ 完整冲击信息：冲击支出占比 %.2f%%\n', shock_max/max(cPolM_fmincon(:))*100);
    end
end

%% 第七步：总结报告

fprintf('\n🎉 最终总结报告\n');
fprintf('========================================================================\n');

fprintf('📊 计算量统计:\n');
fprintf('   生命周期期数:    %d\n', cS.aD_new);
fprintf('   状态空间大小:    %d (k) × %d (e) = %d\n', cS.nk, cS.nw_expanded, cS.nk * cS.nw_expanded);
fprintf('   总决策点数:      %d\n', cS.nk * cS.nw_expanded * cS.aD_new);

fprintf('\n⏱️  性能总结:\n');
fprintf('   原版BGP时间:     %.3f 秒\n', time_bgp);
fprintf('   FMINCON时间:     %.3f 秒\n', time_fmincon);
if success_bgp && success_fmincon
    if time_fmincon < time_bgp
        fprintf('   🚀 性能改进:     %.2fx 加速\n', time_bgp/time_fmincon);
    else
        fprintf('   🐌 性能损失:     %.2fx 减速\n', time_fmincon/time_bgp);
    end
end

if success_bgp && success_fmincon && exist('is_consistent_loose', 'var')
    fprintf('\n🎯 一致性总结:\n');
    if is_consistent_strict
        fprintf('   ✅ 结果高度一致 (数值误差 < %.0e)\n', tolerance_strict);
    elseif is_consistent_loose
        fprintf('   ✅ 结果基本一致 (数值误差 < %.0e)\n', tolerance_loose);
    else
        fprintf('   ⚠️  存在系统性差异，需进一步调查\n');
    end
end

fprintf('\n💡 建议:\n');
if success_bgp && success_fmincon && time_fmincon < time_bgp
    fprintf('   🚀 FMINCON版本在性能上有优势，推荐使用\n');
elseif success_bgp && success_fmincon && time_fmincon > time_bgp * 2
    fprintf('   🔧 FMINCON版本较慢，可能需要优化参数设置\n');
elseif success_bgp && success_fmincon
    fprintf('   ⚖️  两个版本性能相当，可根据需求选择\n');
end

if success_fmincon
    fprintf('   📈 FMINCON版本提供额外的会计信息，有助于调试和分析\n');
end

fprintf('\n========================================================================\n');
fprintf('测试完成！\n');
fprintf('========================================================================\n'); 