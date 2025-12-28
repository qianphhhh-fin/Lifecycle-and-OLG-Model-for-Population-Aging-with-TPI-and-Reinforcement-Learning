% --- TPI.m ---
% =========================================================================
% == CLASSDEF: TPI.m
% == 版本: [v1.3 - 聚合目标与人口对齐修正版]
% ==
% == 目的:
% ==   - 修正 calculate_paths_and_errors 中目标遗赠路径的计算逻辑，确保
% ==     总量到人均的正确转换和时期匹配。
% =========================================================================
classdef TPI
    methods (Static)

% 脚本: TPI.m

function [target_paths, errors, Pol_path, Dist_path, aggr_supply] = calculate_paths_and_errors(r_path, w_hat_path, bequest_gen_raw_pc_path_guess, TR_hat_pc_path, b_hat_path, ss0, ssF, cS, paramSF, valF, polF, dist0)
    % =========================================================================
    % ==          TPI核心计算引擎 (v7.3 - 变量名对齐最终版)
    % ==
    % == 核心修正:
    % ==  - 将第三个输入参数和相关的迭代变量名从 `bequest_gen_hat_pc_path_guess`
    % ==    修改为 `bequest_gen_raw_pc_path_guess`，以准确反映其
    % ==    “未经BGP调整的人均原始遗赠产出”的经济含义。
    % ==  - 在函数末尾，确保访问 `aggr_supply` 和计算误差时，都使用
    % ==    统一的字段名 `Bequest_gen_raw_pc_path`。
    % =========================================================================
    T = cS.T_sim;

    % --- 1. 准备向后求解所需的路径 ---
    pathS = struct();
    pathS.r_path = r_path;
    pathS.w_hat_path = w_hat_path;
    pathS.tr_per_hh_hat_path = TR_hat_pc_path;
    pathS.b_hat_path = b_hat_path;
    pathS = TPI.fill_auxiliary_paths(pathS, cS);

    % --- 2. 向后求解家庭最优策略 ---
    [Pol_path, ~] = household.backward_hh(pathS, cS, paramSF, valF, polF);

    % --- 3. [!!! 核心重构 !!!] 向前模拟：准备新生儿资产并模拟分布 ---
    
    % 3a. 计算【每个新生儿在期初收到的资产】路径 (bequest_per_newborn_hat_path)
    mass_newborns_path = cS.Z_path_raw(1, :);
    total_pop_path = sum(cS.Z_path_raw, 1);
    
    % 获取 t=0 (即ss0) 的总人口
    total_pop_t0 = sum(dist0, 'all');

    % 构造一个 t=0 到 t=T-1 的总人口路径
    total_pop_path_t_minus_1 = [total_pop_t0, total_pop_path(1:T-1)];

    % 构造一个 t=0 到 t=T-1 的【未经BGP调整的人均原始遗赠产出】猜测路径
    % t=0 的值来自 ss0
    bequest_gen_raw_pc_path_guess_t_minus_1 = [ss0.Bequest_gen_hat_raw_ss / total_pop_t0, bequest_gen_raw_pc_path_guess(1:T-1)];

    % 计算 t=0 到 t=T-1 每一期产生的【未经BGP调整的原始遗赠总量】
    bequest_gen_raw_total_path_guess_t_minus_1 = bequest_gen_raw_pc_path_guess_t_minus_1 .* total_pop_path_t_minus_1;
    
    % 将 t-1 期产生的总量，分配给 t 期的新生儿
    bequest_per_newborn_hat_path = zeros(1, T);
    non_zero_newborns_mask = mass_newborns_path > 1e-9;
    bequest_per_newborn_hat_path(non_zero_newborns_mask) = ...
        bequest_gen_raw_total_path_guess_t_minus_1(non_zero_newborns_mask) ./ mass_newborns_path(non_zero_newborns_mask);

    % 3b. 使用正确的初始分布和新生儿资产路径，进行前向模拟
    Dist_path = distribution.simulate_dist_forward(dist0, Pol_path, cS, paramSF, bequest_per_newborn_hat_path);

    % --- 4. 聚合宏观供给量 ---
    aggr_supply = aggregates.get_path_aggregates(Dist_path, Pol_path, cS, ss0, paramSF, dist0);
    L_path_supply = aggr_supply.L_path;
    K_p_hat_total_path_supply = aggr_supply.K_p_hat_total_path;

    % --- 5. 计算引致的新宏观价格和流量路径 (Target Paths) ---
    target_paths = struct();
    target_paths.r_path = zeros(1, T);
    target_paths.w_hat_path = zeros(1, T);
    Y_hat_total_path_calc = zeros(1, T);
    K_g_hat_total_path = (ssF.K_public_hat / sum(cS.Z_path(:,end))) * total_pop_path; % 使用比率来定标

    for t = 1:T
        prices_t = firm.get_prices_at_t(K_p_hat_total_path_supply(t), K_g_hat_total_path(t), L_path_supply(t), cS);
        target_paths.r_path(t) = prices_t.r_mkt_t;
        target_paths.w_hat_path(t) = prices_t.w_hat_t;
        Y_hat_total_path_calc(t) = prices_t.Y_hat_t;
    end
    
    mass_retirees_path = sum(cS.Z_path_raw((cS.aR_new+1):end, :), 1);
    total_pension_pot_hat_path = cS.theta_path .* target_paths.w_hat_path .* L_path_supply;
    target_paths.b_hat_path = zeros(1,T);
    non_zero_retirees_mask = mass_retirees_path > 1e-9;
    target_paths.b_hat_path(non_zero_retirees_mask) = total_pension_pot_hat_path(non_zero_retirees_mask) ./ mass_retirees_path(non_zero_retirees_mask);
    
    TR_total_hat_target_path = TPI.calculate_endogenous_TR_path(target_paths.r_path, target_paths.w_hat_path, K_p_hat_total_path_supply, L_path_supply, Y_hat_total_path_calc, Dist_path, Pol_path, cS, paramSF);
    target_paths.TR_hat_pc_path = TR_total_hat_target_path ./ total_pop_path;
    
    % [!!! 核心修正: 统一变量名 !!!]
    target_paths.bequest_gen_raw_pc_path = aggr_supply.Bequest_gen_raw_pc_path;

    % --- 6. 计算误差 (直接同期比较) ---
    errors = struct();
    update_range = 1:(T-1); % 误差计算应该覆盖到倒数第二期
    errors.r = max(abs(target_paths.r_path(update_range) - r_path(update_range)));
    errors.w = max(abs(target_paths.w_hat_path(update_range) - w_hat_path(update_range)));
    errors.b_hat = max(abs(target_paths.b_hat_path(update_range) - b_hat_path(update_range)));
    errors.tr = max(abs(target_paths.TR_hat_pc_path(update_range) - TR_hat_pc_path(update_range)));
    % [!!! 核心修正: 统一变量名 !!!]
    errors.beq = max(abs(target_paths.bequest_gen_raw_pc_path(update_range) - bequest_gen_raw_pc_path_guess(update_range)));
    errors.total = max([errors.r, errors.w, errors.b_hat, errors.beq, errors.tr]);
end



function pathS = fill_auxiliary_paths(pathS, cS)
            % =========================================================================
            % == 函数: fill_auxiliary_paths
            % == 版本: [v2.0 - 职责简化修正版]
            % ==
            % == 核心修改:
            % ==   - [!!!] 删除了所有内生变量（w, L, b_hat）的估算逻辑。
            % ==     这些变量是TPI的迭代目标，它们的路径由外部的猜测直接给出，
            % ==     不应在此处基于错误的简化假设（如忽略异质性）被重新计算。
            % ==   - 此函数的职责被严格限定为填充真正外生的路径（A, theta）。
            % =========================================================================
            T = cS.T_sim;
            pathS.A_path = cS.A_path;
            pathS.Z_path_raw = cS.Z_path_raw;
            pathS.theta_path = cS.theta_path;

            % 移除了以下错误和冗余的计算:
            % w_path = pathS.w_hat_path .* cS.A_path;
            % L_path_matrix = cS.Z_path_raw(1:cS.aR_new, :) .* cS.ageEffV_new(1:cS.aR_new); % <- 忽略了e维度的异质性
            % L_path = sum(L_path_matrix, 1);
            % pathS.b_hat_path = zeros(1, T);
            % for t = 1:T
            %     mass_retirees_t = sum(cS.Z_path_raw((cS.aR_new+1):end, t));
            %     total_pension_pot = pathS.theta_path(t) * w_path(t) * L_path(t);
            %     b_raw_t = total_pension_pot / max(1e-9, mass_retirees_t);
            %     pathS.b_hat_path(t) = b_raw_t / pathS.A_path(t); % <- b_hat_path是迭代变量，不应在此计算
            % end
        end
        function TR_total_hat_path = calculate_endogenous_TR_path(r_path, w_hat_path, K_p_hat_path, L_path, Y_hat_path, Dist_path, Pol_path, cS, paramS)
            % =========================================================================
            % == 函数: calculate_endogenous_TR_path
            % == 版本: [v1.2 - 总量输出修正版]
            % ==
            % == 核心修正:
            % ==   - 函数的最终输出从人均量 TR_hat_path 修改为【总量】TR_total_hat_path。
            % ==   - 这使得其输出与TPI循环的迭代变量(TR_hat_path)的定义保持一致。
            % =========================================================================
            T = cS.T_sim;
            TR_total_hat_path = zeros(1, T);

            for t=1:T
                % --- 步骤 1: 聚合当期的总税收 (包含PAYG缴费) ---
                total_tax_hat_t_raw = 0;
                dist_t_abs = Dist_path(:,:,:,:,t);
                pol_t = Pol_path{t};
                for ia = 1:cS.aD_new
                    mass_ia_abs = dist_t_abs(:,:,:,ia);
                    total_tax_hat_t_raw = total_tax_hat_t_raw + sum(pol_t(ia).tax_regular .* mass_ia_abs, 'all');
                end

                % --- 步骤 2: 聚合当期的PAYG缴费总额 ---
                payg_contrib_hat_t_raw = 0;
                for ia = 1:cS.aR_new % 只对工作年龄人口计算
                    mass_ia_abs = dist_t_abs(:,:,:,ia);
                    labor_income_slice_hat = w_hat_path(t) * cS.ageEffV_new(ia) .* paramS.leGridV(1:cS.nw_expanded)';
                    payg_contrib_slice_hat = cS.theta_path(t) * labor_income_slice_hat;
                    payg_contrib_grid_hat = repmat(reshape(payg_contrib_slice_hat, [1,1,cS.nw_expanded]), [cS.nk, cS.nkpps, 1]);
                    payg_contrib_hat_t_raw = payg_contrib_hat_t_raw + sum(payg_contrib_grid_hat .* mass_ia_abs, 'all');
                end

                % --- 步骤 3: 计算一般税收收入 (从总税收中剔除PAYG缴费) ---
                general_tax_hat_t_raw = total_tax_hat_t_raw - payg_contrib_hat_t_raw;

                % --- 步骤 4: 计算政府预算的其他项 (使用水平量，level) ---
                General_Tax_Revenue_t_level = general_tax_hat_t_raw * cS.A_path(t);
                Y_t_level = Y_hat_path(t) * cS.A_path(t);
                w_t_level = w_hat_path(t) * cS.A_path(t);
                L_t_level = L_path(t); % 劳动供给已经是水平量
                K_p_t_level = K_p_hat_path(t) * cS.A_path(t);

                K_p_return_gross_t_level = r_path(t) * K_p_t_level + cS.ddk * K_p_t_level;
                Public_Capital_Return_t_level = Y_t_level - w_t_level * L_t_level - K_p_return_gross_t_level;

                G_c_t_level = cS.G_c_to_Y_ratio_ss * Y_t_level;
                I_g_t_level = cS.I_g_to_Y_ratio_ss * Y_t_level;

                Total_General_Gov_Revenue_t_level = General_Tax_Revenue_t_level + Public_Capital_Return_t_level;

                % --- 步骤 5: 求解平衡预算的【总量】TR (level) ---
                TR_total_t_level = Total_General_Gov_Revenue_t_level - G_c_t_level - I_g_t_level;

                % --- 步骤 6: [核心修正] 将TR总量转换为去趋势hat总量 ---
                TR_total_hat_path(t) = TR_total_t_level / cS.A_path(t);
            end
        end



        function check_transition_NIPA(results, cS, ss0, ssF)
    % =========================================================================
    % == 函数: check_transition_NIPA (版本 v3.1 - 深度调试版)
    % == 核心修改:
    % ==   - 增加了详细的 t=1 时刻变量诊断信息，用于追踪投资缺口的根源。
    % ==   - 打印 K(1), K(2), hat(K(1)), hat(K(2)) 等关键存量。
    % ==   - 报告 hat(Y(1)), hat(C(1)) 等关键流量。
    % =========================================================================
    
    fprintf('   正在启动过渡路径国民账户一致性检验 (v3.1 - 深度调试版)...\n');
    T = cS.T_sim;
    
    % --- 1. 获取所有流量和存量路径 (Level) ---
    Y_path = results.Y_path;
    C_path = results.C_path;
    I_g_path = results.I_g_path;
    G_path = cS.G_c_to_Y_ratio_ss * Y_path;
    K_p_path = results.K_p_path;
    
    % --- 2. 获取两种口径的投资 (Level) ---
    I_p_from_NIPA = results.I_p_path; 
    I_p_from_K_law = results.I_p_from_K_law;

    % --- 3. 计算NIPA闭合情况和投资缺口 ---
    Total_Uses_path = C_path + I_p_from_NIPA + I_g_path + G_path;
    nipa_error_path = Y_path - Total_Uses_path;
    investment_gap_path = I_p_from_K_law - I_p_from_NIPA;
    
    safe_Y_path = Y_path;
    safe_Y_path(abs(safe_Y_path) < 1e-9) = 1e-9;
    rel_investment_gap_path = investment_gap_path ./ safe_Y_path;
    
    [max_rel_gap, t_max_gap] = max(abs(rel_investment_gap_path(1:T-1)));
    
    fprintf('   国民账户资源约束误差 (Y - [C+I_nipa+Ig+Gc]): %.4e (应接近于0)\n', max(abs(nipa_error_path(1:T-1))));
    fprintf('   核心诊断: 最大投资缺口 (I_k_law - I_nipa)/Y: %.4e (%.4f%%) at t=%d\n', max_rel_gap, max_rel_gap*100, t_max_gap);
    
    % --- 4. [!!! 深度调试信息 !!!] ---
    fprintf('\n   --- [深度调试] t=1 时刻详细诊断 ---\n');
    % 获取去趋势量 (hat量)
    total_pop_path = sum(cS.Z_path_raw, 1);
    K_p_hat_path = K_p_path ./ cS.A_path;
    Y_hat_path = Y_path ./ cS.A_path;
    C_hat_path = C_path ./ cS.A_path;

    % 打印存量信息
    fprintf('      [存量 Stock]\n');
    fprintf('         初始稳态 ss0.K_private_hat .......... : %15.6f (来自稳态求解器)\n', ss0.K_private_hat);
    fprintf('         t=1 期初资本 K_p_hat(1) ............... : %15.6f (来自TPI聚合dist0)\n', K_p_hat_path(1));
    fprintf('         t=2 期初资本 K_p_hat(2) ............... : %15.6f (来自TPI聚合dist(2))\n', K_p_hat_path(2));
    fprintf('         K_p(1)_level .......................... : %15.6f\n', K_p_path(1));
    fprintf('         K_p(2)_level .......................... : %15.6f\n', K_p_path(2));
    
    % 打印流量信息
    fprintf('      [流量 Flow]\n');
    fprintf('         t=1 产出 Y_hat(1) ...................... : %15.6f\n', Y_hat_path(1));
    fprintf('         t=1 消费 C_hat(1) ...................... : %15.6f\n', C_hat_path(1));
    fprintf('         I_p_from_K_law(1) = K(2)-(1-d)K(1) ..... : %15.6f\n', I_p_from_K_law(1));
    fprintf('         I_p_from_NIPA(1) = Y(1)-C(1)-... ...... : %15.6f\n', I_p_from_NIPA(1));
    
    % 打印关键比率
    fprintf('      [关键比率]\n');
    fprintf('         ss0 K/Y ............................... : %15.6f\n', ss0.K_private_hat / ss0.Y_from_production_hat);
    fprintf('         t=1 K/Y ............................... : %15.6f\n', K_p_hat_path(1) / Y_hat_path(1));

    % 打印诊断结论
    k_diff = K_p_hat_path(1) - ss0.K_private_hat;
    fprintf('      [诊断]\n');
    fprintf('         K_p_hat(1) 与 ss0.K_private_hat 的差异 : %.4e\n', k_diff);
    if abs(k_diff / ss0.K_private_hat) > 1e-5
        fprintf('         >>> 警告: t=1的期初资本总量与初始稳态严重不符！这很可能是问题的根源。\n');
        fprintf('             请检查 test_diffN_ss.m 中 dist0 的构造 和 aggregates.m 中两种聚合方式是否完全一致。\n');
    else
        fprintf('         >>> 信息: t=1的期初资本总量与初始稳态基本一致。问题可能出在 K_p_hat(2) 的计算上。\n');
    end


    % --- 5. 如果投资缺口过大，则打印标准的诊断报告 ---
    if max_rel_gap > 1e-4
        warning('过渡路径存在显著的投资缺口，表明储蓄决策与宏观资源约束不完全一致。');
        fprintf('\n   --- 标准诊断报告 ---\n');
        
        print_gap_at_t = @(t) fprintf(...
            ['\n   [时期 t = %d]\n' ...
             '      (+) 产出 (Y) .................... : %15.6f\n' ...
             '      (-) 私人消费 (C) ................ : %15.6f\n' ...
             '      (-) 公共投资 (I_g) ................ : %15.6f\n' ...
             '      (-) 政府消费 (G_c) ................ : %15.6f\n' ...
             '      ------------------------------------------------------------\n' ...
             '      (=) NIPA隐含投资 (I_nipa) ....... : %15.6f  (%.2f%% of Y)\n' ...
             '      (>) 资本积累要求投资 (I_k_law) .. : %15.6f  (%.2f%% of Y)\n' ...
             '      >>> 投资缺口 (I_k_law - I_nipa) . : %15.6f  (%.4f%% of Y)\n'], ...
            t, Y_path(t), C_path(t), I_g_path(t), G_path(t), ...
            I_p_from_NIPA(t), (I_p_from_NIPA(t)/Y_path(t))*100, ...
            I_p_from_K_law(t), (I_p_from_K_law(t)/Y_path(t))*100, ...
            investment_gap_path(t), rel_investment_gap_path(t)*100);
        
        print_gap_at_t(1);
        if t_max_gap ~= 1 && T > 1, print_gap_at_t(t_max_gap); end
        if T > 2 && t_max_gap ~= T-1 && 1 ~= T-1, print_gap_at_t(T-1); end
        fprintf('\n');
    else
        fprintf('   ✅ 投资缺口在可接受范围内，储蓄决策与宏观资源基本一致。\n');
    end
end        
        

function plot_transition_results(results, cS, ssF, title_suffix)
            % =========================================================================
            % == 函数: plot_transition_results (版本 v1.1 - 稳健性修正版)
            % == 核心修正:
            % ==   - [!!!] 增加了对 cS.start_year 字段是否存在的检查。
            % ==   - 如果该字段不存在，则自动生成一个从 1 到 T 的时间轴，
            % ==     避免程序因 "无法识别的字段" 错误而崩溃。
            % =========================================================================
            if nargin < 4, title_suffix = ''; end
            figure('Name', ['转轨动态路径 ' title_suffix], 'Position', [100 100 1200 800]);
            T = cS.T_sim;

            % [!!! 核心修正: 增加稳健性 !!!]
            if isfield(cS, 'start_year') && ~isempty(cS.start_year)
                time_axis = cS.start_year : cS.time_Step : (cS.start_year + (T-1)*cS.time_Step);
                xlabel_text = '年份';
            else
                time_axis = 1:T;
                xlabel_text = '模型期 (t)';
            end

            ssF_Kp_line = (ones(1, T) * ssF.K_private_hat) .* cS.A_path;
            ssF_Y_line = (ones(1, T) * ssF.Y_from_production_hat) .* cS.A_path;
            ssF_C_line = (ones(1, T) * ssF.C_agg) .* cS.A_path;
            ssF_r_line_annualized = (1 + ssF.r_mkt).^(1/cS.time_Step) - 1;
            ssF_w_line = (ones(1, T) * ssF.w_hat) .* cS.A_path;
            subplot(2,2,1);
            plot(time_axis, results.K_p_path, 'b-o', 'LineWidth', 1.5, 'DisplayName', '私人资本 (K_p)');
            hold on; plot(time_axis, ssF_Kp_line, 'b--', 'LineWidth', 1, 'DisplayName', 'K_p (BGP)');
            title('私人资本路径 (原始水平)'); xlabel(xlabel_text); ylabel('资本');
            legend('show', 'Location', 'best'); grid on;
            subplot(2,2,2);
            plot(time_axis, results.Y_path, 'k-d', 'LineWidth', 1.5, 'DisplayName', '总产出 (Y)');
            hold on; plot(time_axis, results.C_path, 'g-^', 'LineWidth', 1.5, 'DisplayName', '总消费 (C)');
            plot(time_axis, ssF_Y_line, 'k--', 'LineWidth', 1, 'DisplayName', 'Y (BGP)');
            plot(time_axis, ssF_C_line, 'g--', 'LineWidth', 1, 'DisplayName', 'C (BGP)');
            title('产出与消费路径 (原始水平)'); xlabel(xlabel_text); ylabel('总量');
            legend('show', 'Location', 'best'); grid on;
            subplot(2,2,3);
            r_path_annualized = (1 + results.r_path).^(1/cS.time_Step) - 1;
            plot(time_axis, r_path_annualized * 100, 'm-x', 'LineWidth', 1.5, 'DisplayName', '年化市场利率 (r)');
            hold on; yline(ssF_r_line_annualized * 100, 'm--', 'LineWidth', 1, 'DisplayName', 'r (ssF)');
            title('市场利率路径'); xlabel(xlabel_text); ylabel('年化利率 (%)');
            legend('show', 'Location', 'best'); grid on;
            subplot(2,2,4);
            plot(time_axis, results.w_path, 'c-p', 'LineWidth', 1.5, 'DisplayName', '工资 (w)');
            hold on; plot(time_axis, ssF_w_line, 'c--', 'LineWidth', 1, 'DisplayName', 'w (BGP)');
            title('工资路径 (原始水平)'); xlabel(xlabel_text); ylabel('工资');
            legend('show', 'Location', 'best'); grid on;
            sgtitle(['宏观经济转轨动态: ' title_suffix], 'FontSize', 16, 'FontWeight', 'bold');
        end
       
                function plot_population_shock_transition(results, cS, ssF, title_suffix)
            % =========================================================================
            % == 函数: plot_population_shock_transition
            % == 版本: [v1.0 - 人口冲击专用可视化]
            % ==
            % == 目的:
            % ==   1. 详细展示受控人口冲击实验下的宏观经济转轨路径。
            % ==   2. 增加专门的人口动态图，清晰标示冲击的发生。
            % ==   3. 自动处理年度或模型期时间轴。
            % =========================================================================
            if nargin < 4, title_suffix = ''; end
            
            figure('Name', ['转轨动态: ' title_suffix], 'Position', [50 50 1400 900]);
            T = cS.T_sim;

            % --- 1. 设置时间轴 ---
            if isfield(cS, 'start_year') && ~isempty(cS.start_year)
                time_axis = cS.start_year : cS.time_Step : (cS.start_year + (T-1)*cS.time_Step);
                xlabel_text = '年份';
            else
                time_axis = 1:T;
                xlabel_text = '模型期 (t)';
            end

            % --- 2. 准备BGP参考路径 (level) ---
            ssF_Kp_line = (ones(1, T) * ssF.K_private_hat) .* cS.A_path;
            ssF_Y_line = (ones(1, T) * ssF.Y_from_production_hat) .* cS.A_path;
            ssF_C_line = (ones(1, T) * ssF.C_agg) .* cS.A_path;
            ssF_r_line_annualized = (1 + ssF.r_mkt).^(1/cS.time_Step) - 1;
            ssF_w_line = (ones(1, T) * ssF.w_hat) .* cS.A_path;

            % --- 3. 绘制核心宏观变量路径 ---
            % -- 图 A: 私人资本路径 --
            subplot(2, 3, 1);
            plot(time_axis, results.K_p_path, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', '私人资本 (K_p)');
            hold on;
            plot(time_axis, ssF_Kp_line, 'b--', 'LineWidth', 1, 'DisplayName', 'K_p (BGP终态)');
            title('私人资本路径 (原始水平)');
            xlabel(xlabel_text); ylabel('资本总量');
            legend('show', 'Location', 'best'); grid on; xlim([time_axis(1), time_axis(end)]);

            % -- 图 B: 产出与消费路径 --
            subplot(2, 3, 2);
            plot(time_axis, results.Y_path, 'k-d', 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', '总产出 (Y)');
            hold on;
            plot(time_axis, results.C_path, 'g-^', 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', '总消费 (C)');
            plot(time_axis, ssF_Y_line, 'k--', 'LineWidth', 1, 'DisplayName', 'Y (BGP终态)');
            plot(time_axis, ssF_C_line, 'g--', 'LineWidth', 1, 'DisplayName', 'C (BGP终态)');
            title('产出与消费路径 (原始水平)');
            xlabel(xlabel_text); ylabel('总量');
            legend('show', 'Location', 'best'); grid on; xlim([time_axis(1), time_axis(end)]);

            % -- 图 C: 市场利率路径 --
            subplot(2, 3, 4);
            r_path_annualized = (1 + results.r_path).^(1/cS.time_Step) - 1;
            plot(time_axis, r_path_annualized * 100, 'm-x', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '年化市场利率 (r)');
            hold on;
            yline(ssF_r_line_annualized * 100, 'm--', 'LineWidth', 1, 'DisplayName', 'r (ssF终态)');
            title('市场利率路径');
            xlabel(xlabel_text); ylabel('年化利率 (%)');
            legend('show', 'Location', 'best'); grid on; xlim([time_axis(1), time_axis(end)]);

            % -- 图 D: 工资路径 --
            subplot(2, 3, 5);
            plot(time_axis, results.w_path, 'c-p', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '工资 (w)');
            hold on;
            plot(time_axis, ssF_w_line, 'c--', 'LineWidth', 1, 'DisplayName', 'w (BGP终态)');
            title('工资路径 (原始水平)');
            xlabel(xlabel_text); ylabel('工资水平');
            legend('show', 'Location', 'best'); grid on; xlim([time_axis(1), time_axis(end)]);

            % --- 4. 绘制人口动态路径 ---
            % -- 图 E: 总人口与人口增长率 --
            subplot(2, 3, 3);
            yyaxis left
            total_pop_path = sum(cS.Z_path_raw, 1);
            p1 = plot(time_axis, total_pop_path, 'k-o', 'LineWidth', 2, 'MarkerSize', 3, 'MarkerFaceColor', 'k');
            ylabel('总人口 (绝对数量)');
            ax_pop = gca;
            ax_pop.YAxis(1).Color = 'k';

            yyaxis right
            g_n_period = total_pop_path(2:end) ./ total_pop_path(1:end-1) - 1;
            g_n_annual = (1 + g_n_period).^(1/cS.time_Step) - 1;
            time_axis_growth = time_axis(1:end-1) + cS.time_Step/2;
            p2 = plot(time_axis_growth, g_n_annual * 100, 'r-^', 'LineWidth', 1.5, 'MarkerSize', 4);
            ylabel('年化人口增长率 (%)', 'Color', 'r');
            ax_pop.YAxis(2).Color = 'r';
            
            hold on;
            % 标示冲击区域
            shock_start_time = time_axis(3); % 冲击发生在第3期
            shock_end_time = time_axis(4);   % 冲击在第4期结束
            x_patch = [shock_start_time, shock_end_time, shock_end_time, shock_start_time];
            y_limits = ylim;
            y_patch = [y_limits(1), y_limits(1), y_limits(2), y_limits(2)];
            patch(x_patch, y_patch, 'red', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', '出生率冲击期');
            
            title('总人口与增长率');
            xlabel(xlabel_text);
            legend([p1, p2], {'总人口 [左轴]', '年化增长率 [右轴]'}, 'Location', 'best');
            grid on; xlim([time_axis(1), time_axis(end)]);


            % -- 图 F: 人口结构 (老龄化) --
            subplot(2, 3, 6);
            retiree_pop = sum(cS.Z_path_raw((cS.aR_new+1):end, :), 1);
            working_pop = sum(cS.Z_path_raw(1:cS.aR_new, :), 1);
            
            old_age_dependency_ratio = (retiree_pop ./ working_pop) * 100;
            retiree_share = (retiree_pop ./ total_pop_path) * 100;
            
            p3 = plot(time_axis, old_age_dependency_ratio, 'b-d', 'LineWidth', 2, 'DisplayName', '老年抚养比');
            hold on;
            p4 = plot(time_axis, retiree_share, 'g-s', 'LineWidth', 2, 'DisplayName', '退休人口占比');
            
            % 标示冲击区域
            patch(x_patch, [0 0 100 100], 'red', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
            
            title('人口结构动态');
            xlabel(xlabel_text); ylabel('比率 (%)');
            legend([p3, p4], {'老年抚养比', '退休人口占比'}, 'Location', 'best');
            grid on; xlim([time_axis(1), time_axis(end)]); ylim([0, max(old_age_dependency_ratio)*1.1]);

            % --- 5. 添加总标题 ---
            sgtitle(['宏观经济转轨动态: ' title_suffix], 'FontSize', 16, 'FontWeight', 'bold');
        end
        
        function F_error = TPI_error_wrapper(X_guess, T, ss0, ssF, cS, paramSF, valF, polF, dist0)
            len = T - 1;
            r_path_guess_inner = reshape(X_guess(1:len), 1, len);
            w_hat_path_guess_inner = reshape(X_guess(len+1 : 2*len), 1, len);
            Bequest_hat_path_guess_inner = reshape(X_guess(2*len+1 : 3*len), 1, len);
            r_path_guess = [ss0.r_mkt, r_path_guess_inner];
            w_hat_path_guess = [ss0.w_hat, w_hat_path_guess_inner];
            Bequest_hat_path_guess = [ss0.Bequest_distributed_agg, Bequest_hat_path_guess_inner];
            pathS_iter = TPI.calculate_hh_income_paths(r_path_guess, w_hat_path_guess, Bequest_hat_path_guess, cS, ssF);
            [Pol_path, ~] = household.backward_hh(pathS_iter, cS, paramSF, valF, polF);
            [~, aggr_supply] = distribution.simulate_dist_forward(dist0, Pol_path, cS, paramSF, ss0);
            L_supply_path = aggr_supply.L_path;
            K_p_hat_supply_path = aggr_supply.K_p_hat_path;
            r_target_path = zeros(1, T);
            w_hat_target_path = zeros(1, T);
            K_g_hat_fixed = ssF.K_public_hat;
            for t = 1:T
                prices_t = firm.get_prices_at_t(K_p_hat_supply_path(t), K_g_hat_fixed, L_supply_path(t), cS);
                r_target_path(t) = prices_t.r_mkt_t;
                w_hat_target_path(t) = prices_t.w_hat_t;
            end
            Bequest_hat_target_path = zeros(1, T);
            Bequest_hat_target_path(1) = ss0.Bequest_distributed_agg;
            Bequest_hat_target_path(2:T) = aggr_supply.Bequest_gen_hat_path(1:T-1);
            error_r = r_target_path(2:T) - r_path_guess(2:T);
            error_w = w_hat_path_guess(2:T) - w_hat_path_guess(2:T);
            error_beq = Bequest_hat_target_path(2:T) - Bequest_hat_path_guess(2:T);
            F_error = [error_r, error_w, error_beq]';
        end
    end
end
