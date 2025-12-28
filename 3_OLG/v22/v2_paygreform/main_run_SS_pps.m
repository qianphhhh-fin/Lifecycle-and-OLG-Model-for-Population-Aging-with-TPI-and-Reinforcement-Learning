function main_run_SS_pps(pps_path_mode)
    % =========================================================================
    % == SCRIPT: main_run_SS_pps.m
    % == 版本: [v3.2 - 结构性改革模式]
    % ==
    % == 核心修改:
    % ==   - 脚本被封装为一个函数，接收 pps_path_mode 参数。
    % ==   - 在求解终期稳态 (ssF) 时，所有政策参数都设置为改革完成后的值：
    % ==     - theta_path_h 设为最终的低水平。
    % ==     - pps_active 设为 true。
    % ==     - pps_fixed 设为最终的高水平。
    % ==   - 保存的文件名将包含 pps_path_mode，以区分不同改革情景的数据。
    % =========================================================================
    
    if nargin < 1
        pps_path_mode = 'structural_reform';
    end
    fprintf('=== OLG异质性模型稳态求解 (模式: %s) ===\n\n', pps_path_mode);

    % --- 1. 全局设置与初始化 (这部分与之前基本一致) ---
    cS = utils.ParameterValues();
    cS.nk = 30;
    cS.ss0_year = 2023;
    cS.start_year = 2023;
    cS.pop_data_last_year = 2070;
    cS.g_A_ss = 0.015;
    cS.transition_period_years = 100;
    cS.num_hh_types = 4;
    nH = cS.num_hh_types;

    % --- 2. 生成完整的外生路径 ---
    [Z_path, Z_path_raw, cS] = population.generate_Z_path(cS, false);
    T_sim = size(Z_path_raw, 2);
    cS.T_sim = T_sim;
    [cS.A_path, g_path_annual_output, ~] = utils.generate_tfp_path(cS, 'baseline', false);
    [~, ~, ~, cS.nw_expanded] = utils.EarningProcess_AgeDependent(cS);
    
    % [核心] 根据模式生成不同的政策路径
    cS = utils.theta_payg_path(cS, pps_path_mode, false);
    cS = utils.retire_age_path(cS, 'mild', Z_path_raw, false);
    cS = utils.generate_PPS_path(cS, pps_path_mode, false);
    
    cS.Z_path = Z_path;
    cS.Z_path_raw = Z_path_raw;

    % --- 3. 求解初始伪稳态 (ss0, t=0) ---
    % [逻辑不变] ss0 总是代表改革前的世界
    fprintf('\n--- 3. 求解初始伪稳态 (ss0) ---\n');
    cS0 = cS;
    cS0.pps_active = false;
    cS0.nkpps = 1; cS0.n_pps_rate_grid = 1;
    cS0 = utils.generateGrids(cS0);
    cS0.aR_new = cS.aR_new_path(1);
    cS0.n_ss = ((sum(Z_path_raw(:,2))/sum(Z_path_raw(:,1)))^(1/cS.time_Step))-1;
    cS0.s_pathV = cS.s_pathV(:,1);
    % [重要] 初始稳态使用改革前的缴费率
    cS_temp_baseline = utils.theta_payg_path(cS, 'baseline', false);
    cS0.theta_path_h = cS_temp_baseline.theta_path_h(:,1);
    cS0.g_A_ss = g_path_annual_output(1);
    Z0_norm = Z_path_raw(:,1) / sum(Z_path_raw(:,1)) * cS.type_weights';
    paramS0 = struct();
    [paramS0.leGridV, paramS0.TrProbM_by_age, paramS0.leProb1V, ~] = utils.EarningProcess_AgeDependent(cS0);
    params_for_ss0 = struct('Z', Z0_norm, 'A', 1.0, 'g_A_ss', 0, 'n_ss', cS0.n_ss);
    [ss0, dist0_h, ~, ~] = SS.solve_steady_state(cS0, paramS0, params_for_ss0, true, false);
    if isempty(ss0), error('初始稳态(ss0)求解失败！'); end
    fprintf('✅ 初始稳态(ss0)求解成功！\n');

    % --- 4. 求解终期稳态 (ssF, t=T) ---
    % [核心修改] ssF 代表改革完成后的新世界
    fprintf('\n--- 4. 求解终期BGP稳态 (ssF) ---\n');
    cSF = cS;
    cSF.pps_active = true; % 终期PPS是激活的
    cSF.nkpps = 20; cSF.n_pps_rate_grid = 10;
    cSF.pps_fixed = cS.pps_fixed_path(end); % 使用路径的终点值
    cSF.aR_new = cS.aR_new_path(end);
    cSF.g_A_ss = g_path_annual_output(end);
    cSF = utils.generateGrids(cSF);
    cSF.theta_path_h = cS.theta_path_h(:,end); % 使用路径的终点值
    cSF.s_pathV = cS.s_pathV(:,end);
    ZF_norm = Z_path_raw(:,end) / sum(Z_path_raw(:,end)) * cS.type_weights';
    paramSF = struct();
    [paramSF.leGridV, paramSF.TrProbM_by_age, paramSF.leProb1V, ~] = utils.EarningProcess_AgeDependent(cSF);
    params_for_ssF = struct('Z', ZF_norm, 'A', 1.0, 'g_A_ss', cSF.g_A_ss, 'n_ss', cSF.n_ss);
    [ssF, distF_h, polF_h, valF_h] = SS.solve_steady_state(cSF, paramSF, params_for_ssF, true, true);
    if isempty(ssF), error('终期稳态(ssF)求解失败！'); end
    fprintf('✅ 终期稳态(ssF)求解成功！\n');
    
    % --- 5. 保存转轨所需数据 ---
    if ~isempty(fieldnames(ssF))
        fprintf('\n--- 5. 保存转轨分析所需数据 ---\n');
        if ~exist('SS', 'dir'), mkdir('SS'); end
        save_filename = ['SS/data_for_het_transition_pps_', pps_path_mode, '.mat'];

        dist0_abs_h = cell(nH, 1);
        for h = 1:nH
            pop_h_by_age_t1 = cS.Z_path_raw(:, 1) * cS.type_weights(h);
            dist_cond_h = zeros(size(dist0_h{h}));
            for ia = 1:cS.aD_new
                pop_mass_norm_hia = sum(dist0_h{h}(:,:,:,ia), 'all');
                if pop_mass_norm_hia > 1e-12
                    dist_cond_h(:,:,:,ia) = dist0_h{h}(:,:,:,ia) / pop_mass_norm_hia;
                end
            end
            dist0_abs_h{h} = dist_cond_h .* reshape(pop_h_by_age_t1, [1, 1, 1, cS.aD_new]);
        end
        
        cS_for_trans = cSF; 
        cS_for_trans.theta_path_h = cS.theta_path_h;
        cS_for_trans.s_pathV = cS.s_pathV;
        cS_for_trans.Z_path_raw = cS.Z_path_raw;
        cS_for_trans.type_weights = cS.type_weights;
        cS_for_trans.ageEffV_new_h = cS.ageEffV_new_h;
        cS_for_trans.aR_new_path = cS.aR_new_path;

        data_for_transition = struct();
        data_for_transition.ss0 = ss0;
        data_for_transition.ssF = ssF;
        data_for_transition.dist0_h = dist0_abs_h;
        data_for_transition.distF_h = distF_h;
        data_for_transition.polF_h = polF_h;
        data_for_transition.valF_h = valF_h;
        data_for_transition.cS = cS_for_trans;
        data_for_transition.paramSF = paramSF;
        save(save_filename, 'data_for_transition', '-v7.3');
        fprintf('✅ 所有转轨所需数据已成功保存至:\n   %s\n', save_filename);
    end
end