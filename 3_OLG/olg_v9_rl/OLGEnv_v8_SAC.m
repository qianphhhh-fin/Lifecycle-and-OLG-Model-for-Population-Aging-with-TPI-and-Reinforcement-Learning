% --- START OF FILE OLGEnv_v8_SAC.m ---
classdef OLGEnv_v8_SAC < rl.env.MATLABEnvironment
    %OLGENV_V8_SAC OLG环境，用于SAC Agent训练 (期望化死亡率版本)
    %   家庭在给定的宏观经济状态M下进行生命周期决策
    %   观测: [k, k_pps, age_idx, eps_idx, M_vars(6)] (归一化)
    %   行动: [prop_pps_contrib, prop_non_pps_saving] (0-1)
    %   特性: 死亡率通过期望值纳入，确保固定长度的episodes

    properties
        cS % OLG模型参数 (来自ParameterValues_HuggettStyle)
        paramS_rl % RL相关的派生参数 (如劳动效率过程)
        rng_M % 宏观变量的采样范围

        % 当前episode的宏观经济状态 (在reset时采样)
        current_M struct
        current_bV_payg % 当前episode的PAYG福利向量 (由current_M.b_payg_avg_retiree生成)

        % 当前个体状态
        current_age_idx % 模型年龄组索引 (1 to aD_new)
        current_k_val
        current_k_pps_val
        current_eps_idx % 劳动效率状态索引 (1 to nw)

        % RL规范
        ObservationInfo_
        ActionInfo_
        
        % 辅助归一化参数 (在reset中计算一次，或预定义)
        obs_norm_min
        obs_norm_max
        obs_norm_mean
        obs_norm_range
    end

    methods
        % 构造函数
        function this = OLGEnv_v8_SAC(cS_input, paramS_rl_input, rng_M_input, obsInfo, actInfo)
            % --- 调用父类构造函数 ---
            this = this@rl.env.MATLABEnvironment(obsInfo, actInfo); % <--- 修改点

            % --- 然后再进行子类的初始化 ---
            % this.ObservationInfo_ = obsInfo; % 父类构造函数已经处理了这个
            % this.ActionInfo_ = actInfo;   % 父类构造函数已经处理了这个
            
            this.cS = cS_input;
            this.paramS_rl = paramS_rl_input;
            this.rng_M = rng_M_input;
            
            % 初始化 current_M 结构体，避免reset时的空结构体错误
            this.current_M = struct();
            this.current_M.R_k_net_factor = 0;
            this.current_M.w_gross = 0;
            this.current_M.TR_total = 0;
            this.current_M.tau_l = 0;
            this.current_M.theta_payg_actual = 0;
            this.current_M.b_payg_avg_for_obs = 0;
            
            % 初始化其他状态变量
            this.current_bV_payg = zeros(1, this.cS.aD_new);
            this.current_age_idx = 1;
            this.current_k_val = this.cS.kMin;
            this.current_k_pps_val = this.cS.kppsMin;
            this.current_eps_idx = 1;

            % 初始化归一化参数 (更稳健的做法是基于数据或合理范围)
            % [k, k_pps, age_idx, eps_idx, R_k_net, w_g, TR, b_payg_avg, tau_l, theta_payg]
            this.obs_norm_min = [this.cS.kMin, this.cS.kppsMin, 1, 1, ...
                                 this.rng_M.R_k_net_factor(1), this.rng_M.w_gross(1), ...
                                 this.rng_M.TR_total(1), this.rng_M.b_payg_avg_retiree(1), ...
                                 this.rng_M.tau_l(1), this.rng_M.theta_payg_actual(1)];
            this.obs_norm_max = [this.cS.kMax, this.cS.kppsMax, this.cS.aD_new, this.cS.nw, ...
                                 this.rng_M.R_k_net_factor(2), this.rng_M.w_gross(2), ...
                                 this.rng_M.TR_total(2), this.rng_M.b_payg_avg_retiree(2), ...
                                 this.rng_M.tau_l(2), this.rng_M.theta_payg_actual(2)];
            this.obs_norm_range = this.obs_norm_max - this.obs_norm_min;
            this.obs_norm_range(this.obs_norm_range < 1e-6) = 1; % 避免除以零
            this.obs_norm_mean = (this.obs_norm_min + this.obs_norm_max) / 2;
        end

        

        % 重置环境 (每个episode开始时调用)
        function initialObservation = reset(this)
            % 1. 采样新的宏观经济状态 M for this episode
            this.current_M.R_k_net_factor = this.rng_M.R_k_net_factor(1) + rand() * (this.rng_M.R_k_net_factor(2) - this.rng_M.R_k_net_factor(1));
            this.current_M.w_gross = this.rng_M.w_gross(1) + rand() * (this.rng_M.w_gross(2) - this.rng_M.w_gross(1));
            this.current_M.TR_total = this.rng_M.TR_total(1) + rand() * (this.rng_M.TR_total(2) - this.rng_M.TR_total(1));
            b_payg_avg = this.rng_M.b_payg_avg_retiree(1) + rand() * (this.rng_M.b_payg_avg_retiree(2) - this.rng_M.b_payg_avg_retiree(1));
            this.current_M.tau_l = this.rng_M.tau_l(1) + rand() * (this.rng_M.tau_l(2) - this.rng_M.tau_l(1));
            this.current_M.theta_payg_actual = this.rng_M.theta_payg_actual(1) + rand() * (this.rng_M.theta_payg_actual(2) - this.rng_M.theta_payg_actual(1));
            
            % 从 b_payg_avg 生成 bV_payg (与VFI中类似)
            this.current_bV_payg = zeros(1, this.cS.aD_new);
            if this.cS.aR_new < this.cS.aD_new
                this.current_bV_payg(this.cS.aR_new + 1 : this.cS.aD_new) = b_payg_avg;
            end
            this.current_M.b_payg_avg_for_obs = b_payg_avg; % 用于观测向量

            % 2. 初始化个体状态 (例如，年龄为1，资产为0，效率随机抽取)
            this.current_age_idx = 1;
            this.current_k_val = this.cS.kMin; % 或者从一个小的初始分布中抽样
            this.current_k_pps_val = this.cS.kppsMin;
            % 初始效率状态 (从paramS_rl.leProb1V中抽取)
            this.current_eps_idx = find(rand() <= cumsum(this.paramS_rl.leProb1V), 1, 'first');
            if isempty(this.current_eps_idx), this.current_eps_idx = 1; end

            initialObservation = getObservation(this);
        end

        % 执行一步 (agent采取行动后环境的响应)
        function [nextObservation, reward, isDone, loggedSignals] = step(this, action)
            loggedSignals = []; % 可用于记录额外信息

            % 0. 解析行动 (比例变量)
            prop_pps_contrib = action(1);
            prop_non_pps_saving = action(2);
            
            prop_pps_contrib = max(0, min(1, prop_pps_contrib)); % 确保在[0,1]
            prop_non_pps_saving = max(0, min(1, prop_non_pps_saving));

            % 1. 计算实际的PPS缴费 c_pps
            actual_c_pps = 0;
            max_permissible_cpps_for_action = 0;
            current_epsilon_val = this.paramS_rl.leGridV(this.current_eps_idx);
            
            if this.current_age_idx <= this.cS.aR_new && ...
               this.cS.physAgeMap{this.current_age_idx}(1) <= this.cS.pps_contribution_age_max_idx && ...
               this.cS.pps_active
                
                age_efficiency = this.cS.ageEffV_new(this.current_age_idx);
                current_gross_labor_income = this.current_M.w_gross * age_efficiency * current_epsilon_val;
                
                if current_gross_labor_income > 1e-6
                    max_cpps_by_frac = current_gross_labor_income * this.cS.pps_max_contrib_frac;
                    max_permissible_cpps_for_action = min(this.cS.pps_contrib_limit, max_cpps_by_frac);
                    max_permissible_cpps_for_action = max(0, max_permissible_cpps_for_action);
                end
                actual_c_pps = prop_pps_contrib * max_permissible_cpps_for_action;
                actual_c_pps = max(0, min(actual_c_pps, max_permissible_cpps_for_action)); % 再次确保
            end

            % 2. 计算可用于非PPS储蓄和消费的资源
            % 构建临时的paramS_hh用于HHIncome_Huggett
            paramS_hh_step.tau_l = this.current_M.tau_l;
            paramS_hh_step.theta_payg_actual_for_hh = this.current_M.theta_payg_actual;
            paramS_hh_step.pps_tax_deferral_active = this.cS.pps_active; % 或从 M 获取
            
            b_payg_this_age = this.current_bV_payg(this.current_age_idx);

            [resources_after_pps, ~, ~] = main_olg_v8_utils.HHIncome_Huggett(...
                this.current_k_val, this.current_M.R_k_net_factor, this.current_M.w_gross, ...
                this.current_M.TR_total, b_payg_this_age, actual_c_pps, ...
                this.current_age_idx, paramS_hh_step, this.cS, current_epsilon_val);

            % 3. 计算实际的非PPS储蓄 k_prime 和消费 c
            consumption_floor_spending = this.cS.cFloor * (1 + this.cS.tau_c);
            resources_for_kprime_and_c_above_floor = resources_after_pps - consumption_floor_spending;
            
            actual_k_prime = 0;
            current_c = this.cS.cFloor;

            if resources_for_kprime_and_c_above_floor >= 0
                % prop_non_pps_saving 应用于 (总资源 - 最低消费开销)
                actual_k_prime = prop_non_pps_saving * resources_for_kprime_and_c_above_floor;
                % 确保k_prime在可行范围内
                actual_k_prime = max(this.cS.kMin, min(actual_k_prime, resources_for_kprime_and_c_above_floor));
            else
                % 如果资源不足以支付最低消费，则k_prime为kMin
                actual_k_prime = this.cS.kMin;
            end
            actual_k_prime = max(this.cS.kMin, min(actual_k_prime, this.cS.kMax)); % 最终约束

            consumption_expenditure = resources_after_pps - actual_k_prime;
            current_c = max(this.cS.cFloor, consumption_expenditure / (1 + this.cS.tau_c));
            
            % 4. 计算奖励 (当期效用，期望化死亡率处理)
            [~, utility] = main_olg_v8_utils.CES_utility(current_c, this.cS.sigma, this.cS);
            
            % 获取存活概率用于期望化处理
            survival_prob = 1.0;  % 默认值
            if this.current_age_idx <= length(this.cS.s_1yr_transitionV)
                survival_prob = this.cS.s_1yr_transitionV(this.current_age_idx);
            end
            
            if ~isfinite(utility) || utility < -1e9 % 惩罚无效消费
                reward = -1000 - abs(current_c - this.cS.cFloor)*100; % 重罚
            else
                % 期望化处理：死亡率通过RL的价值函数学习自动纳入期望
                % 保持原始效用，让SAC的折现机制和价值函数学习处理死亡风险
                reward = utility;
            end
            
            % 5. 更新个体状态到下一期
            % 5.1 PPS资产演化
            k_pps_next = this.current_k_pps_val; % 默认
            if this.cS.pps_active
                pps_withdrawal = 0;
                annual_age_phys = this.cS.physAgeV_new(this.current_age_idx); % 代表性年度年龄
                % 注意: cS.pps_withdrawal_age_min_idx 是年度年龄索引, 不是模型年龄组索引
                % 为了简化，我们这里也基于模型年龄组判断是否退休
                is_retired_model_group = (this.current_age_idx > this.cS.aR_new);
                if is_retired_model_group % && annual_age_phys >= this.cS.physAgeV_orig(this.cS.pps_withdrawal_age_min_idx)
                    pps_withdrawal = this.current_k_pps_val * this.cS.pps_withdrawal_rate;
                end
                
                pps_return_factor = 1 + ((this.current_M.R_k_net_factor - 1) + this.cS.pps_return_rate_premium);
                k_pps_next_unclamped = (this.current_k_pps_val + actual_c_pps - pps_withdrawal) * pps_return_factor;
                k_pps_next = max(this.cS.kppsMin, min(this.cS.kppsMax, k_pps_next_unclamped));
            end
            
            % 5.2 非PPS资产由 actual_k_prime 决定
            k_next = actual_k_prime;

            % 5.3 年龄演化
            age_next_idx = this.current_age_idx + 1;
            
            % 5.4 效率冲击演化 (基于 paramS_rl.leTrProbM)
            eps_next_idx = this.current_eps_idx; % 默认
            if this.current_age_idx < this.cS.aD_new % 只有在不是最后年龄时才转移
                 trans_probs_eps = this.paramS_rl.leTrProbM(this.current_eps_idx, :);
                 eps_next_idx = find(rand() <= cumsum(trans_probs_eps), 1, 'first');
                 if isempty(eps_next_idx), eps_next_idx = this.current_eps_idx; end % 保持原状以防万一
            end

            % 6. 判断是否终止 (只在到达最大年龄时终止，期望化死亡率处理)
            isDone = false;
            if age_next_idx > this.cS.aD_new
                isDone = true;
                % 可以考虑添加终期奖励，如遗赠效用
                % 当前模型中，自然死亡即结束，无额外遗赠效用
            end
            
            % 移除随机死亡：死亡率现在通过期望值纳入RL的价值函数学习
            % 这确保每个episode都是完整的生命周期，提高训练稳定性
            % RL agent会自动学习在高死亡率年龄段调整决策
            
            % 7. 更新内部状态 (如果未终止)
            if ~isDone
                this.current_age_idx = age_next_idx;
                this.current_k_val = k_next;
                this.current_k_pps_val = k_pps_next;
                this.current_eps_idx = eps_next_idx;
            end

            nextObservation = getObservation(this);
        end

        % 获取当前观测 (并归一化)
        function obs = getObservation(this)
            % 原始观测向量
            raw_obs_vec = [this.current_k_val, this.current_k_pps_val, ...
                           this.current_age_idx, this.current_eps_idx, ...
                           this.current_M.R_k_net_factor, this.current_M.w_gross, ...
                           this.current_M.TR_total, this.current_M.b_payg_avg_for_obs, ...
                           this.current_M.tau_l, this.current_M.theta_payg_actual];
            
            % 归一化 (min-max scaling to [0,1] or [-1,1] or z-score)
            % 这里使用简单的min-max到[0,1]
            obs = (raw_obs_vec - this.obs_norm_min) ./ this.obs_norm_range;
            obs = max(0, min(1, obs)); % 确保在[0,1]范围内
            
            % 确保是行向量
            obs = obs(:)';
        end
        
        % (可选) 设置宏观参数的方法，用于评估
        function setMacroParameters(this, M_fixed)
            this.current_M = M_fixed;
            % 并据此更新 current_bV_payg
            this.current_bV_payg = zeros(1, this.cS.aD_new);
            if this.cS.aR_new < this.cS.aD_new
                this.current_bV_payg(this.cS.aR_new + 1 : this.cS.aD_new) = M_fixed.b_payg_avg_retiree; % 假设M_fixed有此字段
            end
             this.current_M.b_payg_avg_for_obs = M_fixed.b_payg_avg_retiree;
            fprintf('Environment macro parameters set externally for evaluation.\n');
        end
    end
end
% --- END OF FILE OLGEnv_v8_SAC.m ---