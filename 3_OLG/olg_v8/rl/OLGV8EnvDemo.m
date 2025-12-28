classdef OLGV8EnvDemo < rl.env.MATLABEnvironment
    % OLG V8 强化学习环境类
    
    properties
        % 模型参数
        cS
        fixedMacro
        
        % 环境状态
        State  % [k, k_pps, eps_idx, age_idx]
        EpisodeCount
    end
    
    methods
        function this = OLGV8EnvDemo(envConstantParams)
            % 构造函数
            
            % 定义观察空间
            obsInfo = rlNumericSpec([4 1]);
            obsInfo.Name = 'OLG V8 States';
            obsInfo.Description = 'Normalized k, k_pps, epsilon_index, normalized_age_group_index';
            
            % 定义动作空间
            actInfo = rlNumericSpec([2 1]);
            actInfo.Name = 'OLG V8 Actions';
            actInfo.LowerLimit = [-1; -1];
            actInfo.UpperLimit = [1; 1];
            
            % 调用父类构造函数
            this = this@rl.env.MATLABEnvironment(obsInfo, actInfo);
            
            % 初始化环境参数
            this.cS = envConstantParams.cS;
            this.fixedMacro = envConstantParams.fixedMacro;
            this.State = zeros(4,1);
            this.EpisodeCount = 0;
        end
        
        function [Observation, Reward, IsDone, LoggedSignals] = step(this, Action)
            % 执行一步环境交互
            
            LoggedSignals = struct();
            cS = this.cS;
            fixedMacro = this.fixedMacro;
            
            % 获取当前状态
            k_now = this.State(1);
            k_pps_now = this.State(2);
            eps_idx_now = this.State(3);
            age_idx_now = this.State(4);
            
            % 安全地获取epsilon值
            if isfield(cS, 'leGridV') && length(cS.leGridV) >= eps_idx_now
                epsilon_val_now = cS.leGridV(eps_idx_now);
            else
                epsilon_val_now = 1.0; % 默认值
            end
            
            % 反归一化动作
            action_denorm = this.denormalizeAndConstrainAction(...
                Action, k_now, k_pps_now, epsilon_val_now, age_idx_now);
            k_prime_chosen = action_denorm(1);
            c_pps_chosen_effective = action_denorm(2);
            
            % 计算家庭收入和资源
            paramS_hh_step = struct('tau_l', fixedMacro.tau_l, ...
                                    'theta_payg_actual_for_hh', fixedMacro.theta_payg_actual_for_hh, ...
                                    'pps_tax_deferral_active', cS.pps_active);
            
            % 安全地添加可选字段
            if isfield(cS, 'leGridV')
                paramS_hh_step.leGridV = cS.leGridV;
            end
            if isfield(cS, 'leTrProbM')
                paramS_hh_step.leTrProbM = cS.leTrProbM;
            end
            
            % 简化的收入计算（如果main_olg_v8_utils不可用）
            try
                [resources_for_c_k_prime, ~, ~] = main_olg_v8_utils.HHIncome_Huggett(...
                    k_now, fixedMacro.R_k_net_factor_hh, fixedMacro.MPL_gross, ...
                    fixedMacro.TR_total, fixedMacro.bV_payg(age_idx_now), c_pps_chosen_effective, ...
                    age_idx_now, paramS_hh_step, cS, epsilon_val_now);
            catch
                % 简化的收入计算
                if isfield(cS, 'ageEffV_new') && length(cS.ageEffV_new) >= age_idx_now
                    age_efficiency = cS.ageEffV_new(age_idx_now);
                else
                    age_efficiency = 1.0;
                end
                labor_income = fixedMacro.MPL_gross * age_efficiency * epsilon_val_now;
                capital_income = k_now * (fixedMacro.R_k_net_factor_hh - 1);
                resources_for_c_k_prime = labor_income + capital_income + k_now + fixedMacro.TR_total;
            end
            
            % 计算消费
            consumption_expenditure = resources_for_c_k_prime - k_prime_chosen;
            
            % 安全地获取参数
            if isfield(cS, 'cFloor')
                cFloor_val = cS.cFloor;
            else
                cFloor_val = 0.01;
            end
            
            if isfield(cS, 'tau_c')
                tau_c_val = cS.tau_c;
            else
                tau_c_val = 0.05;
            end
            
            current_c_quantity = max(cFloor_val, consumption_expenditure / (1 + tau_c_val));
            
            % 计算奖励（当前效用）
            try
                [~, utility_now] = main_olg_v8_utils.CES_utility(current_c_quantity, cS.sigma, cS);
            catch
                % 简化的效用函数
                if isfield(cS, 'sigma')
                    sigma_val = cS.sigma;
                else
                    sigma_val = 2;
                end
                
                if sigma_val == 1
                    utility_now = log(current_c_quantity);
                else
                    utility_now = (current_c_quantity^(1-sigma_val) - 1) / (1-sigma_val);
                end
            end
            
            Reward = utility_now;
            
            % 惩罚过低消费
            if current_c_quantity < cFloor_val + 1e-3
                Reward = Reward - 10;
            end
            
            if ~isfinite(Reward)
                Reward = -1e7;
            end
            
            % 确定下一状态
            k_next = k_prime_chosen;
            
            % 下一期PPS资产
            pps_withdrawal_pretax = 0;
            if isfield(cS, 'physAgeMap') && length(cS.physAgeMap) >= age_idx_now
                model_age_group_start_year_idx = cS.physAgeMap{age_idx_now}(1);
            else
                model_age_group_start_year_idx = age_idx_now;
            end
            
            if isfield(cS, 'aR_new')
                is_retired_group_now = (age_idx_now > cS.aR_new);
            else
                is_retired_group_now = (age_idx_now > 7); % 默认工作年龄组数
            end
            
            % 安全地获取PPS参数
            if isfield(cS, 'pps_withdrawal_age_min_idx')
                pps_withdrawal_age_min = cS.pps_withdrawal_age_min_idx;
            else
                pps_withdrawal_age_min = 65;
            end
            
            if isfield(cS, 'pps_active')
                pps_active_val = cS.pps_active;
            else
                pps_active_val = 1;
            end
            
            if isfield(cS, 'pps_withdrawal_rate')
                pps_withdrawal_rate_val = cS.pps_withdrawal_rate;
            else
                pps_withdrawal_rate_val = 0.04;
            end
            
            if isfield(cS, 'pps_return_rate_premium')
                pps_return_premium = cS.pps_return_rate_premium;
            else
                pps_return_premium = 0.01;
            end
            
            if is_retired_group_now && model_age_group_start_year_idx >= pps_withdrawal_age_min && pps_active_val
                pps_withdrawal_pretax = k_pps_now * pps_withdrawal_rate_val;
            end
            
            pps_return_factor = 1 + ((fixedMacro.R_k_net_factor_hh - 1) + pps_return_premium);
            k_pps_next = (k_pps_now + c_pps_chosen_effective - pps_withdrawal_pretax) * pps_return_factor;
            
            % 安全地获取边界值
            if isfield(cS, 'kppsMin')
                kppsMin_val = cS.kppsMin;
            else
                kppsMin_val = 0;
            end
            
            if isfield(cS, 'kppsMax')
                kppsMax_val = cS.kppsMax;
            else
                kppsMax_val = 5;
            end
            k_pps_next = max(kppsMin_val, min(kppsMax_val, k_pps_next));
            
            % 下一期epsilon冲击
            if age_idx_now < cS.aD_new
                if isfield(cS, 'leTrProbM') && size(cS.leTrProbM, 1) >= eps_idx_now
                    trans_probs = cS.leTrProbM(eps_idx_now, :);
                    eps_idx_next = find(rand <= cumsum(trans_probs), 1, 'first');
                else
                    eps_idx_next = randi([1, cS.nw]); % 随机选择
                end
            else
                eps_idx_next = eps_idx_now;
            end
            
            % 下一期年龄
            age_idx_next = age_idx_now + 1;
            
            % 检查是否结束
            IsDone = false;
            if age_idx_next > cS.aD_new
                IsDone = true;
                age_idx_next = cS.aD_new;
            end
            
            % 更新状态
            this.State(1) = k_next;
            this.State(2) = k_pps_next;
            this.State(3) = eps_idx_next;
            this.State(4) = age_idx_next;
            
            % 归一化观察
            Observation = this.normalizeObservation(this.State);
            
            % 记录信号
            LoggedSignals.k_now = k_now;
            LoggedSignals.k_pps_now = k_pps_now;
            LoggedSignals.c_pps_chosen = c_pps_chosen_effective;
            LoggedSignals.k_prime_chosen = k_prime_chosen;
            LoggedSignals.consumption = current_c_quantity;
            LoggedSignals.raw_reward = utility_now;
        end
        
        function InitialObservation = reset(this)
            % 重置环境
            
            this.EpisodeCount = this.EpisodeCount + 1;
            cS = this.cS;
            
            % 初始年龄
            this.State(4) = 1;
            
            % 初始资产
            if isfield(cS, 'kMin')
                kMin_val = cS.kMin;
            else
                kMin_val = 0;
            end
            
            if isfield(cS, 'kppsMin')
                kppsMin_val = cS.kppsMin;
            else
                kppsMin_val = 0;
            end
            
            if isfield(cS, 'nw')
                nw_val = cS.nw;
            else
                nw_val = 3;
            end
            
            this.State(1) = kMin_val + 1e-3;
            this.State(2) = kppsMin_val + 1e-3;
            
            % 初始epsilon冲击
            if isfield(cS,'leProb1V') && ~isempty(cS.leProb1V)
                P0 = cS.leProb1V;
                this.State(3) = find(rand <= cumsum(P0), 1, 'first');
            else
                this.State(3) = randi([1, nw_val]);
            end
            
            % 归一化观察
            InitialObservation = this.normalizeObservation(this.State);
        end
    end
    
    methods (Access = public)
        function obs_norm = normalizeObservation(this, state_unnorm)
            % 归一化状态
            cS = this.cS;
            obs_norm = zeros(4,1);
            
            % 安全地获取边界值
            if isfield(cS, 'kMin')
                kMin_val = cS.kMin;
            else
                kMin_val = 0;
            end
            
            if isfield(cS, 'kMax')
                kMax_val = cS.kMax;
            else
                kMax_val = 10;
            end
            
            if isfield(cS, 'kppsMin')
                kppsMin_val = cS.kppsMin;
            else
                kppsMin_val = 0;
            end
            
            if isfield(cS, 'kppsMax')
                kppsMax_val = cS.kppsMax;
            else
                kppsMax_val = 5;
            end
            
            if isfield(cS, 'nw')
                nw_val = cS.nw;
            else
                nw_val = 3;
            end
            
            if isfield(cS, 'aD_new')
                aD_new_val = cS.aD_new;
            else
                aD_new_val = 10;
            end
            
            obs_norm(1) = (state_unnorm(1) - kMin_val) / (kMax_val - kMin_val + 1e-6);
            obs_norm(2) = (state_unnorm(2) - kppsMin_val) / (kppsMax_val - kppsMin_val + 1e-6);
            obs_norm(3) = (state_unnorm(3)-1) / (nw_val -1 + 1e-6);
            if nw_val == 1, obs_norm(3) = 0.5; end
            obs_norm(4) = (state_unnorm(4)-1) / (aD_new_val -1 + 1e-6);
            obs_norm = max(0, min(1, obs_norm));
        end
        
        function action_denorm_constrained = denormalizeAndConstrainAction(this, action_norm, k_now, k_pps_now, epsilon_val_now, age_idx_now)
            % 反归一化并约束动作
            cS = this.cS;
            fixedMacro = this.fixedMacro;
            
            % 动作1: k_prime
            % 安全地获取边界值
            if isfield(cS, 'kMin')
                kMin_val = cS.kMin;
            else
                kMin_val = 0;
            end
            
            if isfield(cS, 'kMax')
                kMax_val = cS.kMax;
            else
                kMax_val = 10;
            end
            
            k_prime_denorm = kMin_val + (action_norm(1) + 1)/2 * (kMax_val - kMin_val);
            k_prime_denorm = max(kMin_val, min(kMax_val, k_prime_denorm));
            
            % 动作2: c_pps
            c_pps_denorm_effective = 0;
            
            % 安全地获取年龄组信息
            if isfield(cS, 'physAgeMap') && length(cS.physAgeMap) >= age_idx_now
                model_age_group_start_year_idx = cS.physAgeMap{age_idx_now}(1);
            else
                model_age_group_start_year_idx = age_idx_now;
            end
            
            % 检查PPS贡献资格
            if isfield(cS, 'aR_new')
                aR_new_val = cS.aR_new;
            else
                aR_new_val = 7;
            end
            
            if isfield(cS, 'pps_contribution_age_max_idx')
                pps_contrib_age_max = cS.pps_contribution_age_max_idx;
            else
                pps_contrib_age_max = 60;
            end
            
            if isfield(cS, 'pps_active')
                pps_active_val = cS.pps_active;
            else
                pps_active_val = 1;
            end
            
            is_pps_contrib_eligible = (age_idx_now <= aR_new_val && ...
                                    model_age_group_start_year_idx <= pps_contrib_age_max && ...
                                    pps_active_val);
            
            if is_pps_contrib_eligible
                % 安全地获取年龄效率
                if isfield(cS, 'ageEffV_new') && length(cS.ageEffV_new) >= age_idx_now
                    age_efficiency = cS.ageEffV_new(age_idx_now);
                else
                    age_efficiency = 1.0; % 默认效率
                end
                
                current_gross_labor_income = fixedMacro.MPL_gross * age_efficiency * epsilon_val_now;
                
                if current_gross_labor_income > 1e-6
                    if isfield(cS, 'pps_max_contrib_frac')
                        pps_max_frac = cS.pps_max_contrib_frac;
                    else
                        pps_max_frac = 0.1;
                    end
                    
                    if isfield(cS, 'pps_annual_contrib_limit')
                        pps_annual_limit = cS.pps_annual_contrib_limit;
                    else
                        pps_annual_limit = 1000;
                    end
                    
                    max_cpps_by_frac = current_gross_labor_income * pps_max_frac;
                    max_permissible_cpps = min(pps_annual_limit, max_cpps_by_frac);
                    max_permissible_cpps = max(0, max_permissible_cpps);
                    
                    c_pps_denorm_effective = (action_norm(2) + 1)/2 * max_permissible_cpps;
                    c_pps_denorm_effective = max(0, min(max_permissible_cpps, c_pps_denorm_effective));
                end
            end
            
            action_denorm_constrained = [k_prime_denorm; c_pps_denorm_effective];
        end
    end
end 