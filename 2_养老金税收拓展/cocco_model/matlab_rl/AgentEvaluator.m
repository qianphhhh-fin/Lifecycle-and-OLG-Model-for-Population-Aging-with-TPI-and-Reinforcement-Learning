classdef AgentEvaluator < handle
    % 描述:
    % 一个用于严谨评估和对比 RL 与 VFI 策略的类 (重构版)。
    % - 直接使用 CoccoEnv 实例进行所有模拟。
    % - 通过重置随机数生成器种子，确保 RL 和 VFI 在相同的冲击序列下评估。
    
    properties
        % --- 输入 ---
        agent
        stats_path
        vfi_path
        results_path
        
        % --- 模拟参数 ---
        nsim = 1000
        env % normalizedcoccoenv包括的env
        raw_env % 原始的env
        
        % --- 加载的数据 ---
        obs_rms_data
        vfi_policy_interp
        
        % --- 结果 ---
        rl_results
        vfi_results_sim
        
        % --- 状态标志 ---
        is_rl_eval (1,1) logical
        is_vfi_eval (1,1) logical
        is_normalizedEnv (1,1) logical
    end
    
    methods
        % --- 构造函数 ---
        function this = AgentEvaluator(agent, stats_path, vfi_path, results_path)
            this.agent = agent;
            this.agent.UseExplorationPolicy = 0; % 使用确定性策略
            this.stats_path = stats_path;
            this.vfi_path = vfi_path;
            this.results_path = results_path;
            
            this.is_rl_eval = ~isempty(this.agent);
            this.is_vfi_eval = ~isempty(this.vfi_path) && exist(this.vfi_path, 'file');
            
            if ~this.is_rl_eval && ~this.is_vfi_eval
                error('AgentEvaluator: 无效输入', '必须提供有效的 agent 或 vfi_path。');
            end
            
            % 创建一个用于评估的环境实例
            
          
            % 预加载所需数据
            if this.is_rl_eval 
                this.load_rl_eval_env(); 
            end
            if this.is_vfi_eval, this.load_vfi_policies(); end
        end
        
        % --- 主运行函数 ---
        function run(this)
            fprintf('开始评估流程...\n');
            if this.is_rl_eval
                this.run_rl_simulation();
            end
            if this.is_vfi_eval
                this.run_vfi_simulation();
            end
            this.plot_results();
        end

        % --- RL 智能体模拟方法 ---
% AgentEvaluator.m

        % --- RL 智能体模拟方法 ---
        function run_rl_simulation(this)
            fprintf('  [RL] 正在进行模拟...\n');
            tn = this.raw_env.tn;
            simC = zeros(tn, this.nsim); simY = zeros(tn, this.nsim); 
            simA = zeros(tn, this.nsim); simW = zeros(tn, this.nsim);
            all_utilities = zeros(1, this.nsim);
            
            rng(123, 'twister');

 % --- 关键修正: 从智能体中提取行动者(actor) ---


            for i = 1:this.nsim
                [normObs, logged_reset] = this.env.reset();
                
                % 初始化 t=1 的财富和收入
                W_current = logged_reset.AbsoluteWealth;
                Y_current = logged_reset.AbsoluteIncome;

                for t = 1:tn
                    simW(t,i) = W_current;
                    simY(t,i) = Y_current;

                    action = getAction(this.agent, {normObs});
                    [normObs, ~, is_done, logged_step] = this.env.step(action{1});
                    
                    simC(t, i) = logged_step.AbsoluteConsumption;
                    simC_prop(t,i) = logged_step.C_prop;
                    simA(t, i) = logged_step.RiskyShare;
                    
                    % 从 logged_step 更新下一期的 W 和 Y
                    W_current = logged_step.AbsoluteWealth_next;
                    Y_current = logged_step.AbsoluteIncome_next;

                    if is_done
                        [simC, simY, simA, simW] = this.handle_death(t, i, simC, simY, simA, simW);
                        break;
                    end
                end
                all_utilities(i) = this.calculate_lifetime_utility(simC(:, i));
            end
            this.rl_results.W = nanmean(simW, 2); this.rl_results.C = nanmean(simC, 2);this.rl_results.C_prop = nanmean(simC_prop, 2);
            this.rl_results.Y = nanmean(simY, 2); this.rl_results.A = nanmean(simA, 2);
            this.rl_results.MeanUtility = mean(all_utilities);
            fprintf('  [RL] 模拟完成。平均终身效用: %.4f\n', this.rl_results.MeanUtility);

            plot(this.rl_results.C_prop)
            hold on
            plot(this.rl_results.A)
            legend({'c-prop','alpha'})
        end        
        % --- VFI 策略模拟方法 ---
        function run_vfi_simulation(this)
            fprintf('  [VFI] 正在进行模拟...\n');
            tn = this.env.tn;
            simC = zeros(tn, this.nsim); simY = zeros(tn, this.nsim); 
            simA = zeros(tn, this.nsim); simW = zeros(tn, this.nsim);
            all_utilities = zeros(1, this.nsim);

            % **核心**: 重置为与RL模拟完全相同的随机种子
            rng(123, 'twister');

            for i = 1:this.nsim
                this.env.reset();
                for t = 1:tn
                    simW(t,i) = this.env.W;

                    % a. 获取当前状态以计算收入和归一化变量
                    P = this.env.P;
                    P_retire = this.env.P_retirement;
                    age = this.env.age;
                    det_Y = this.env.f_y(age - this.env.tb + 1);

                    % b. 使用VFI策略计算动作
                    if age < this.env.tr
                        P_total = det_Y * P;
                        Y_gross = P_total * exp(this.env.shock_grid(randsample(this.env.n_shocks,1,true,this.env.shock_weig))*this.env.smay); % 临时冲击只影响当期收入，不影响状态
                        Y = Y_gross * (1-this.env.tau_y);
                        X = simW(t,i) + Y;
                        normalized_cash = X / P_total;
                        c_norm = this.vfi_policy_interp.c{t}(normalized_cash);
                        alpha = this.vfi_policy_interp.a{t}(normalized_cash);
                        C = min(c_norm * P_total, X * 0.9999);
                    else
                        Y = this.env.ret_fac * P_retire;
                        X = simW(t,i) + Y;
                        C = min(this.vfi_policy_interp.c{t}(X), X * 0.9999);
                        alpha = this.vfi_policy_interp.a{t}(X);
                    end
                    if C<1e-6, C=1e-6; end

                    % c. 执行 "dummy" step 来驱动环境的随机过程
                    % 我们只关心环境内部状态的转移 (P, age)
                    % 这个step的动作和奖励将被我们手动计算的值覆盖
                    [~, ~, is_done, ~] = this.env.step([C/X, alpha]); % 传入计算出的动作比例
                    
                    simC(t, i) = C;
                    simY(t, i) = Y;
                    simA(t, i) = alpha;

                    if is_done
                        [simC, simY, simA, simW] = this.handle_death(t, i, simC, simY, simA, simW);
                        break;
                    end
                end
                 all_utilities(i) = this.calculate_lifetime_utility(simC(:, i));
            end
            this.vfi_results_sim.W = nanmean(simW, 2); this.vfi_results_sim.C = nanmean(simC, 2);
            this.vfi_results_sim.Y = nanmean(simY, 2); this.vfi_results_sim.A = nanmean(simA, 2);
            this.vfi_results_sim.MeanUtility = mean(all_utilities);
            fprintf('  [VFI] 模拟完成。平均终身效用: %.4f\n', this.vfi_results_sim.MeanUtility);
        end
        
        % --- 绘图方法 (与之前版本相同) ---
        function plot_results(this)
            fprintf('  正在生成图表...\n');
            figure('Position', [100, 100, 1200, 900]);
            
            subplot(2,1,1); hold on;
            if this.is_rl_eval
                plot(this.raw_env.tb:this.raw_env.td, this.rl_results.W, 'b--', 'LineWidth', 2, 'DisplayName', 'Wealth (RL)');
                plot(this.raw_env.tb:this.raw_env.td, this.rl_results.C, 'r--', 'LineWidth', 2, 'DisplayName', 'Consumption (RL)');
            end
            if this.is_vfi_eval
                plot(this.raw_env.tb:this.raw_env.td, this.vfi_results_sim.W, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Wealth (VFI)');
                plot(this.raw_env.tb:this.raw_env.td, this.vfi_results_sim.C, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Consumption (VFI)');
            end
            if this.is_rl_eval || this.is_vfi_eval
                income_to_plot = this.is_rl_eval*this.rl_results.Y + (1-this.is_rl_eval)*this.vfi_results_sim.Y;
                plot(this.raw_env.tb:this.raw_env.td, income_to_plot, 'k:', 'LineWidth', 2, 'DisplayName', 'Income');
                xline(this.raw_env.tr-1, '--', 'Retirement', 'LineWidth', 1.5);
            end
            hold off; title('Mean Lifecycle Profiles'); xlabel('Age'); ylabel('Value (Absolute Units)');
            legend('show', 'Location', 'northwest'); grid on; xlim([this.raw_env.tb, this.raw_env.td]);

            subplot(2,1,2); hold on;
            if this.is_rl_eval, plot(this.raw_env.tb:this.raw_env.td, this.rl_results.A, 'm--', 'LineWidth', 2, 'DisplayName', 'Risky Share (RL)'); end
            if this.is_vfi_eval, plot(this.raw_env.tb:this.raw_env.td, this.vfi_results_sim.A, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Risky Share (VFI)'); end
            if this.is_rl_eval || this.is_vfi_eval, xline(this.raw_env.tr-1, '--', 'Retirement', 'LineWidth', 1.5); end
            hold off; title('Mean Portfolio Choice'); xlabel('Age'); ylabel('Share');
            legend('show', 'Location', 'best'); grid on; xlim([this.raw_env.tb, this.raw_env.td]); ylim([-0.05, 1.05]);

            if this.is_rl_eval && this.is_vfi_eval, base_title = 'RL vs. VFI Solution Comparison (Identical Shocks)'; file_name = 'rl_vs_vfi_comparison.png';
            elseif this.is_rl_eval, base_title = 'RL Agent Evaluation'; file_name = 'rl_only_evaluation.png';
            else, base_title = 'VFI Policy Simulation'; file_name = 'vfi_only_simulation.png'; end
            
            utility_title = '';
            if this.is_rl_eval, utility_title = [utility_title, sprintf('RL Utility: %.4f', this.rl_results.MeanUtility)]; end
            if this.is_vfi_eval
                if ~isempty(utility_title), utility_title = [utility_title, '  |  ']; end
                utility_title = [utility_title, sprintf('VFI Utility: %.4f', this.vfi_results_sim.MeanUtility)];
            end
            sgtitle({base_title; utility_title}, 'FontSize', 14);
            
            output_filename = fullfile(this.results_path, file_name);
            saveas(gcf, output_filename);
            fprintf('评估图表已保存到: %s\n', output_filename);
        end
    end
    
    methods (Access = private)
        % --- 辅助方法 ---
        function load_rl_eval_env(this)
            if this.is_normalizedEnv
            a = load(this.stats_path); 
            this.env = a.eval_env; 
                this.raw_env = a.env;
            else
                this.raw_env = CoccoEnv('disable_random_death', false);
                this.env = CoccoEnv('disable_random_death', false);
            end
        end
        function load_vfi_policies(this)
            load(this.vfi_path, 'vfi_results', 'cS');
            tn = this.env.tn;
            this.vfi_policy_interp.c = cell(tn, 1); this.vfi_policy_interp.a = cell(tn, 1);
            for t = 1:tn, this.vfi_policy_interp.c{t} = griddedInterpolant(cS.gcash, vfi_results.C_policy(:,t), 'linear', 'nearest'); this.vfi_policy_interp.a{t} = griddedInterpolant(cS.gcash, vfi_results.A_policy(:,t), 'linear', 'nearest'); end
        end
        % function normObs = normalize_obs(this, rawObs),
        %     normObs = (rawObs - this.obs_rms_data.Mean) ./ sqrt(this.obs_rms_data.Variance + 1e-8); 
        %     normObs(3) = rawObs(3); 
        % end
        function lifetime_utility = calculate_lifetime_utility(this, consumption_path)
            utility_path = (consumption_path.^(1 - this.raw_env.gamma)) / (1 - this.raw_env.gamma);
            if this.raw_env.gamma == 1, 
                utility_path = log(consumption_path); 
            end
            discount_factors = (this.raw_env.beta .^ (0:this.raw_env.tn-1))';
            lifetime_utility = nansum(utility_path .* discount_factors);
        end
        function [C, Y, A, W] = handle_death(~, t, i, C, Y, A, W), C(t+1:end, i) = NaN; Y(t+1:end, i) = NaN; A(t+1:end, i) = NaN; W(t+1:end, i) = NaN; end
    end
end