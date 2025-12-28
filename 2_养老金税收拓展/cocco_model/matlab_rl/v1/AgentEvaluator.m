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
        nsim = 50
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
% AgentEvaluator.m

        function run_rl_simulation(this)
            fprintf('  [RL] 正在进行模拟 (与训练时评估对齐)...\n');
            tn = this.raw_env.tn;
            simC = zeros(tn, this.nsim); simY = zeros(tn, this.nsim);
            simA = zeros(tn, this.nsim); simW = zeros(tn, this.nsim);
            simC_prop = zeros(tn, this.nsim); % 增加 C_prop 记录
            all_utilities = zeros(1, this.nsim);
            
            % 使用确定性策略
            this.agent.UseExplorationPolicy = false;
            
            % **核心**: 重置为与VFI模拟完全相同的随机种子
            rng(123, 'twister');

            for i = 1:this.nsim
                % **核心修改**: 模仿 trainAgentWithLivePlots 的评估循环
                % 1. 从底层原始环境 reset，获取真实的初始状态
                [rawObs, logged_reset] = this.raw_env.reset();
                
                % 2. 手动归一化初始观测值
                normObs = this.env.normalize_obs(rawObs);

                % 初始化 t=1 的财富和收入
                W_current = logged_reset.AbsoluteWealth;
                Y_current = logged_reset.AbsoluteIncome;

                for t = 1:tn
                    simW(t,i) = W_current;
                    simY(t,i) = Y_current;

                    % 3. 智能体使用归一化的观测值进行决策
                    action = getAction(this.agent, {normObs});
                    
                    % 4. 在底层原始环境中执行动作，以获取真实的 logged_step
                    [rawNextObs, ~, is_done, logged_step] = this.raw_env.step(action{1});
                    
                    % 5. 手动归一化下一个观测值，为循环的下一步做准备
                    normObs = this.env.normalize_obs(rawNextObs);
                    
                    % 从 logged_step 记录结果
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
            
            % 汇总结果
            this.rl_results.W = nanmean(simW, 2); 
            this.rl_results.C = nanmean(simC, 2);
            this.rl_results.C_prop = nanmean(simC_prop, 2); % 保存 C_prop 的均值
            this.rl_results.Y = nanmean(simY, 2); 
            this.rl_results.A = nanmean(simA, 2);
            this.rl_results.MeanUtility = mean(all_utilities);
            
            fprintf('  [RL] 模拟完成。平均终身效用: %.4f\n', this.rl_results.MeanUtility);

            % 绘图部分可以暂时移到主 plot_results 函数中
            % figure; % 创建新图窗以避免覆盖
            % plot(this.raw_env.tb:this.raw_env.td, this.rl_results.C_prop);
            % hold on;
            % plot(this.raw_env.tb:this.raw_env.td, this.rl_results.A);
            % legend({'c-prop','alpha'});
            % title('RL Agent Policy Profile');
            % xlabel('Age');
            % ylabel('Action Value');
        end       
        
        % --- VFI 策略模拟方法 ---

        % --- VFI 策略模拟方法 ---
        function run_vfi_simulation(this)
            fprintf('  [VFI] 正在进行模拟 (与RL环境对齐)...\n');
            tn = this.raw_env.tn;
            simC = zeros(tn, this.nsim); simY = zeros(tn, this.nsim);
            simA = zeros(tn, this.nsim); simW = zeros(tn, this.nsim);
            simC_prop = zeros(tn, this.nsim);
            all_utilities = zeros(1, this.nsim);

            % **核心**: 重置为与RL模拟完全相同的随机种子
            rng(123, 'twister');

            for i = 1:this.nsim
                % **核心修改**: VFI 也使用 env.reset() 来初始化
                [~, logged_reset] = this.env.reset(); % NormObs 被忽略

                % 初始化 t=1 的财富和收入
                
                W_current = logged_reset.AbsoluteWealth;
                Y_current = logged_reset.AbsoluteIncome;
                X_current = W_current + Y_current;
                for t = 1:tn
                    simW(t,i) = W_current;
                    simY(t,i) = Y_current;

                    % a. 使用VFI策略插值计算当前状态下的最优动作
                    % 注意：VFI的策略函数是基于绝对现金持有量(X_t)
                    C_vfi = this.vfi_policy_interp.c{t}(X_current);
                    alpha_vfi = this.vfi_policy_interp.a{t}(X_current);
                    
                    % 将绝对消费转换为消费比例 c_prop，以作为 env.step 的输入
                    c_prop_vfi = C_vfi / X_current;
                    
                    % 确保动作在 [0, 1] 区间内
                    c_prop_vfi = max(0, min(1, c_prop_vfi));
                    alpha_vfi  = max(0, min(1, alpha_vfi));

                    % b. 将计算出的VFI动作喂给环境的step函数
                    % 我们不关心step返回的obs和reward，只关心它如何驱动状态演化
                    % 并记录 logged_step 中的真实数值
                    [~, ~, is_done, logged_step] = this.env.step([c_prop_vfi, alpha_vfi]);
                    
                    % c. 从 logged_step 中记录结果 (这保证了与RL的绝对一致性)
                    simC(t, i) = logged_step.AbsoluteConsumption;
                    simC_prop(t,i) = logged_step.C_prop;
                    simA(t, i) = logged_step.RiskyShare;
                    
                    % d. 从 logged_step 更新下一期的 W, Y 和 X
                    Y_current = logged_step.AbsoluteIncome_next;
                    W_current = logged_step.AbsoluteWealth_next;
                    X_current = Y_current + W_current;

                    if is_done
                        [simC, simY, simA, simW] = this.handle_death(t, i, simC, simY, simA, simW);
                        break;
                    end
                end
                 all_utilities(i) = this.calculate_lifetime_utility(simC(:, i));
            end
            
            % 汇总并保存结果
            this.vfi_results_sim.W = nanmean(simW, 2); 
            this.vfi_results_sim.C = nanmean(simC, 2);
            this.vfi_results_sim.C_prop = nanmean(simC_prop, 2);
            this.vfi_results_sim.Y = nanmean(simY, 2); 
            this.vfi_results_sim.A = nanmean(simA, 2);
            this.vfi_results_sim.MeanUtility = mean(all_utilities);
            
            fprintf('  [VFI] 模拟完成。平均终身效用: %.4f\n', this.vfi_results_sim.MeanUtility);
        end        
       
        
        % --- 绘图方法 (与之前版本相同) ---
% AgentEvaluator.m

        % --- 绘图方法 (重构版) ---
        function plot_results(this)
            fprintf('  正在生成详细的比较图表...\n');
            
            % --- 1. 初始化 ---
            fig = figure('Position', [100, 100, 1400, 1000], 'Name', 'RL vs. VFI Detailed Comparison');
            ages = this.raw_env.tb:this.raw_env.td;
            retirement_age_line = this.raw_env.tr - 1;

            % --- 3. 子图 2: 消费策略 (c_prop) ---
            ax2 = subplot(2, 2, 1);
            hold(ax2, 'on');
            if this.is_rl_eval && isfield(this.rl_results, 'C_prop')
                plot(ax2, ages, this.rl_results.C_prop, 'r--', 'LineWidth', 2, 'DisplayName', 'c-prop (RL)');
            end
            % VFI的c_prop需要根据其模拟结果计算: C / (W+Y)
            if this.is_vfi_eval
                plot(ax2, ages, this.vfi_results_sim.C_prop, 'r-', 'LineWidth', 1.5, 'DisplayName', 'c-prop (VFI)');
            end
            if this.is_rl_eval || this.is_vfi_eval
                xline(ax2, retirement_age_line, '--', 'Retirement', 'LineWidth', 1.5);
            end
            hold(ax2, 'off');
            title(ax2, 'c-prop');
            xlabel(ax2, 'Age');
            ylabel(ax2, 'Proportion of Cash-on-Hand');
            legend(ax2, 'show', 'Location', 'best');
            grid(ax2, 'on');
            xlim(ax2, [ages(1), ages(end)]);
            ylim(ax2, [-0.05, 1.05]);
            
            % --- 4. 子图 3: 投资策略 (alpha) ---
            ax3 = subplot(2, 2, 2);
            hold(ax3, 'on');
            if this.is_rl_eval
                plot(ax3, ages, this.rl_results.A, 'm--', 'LineWidth', 2, 'DisplayName', 'Risky Share (RL)');
            end
            if this.is_vfi_eval
                plot(ax3, ages, this.vfi_results_sim.A, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Risky Share (VFI)');
            end
            if this.is_rl_eval || this.is_vfi_eval
                xline(ax3, retirement_age_line, '--', 'Retirement', 'LineWidth', 1.5);
            end
            hold(ax3, 'off');
            title(ax3, 'Investment Policy (\alpha)');
            xlabel(ax3, 'Age');
            ylabel(ax3, 'Share in Risky Asset');
            legend(ax3, 'show', 'Location', 'best');
            grid(ax3, 'on');
            xlim(ax3, [ages(1), ages(end)]);
            ylim(ax3, [-0.05, 1.05]);

                        % --- 2. 子图 3: 宏观变量 (财富) ---
            ax1 = subplot(2, 2, 3);
            hold(ax1, 'on');
            if this.is_rl_eval
                plot(ax1, ages, this.rl_results.W, 'b--', 'LineWidth', 2, 'DisplayName', 'Wealth (RL)');
            end
            if this.is_vfi_eval
                plot(ax1, ages, this.vfi_results_sim.W, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Wealth (VFI)');
            end
            % if this.is_rl_eval || this.is_vfi_eval
            %     % 绘制收入曲线 (两者应几乎一致)
            %     income_to_plot = this.rl_results.Y; % 默认使用RL的，因为两者应相同
            %     if isempty(income_to_plot) && isfield(this.vfi_results_sim, 'Y')
            %         income_to_plot = this.vfi_results_sim.Y;
            %     end
            %     plot(ax1, ages, income_to_plot, 'k:', 'LineWidth', 2, 'DisplayName', 'Income');
            %     xline(ax1, retirement_age_line, '--', 'Retirement', 'LineWidth', 1.5);
            % end
            hold(ax1, 'off');
            title(ax1, 'Mean Lifecycle Profiles (Absolute Units)');
            xlabel(ax1, 'Age');
            ylabel(ax1, 'Value');
            legend(ax1, 'show', 'Location', 'northwest');
            grid(ax1, 'on');
            xlim(ax1, [ages(1), ages(end)]);

            
            % --- 5. 设置总标题和保存 ---
            if this.is_rl_eval && this.is_vfi_eval, base_title = 'RL vs. VFI Solution Comparison (Identical Shocks)'; file_name = 'rl_vs_vfi_comparison_detailed.png';
            elseif this.is_rl_eval, base_title = 'RL Agent Evaluation'; file_name = 'rl_only_evaluation_detailed.png';
            else, base_title = 'VFI Policy Simulation'; file_name = 'vfi_only_simulation_detailed.png'; end
            
            utility_title = '';
            if this.is_rl_eval, utility_title = [utility_title, sprintf('RL Utility: %.4f', this.rl_results.MeanUtility)]; end
            if this.is_vfi_eval
                if ~isempty(utility_title), utility_title = [utility_title, '  |  ']; end
                utility_title = [utility_title, sprintf('VFI Utility: %.4f', this.vfi_results_sim.MeanUtility)];
            end
            sgtitle({base_title; utility_title}, 'FontSize', 14);
            
            output_filename = fullfile(this.results_path, file_name);
            saveas(fig, output_filename);
            fprintf('评估图表已保存到: %s\n', output_filename);
        end
    
    
    end
    
    methods (Access = private)
        % --- 辅助方法 ---
% AgentEvaluator.m

        function load_rl_eval_env(this)
            % 描述:
            % 直接加载训练阶段结束时保存的、包含了完整归一化统计信息的
            % NormalizedCoccoEnv 对象。
            % **核心修改**: 不再分别加载 agent 和 stats，而是加载包含 env 的文件。
            %              我们假设 STATS_PATH 指向的是由 env.saveStats 保存的文件。
            
            fprintf('  正在加载用于评估的归一化环境从: %s\n', this.stats_path);
            
            % STATS_PATH 文件中保存了一个名为 'eval_env' 的变量，
            % 这个变量就是包含了统计数据的 NormalizedCoccoEnv 对象。
            load_data = load(this.stats_path, 'eval_env');
            this.env = load_data.eval_env;
            this.raw_env = this.env.env; % 从封装器中获取底层的 CoccoEnv
            
            % 确保在评估时，环境不会再更新其统计数据
            this.env.is_training = false;
            
            fprintf('  评估环境加载完成。\n');
        end
% AgentEvaluator.m

        function load_vfi_policies(this)
            % 描述:
            % 加载与 RL 环境匹配的 VFI 基准策略。
            % **核心修改**: 从 vfi_benchmark_results.mat 文件中加载策略，
            % 并使用 cS.x_grid (绝对现金网格) 来创建插值对象。
            
            fprintf('  正在加载 VFI 基准策略从: %s\n', this.vfi_path);
            
            % 加载包含 vfi_results 和 cS 的 .mat 文件
            data = load(this.vfi_path, 'vfi_results', 'cS');
            vfi_results = data.vfi_results;
            cS = data.cS;
            
            tn = this.raw_env.tn;
            
            % 确保加载的 VFI 结果与环境的时间维度匹配
            if size(vfi_results.C_policy, 2) ~= tn
                error('AgentEvaluator: VFI维度不匹配', ...
                    '加载的 VFI 策略的时间维度 (%d) 与环境的时间维度 (%d) 不符。', ...
                    size(vfi_results.C_policy, 2), tn);
            end
            
            this.vfi_policy_interp.c = cell(tn, 1);
            this.vfi_policy_interp.a = cell(tn, 1);
            
            % 使用绝对现金网格 cS.x_grid 创建插值
            for t = 1:tn
                % 创建消费策略的插值对象
                this.vfi_policy_interp.c{t} = griddedInterpolant(cS.x_grid, vfi_results.C_policy(:,t), 'spline', 'linear');
                
                % 创建风险资产配置策略的插值对象
                this.vfi_policy_interp.a{t} = griddedInterpolant(cS.x_grid, vfi_results.A_policy(:,t), 'spline', 'linear');
            end
            
            fprintf('  VFI 基准策略加载并插值完成。\n');
        end
        
        
        % function normObs = normalize_obs(this, rawObs),
        %     normObs = (rawObs - this.obs_rms_data.Mean) ./ sqrt(this.obs_rms_data.Variance + 1e-8); 
        %     normObs(3) = rawObs(3); 
        % end
        function lifetime_utility = calculate_lifetime_utility(this, consumption_path)
            utility_path = (consumption_path.^(1 - this.raw_env.gamma)) / (1 - this.raw_env.gamma)+1;
            if this.raw_env.gamma == 1, 
                utility_path = log(consumption_path); 
            end
            discount_factors = (this.raw_env.beta .^ (0:this.raw_env.tn-1))';
            lifetime_utility = nansum(utility_path .* discount_factors);
        end

        function [C, Y, A, W] = handle_death(~, t, i, C, Y, A, W), 
            C(t+1:end, i) = NaN; 
            Y(t+1:end, i) = NaN; 
            A(t+1:end, i) = NaN; 
            W(t+1:end, i) = NaN; 
        end
    end
end