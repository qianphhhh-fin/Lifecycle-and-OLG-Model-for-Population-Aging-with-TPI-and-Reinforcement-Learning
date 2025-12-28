%% main_rl_cocco.m (修正版 V4 - 带周期性评估)
% 描述:
% 使用MATLAB强化学习工具箱 (SAC算法) 求解Cocco(2005)生命周期模型。
% **V4修正**: 修正了评估时获取确定性动作的函数调用。
%            错误: getActionForExploitation (不存在)
%            正确: 使用 getActor() 和 evaluate() 来获取动作均值。

% --- 初始化 ---
clear;
close all;
clc;
rng(42); % 设置随机种子以保证可复现性

fprintf('开始执行 Cocco 生命周期模型的 RL 求解 (带周期性评估)...\n');

% --- A. 路径和参数设置 ---
results_path = 'results_rl';
if ~exist(results_path, 'dir'), mkdir(results_path); end

VFI_MAT_FILE = 'result/vfi_results_crra_nopps_solo.mat';
if ~exist(VFI_MAT_FILE, 'file')
    warning('VFI结果文件 "%s" 未找到。最终评估将无法对比VFI。', VFI_MAT_FILE);
    vfi_path_for_eval = [];
else
    vfi_path_for_eval = VFI_MAT_FILE;
end

% 训练参数
TOTAL_EPISODES_PHASE1 = 3000;
EVAL_FREQUENCY_EPISODES = 200; % 每隔多少回合评估一次
NUM_EVAL_SIMULATIONS = 400; % 每次评估使用的模拟次数
AGENT_SAMPLE_TIME = 1.0;
DISCOUNT_FACTOR = 0.95;

% 文件路径
best_agent_phase1_path = fullfile(results_path, 'best_agent_phase1.mat');
best_stats_phase1_path = fullfile(results_path, 'best_stats_phase1.mat');
final_agent_path = fullfile(results_path, 'final_agent.mat');


% --- B. 阶段一: 预训练 (无随机死亡) 并周期性评估 ---
fprintf('\n===== 阶段一: 预训练 (带周期性评估) =====\n');

% 1. 创建简化环境 (禁用随机死亡)
env_p1_base = CoccoEnv('disable_random_death', true);
env_p1 = NormalizedCoccoEnv(env_p1_base, 'norm_obs_idx', [1, 2]);

% 2. 创建SAC智能体
obsInfo = getObservationInfo(env_p1);
actInfo = getActionInfo(env_p1);
agent_p1 = createAgentCocco(obsInfo, actInfo, DISCOUNT_FACTOR, AGENT_SAMPLE_TIME);

% 3. 初始化评估和绘图
best_mean_utility = -inf;
eval_episodes = 0:EVAL_FREQUENCY_EPISODES:TOTAL_EPISODES_PHASE1;
utility_history = NaN(1, length(eval_episodes));
h_fig = figure('Name', 'Live Evaluation During Training', 'Position', [100, 100, 1000, 800]);
ax1 = subplot(2,1,1);
ax2 = subplot(2,1,2);
sgtitle(h_fig, 'Training Phase 1: Evaluating Agent Policy...', 'FontSize', 14);

% 4. 手动训练循环
for i = 1:TOTAL_EPISODES_PHASE1
    
    % --- 运行一个回合以采集经验 (经验自动存入智能体内部) ---
    sim(agent_p1, env_p1, rlSimulationOptions('MaxSteps', env_p1_base.tn));
    
    % --- 从经验缓冲区中采样并训练 ---
    if agent_p1.ExperienceBuffer.Length > agent_p1.AgentOptions.MiniBatchSize
        agent_p1 = train(agent_p1);
    end
    
    % --- 检查是否到达评估点 ---
    if mod(i, EVAL_FREQUENCY_EPISODES) == 0 || i == TOTAL_EPISODES_PHASE1
        
        fprintf('\n--- 评估点: 回合 %d / %d ---\n', i, TOTAL_EPISODES_PHASE1);
        
        % a. 准备评估
        env_p1.is_training = false;
        temp_agent_for_eval = copy(agent_p1);
        actor_for_eval = getActor(temp_agent_for_eval); % 获取actor用于评估
        
        % b. 运行模拟评估
        all_utilities = zeros(1, NUM_EVAL_SIMULATIONS);
        simC = zeros(env_p1_base.tn, NUM_EVAL_SIMULATIONS);
        simA = zeros(env_p1_base.tn, NUM_EVAL_SIMULATIONS);
        
        for sim_idx = 1:NUM_EVAL_SIMULATIONS
            rawObs = env_p1_base.reset();
            is_done = false;
            t = 1;
            while ~is_done && t <= env_p1_base.tn
                normObs = env_p1.normalize_obs(rawObs);
                
                % --- MODIFICATION START ---
                % 使用 actor 的 evaluate 方法获取确定性动作 (分布的均值)
                [actionMean, ~] = evaluate(actor_for_eval, {normObs});
                action = actionMean{1}; % 提取动作
                % --- MODIFICATION END ---
                
                [rawObs, ~, is_done, loggedSignals] = env_p1_base.step(action);
                simC(t, sim_idx) = loggedSignals.C_prop;
                simA(t, sim_idx) = loggedSignals.RiskyShare;
                t = t + 1;
            end
            
            % 计算该次模拟的总效用
            consumption_path = simC(:, sim_idx);
            utility_path = (consumption_path.^(1 - env_p1_base.gamma)) / (1 - env_p1_base.gamma);
            if env_p1_base.gamma == 1, utility_path = log(consumption_path); end
            discount_factors = (env_p1_base.beta .^ (0:env_p1_base.tn-1))';
            all_utilities(sim_idx) = nansum(utility_path .* discount_factors);
        end
        
        mean_utility_current = mean(all_utilities);
        
        % c. 更新历史记录和动态图表
        eval_idx = find(eval_episodes == i);
        if ~isempty(eval_idx), utility_history(eval_idx) = mean_utility_current; end
        
        meanC = nanmean(simC, 2);
        meanA = nanmean(simA, 2);
        
        plot(ax1, env_p1_base.tb:env_p1_base.td, meanC, 'LineWidth', 1.5);
        title(ax1, 'Mean Consumption Over Lifecycle');
        xlabel(ax1, 'Age'); ylabel(ax1, 'Consumption'); grid(ax1, 'on'); xlim(ax1, [env_p1_base.tb, env_p1_base.td]);
        
        plot(ax2, env_p1_base.tb:env_p1_base.td, meanA, 'LineWidth', 1.5);
        title(ax2, 'Mean Risky Share (Alpha) Over Lifecycle');
        xlabel(ax2, 'Age'); ylabel(ax2, 'Alpha'); grid(ax2, 'on'); xlim(ax2, [env_p1_base.tb, env_p1_base.td]); ylim(ax2, [-0.05, 1.05]);
        
        sgtitle(h_fig, sprintf('Phase 1 Eval @ Episode %d | Mean Utility: %.4f', i, mean_utility_current), 'FontSize', 14);
        drawnow;
        
        % d. 检查并保存最优模型
        if mean_utility_current > best_mean_utility
            best_mean_utility = mean_utility_current;
            fprintf('>>> 新的最优模型! 平均效用: %.4f. 正在保存...\n', best_mean_utility);
            
            best_agent_p1 = agent_p1;
            save(best_agent_phase1_path, 'best_agent_p1');
            env_p1.saveStats(best_stats_phase1_path);
        end
        
        % e. 恢复训练模式
        env_p1.is_training = true;
    end
end

fprintf('\n阶段一训练和评估完成。最优智能体已保存到 %s\n', best_agent_phase1_path);
fprintf('最优归一化统计数据已保存到 %s\n', best_stats_phase1_path);
saveas(h_fig, fullfile(results_path, 'phase1_live_evaluation_final.png'));

% --- 加载最优模型以进行后续步骤 ---
try
    load(best_agent_phase1_path, 'best_agent_p1');
    agentToLoad = best_agent_p1;
    fprintf('已成功加载阶段一最优智能体。\n');
catch
    warning('无法加载最优智能体，将使用训练结束时的最后一个智能体。');
    agentToLoad = agent_p1;
end


% --- C. 阶段二: 微调 (有随机死亡) ---
fprintf('\n===== 阶段二: 微调 (有随机死亡) =====\n');
% 1. 创建完整环境
env_p2_base = CoccoEnv('disable_random_death', false);

% 2. 创建归一化封装器并 **加载** 阶段一的统计数据
env_p2 = NormalizedCoccoEnv(env_p2_base, 'norm_obs_idx', [1, 2]);
env_p2.loadStats(best_stats_phase1_path);
fprintf('已成功加载阶段一最优的归一化统计数据。\n');

% 3. 加载阶段一训练好的智能体
agent_p2 = agentToLoad;

% 4. 设置微调的训练选项
trainOpts_p2 = rlTrainingOptions(...
    'MaxEpisodes', 5000, ...
    'MaxStepsPerEpisode', env_p2_base.tn, ...
    'ScoreAveragingWindowLength', 50, ...
    'Verbose', true, ...
    'Plots', 'training-progress', ...
    'StopTrainingCriteria', 'AverageReward', ...
    'StopTrainingValue', 1000, ...
    'SaveAgentCriteria', 'EpisodeFrequency',...
    'SaveAgentValue', 200, ...
    'SaveAgentDirectory', fullfile(results_path, 'p2_checkpoints'));

% 5. 继续训练 (微调)
if trainOpts_p2.MaxEpisodes > 0
    trainingStats_p2 = train(agent_p2, env_p2, trainOpts_p2);
    final_agent = trainingStats_p2.Agent;
else
    final_agent = agent_p2;
end
    
% 6. 保存最终的智能体
save(final_agent_path, 'final_agent');
fprintf('阶段二微调完成。最终智能体已保存到 %s\n', final_agent_path);


% --- D. 最终评估 ---
fprintf('\n===== 最终评估 =====\n');
fprintf('正在使用最终微调过的智能体进行评估...\n');
evaluator = AgentEvaluator(final_agent, best_stats_phase1_path, vfi_path_for_eval, results_path);
evaluator.nsim = 1000;
evaluator.run();


fprintf('\n所有流程执行完毕。\n');