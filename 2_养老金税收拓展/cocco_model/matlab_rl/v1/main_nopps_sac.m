%% main_rl_cocco.m
% 描述:
% 使用MATLAB强化学习工具箱 (SAC算法) 求解Cocco(2005)生命周期模型。
% 采用两阶段课程学习策略，并实现了关键的观测/奖励归一化。
%
% 流程:
% 1. 设置参数。
% 2. 阶段一: 在无随机死亡的简化环境中进行预训练。
% 3. 阶段二: 加载预训练模型，在有随机死亡的完整环境中进行微调。
% 4. 评估最终模型，并与VFI结果进行对比。

% --- 初始化 ---
clear;
close all;
% clc;
rng(42); % 设置随机种子以保证可复现性

fprintf('开始执行 Cocco 生命周期模型的 RL 求解...\n');

% --- A. 路径和参数设置 ---
results_path = 'results_rl';
if ~exist(results_path, 'dir'), mkdir(results_path); end

% VFI_MAT_FILE = 'result/vfi_results_crra_nopps_solo.mat';
% if ~exist(VFI_MAT_FILE, 'file')
%     error('VFI结果文件 "%s" 未找到。请先运行 main_nopps.m。', VFI_MAT_FILE);
% end

% 训练参数
TIMESTEPS_PHASE1 = 3000;  % 阶段一训练步数
TIMESTEPS_PHASE2 = 1000;  % 阶段二训练步数
AGENT_SAMPLE_TIME = 1.0;
DISCOUNT_FACTOR = 0.95;    % 必须与环境中的 beta 匹配

% 文件路径
agent_phase1_path = fullfile(results_path, 'agent_phase1.mat');
stats_phase1_path = fullfile(results_path, 'stats_phase1.mat');
final_agent_path = fullfile(results_path, 'final_agent.mat');

% --- B. 阶段一: 预训练 (无随机死亡) ---
fprintf('\n===== 阶段一: 预训练 (无随机死亡) =====\n');

% 1. 创建简化环境 (禁用随机死亡)
env_p1_base = CoccoEnv('disable_random_death', true);

% 2. 使用归一化封装器
env_p1 = NormalizedCoccoEnv(env_p1_base, 'norm_obs_idx', [1]); % 只归一化财富和永久收入
% env_p1 = env_p1_base;

% 3. 创建SAC智能体
obsInfo = getObservationInfo(env_p1);
actInfo = getActionInfo(env_p1);
% agent_p1 = createAgentCocco(obsInfo, actInfo, DISCOUNT_FACTOR, AGENT_SAMPLE_TIME);
agent_p1 = createAgentCocco_sep(obsInfo, actInfo, DISCOUNT_FACTOR, AGENT_SAMPLE_TIME);

% 4. 设置训练选项
trainOpts_p1 = rlTrainingOptions(...
    'MaxEpisodes', 30000, ...
    'MaxStepsPerEpisode', env_p1_base.tn, ...
    'ScoreAveragingWindowLength', 50, ...
    'Verbose', true, ...
    'Plots', 'training-progress', ...
    'StopTrainingCriteria', 'AverageReward', ...
    'StopTrainingValue', 1000, ... % 设置一个较高的值，让训练由步数决定
    'SaveAgentCriteria', 'EpisodeFrequency',...
    'SaveAgentValue', 50, ...
    'SaveAgentDirectory', fullfile(results_path, 'p1_checkpoints'));

% if TIMESTEPS_PHASE1 > 0
%     % 5. 训练智能体
%     trainingStats_p1 = train(agent_p1, env_p1, trainOpts_p1);
% 
%     % 6. 保存智能体和归一化统计数据
%     save(agent_phase1_path, 'agent_p1');
%     env_p1.saveStats(stats_phase1_path);
% 
%     fprintf('阶段一训练完成。智能体已保存到 %s\n', agent_phase1_path);
%     fprintf('归一化统计数据已保存到 %s\n', stats_phase1_path);
% 
%     % 加载最优的智能体以进入下一阶段
%     agentToLoad = trained_agent_p1; 
% else
%     fprintf('跳过阶段一训练，尝试加载已有文件...\n');
%     load(agent_phase1_path, 'trained_agent_p1');
%     agentToLoad = trained_agent_p1;
% end

if TIMESTEPS_PHASE1 > 0
        % **新增**: 定义阶段一最佳模型的保存路径
    bestModelDir = fullfile(results_path, 'p1_best_model');
    bestAgentPath_p1 = fullfile(bestModelDir, 'agent_phase1_best.mat');
    bestStatsPath_p1 = fullfile(bestModelDir, 'stats_phase1_best.mat');

    % 5. 训练智能体 (使用我们的自定义函数)
    trainingStats_p1 = trainAgentWithLivePlots(agent_p1, env_p1, trainOpts_p1, ...
        bestAgentPath_p1, bestStatsPath_p1);

    % 6. 保存智能体和归一化统计数据
    % 注意：自定义循环不返回训练好的智能体，agent_p1 本身就是训练后的版本
    save(agent_phase1_path, 'agent_p1');
    env_p1.saveStats(stats_phase1_path);

    fprintf('阶段一训练完成。智能体已保存到 %s\n', agent_phase1_path);
    fprintf('归一化统计数据已保存到 %s\n', stats_phase1_path);
    aaa
    % 加载最优的智能体以进入下一阶段
    agentToLoad = agent_p1; % 直接使用训练后的 agent_p1
else
    fprintf('跳过阶段一训练，尝试加载已有文件...\n');
    load(agent_phase1_path, 'agent_p1'); % <-- 确保加载的变量名是 agent_p1
    agentToLoad = agent_p1;
end

aaa
% --- C. 阶段二: 微调 (有随机死亡) ---
fprintf('\n===== 阶段二: 微调 (有随机死亡) =====\n');

% 1. 创建完整环境
env_p2_base = CoccoEnv('disable_random_death', false);

% 2. 创建归一化封装器并 **加载** 阶段一的统计数据
env_p2 = NormalizedCoccoEnv(env_p2_base, 'norm_obs_idx', [1, 2]);
env_p2.loadStats(stats_phase1_path);
fprintf('已成功加载阶段一的归一化统计数据。\n');

% 3. 加载阶段一训练好的智能体
agent_p2 = agentToLoad;

% 4. 设置微调的训练选项 (通常使用更小的学习率，但这里为了简化保持一致)
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

if TIMESTEPS_PHASE2 > 0
    % 5. 继续训练 (微调)
    trainingStats_p2 = train(agent_p2, env_p2, trainOpts_p2);
    
    % 6. 保存最终的智能体
    final_agent = trainingStats_p2.Agent;
    save(final_agent_path, 'final_agent');
    fprintf('阶段二微调完成。最终智能体已保存到 %s\n', final_agent_path);
else
    fprintf('跳过阶段二训练。\n');
    final_agent = agentToLoad; % 如果不进行阶段二, 则阶段一的结果就是最终结果
    save(final_agent_path, 'final_agent');
end

% --- D. 最终评估 ---
fprintf('\n===== 最终评估 =====\n');
evaluate_agent(final_agent, stats_phase1_path, VFI_MAT_FILE, results_path);

fprintf('\n所有流程执行完毕。\n');


% main_nopps_sac.m (在文件末尾的函数)

% main_nopps_sac.m (在文件末尾的函数)

function trainingStats = trainAgentWithLivePlots(agent, env, trainOpts, bestAgentPath, bestStatsPath)
    % 描述:
    % 一个自定义的训练循环，用于在训练过程中实时可视化智能体的策略演变。
    % V5: 如果原始奖励有提升，则自动保存最佳模型和统计数据。
    %
    % 输入:
    %   agent - 要训练的 RL 智能体
    %   env - 强化学习环境
    %   trainOpts - 包含训练参数的 rlTrainingOptions 对象
    %   bestAgentPath - 最佳智能体模型的保存路径 (string)
    %   bestStatsPath - 最佳统计数据的保存路径 (string)
    %
    % 输出:
    %   trainingStats - 包含训练统计信息（如奖励）的结构体

    fprintf('开始使用自定义循环进行训练，并实时显示策略...\n');

    % --- 1. 从 trainOpts 中提取参数和新配置 ---
    maxEpisodes = trainOpts.MaxEpisodes;
    maxStepsPerEpisode = trainOpts.MaxStepsPerEpisode;
    scoreAvgWindow = trainOpts.ScoreAveragingWindowLength;
    plotFrequency = 50;
    rewardPlotWindow = 10000;

    % **新增**: 创建保存最佳模型的目录 (如果不存在)
    [bestModelDir, ~, ~] = fileparts(bestAgentPath);
    if ~exist(bestModelDir, 'dir')
        mkdir(bestModelDir);
        fprintf('创建目录: %s\n', bestModelDir);
    end

    % --- 2. 初始化绘图窗口 ---
    fig = figure('Position', [200, 200, 1200, 600], 'Name', 'RL Agent Training Monitor');
    rewardPlotAx = subplot(1,2,1);
    rewardPlotLine = plot(rewardPlotAx, NaN, NaN);
    title(rewardPlotAx, 'Episode Discounted Reward (Recent 500)');
    xlabel(rewardPlotAx, 'Episode'); ylabel(rewardPlotAx, 'Total Normalized Reward');
    grid(rewardPlotAx, 'on');
    policyPlotAx = subplot(1,2,2);
    hold(policyPlotAx, 'on');
    ages = env.env.tb : env.env.td;
    policyPlotC = plot(policyPlotAx, ages, nan(1, maxStepsPerEpisode), 'b-', 'LineWidth', 2, 'DisplayName', 'c-prop');
    policyPlotA = plot(policyPlotAx, ages, nan(1, maxStepsPerEpisode), 'r-', 'LineWidth', 2, 'DisplayName', 'alpha');
    title(policyPlotAx, 'Agent Policy (Deterministic)');
    xlabel(policyPlotAx, 'Age'); ylabel(policyPlotAx, 'Action Value');
    legend(policyPlotAx, 'show', 'Location', 'best');
    grid(policyPlotAx, 'on');
    ylim(policyPlotAx, [-0.1, 1.1]); xlim(policyPlotAx, [env.env.tb, env.env.td]);
    hold(policyPlotAx, 'off');

    % --- 3. 初始化训练变量 ---
    episodeRewards = zeros(1, maxEpisodes);
    % **核心修改**: 新增一个数组用于存储每个 episode 的原始折现奖励 (终身效用)
    episodeRawDiscountedRewards = zeros(1, maxEpisodes);
    best_eval_reward_so_far = -inf; 
    trainingStats.EpisodeReward = [];
    totalSteps = 0;
    tic; 

    % --- 4. 主训练循环 ---
    for iEpisode = 1:maxEpisodes
        
        obs = env.reset();
        isDone = false;
        currentEpisodeReward = 0;
        currentEpisodeRawDiscountedReward = 0; % **核心修改**: 初始化当前 episode 的原始折现奖励
        
        for iStep = 1:maxStepsPerEpisode
            actionCell = getAction(agent, {obs}); 
            % **核心修改**: step 函数现在只返回原始奖励
            [nextObs, rawReward, isDone, ~] = env.step(actionCell{1}); 
            
            experience.Observation = {obs};
            experience.Action = actionCell;
            experience.Reward = rawReward; % **核心修改**: 传给 agent.learn 的也是原始奖励
            experience.NextObservation = {nextObs};
            experience.IsDone = double(isDone);
            
            learn(agent, experience); 
            
            obs = nextObs;
            currentEpisodeReward = currentEpisodeReward + rawReward;
            % **核心修改**: 累加计算当前 episode 的原始折现奖励
            currentEpisodeRawDiscountedReward = currentEpisodeRawDiscountedReward + (env.env.beta^(iStep-1)) * rawReward;

            totalSteps = totalSteps + 1;
            
            if isDone, break; end
        end
        % **核心修改**: 保存两种不同的奖励值
        episodeRewards(iEpisode) = currentEpisodeReward; % 这个是用于绘图的未折现奖励总和
        episodeRawDiscountedRewards(iEpisode) = currentEpisodeRawDiscountedReward; % 这个是用于日志汇报的折现奖励总和

        
        % 更新绘图 (使用折现的奖励总和)
        plotStartIdx = max(1, iEpisode - rewardPlotWindow + 1);
        set(rewardPlotLine, 'XData', plotStartIdx:iEpisode, 'YData', episodeRawDiscountedRewards(plotStartIdx:iEpisode));
        xlim(rewardPlotAx, [plotStartIdx-1, iEpisode+1]); 
        
        if mod(iEpisode, plotFrequency) == 0 || iEpisode == 1
            agent.UseExplorationPolicy = false;
            evalObsRaw = env.env.reset(); 
            evalObs = env.normalize_obs(evalObsRaw);
            isDoneEval = false;
            c_prop_path = nan(1, maxStepsPerEpisode);
            alpha_path = nan(1, maxStepsPerEpisode);
            eval_reward = 0;
            for t = 1:maxStepsPerEpisode
                evalActionCell = getAction(agent, {evalObs});
                [evalNextObsRaw, reward, isDoneEval, loggedStep] = env.env.step(evalActionCell{1}); % 
                c_prop_path(t) =loggedStep.C_prop;
                alpha_path(t) = loggedStep.RiskyShare;
                evalObs = env.normalize_obs(evalNextObsRaw);
                % **核心修改**: 修正了 beta 的指数
                eval_reward = eval_reward + (env.env.beta^(t-1)) * reward; 
                if isDoneEval, 
                    break; 
                end

            end
            


            % 检查并保存最佳模型
            if eval_reward > best_eval_reward_so_far
                best_eval_reward_so_far = eval_reward;
                
                % 保存智能体
                save(bestAgentPath, 'agent');
                % 保存归一化统计数据
                env.saveStats(bestStatsPath);
                
                fprintf('*** Episode %d: New best raw reward: %.4f! Saving model\n', ...
                    iEpisode, best_eval_reward_so_far);

                            agent.UseExplorationPolicy = true;
            
            set(policyPlotC, 'YData', c_prop_path);
            set(policyPlotA, 'YData', alpha_path);
            title(policyPlotAx, sprintf('Agent Policy (Ep: %d) | Best Raw Reward: %.2f', ...
                iEpisode, best_eval_reward_so_far));
            

            end
            drawnow;

        end
        
        % 打印日志
        if mod(iEpisode, scoreAvgWindow) == 0 && iEpisode > 0
             % **核心修改**: 计算原始折现奖励的平均值并用于汇报
             avgDiscountedReward = mean(episodeRawDiscountedRewards(max(1, iEpisode-scoreAvgWindow+1) : iEpisode));
             elapsedTimeStr = datestr(toc/86400, 'HH:MM:SS');
             fprintf('Episode: %d | Total Steps: %d | Elapsed: %s | Avg Raw Discounted Reward (last %d): %.4f\n', ...
                 iEpisode, totalSteps, elapsedTimeStr, scoreAvgWindow, avgDiscountedReward);
        end
    end
    
    trainingStats.EpisodeReward = episodeRewards;
    fprintf('自定义训练循环完成。总用时: %s\n', datestr(toc/86400, 'HH:MM:SS'));
end