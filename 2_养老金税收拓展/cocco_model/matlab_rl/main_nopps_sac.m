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
clc;
rng(42); % 设置随机种子以保证可复现性

fprintf('开始执行 Cocco 生命周期模型的 RL 求解...\n');

% --- A. 路径和参数设置 ---
results_path = 'results_rl';
if ~exist(results_path, 'dir'), mkdir(results_path); end

VFI_MAT_FILE = 'result/vfi_results_crra_nopps_solo.mat';
if ~exist(VFI_MAT_FILE, 'file')
    error('VFI结果文件 "%s" 未找到。请先运行 main_nopps.m。', VFI_MAT_FILE);
end

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
agent_p1 = createAgentCocco(obsInfo, actInfo, DISCOUNT_FACTOR, AGENT_SAMPLE_TIME);

% 4. 设置训练选项
trainOpts_p1 = rlTrainingOptions(...
    'MaxEpisodes', 3000, ...
    'MaxStepsPerEpisode', env_p1_base.tn, ...
    'ScoreAveragingWindowLength', 50, ...
    'Verbose', true, ...
    'Plots', 'training-progress', ...
    'StopTrainingCriteria', 'AverageReward', ...
    'StopTrainingValue', 1000, ... % 设置一个较高的值，让训练由步数决定
    'SaveAgentCriteria', 'EpisodeFrequency',...
    'SaveAgentValue', 200, ...
    'SaveAgentDirectory', fullfile(results_path, 'p1_checkpoints'));

if TIMESTEPS_PHASE1 > 0
    % 5. 训练智能体
    trainingStats_p1 = train(agent_p1, env_p1, trainOpts_p1);

    % 6. 保存智能体和归一化统计数据
    save(agent_phase1_path, 'agent_p1');
    env_p1.saveStats(stats_phase1_path);

    fprintf('阶段一训练完成。智能体已保存到 %s\n', agent_phase1_path);
    fprintf('归一化统计数据已保存到 %s\n', stats_phase1_path);
    
    % 加载最优的智能体以进入下一阶段
    agentToLoad = trained_agent_p1; 
else
    fprintf('跳过阶段一训练，尝试加载已有文件...\n');
    load(agent_phase1_path, 'trained_agent_p1');
    agentToLoad = trained_agent_p1;
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