%% evaluate_checkpoint.m
% 描述:
% 这是一个独立的评估脚本，用于加载训练阶段一过程中保存的
% 某个特定检查点 (checkpoint)，并对其进行详细的性能评估。

% --- 初始化 ---
clear;
close all;


fprintf('开始评估指定的 RL 检查点...\n');

% --- 1. 参数设置 ---

% ** 请在这里修改为您想要评估的检查点文件路径 **
% 路径通常是 'results_rl/p1_checkpoints/AgentK.mat'
% 其中 K 是您想评估的 episode 编号 (例如 200, 400, ...)
AGENT_CHECKPOINT_PATH = 'results_rl/agent_phase1.mat'; 

% 其他必要文件的路径 (通常不需要修改)
STATS_PATH = 'results_rl/stats_phase1.mat';
VFI_MAT_FILE = 'result/vfi_results_crra_nopps_solo.mat';
RESULTS_SAVE_PATH = 'results_rl'; % 评估结果图的保存位置


% --- 检查文件是否存在 ---
if ~exist(AGENT_CHECKPOINT_PATH, 'file')
    error('指定的智能体检查点文件 "%s" 未找到。请检查路径和文件名是否正确。', AGENT_CHECKPOINT_PATH);
end
if ~exist(STATS_PATH, 'file')
    error('归一化统计数据文件 "%s" 未找到。请确保已至少成功运行一部分阶段一训练。', STATS_PATH);
end
if ~exist(VFI_MAT_FILE, 'file')
    error('VFI结果文件 "%s" 未找到。', VFI_MAT_FILE);
end

% --- 2. 加载智能体 ---
fprintf('正在加载智能体检查点从: %s\n', AGENT_CHECKPOINT_PATH);
load(AGENT_CHECKPOINT_PATH); % 检查点文件通常将智能体保存在名为 'agent' 的变量中

% --- 3. 调用主评估函数 ---
fprintf('智能体已加载，开始调用评估函数...\n');
evaluator = AgentEvaluator(agent_p1, STATS_PATH,[] , RESULTS_SAVE_PATH); % VFI_MAT_FILE
evaluator.run(); % 执行模拟和绘图


% 调用我们已有的 evaluate_agent 函数进行评估
% 它会处理所有的模拟和绘图工作

fprintf('\n检查点评估完成。\n');