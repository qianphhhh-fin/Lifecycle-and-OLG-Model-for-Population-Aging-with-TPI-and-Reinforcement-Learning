% --- START OF FILE main_olg_v8_TD3.m ---
clc; % 清除命令行窗口
clear; % 清除工作区变量
close all; % 关闭所有图形窗口
fprintf('=== OLG 模型 V8 - TD3 智能体训练 (连续动作原生支持) ===\n'); % 打印模型信息
fprintf('    (在线RL，宏观变量M作为环境参数)\n'); % 打印环境参数信息
fprintf('    (决策变量：PPS缴费比例, 非PPS储蓄比例)\n'); % 打印决策变量信息
fprintf('    (网络架构：Critic双网络[256,256], Actor确定性[256,256], ReLU激活)\n'); % 打印网络架构信息

%% 1. 初始化参数
fprintf('\n--- 1. 初始化参数 ---\n'); % 打印初始化参数信息
cS = main_olg_v8_utils.ParameterValues_HuggettStyle(); % 获取参数值
paramS_for_rl = struct(); % 创建一个空结构体用于强化学习参数
if ~isfield(cS, 'leGridV') || ~isfield(cS, 'leTrProbM') || ~isfield(cS, 'leProb1V') % 检查参数是否存在
    [paramS_for_rl.leLogGridV, paramS_for_rl.leTrProbM, paramS_for_rl.leProb1V] = main_olg_v8_utils.EarningProcess_olgm(cS); % 计算收入过程
    paramS_for_rl.leGridV = exp(paramS_for_rl.leLogGridV(:)); % 计算收入网格
    cS.leGridV = paramS_for_rl.leGridV; % 更新参数
    cS.leTrProbM = paramS_for_rl.leTrProbM; % 更新参数
    cS.leProb1V = paramS_for_rl.leProb1V; % 更新参数
else
    paramS_for_rl.leGridV = cS.leGridV; % 使用现有参数
    paramS_for_rl.leTrProbM = cS.leTrProbM; % 使用现有参数
    paramS_for_rl.leProb1V = cS.leProb1V; % 使用现有参数
end
paramS_for_rl.ageEffV_new = cS.ageEffV_new; % 更新年龄效率参数

%% 2. 定义宏观状态 M 的采样范围/分布
fprintf('\n--- 2. 定义宏观状态 M 的采样范围 ---\n'); % 打印采样范围信息
rng_M.R_k_net_factor = [1.01, 1.05]; % 定义净资本回报率范围
rng_M.w_gross = [1.5, 2.5]; % 定义总工资率范围
rng_M.TR_total = [0.0, 0.2]; % 定义总转移支付范围
rng_M.b_payg_avg_retiree = [0.1, 0.8]; % 定义平均PAYG福利范围
rng_M.tau_l = [0.05, 0.25]; % 定义劳动所得税率范围
rng_M.theta_payg_actual = [0.05, 0.20]; % 定义实际PAYG税率范围
fprintf('宏观参数采样范围已定义。\n'); % 打印采样范围定义完成信息

%% 3. 创建强化学习环境
fprintf('\n--- 3. 创建强化学习环境 ---\n'); % 打印创建环境信息
obsDim = 1 + 1 + 1 + 1 + 6; % 定义观察维度
obsInfo = rlNumericSpec([1, obsDim], 'Name', 'OLG States and Macro Conditions'); % 定义观察信息
actDim = 2; % 定义动作维度
actInfo = rlNumericSpec([actDim,1], 'LowerLimit', zeros(actDim,1), 'UpperLimit', ones(actDim,1), 'Name', 'Proportional Choices'); % 定义动作信息
env = OLGEnv_v8_SAC(cS, paramS_for_rl, rng_M, obsInfo, actInfo); % 创建强化学习环境
fprintf('RL环境已创建。\n'); % 打印环境创建完成信息

%% 4. 创建 TD3 Agent
fprintf('\n--- 4. 创建 TD3 Agent ---\n'); % 打印创建TD3 Agent信息
numObs = obsInfo.Dimension(2); % 获取观察维度
numAct = actInfo.Dimension(1); % 获取动作维度

% --- Critic 网络 (TD3双Critic架构) ---
% TD3使用Twin Critic，每个Critic接收状态和动作作为输入
statePathCritic = [ % 定义状态路径的Critic网络
    featureInputLayer(numObs, 'Normalization', 'none', 'Name', 'state_input_critic') % 特征输入层
    fullyConnectedLayer(256, 'Name', 'CriticStateFC1') % 第一层：256单元
    reluLayer('Name', 'CriticStateRelu1') % ReLU激活层
    fullyConnectedLayer(256, 'Name', 'CriticStateFC2') % 第二层：256单元
    reluLayer('Name', 'CriticStateRelu2') % ReLU激活层
];

actionPathCritic = [ % 定义动作路径的Critic网络
    featureInputLayer(numAct, 'Normalization', 'none', 'Name', 'action_input_critic') % 特征输入层
    fullyConnectedLayer(64, 'Name', 'CriticActionFC1') % 动作处理层：64单元
    reluLayer('Name', 'CriticActionRelu1') % ReLU激活层
];

commonPathCritic = [ % 定义公共路径的Critic网络
    concatenationLayer(1, 2, 'Name', 'critic_concat') % 拼接层
    fullyConnectedLayer(256, 'Name', 'CriticCombinedFC1') % 合并层：256单元
    reluLayer('Name', 'CriticCombinedRelu1') % ReLU激活层
    fullyConnectedLayer(1, 'Name', 'qValue') % 输出层：Q值
];

% 创建第一个Critic网络
criticLGraph1 = layerGraph(); % 创建第一个Critic图层图
criticLGraph1 = addLayers(criticLGraph1, statePathCritic); % 添加状态路径层
criticLGraph1 = addLayers(criticLGraph1, actionPathCritic); % 添加动作路径层
criticLGraph1 = addLayers(criticLGraph1, commonPathCritic); % 添加公共路径层
criticLGraph1 = connectLayers(criticLGraph1, 'CriticStateRelu2', 'critic_concat/in1'); % 连接状态路径
criticLGraph1 = connectLayers(criticLGraph1, 'CriticActionRelu1', 'critic_concat/in2'); % 连接动作路径

% 创建第二个Critic网络
criticLGraph2 = layerGraph(); % 创建第二个Critic图层图
criticLGraph2 = addLayers(criticLGraph2, statePathCritic); % 添加状态路径层
criticLGraph2 = addLayers(criticLGraph2, actionPathCritic); % 添加动作路径层
criticLGraph2 = addLayers(criticLGraph2, commonPathCritic); % 添加公共路径层
criticLGraph2 = connectLayers(criticLGraph2, 'CriticStateRelu2', 'critic_concat/in1'); % 连接状态路径
criticLGraph2 = connectLayers(criticLGraph2, 'CriticActionRelu1', 'critic_concat/in2'); % 连接动作路径

exampleStateInputCritic = dlarray(zeros(numObs, 1), 'CB'); % 创建示例状态输入
exampleActionInputCritic = dlarray(zeros(numAct, 1), 'CB'); % 创建示例动作输入

criticDlNetwork1 = dlnetwork(criticLGraph1, exampleStateInputCritic, exampleActionInputCritic, 'OutputNames', {'qValue'}); % 创建第一个Critic深度学习网络
criticDlNetwork2 = dlnetwork(criticLGraph2, exampleStateInputCritic, exampleActionInputCritic, 'OutputNames', {'qValue'}); % 创建第二个Critic深度学习网络

fprintf("Successfully converted critic layerGraphs to dlnetworks.\n"); % 打印转换成功信息

% --- Critic网络选项 ---
criticOptions = rlRepresentationOptions(); % 创建一个空的选项对象
criticOptions.LearnRate = 1e-3; % TD3推荐的Critic学习率
criticOptions.GradientThreshold = 1; % 设置梯度阈值
criticOptions.UseDevice = "cpu"; % 设置使用设备为CPU

% 创建Critic表示
critic1 = rlQValueRepresentation(criticDlNetwork1, obsInfo, actInfo, criticOptions); % 创建第一个Q值表示
critic2 = rlQValueRepresentation(criticDlNetwork2, obsInfo, actInfo, criticOptions); % 创建第二个Q值表示
fprintf("Successfully created critic representations.\n"); % 打印创建成功信息

% --- Actor 网络 (TD3确定性Actor架构) ---
% TD3使用确定性策略，输出确定性动作
actorLayers = [ % 定义TD3 Actor网络
    featureInputLayer(numObs, 'Normalization', 'none', 'Name', 'actor_input') % 状态输入层
    fullyConnectedLayer(256, 'Name', 'ActorFC1') % 第一层：256单元
    reluLayer('Name', 'ActorRelu1') % ReLU激活层
    fullyConnectedLayer(256, 'Name', 'ActorFC2') % 第二层：256单元
    reluLayer('Name', 'ActorRelu2') % ReLU激活层
    fullyConnectedLayer(numAct, 'Name', 'action_output') % 输出层：确定性动作
    sigmoidLayer('Name', 'action_sigmoid') % Sigmoid激活确保输出在[0,1]范围
];

actorLgraph = layerGraph(); % 创建Actor图层图
actorLgraph = addLayers(actorLgraph, actorLayers); % 添加Actor层
exampleInputActor = dlarray(zeros(numObs, 1), 'CB'); % 创建示例Actor输入
actorNetwork = dlnetwork(actorLgraph, exampleInputActor, 'OutputNames', {'action_sigmoid'}); % 创建Actor深度学习网络

fprintf('Actor Network - Actual InputNames:'); disp(actorNetwork.InputNames); % 打印Actor输入名称
fprintf('Actor Network - Actual OutputNames:'); disp(actorNetwork.OutputNames); % 打印Actor输出名称
fprintf("Successfully converted actorLgraph to dlnetwork.\n"); % 打印转换成功信息

% --- Actor网络选项 ---  
actorOptionsRep = rlRepresentationOptions(); % 创建Actor选项对象
actorOptionsRep.LearnRate = 1e-3; % TD3推荐的Actor学习率
actorOptionsRep.GradientThreshold = 1; % 设置梯度阈值
actorOptionsRep.UseDevice = "cpu"; % 设置使用设备为CPU

% 创建确定性Actor表示
actor = rlDeterministicActorRepresentation(actorNetwork, obsInfo, actInfo, actorOptionsRep); % 创建确定性Actor表示
fprintf("Successfully created deterministic actor representation.\n"); % 打印创建成功信息

% --- TD3 Agent 选项 ---
% 设置TD3特有的参数
discountFactor = 0.97; % 设置折扣因子

% 创建噪声模型
explorationNoise = rl.option.GaussianActionNoise; % 创建高斯探索噪声
explorationNoise.Variance = 0.1^2 * ones(numAct, 1); % 设置探索噪声方差

targetSmoothNoise = rl.option.GaussianActionNoise; % 创建目标策略平滑噪声  
targetSmoothNoise.Variance = 0.2^2 * ones(numAct, 1); % 设置目标噪声方差
% 注意：MATLAB的GaussianActionNoise不支持ClipLimit属性，噪声剪切由TD3算法内部处理

agentOptions = rlTD3AgentOptions(... % 创建TD3 Agent选项
    'SampleTime', 1, ... % 设置采样时间
    'TargetSmoothFactor', 5e-3, ... % 目标网络软更新系数
    'ExperienceBufferLength', 1e6, ... % 经验缓冲区大小
    'MiniBatchSize', 256, ... % 小批量大小
    'NumStepsToLookAhead', 1, ... % 前瞻步数
    'DiscountFactor', discountFactor, ... % 折扣因子
    'PolicyUpdateFrequency', 2, ... % TD3特色：延迟策略更新，每2次Critic更新才更新1次Actor
    'TargetPolicySmoothModel', targetSmoothNoise, ... % 目标策略平滑噪声模型
    'ExplorationModel', explorationNoise); % 探索噪声模型

fprintf('TD3噪声模型配置完成：\n');
fprintf('  探索噪声方差=%.3f\n', explorationNoise.Variance(1));
fprintf('  目标平滑噪声方差=%.3f\n', targetSmoothNoise.Variance(1));
fprintf('  噪声剪切由TD3算法内部处理\n');

% 创建TD3 Agent
agent = rlTD3Agent(actor, [critic1, critic2], agentOptions); % 创建TD3 Agent
fprintf('TD3 Agent 已创建。\n'); % 打印Agent创建成功信息

%% 5. 训练 Agent
fprintf('\n--- 5. 训练 Agent ---\n'); % 打印训练Agent信息
maxEpisodes = 15000; % 设置最大训练回合数
maxStepsPerEpisode = cS.aD_new; % 设置每回合最大步数 (完整生命周期)
stopTrainingValue = -25; % 设置停止训练阈值
scoreAveragingWindow = 50; % 设置得分平均窗口
verbose = true; % 设置详细输出

trainOpts = rlTrainingOptions(... % 创建训练选项
    'MaxEpisodes', maxEpisodes, ... % 设置最大回合数
    'MaxStepsPerEpisode', maxStepsPerEpisode, ... % 设置每回合最大步数
    'ScoreAveragingWindowLength', scoreAveragingWindow, ... % 设置得分平均窗口长度
    'Verbose', verbose, ... % 设置详细输出
    'Plots', "training-progress", ... % 设置绘图选项
    'StopTrainingCriteria', "AverageReward", ... % 设置停止训练标准
    'StopTrainingValue', stopTrainingValue); % 设置停止训练值

% 在训练前测试环境和智能体
fprintf('测试环境和智能体初始化...\n'); % 打印测试初始化信息
obs0 = reset(env); % 重置环境
fprintf('环境重置成功，观察维度: %s\n', mat2str(size(obs0))); % 打印环境重置成功信息

% 测试智能体的初始动作
act0 = getAction(agent, obs0); % 获取智能体的初始动作

% 安全地获取动作信息
try
    if iscell(act0)
        act0_values = act0{1}; % 如果是cell，取第一个元素
    elseif isstruct(act0)
        if isfield(act0, 'Data')
            act0_values = act0.Data; % 如果是结构体，取Data字段
        else
            act0_values = struct2array(act0); % 转换结构体为数组
        end
    else
        act0_values = act0; % 已经是数值类型
    end
    
    fprintf('智能体初始动作生成成功，动作维度: %s，动作值: [%.4f, %.4f]\n', ...
        mat2str(size(act0_values)), act0_values(1), act0_values(2)); % 打印动作信息
catch ME
    fprintf('智能体初始动作生成成功，但解析动作值时出错: %s\n', ME.message); % 打印错误信息
end

fprintf('开始TD3训练...\n'); % 打印训练开始信息
finalAgent = agent; % 预设finalAgent，以防训练失败
trainingStats = train(agent, env, trainOpts); % 训练智能体
fprintf('TD3训练完成。\n'); % 打印训练完成信息
finalAgent = agent; % 更新finalAgent

% 绘制训练统计
figure; plot(trainingStats.EpisodeReward); hold on; % 绘制回合奖励图
plot(trainingStats.AverageReward); legend('Episode Reward','Average Reward'); % 绘制平均奖励图
title('TD3 Training Stats'); xlabel('Episode'); ylabel('Reward'); grid on; % 设置图表标题和标签

save('final_TD3_Agent_OLG.mat','finalAgent', 'cS', 'paramS_for_rl', 'rng_M', 'trainOpts', 'obsInfo', 'actInfo'); % 保存最终Agent
fprintf('最终TD3 Agent已保存到 final_TD3_Agent_OLG.mat\n'); % 打印保存成功信息

%% 6. (可选) 评估训练好的 Agent
fprintf('\n--- 6. 评估训练好的 Agent ---\n'); % 打印评估Agent信息
if exist('finalAgent','var') && isa(finalAgent, 'rl.agent.AbstractAgent') % 检查Agent是否存在
    simOptions = rlSimulationOptions('MaxSteps',maxStepsPerEpisode,'NumSimulations',100); % 创建模拟选项
    
    % 使用与compare_rl_and_vfi.m完全一致的测试参数
    M_eval.R_k_net_factor = 1.03;    % 净资本回报率
    M_eval.w_gross = 2.0;            % 总工资率
    M_eval.TR_total = 0.1;           % 总转移支付
    M_eval.b_payg_avg_retiree = 0.4; % 平均PAYG福利
    M_eval.tau_l = 0.15;             % 劳动所得税率
    M_eval.theta_payg_actual = 0.12; % 实际PAYG税率
    M_eval.b_payg_avg_for_obs = M_eval.b_payg_avg_retiree; % 设置观察用平均PAYG福利
    
    fprintf('使用固定测试参数（与compare_rl_and_vfi.m一致）:\n'); % 打印测试参数信息
    fprintf('  R_k_net_factor = %.3f\n', M_eval.R_k_net_factor); % 打印净资本回报率
    fprintf('  w_gross = %.3f\n', M_eval.w_gross); % 打印总工资率
    fprintf('  tau_l = %.3f\n', M_eval.tau_l); % 打印劳动所得税率
    fprintf('  theta_payg_actual = %.3f\n', M_eval.theta_payg_actual); % 打印实际PAYG税率

    if ismethod(env, 'setMacroParameters') % 检查环境方法
        setMacroParameters(env, M_eval); % 设置宏观参数
        fprintf('Simulating trained TD3 agent with fixed M...\n'); % 打印模拟信息
        
        % 测试单步评估
        obs_test = reset(env); % 重置环境
        fprintf('评估环境重置成功，观察维度: %s\n', mat2str(size(obs_test))); % 打印重置成功信息
        
        act_test = getAction(finalAgent, obs_test); % 获取动作
        
        % 处理可能的RL数据结构
        if iscell(act_test)
            act_values = act_test{1}; % 如果是cell，取第一个元素
        elseif isstruct(act_test)
            if isfield(act_test, 'Data')
                act_values = act_test.Data; % 如果是结构体，取Data字段
            else
                act_values = struct2array(act_test); % 转换结构体为数组
            end
        elseif isa(act_test, 'double') || isa(act_test, 'single')
            act_values = act_test; % 已经是数值类型
        else
            try
                act_values = double(act_test); % 尝试转换为double
            catch
                act_values = [NaN, NaN]; % 如果都失败，使用NaN
                fprintf('警告：无法解析动作数据类型\n'); % 打印警告信息
            end
        end
        
        fprintf('智能体评估动作生成成功，动作值: [%.4f, %.4f]\n', act_values(1), act_values(2)); % 打印动作信息
        
        % 进行模拟
        experience_eval = sim(finalAgent, env, simOptions); % 进行模拟
        totalRewards_eval = zeros(simOptions.NumSimulations,1); % 初始化总回报
        for i = 1:simOptions.NumSimulations, totalRewards_eval(i) = sum(experience_eval(i).Reward.Data); end % 计算总回报
        avgReward_eval = mean(totalRewards_eval); % 计算平均回报
        fprintf('评估完成。在固定M下的平均回报 (来自 %d 次模拟): %.2f\n', simOptions.NumSimulations, avgReward_eval); % 打印评估结果
    else
        fprintf('环境缺少 setMacroParameters 方法，无法进行固定M的评估。\n'); % 打印缺少方法信息
    end
else
    fprintf('finalAgent 未正确定义，跳过评估。\n'); % 打印跳过评估信息
end

%% 7. 将训练好的Agent转换为与VFI兼容的策略函数格式
fprintf('\n--- 7. 转换Agent为策略函数 ---\n'); % 打印转换策略函数信息
if exist('finalAgent','var') && isa(finalAgent, 'rl.agent.AbstractAgent') && exist('M_eval','var') % 检查Agent和参数是否存在
    % 需要为TD3创建专门的转换函数，因为TD3是确定性策略
    if exist('convert_TD3_to_PolicyFunctions', 'file') == 2
        [cPolM_td3, kPolM_td3, cPpsPolM_choice_td3, ~] = convert_TD3_to_PolicyFunctions(... % 转换为策略函数
            finalAgent, cS, paramS_for_rl, M_eval); % 调用转换函数
        save('converted_TD3_policies.mat', 'cPolM_td3', 'kPolM_td3', 'cPpsPolM_choice_td3', 'M_eval'); % 保存策略函数
        fprintf('TD3策略已转换为矩阵并保存到 converted_TD3_policies.mat (针对M_eval)。\n'); % 打印保存成功信息
    else
        fprintf('convert_TD3_to_PolicyFunctions 函数不存在，跳过策略转换。\n'); % 打印跳过转换信息
        fprintf('注意：可以尝试使用 convert_SAC_to_PolicyFunctions 作为基础进行修改。\n'); % 打印提示信息
    end
else
    fprintf('finalAgent 或 M_eval 未定义，跳过策略转换。\n'); % 打印跳过转换信息
end

fprintf('TD3 Agent 训练和处理框架完成。\n'); % 打印框架完成信息
% --- END OF FILE main_olg_v8_TD3.m --- 