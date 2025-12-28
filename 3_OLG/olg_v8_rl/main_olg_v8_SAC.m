% --- START OF FILE main_olg_v8_SAC.m ---
clc; % 清除命令行窗口
clear; % 清除工作区变量
close all; % 关闭所有图形窗口
fprintf('=== OLG 模型 V8 - SAC 智能体训练 (严格官方标准架构) ===\n'); % 打印模型信息
fprintf('    (在线RL，宏观变量M作为环境参数)\n'); % 打印环境参数信息
fprintf('    (决策变量：PPS缴费比例, 非PPS储蓄比例)\n'); % 打印决策变量信息
fprintf('    (网络架构：官方标准SAC - Critic[256,256], Actor[256,256], 官方超参数)\n'); % 打印网络架构信息

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

%% 4. 创建 SAC Agent
fprintf('\n--- 4. 创建 SAC Agent ---\n'); % 打印创建SAC Agent信息
numObs = obsInfo.Dimension(2); % 获取观察维度
numAct = actInfo.Dimension(1); % 获取动作维度

% --- Critic网络 (MATLAB兼容的分离式架构) ---
% 参考MATLAB官方SAC文档的Continuous Action Generation部分
% 对于连续动作空间，使用分离的状态和动作输入更兼容

% 第一个Critic网络 - 分离式架构
% 状态输入路径
statePath1 = [
    featureInputLayer(numObs, 'Normalization', 'none', 'Name', 'state_input1')
    fullyConnectedLayer(256, 'Name', 'state_fc1_1')
    reluLayer('Name', 'state_relu1_1')
    fullyConnectedLayer(256, 'Name', 'state_fc2_1')
    reluLayer('Name', 'state_relu2_1')
];

% 动作输入路径
actionPath1 = [
    featureInputLayer(numAct, 'Normalization', 'none', 'Name', 'action_input1')
    fullyConnectedLayer(64, 'Name', 'action_fc1_1')
    reluLayer('Name', 'action_relu1_1')
];

% 组合路径 - 使用拼接层正确组合状态和动作
combinedPath1 = [
    concatenationLayer(1, 2, 'Name', 'combined_concat1') % 在特征维度拼接 (256+64=320)
    fullyConnectedLayer(256, 'Name', 'combined_fc1_1')
    reluLayer('Name', 'combined_relu1_1')
    fullyConnectedLayer(1, 'Name', 'critic_output1')
];

% 构建第一个Critic的layerGraph
criticLGraph1 = layerGraph();
criticLGraph1 = addLayers(criticLGraph1, statePath1);
criticLGraph1 = addLayers(criticLGraph1, actionPath1);
criticLGraph1 = addLayers(criticLGraph1, combinedPath1);
criticLGraph1 = connectLayers(criticLGraph1, 'state_relu2_1', 'combined_concat1/in1');
criticLGraph1 = connectLayers(criticLGraph1, 'action_relu1_1', 'combined_concat1/in2');

% 第二个Critic网络 - 相同架构，不同名称
statePath2 = [
    featureInputLayer(numObs, 'Normalization', 'none', 'Name', 'state_input2')
    fullyConnectedLayer(256, 'Name', 'state_fc1_2')
    reluLayer('Name', 'state_relu1_2')
    fullyConnectedLayer(256, 'Name', 'state_fc2_2')
    reluLayer('Name', 'state_relu2_2')
];

actionPath2 = [
    featureInputLayer(numAct, 'Normalization', 'none', 'Name', 'action_input2')
    fullyConnectedLayer(64, 'Name', 'action_fc1_2')
    reluLayer('Name', 'action_relu1_2')
];

combinedPath2 = [
    concatenationLayer(1, 2, 'Name', 'combined_concat2') % 在特征维度拼接 (256+64=320)
    fullyConnectedLayer(256, 'Name', 'combined_fc1_2')
    reluLayer('Name', 'combined_relu1_2')
    fullyConnectedLayer(1, 'Name', 'critic_output2')
];

criticLGraph2 = layerGraph();
criticLGraph2 = addLayers(criticLGraph2, statePath2);
criticLGraph2 = addLayers(criticLGraph2, actionPath2);
criticLGraph2 = addLayers(criticLGraph2, combinedPath2);
criticLGraph2 = connectLayers(criticLGraph2, 'state_relu2_2', 'combined_concat2/in1');
criticLGraph2 = connectLayers(criticLGraph2, 'action_relu1_2', 'combined_concat2/in2');

% 为分离式架构创建dlnetwork (MATLAB标准方式)
% 对于多输入网络，使用结构体格式或让MATLAB自动推断
criticDlNetwork1 = dlnetwork(criticLGraph1);
criticDlNetwork2 = dlnetwork(criticLGraph2);

fprintf("Successfully created separated critic dlnetworks (MATLAB官方兼容架构).\n");
fprintf("Critic1 InputNames: "); disp(criticDlNetwork1.InputNames);
fprintf("Critic1 OutputNames: "); disp(criticDlNetwork1.OutputNames);

% --- Critic网络选项 (严格官方标准) ---
criticOptions = rlRepresentationOptions();
criticOptions.LearnRate = 1e-3; % 使用MATLAB官方推荐学习率
criticOptions.GradientThreshold = 1;
criticOptions.UseDevice = "cpu";

% 使用MATLAB标准方式创建Critic表示
critic1 = rlQValueRepresentation(criticDlNetwork1, obsInfo, actInfo, criticOptions);
critic2 = rlQValueRepresentation(criticDlNetwork2, obsInfo, actInfo, criticOptions);
fprintf("Successfully created critic representations with MATLAB官方兼容的分离式架构.\n");

% --- Actor 网络 (MATLAB官方rlSACAgent标准架构) ---
% 参考 https://ww2.mathworks.cn/help/reinforcement-learning/ref/rl.agent.rlsacagent.html
% 使用MATLAB官方推荐的连续动作SAC Actor架构 - 双输出设计

% 创建主要特征层
mainLayers = [
    featureInputLayer(numObs, 'Normalization', 'none', 'Name', 'observation')
    fullyConnectedLayer(256, 'Name', 'fc1')
    reluLayer('Name', 'relu1')
    fullyConnectedLayer(256, 'Name', 'fc2') 
    reluLayer('Name', 'relu2')
];

% 创建均值输出分支
meanLayers = [
    fullyConnectedLayer(numAct, 'Name', 'mean_fc')
];

% 创建标准差输出分支 (确保非负输出)
stdLayers = [
    fullyConnectedLayer(numAct, 'Name', 'std_fc')
    softplusLayer('Name', 'std_softplus') % 确保标准差为正数
];

% 构建layerGraph
actorLgraph = layerGraph();
actorLgraph = addLayers(actorLgraph, mainLayers);
actorLgraph = addLayers(actorLgraph, meanLayers);
actorLgraph = addLayers(actorLgraph, stdLayers);

% 连接主要特征到两个输出分支
actorLgraph = connectLayers(actorLgraph, 'relu2', 'mean_fc');
actorLgraph = connectLayers(actorLgraph, 'relu2', 'std_fc');

% 创建dlnetwork
exampleInputActor = dlarray(zeros(numObs, 1), 'CB');
actorNetwork = dlnetwork(actorLgraph, exampleInputActor);

fprintf('MATLAB官方标准Actor网络已创建。\n');
fprintf('Actor InputNames: '); disp(actorNetwork.InputNames);
fprintf('Actor OutputNames: '); disp(actorNetwork.OutputNames);

% --- Actor网络选项 (MATLAB官方标准) ---
actorOptions = rlRepresentationOptions();
actorOptions.LearnRate = 1e-3; % MATLAB官方推荐学习率
actorOptions.GradientThreshold = 1;
actorOptions.UseDevice = "cpu";

% 根据MATLAB官方rlSACAgent文档创建Actor表示
actor = rlStochasticActorRepresentation(actorNetwork, obsInfo, actInfo, ...
    'Options', actorOptions);
fprintf("Successfully created actor using MATLAB官方标准方法.\n");

% --- SAC Agent 创建 (MATLAB官方rlSACAgent标准) ---
% 参考 https://ww2.mathworks.cn/help/reinforcement-learning/ref/rl.agent.rlsacagent.html

% 创建MATLAB官方标准的rlSACAgentOptions
agentOptions = rlSACAgentOptions();
agentOptions.SampleTime = 1;
agentOptions.DiscountFactor = 0.97; % MATLAB官方推荐值
agentOptions.TargetSmoothFactor = 5e-3; % MATLAB官方推荐的软更新系数
agentOptions.ExperienceBufferLength = 1e6; % MATLAB官方推荐的经验缓冲区大小
agentOptions.MiniBatchSize = 64; % MATLAB官方推荐的批量大小

% 设置熵系数 (如果支持的话)
if isprop(agentOptions, 'EntropyWeightOptions')
    agentOptions.EntropyWeightOptions.LearnRate = 1e-3;
    agentOptions.EntropyWeightOptions.GradientThreshold = 1;
    fprintf("设置熵权重学习选项。\n");
elseif isprop(agentOptions, 'Alpha')
    agentOptions.Alpha = 0.2; % MATLAB官方推荐的熵系数
    fprintf("设置固定熵系数 Alpha = %.2f\n", agentOptions.Alpha);
else
    fprintf("使用默认熵设置。\n");
end

% 根据MATLAB官方文档创建rlSACAgent
agent = rlSACAgent(actor, [critic1, critic2], agentOptions);
fprintf('MATLAB官方标准rlSACAgent已创建。\n');

%% 5. 训练 Agent (官方标准设置)
fprintf('\n--- 5. 训练 Agent (官方标准配置) ---\n'); % 打印训练Agent信息
maxEpisodes = 15000; % 官方推荐的训练回合数
maxStepsPerEpisode = cS.aD_new; % 设置每回合最大步数 (完整生命周期，期望化死亡率)
stopTrainingValue = -20; % 设置更积极的训练目标（官方标准通常不过于保守）
scoreAveragingWindow = 100; % 官方推荐的得分平均窗口
% saveAgentDirectory = "saved_SAC_Agents_OLG_R2024b"; % 已注释：训练中不保存Agent
% if ~exist(saveAgentDirectory, 'dir'), mkdir(saveAgentDirectory); end % 已注释：训练中不创建保存目录
verbose = true; % 设置详细输出

trainOpts = rlTrainingOptions(... % 创建训练选项
    'MaxEpisodes', maxEpisodes,... % 设置最大回合数
    'MaxStepsPerEpisode', maxStepsPerEpisode,... % 设置每回合最大步数
    'ScoreAveragingWindowLength', scoreAveragingWindow,... % 设置得分平均窗口长度
    'Verbose', verbose,... % 设置详细输出
    'Plots', "training-progress",... % 设置绘图选项
    'StopTrainingCriteria', "AverageReward",... % 设置停止训练标准
    'StopTrainingValue', stopTrainingValue); % 设置停止训练值
    % 'SaveAgentCriteria', "AverageReward", ... % Save if average reward improves % 已注释：训练中不保存Agent
    % 'SaveAgentValue', -inf, ... % Save any agent that gives non -inf average reward initially % 已注释：训练中不保存Agent
    % 'SaveAgentDirectory', saveAgentDirectory); % 已注释：训练中不保存Agent目录

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

fprintf('开始标准 train 函数训练...\n'); % 打印训练开始信息
finalAgent = agent; % 预设finalAgent，以防训练失败
trainingStats = train(agent,env,trainOpts); % 训练智能体
fprintf('标准 train 函数训练完成。\n'); % 打印训练完成信息
finalAgent = agent; % 更新finalAgent
figure; plot(trainingStats.EpisodeReward); hold on; % 绘制回合奖励图
plot(trainingStats.AverageReward); legend('Episode Reward','Average Reward'); % 绘制平均奖励图
title('Training Stats from train function (R2024b)'); xlabel('Episode'); ylabel('Reward'); grid on; % 设置图表标题和标签

save('final_SAC_Agent_OLG_R2024b.mat','finalAgent', 'cS', 'paramS_for_rl', 'rng_M', 'trainOpts', 'obsInfo', 'actInfo'); % 保存最终Agent
fprintf('最终Agent已保存到 final_SAC_Agent_OLG_R2024b.mat\n'); % 打印保存成功信息

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
        fprintf('Simulating trained agent with fixed M...\n'); % 打印模拟信息
        
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
        experience_eval = sim(finalAgent,env,simOptions); % 进行模拟
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
    % convert_SAC_to_PolicyFunctions 只需要4个参数：agent, cS, paramS_rl, M_fixed
    [cPolM_sac, kPolM_sac, cPpsPolM_choice_sac, ~] = convert_SAC_to_PolicyFunctions(... % 转换为策略函数
        finalAgent, cS, paramS_for_rl, M_eval); % 调用转换函数
    save('converted_SAC_policies_R2024b.mat', 'cPolM_sac', 'kPolM_sac', 'cPpsPolM_choice_sac', 'M_eval'); % 保存策略函数
    fprintf('SAC策略已转换为矩阵并保存到 converted_SAC_policies_R2024b.mat (针对M_eval)。\n'); % 打印保存成功信息
else
    fprintf('finalAgent 或 M_eval 未定义，跳过策略转换。\n'); % 打印跳过转换信息
end

fprintf('SAC Agent 训练和处理框架完成 (R2024b attempt)。\n'); % 打印框架完成信息
% --- END OF FILE main_olg_v8_SAC.m ---