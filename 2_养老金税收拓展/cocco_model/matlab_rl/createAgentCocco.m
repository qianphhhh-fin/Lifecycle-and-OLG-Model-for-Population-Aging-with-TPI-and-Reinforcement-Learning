% createAgentCocco.m

function agent = createAgentCocco(obsInfo, actInfo, discountFactor, sampleTime)
    % 描述:
    % 创建用于Cocco环境的SAC智能体。
    
    % --- 通用网络设置 ---
    statePath = [
        featureInputLayer(obsInfo.Dimension(1), 'Normalization', 'none', 'Name', 'state')
        fullyConnectedLayer(32, 'Name', 's_fc1')
        reluLayer('Name', 's_relu1')
        fullyConnectedLayer(32, 'Name', 's_fc2')
        reluLayer('Name', 's_relu2')
        ];
    
    actionPath = [
        featureInputLayer(actInfo.Dimension(1), 'Normalization', 'none', 'Name', 'action')
        fullyConnectedLayer(32, 'Name', 'a_fc1')
        reluLayer('Name', 'a_relu1')
        ];

    commonPath = [
        concatenationLayer(1, 2, 'Name', 'concat')
        reluLayer('Name', 'c_relu1')
        fullyConnectedLayer(1, 'Name', 'c_fc_out')
        ];
    
    % --- Critic 网络 ---
    criticNetwork = layerGraph(statePath);
    criticNetwork = addLayers(criticNetwork, actionPath);
    criticNetwork = addLayers(criticNetwork, commonPath);
    criticNetwork = connectLayers(criticNetwork, 's_relu2', 'concat/in1');
    criticNetwork = connectLayers(criticNetwork, 'a_relu1', 'concat/in2');

    optimizerOpts = rlOptimizerOptions('Optimizer', 'adam', 'LearnRate', 1e-4, 'GradientThreshold', 10);
    
    critic1 = rlQValueFunction(criticNetwork, obsInfo, actInfo, 'ObservationInputNames', 'state', 'ActionInputNames', 'action');
    critic2 = rlQValueFunction(criticNetwork, obsInfo, actInfo, 'ObservationInputNames', 'state', 'ActionInputNames', 'action');

    % --- Actor 网络 ---
    actorNetwork = layerGraph(statePath);
    
    % --- 关键修正: 移除 Sigmoid 层 ---
    % Sigmoid 层将输出限制在 (0, 1) 开区间，无法达到边界值 0 或 1。
    % 此外，当输出接近边界时，会引发梯度消失问题，阻碍学习。
    % 我们让网络直接输出动作分布的均值，不加激活函数。
    % 最终的动作由环境 step 函数中的 max/min 操作来确保在 [0, 1] 闭区间内。
    meanPath = [
        fullyConnectedLayer(actInfo.Dimension(1), 'Name', 'mean_fc')
        % <--- 移除了这里的 sigmoidLayer
        ]; 
    
    stdPath = [
        fullyConnectedLayer(actInfo.Dimension(1), 'Name', 'std_fc')
        softplusLayer('Name','softplus')
        ];
        
    actorNetwork = addLayers(actorNetwork, meanPath);
    actorNetwork = addLayers(actorNetwork, stdPath);
    
    actorNetwork = connectLayers(actorNetwork,'s_relu2','mean_fc');
    actorNetwork = connectLayers(actorNetwork,'s_relu2','std_fc');
    
    actor = rlContinuousGaussianActor(actorNetwork, obsInfo, actInfo, ...
        'ActionMeanOutputNames', 'mean_fc',... % <--- 输出层名字现在是 'mean_fc'
        'ActionStandardDeviationOutputNames', 'softplus',...
        'ObservationInputNames', 'state');

    % --- 智能体选项 ---
    agentOptions = rlSACAgentOptions(...
        'SampleTime', sampleTime,...
        'DiscountFactor', discountFactor,...
        'ExperienceBufferLength', 500000,...
        'MiniBatchSize', 256,...
        'NumStepsToLookAhead', 1, ...
        'TargetSmoothFactor', 5e-3);
    
    % 确保熵调节是开启的 (默认行为)
    agentOptions.ActorOptimizerOptions = optimizerOpts;
    agentOptions.CriticOptimizerOptions = optimizerOpts;
    
    % --- 创建智能体 ---
    agent = rlSACAgent(actor, [critic1; critic2], agentOptions);
end