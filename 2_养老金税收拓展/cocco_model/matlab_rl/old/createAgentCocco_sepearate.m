function agent = createAgentCocco(obsInfo, actInfo, discountFactor, sampleTime)
    % 描述:
    % 创建用于Cocco环境的SAC智能体。
    % **V3修正**: 在Actor网络图内部直接加入 tanh 和 scaling 层，以确保
    %            输出的动作始终在 actInfo 定义的范围内。
    %            删除了不存在的 setActionMeanScalingLayer 函数调用。
    
    % --- 0. 从动作空间信息中获取上下限 ---
    actionLowerLimit = actInfo.LowerLimit;
    actionUpperLimit = actInfo.UpperLimit;
    
    % --- 1. 通用网络设置 ---
    statePath = [
        featureInputLayer(obsInfo.Dimension(1), 'Normalization', 'none', 'Name', 'state')
        fullyConnectedLayer(64, 'Name', 's_fc1')
        reluLayer('Name', 's_relu1')
        fullyConnectedLayer(64, 'Name', 's_fc2')
        reluLayer('Name', 's_relu2')
        ];
    
    actionPath = [
        featureInputLayer(actInfo.Dimension(1), 'Normalization', 'none', 'Name', 'action')
        fullyConnectedLayer(64, 'Name', 'a_fc1')
        reluLayer('Name', 'a_relu1')
        ];

    commonPath = [
        concatenationLayer(1, 2, 'Name', 'concat_critic')
        reluLayer('Name', 'c_relu1')
        fullyConnectedLayer(1, 'Name', 'c_fc_out')
        ];
    
    % --- 2. Critic 网络 (保持不变) ---
    criticNetwork = layerGraph(statePath);
    criticNetwork = addLayers(criticNetwork, actionPath);
    criticNetwork = addLayers(criticNetwork, commonPath);
    criticNetwork = connectLayers(criticNetwork, 's_relu2', 'concat_critic/in1');
    criticNetwork = connectLayers(criticNetwork, 'a_relu1', 'concat_critic/in2');

    optimizerOpts = rlOptimizerOptions('Optimizer', 'adam', 'LearnRate', 5e-3, 'GradientThreshold', 10);
    
    critic1 = rlQValueFunction(criticNetwork, obsInfo, actInfo, 'ObservationInputNames', 'state', 'ActionInputNames', 'action');
    critic2 = rlQValueFunction(criticNetwork, obsInfo, actInfo, 'ObservationInputNames', 'state', 'ActionInputNames', 'action');

    % --- 3. Actor 网络 (多头结构) ---
    actorNetwork = layerGraph(statePath);
    
    % --- 分离的消费决策头 (Consumption Head) ---
    consumptionHead = [
        fullyConnectedLayer(32, 'Name', 'c_head_fc1')
        reluLayer('Name', 'c_head_relu1')
    ];
    actorNetwork = addLayers(actorNetwork, consumptionHead);
    actorNetwork = connectLayers(actorNetwork, 's_relu2', 'c_head_fc1');
    
    % --- 消费均值路径 (带缩放) ---
    scale_c = (actionUpperLimit(1) - actionLowerLimit(1)) / 2;
    bias_c = (actionUpperLimit(1) + actionLowerLimit(1)) / 2;
    consumptionMeanPath = [
        fullyConnectedLayer(1, 'Name', 'c_mean_fc_raw')
        tanhLayer('Name', 'c_mean_tanh')
        scalingLayer('Name', 'c_scaling', 'Scale', scale_c, 'Bias', bias_c) % 直接在图中缩放
    ];
    consumptionStdPath = [
        fullyConnectedLayer(1, 'Name', 'c_std_fc')
        softplusLayer('Name', 'c_softplus')
    ];
    actorNetwork = addLayers(actorNetwork, consumptionMeanPath);
    actorNetwork = addLayers(actorNetwork, consumptionStdPath);
    actorNetwork = connectLayers(actorNetwork, 'c_head_relu1', 'c_mean_fc_raw');
    actorNetwork = connectLayers(actorNetwork, 'c_head_relu1', 'c_std_fc');
    
    % --- 分离的投资决策头 (Alpha/Portfolio Head) ---
    alphaHead = [
        fullyConnectedLayer(32, 'Name', 'alpha_head_fc1')
        reluLayer('Name', 'alpha_head_relu1')
    ];
    actorNetwork = addLayers(actorNetwork, alphaHead);
    actorNetwork = connectLayers(actorNetwork, 's_relu2', 'alpha_head_fc1');
    
    % --- 投资均值路径 (带缩放) ---
    scale_a = (actionUpperLimit(2) - actionLowerLimit(2)) / 2;
    bias_a = (actionUpperLimit(2) + actionLowerLimit(2)) / 2;
    alphaMeanPath = [
        fullyConnectedLayer(1, 'Name', 'alpha_mean_fc_raw')
        tanhLayer('Name', 'alpha_mean_tanh')
        scalingLayer('Name', 'alpha_scaling', 'Scale', scale_a, 'Bias', bias_a) % 直接在图中缩放
    ];
    alphaStdPath = [
        fullyConnectedLayer(1, 'Name', 'alpha_std_fc')
        softplusLayer('Name', 'alpha_softplus')
    ];
    actorNetwork = addLayers(actorNetwork, alphaMeanPath);
    actorNetwork = addLayers(actorNetwork, alphaStdPath);
    actorNetwork = connectLayers(actorNetwork, 'alpha_head_relu1', 'alpha_mean_fc_raw');
    actorNetwork = connectLayers(actorNetwork, 'alpha_head_relu1', 'alpha_std_fc');

    % --- 组合最终输出 ---
    concatMeanLayer = concatenationLayer(1, 2, 'Name', 'mean_concat');
    concatStdLayer = concatenationLayer(1, 2, 'Name', 'std_concat');
    
    actorNetwork = addLayers(actorNetwork, concatMeanLayer);
    actorNetwork = addLayers(actorNetwork, concatStdLayer);
    
    % 连接到缩放层的输出
    actorNetwork = connectLayers(actorNetwork, 'c_scaling', 'mean_concat/in1');
    actorNetwork = connectLayers(actorNetwork, 'alpha_scaling', 'mean_concat/in2');
    
    actorNetwork = connectLayers(actorNetwork, 'c_softplus', 'std_concat/in1');
    actorNetwork = connectLayers(actorNetwork, 'alpha_softplus', 'std_concat/in2');
    
    % --- 创建 Actor 对象 ---
    actor = rlContinuousGaussianActor(actorNetwork, obsInfo, actInfo, ...
        'ActionMeanOutputNames', 'mean_concat',...
        'ActionStandardDeviationOutputNames', 'std_concat',...
        'ObservationInputNames', 'state');

    % --- 4. 智能体选项与创建 (保持不变) ---
    agentOptions = rlSACAgentOptions(...
        'SampleTime', sampleTime,...
        'DiscountFactor', discountFactor,...
        'ExperienceBufferLength', 500000,...
        'MiniBatchSize', 256,...
        'NumStepsToLookAhead', 1, ...
        'TargetSmoothFactor', 5e-3);

    agentOptions.ActorOptimizerOptions = optimizerOpts;
    agentOptions.CriticOptimizerOptions = optimizerOpts;
    
    agent = rlSACAgent(actor, [critic1; critic2], agentOptions);
end