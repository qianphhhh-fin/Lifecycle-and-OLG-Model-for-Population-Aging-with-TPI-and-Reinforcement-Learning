function agent = createAgentCocco_sep(obsInfo, actInfo, discountFactor, sampleTime)
    % 描述:
    % 创建用于Cocco环境的SAC智能体。
    % **核心修改**: 此版本使用一个分离的多头Actor网络结构。
    % 网络有一个处理状态的共享主干，然后分为两个独立的“头”，
    % 一个头专门用于预测 c_prop，另一个专门用于预测 alpha。
    % 这种结构可能有助于学习两个动作之间解耦的策略。
    
    % --- 通用网络设置 ---
    layer_size = 128; 
    
    % 1. 共享主干 (State Path) - 处理状态观测
    statePath = [
        featureInputLayer(obsInfo.Dimension(1), 'Normalization', 'none', 'Name', 'state')
        fullyConnectedLayer(layer_size, 'Name', 's_fc1')
        reluLayer('Name', 's_relu1')
        fullyConnectedLayer(layer_size, 'Name', 's_fc2')
        reluLayer('Name', 's_relu2', 'Name', 'shared_features') % 命名主干的输出层
        ];
    
    % --- Critic 网络 (保持不变) ---
    % Critic网络评估的是 (状态, 动作) 对的价值，其结构不需要改变。
    actionPath = [
        featureInputLayer(actInfo.Dimension(1), 'Normalization', 'none', 'Name', 'action')
        fullyConnectedLayer(layer_size, 'Name', 'a_fc1')
        reluLayer('Name', 'a_relu1')
        ];

    commonPath = [
        concatenationLayer(1, 2, 'Name', 'concat_critic')
        reluLayer('Name', 'c_relu1')
        fullyConnectedLayer(1, 'Name', 'c_fc_out')
        ];
    
    criticNetwork = layerGraph(statePath);
    criticNetwork = addLayers(criticNetwork, actionPath);
    criticNetwork = addLayers(criticNetwork, commonPath);
    criticNetwork = connectLayers(criticNetwork, 'shared_features', 'concat_critic/in1');
    criticNetwork = connectLayers(criticNetwork, 'a_relu1', 'concat_critic/in2');

    criticOptimizerOpts = rlOptimizerOptions(...
        'Optimizer', 'adam', 'LearnRate', 3e-5, 'GradientThreshold', 10);
    actorOptimizerOpts = rlOptimizerOptions(...
        'Optimizer', 'adam', 'LearnRate', 1e-5, 'GradientThreshold', 10);

    critic1 = rlQValueFunction(criticNetwork, obsInfo, actInfo, 'ObservationInputNames', 'state', 'ActionInputNames', 'action');
    critic2 = rlQValueFunction(criticNetwork, obsInfo, actInfo, 'ObservationInputNames', 'state', 'ActionInputNames', 'action');

    % --- 多头 Actor 网络 ---
    actorNetwork = layerGraph(statePath);
    
    % **核心修改**: 为了促进持续探索，尤其是在alpha维度上，
    % 我们给标准差的输出层增加一个小的正偏置（bias）。
    % 这确保了即使网络权重收敛，标准差也不会趋近于0，从而保证了探索的持续进行。
    log_std_min = -5;
    log_std_max = 2;
    % 初始化一个小的正偏置，鼓励标准差的初始值不为0
    initial_bias = 1* (log_std_min + log_std_max);

    % 2a. c_prop 动作头 (Head 1)
    c_prop_mean_path = fullyConnectedLayer(1, 'Name', 'c_prop_mean');
    c_prop_std_path = [
        fullyConnectedLayer(1, 'Name', 'c_prop_std_fc', 'BiasLearnRateFactor', 1, 'Bias', initial_bias)
        softplusLayer('Name', 'c_prop_softplus')];
    
    % 2b. alpha 动作头 (Head 2)
    alpha_mean_path = fullyConnectedLayer(1, 'Name', 'alpha_mean');
    alpha_std_path = [
        fullyConnectedLayer(1, 'Name', 'alpha_std_fc', 'BiasLearnRateFactor', 2, 'Bias', initial_bias)
        softplusLayer('Name', 'alpha_softplus')];

    % 3. 将所有头添加到网络图中
    actorNetwork = addLayers(actorNetwork, c_prop_mean_path);
    actorNetwork = addLayers(actorNetwork, c_prop_std_path);
    actorNetwork = addLayers(actorNetwork, alpha_mean_path);
    actorNetwork = addLayers(actorNetwork, alpha_std_path);
    
    % 4. 连接主干和各个头
    actorNetwork = connectLayers(actorNetwork, 'shared_features', 'c_prop_mean');
    actorNetwork = connectLayers(actorNetwork, 'shared_features', 'c_prop_std_fc');
    actorNetwork = connectLayers(actorNetwork, 'shared_features', 'alpha_mean');
    actorNetwork = connectLayers(actorNetwork, 'shared_features', 'alpha_std_fc');
    
    % 5. 将分离的输出合并，以匹配Actor对象的输入要求
    concatMeanLayer = concatenationLayer(1, 2, 'Name', 'concat_mean');
    concatStdLayer = concatenationLayer(1, 2, 'Name', 'concat_std');
    
    actorNetwork = addLayers(actorNetwork, concatMeanLayer);
    actorNetwork = addLayers(actorNetwork, concatStdLayer);
    
    actorNetwork = connectLayers(actorNetwork, 'c_prop_mean', 'concat_mean/in1');
    actorNetwork = connectLayers(actorNetwork, 'alpha_mean', 'concat_mean/in2');
    actorNetwork = connectLayers(actorNetwork, 'c_prop_softplus', 'concat_std/in1');
    actorNetwork = connectLayers(actorNetwork, 'alpha_softplus', 'concat_std/in2');
    
    % 6. 创建 Actor 对象
    actor = rlContinuousGaussianActor(actorNetwork, obsInfo, actInfo, ...
        'ActionMeanOutputNames', 'concat_mean',...
        'ActionStandardDeviationOutputNames', 'concat_std',...
        'ObservationInputNames', 'state');

    % --- 智能体选项 ---
    agentOptions = rlSACAgentOptions(...
        'SampleTime', sampleTime,...
        'DiscountFactor', discountFactor,...
        'ExperienceBufferLength', 500000,...
        'MiniBatchSize', 256,...
        'NumStepsToLookAhead', 1, ...
        'TargetSmoothFactor', 5e-3);
        
    % agentOptions.EntropyWeight = 0.005;

    agentOptions.ActorOptimizerOptions = actorOptimizerOpts;
    agentOptions.CriticOptimizerOptions = criticOptimizerOpts;
    % agentOptions.EntropyWeightOptions.TargetEntropy=2;
    agentOptions.EntropyWeightOptions.EntropyWeight = 0.1;
    
    % --- 创建智能体 ---
    agent = rlSACAgent(actor, [critic1; critic2], agentOptions);
end