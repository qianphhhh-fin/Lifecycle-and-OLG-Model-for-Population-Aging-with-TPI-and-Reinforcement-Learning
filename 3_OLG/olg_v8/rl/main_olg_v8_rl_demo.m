% --- main_olg_v8_rl_demo.m ---
% 演示：在给定宏观变量和价格下，使用DRL (SAC算法) 求解OLG模型V8的家庭优化问题

clc;
clear;
close all;
rng(0); % For reproducibility

fprintf('=== OLG V8 - DRL Demo (Fixed Macro Environment) ===\n');

%% 1. 初始化V8模型参数 (SAC训练必要参数)
fprintf('\n--- 1. Loading Model Parameters (cS) ---\n');
% 尝试使用现有的参数函数，如果失败则创建简化版本
try
    cS_v8 = main_olg_v8_utils.ParameterValues_HuggettStyle();
    fprintf('  Successfully loaded parameters from main_olg_v8_utils\n');
catch
    fprintf('  main_olg_v8_utils not available, creating simplified parameters\n');
    % 创建简化的参数结构 - 只包含SAC训练必要的参数
    cS_v8 = struct();
    
    % === 状态空间维度参数 (必要) ===
    cS_v8.nw = 3;           % 收入冲击状态数
    cS_v8.aD_new = 10;      % 年龄组数
    
    % === 状态空间边界参数 (必要) ===
    cS_v8.kMin = 0;         % 普通资产下界
    cS_v8.kMax = 10;        % 普通资产上界
    cS_v8.kppsMin = 0;      % PPS资产下界
    cS_v8.kppsMax = 5;      % PPS资产上界
    
    % === 收入过程参数 (必要) ===
    cS_v8.leGridV = [0.8, 1.0, 1.2]; % epsilon网格
    cS_v8.leTrProbM = [0.7 0.2 0.1; 0.2 0.6 0.2; 0.1 0.2 0.7]; % 转移概率矩阵
    cS_v8.leProb1V = [0.33, 0.34, 0.33]; % 初始分布
    
    % === 年龄效率参数 (必要) ===
    cS_v8.ageEffV_new = ones(10,1); % 年龄效率
    cS_v8.physAgeMap = cell(10,1);
    for i = 1:10
        cS_v8.physAgeMap{i} = i;
    end
    cS_v8.aR_new = 7; % 工作年龄组数
    
    % === PPS系统参数 (必要) ===
    cS_v8.pps_active = 1;
    cS_v8.pps_withdrawal_age_min_idx = 65;
    cS_v8.pps_withdrawal_rate = 0.04;
    cS_v8.pps_return_rate_premium = 0.01;
    cS_v8.pps_contribution_age_max_idx = 60;
    cS_v8.pps_max_contrib_frac = 0.1;
    cS_v8.pps_annual_contrib_limit = 1000;
    
    % === 效用函数参数 (必要) ===
    cS_v8.cFloor = 0.01;    % 消费下界
    cS_v8.tau_c = 0.05;     % 消费税率
    cS_v8.sigma = 2;        % 风险厌恶系数
    
    % === 贴现因子参数 (必要) ===
    cS_v8.beta = 0.96;      % 贴现因子
    cS_v8.yearStep = 7.5;   % 年龄组跨度
    
    % === 年龄映射参数 (必要) ===
    cS_v8.age1_orig = 25;
    cS_v8.ageLast_orig = 100;
    cS_v8.ageRetire_orig = 65;
end

% 可以按需修改cS_v8中的某些参数用于RL演示，例如减少网格点数以加速初始测试
cS_v8.nw = 3;
cS_v8.aD_new = 10; % 减少年龄组数量可以显著加快RL初步测试
if isfield(cS_v8, 'ageLast_orig') && isfield(cS_v8, 'age1_orig')
    cS_v8.yearStep = (cS_v8.ageLast_orig - cS_v8.age1_orig + 1) / cS_v8.aD_new; % 重新计算
else
    cS_v8.yearStep = 7.5; % 默认值
end

fprintf('  Updated: nw=%d, aD_new=%d, yearStep=%.1f\n', cS_v8.nw, cS_v8.aD_new, cS_v8.yearStep);

%% 2. 定义固定的宏观环境参数 (用于本次RL演示)
% 这些参数在V8的完整一般均衡求解中是内生的或迭代更新的
% 但在这里，为了演示DRL求解家庭问题，我们将其固定
fprintf('\n--- 2. Defining Fixed Macro Environment for RL Agent ---\n');
fixedMacro = struct();
fixedMacro.R_k_net_factor_hh = 1.02; % 示例: 家庭税后资本净回报因子 (1+r_net)
fixedMacro.MPL_gross = 1.5;          % 示例: 市场毛工资率
fixedMacro.TR_total = 0.05;          % 示例: 总的政府 lump-sum 转移支付 (人均)
fixedMacro.tau_l = 0.10;             % 示例: 固定劳动所得税率
fixedMacro.theta_payg_actual_for_hh = 0.10; % 示例: 固定实际PAYG税率
fixedMacro.bV_payg = zeros(1, cS_v8.aD_new); % PAYG福利向量
% 假设退休年龄组（例如，最后1/3的年龄组）获得固定福利
% retirement_start_group = ceil(cS_v8.aR_new) + 1; % V8中aR_new是工作年龄组数量
retirement_start_group = ceil( (cS_v8.ageRetire_orig - cS_v8.age1_orig +1) / cS_v8.yearStep ) +1;

if retirement_start_group <= cS_v8.aD_new
    fixedMacro.bV_payg(retirement_start_group:cS_v8.aD_new) = 0.2; % 示例: 固定PAYG福利额
end

% 将cS_v8和fixedMacro传递给环境
envConstantParams = struct('cS', cS_v8, 'fixedMacro', fixedMacro);
fprintf('  Fixed r_k_net_hh: %.4f, MPL_gross: %.2f, tau_l: %.2f, theta_payg: %.2f\n', ...
    fixedMacro.R_k_net_factor_hh-1, fixedMacro.MPL_gross, fixedMacro.tau_l, fixedMacro.theta_payg_actual_for_hh);

%% 3. 创建强化学习环境
fprintf('\n--- 3. Creating Reinforcement Learning Environment ---\n');

% 创建环境实例
env = OLGV8EnvDemo(envConstantParams);

% 验证环境
fprintf('  Environment created successfully!\n');

% (可选) 验证环境
% validateEnvironment(env);

%% 4. 创建行动者和评论家网络
fprintf('\n--- 4. Creating Actor and Critic Networks ---\n');

% 神经网络的超参数
obsInfo = getObservationInfo(env);
actInfo = getActionInfo(env);
stateDim = obsInfo.Dimension(1);
actionDim = actInfo.Dimension(1);
hiddenUnits = [128 128]; % 示例，可以调整

% 创建行动者网络 - SAC需要分别的均值和标准差输出
% 共享层
sharedLayers = [
    featureInputLayer(stateDim,'Normalization','none','Name','observation')
    fullyConnectedLayer(hiddenUnits(1),'Name','actor_fc1')
    reluLayer('Name','actor_relu1')
    fullyConnectedLayer(hiddenUnits(2),'Name','actor_fc2')
    reluLayer('Name','actor_relu2')
    ];

% 均值输出分支
meanLayers = [
    fullyConnectedLayer(actionDim,'Name','action_mean')
    tanhLayer('Name','action_mean_tanh')
    ];

% 标准差输出分支
stdLayers = [
    fullyConnectedLayer(actionDim,'Name','action_std_raw')
    softplusLayer('Name','action_std')
    ];

% 创建层图
actorLayers = layerGraph();
actorLayers = addLayers(actorLayers, sharedLayers);
actorLayers = addLayers(actorLayers, meanLayers);
actorLayers = addLayers(actorLayers, stdLayers);

% 连接层
actorLayers = connectLayers(actorLayers, 'actor_relu2', 'action_mean');
actorLayers = connectLayers(actorLayers, 'actor_relu2', 'action_std_raw');

actorDlNetwork = dlnetwork(actorLayers);
% actorDlNetwork = initialize(actorDlNetwork); % 可能需要

% 创建评论家网络 (输入状态和行动，输出Q值)
% 状态路径
statePath = [
    featureInputLayer(stateDim,'Normalization','none','Name','observation')
    fullyConnectedLayer(hiddenUnits(1),'Name','critic_fc1_state')
    reluLayer('Name','critic_relu1_state')
    ];

% 行动路径
actionPath = [
    featureInputLayer(actionDim,'Normalization','none','Name','action')
    fullyConnectedLayer(hiddenUnits(1),'Name','critic_fc1_action')
    reluLayer('Name','critic_relu1_action')
    ];

% 合并路径
combinedPath = [
    concatenationLayer(1,2,'Name','concat')
    fullyConnectedLayer(hiddenUnits(2),'Name','critic_fc2_concat')
    reluLayer('Name','critic_relu2_concat')
    fullyConnectedLayer(1,'Name','value')
    ];

% 创建层图
criticLayers = layerGraph();
criticLayers = addLayers(criticLayers, statePath);
criticLayers = addLayers(criticLayers, actionPath);
criticLayers = addLayers(criticLayers, combinedPath);

% 连接层
criticLayers = connectLayers(criticLayers, 'critic_relu1_state', 'concat/in1');
criticLayers = connectLayers(criticLayers, 'critic_relu1_action', 'concat/in2');

% 创建第一个评论家网络
criticDlNetwork1 = dlnetwork(criticLayers);

% 创建第二个评论家网络（独立的网络结构）
criticDlNetwork2 = dlnetwork(criticLayers);

% 创建表示 (Representation)
actorOptions = rlRepresentationOptions('Optimizer','adam','LearnRate',5e-4, 'GradientThreshold',1);
criticOptions = rlRepresentationOptions('Optimizer','adam','LearnRate',1e-3, 'GradientThreshold',1);

% *** 指定使用GPU ***
if gpuDeviceCount > 0
    fprintf('  Using GPU for training.\n');
    actorOptions.UseDevice = 'gpu';
    criticOptions.UseDevice = 'gpu';
else
    fprintf('  No GPU available, using CPU for training.\n');
    actorOptions.UseDevice = 'cpu';
    criticOptions.UseDevice = 'cpu';
end

actor = rlContinuousGaussianActor(actorDlNetwork,obsInfo,actInfo,...
    'ObservationInputNames','observation',...
    'ActionMeanOutputNames','action_mean_tanh',...
    'ActionStandardDeviationOutputNames','action_std');

% 创建第一个评论家
critic1 = rlQValueFunction(criticDlNetwork1,obsInfo,actInfo,...
    'ObservationInputNames','observation',...
    'ActionInputNames','action');

% 创建第二个评论家（使用独立的网络）
critic2 = rlQValueFunction(criticDlNetwork2,obsInfo,actInfo,...
    'ObservationInputNames','observation',...
    'ActionInputNames','action');

%% 5. 创建SAC代理人
fprintf('\n--- 5. Creating SAC Agent ---\n');
% 确保贴现因子小于等于1
discountFactor = min(0.99, cS_v8.beta^cS_v8.yearStep);
fprintf('  Beta: %.4f, YearStep: %.2f, DiscountFactor: %.4f\n', cS_v8.beta, cS_v8.yearStep, discountFactor);

agentOptions = rlSACAgentOptions(...
    'SampleTime', 1, ...
    'TargetSmoothFactor',1e-3,...
    'ExperienceBufferLength',1e6,...
    'MiniBatchSize',256, ...
    'NumWarmStartSteps', 1000, ...
    'DiscountFactor', discountFactor);

% 创建SAC代理
agent = rlSACAgent(actor,[critic1,critic2],agentOptions);

%% 6. 训练代理人
fprintf('\n--- 6. Training the Agent ---\n');
trainingOptions = rlTrainingOptions(...
    'MaxEpisodes', 5000, ...
    'MaxStepsPerEpisode', cS_v8.aD_new, ...
    'ScoreAveragingWindowLength',50,...
    'Verbose',true,...
    'Plots','training-progress',...
    'StopTrainingCriteria','AverageReward',...
    'StopTrainingValue', -5, ...
    'SaveAgentCriteria',"EpisodeCount",...
    'SaveAgentValue',1000,...
    'SaveAgentDirectory', "savedAgents");

% 开始训练
doTraining = true;
if doTraining
    trainingStats = train(agent,env,trainingOptions);
else
    % (可选) 加载预训练的代理人
    % load('savedAgents/Agent1000.mat','agent'); % 修改为实际保存的文件名
    fprintf('Skipping training, attempting to load a pre-trained agent (if available).\n');
end

%% 7. 评估学习到的策略并与VFI比较 (概念性)
fprintf('\n--- 7. Evaluating Learned Policy and Comparing with VFI (Conceptual) ---\n');

if exist('agent', 'var')
    % 提取DRL学到的策略 (从行动者网络)
    % 这需要一个函数来从actor对象中获取策略，或者直接使用getAction(actor, observation)
    % DRL策略是在归一化状态和行动空间中定义的

    % 示例：在一个特定的状态下获取DRL的行动
    testObs = reset(env); % 获取一个初始观察
    fprintf('  Initial observation size: %s\n', mat2str(size(testObs)));
    
    testActionDRL_norm = getAction(agent, testObs); % 获取归一化的行动
    fprintf('  Action type: %s\n', class(testActionDRL_norm));
    
    % 检查返回的动作格式并正确提取
    if iscell(testActionDRL_norm)
        action_vec = testActionDRL_norm{1};
        fprintf('  Action is cell, extracted size: %s\n', mat2str(size(action_vec)));
    else
        action_vec = testActionDRL_norm;
        fprintf('  Action is not cell, size: %s\n', mat2str(size(action_vec)));
    end
    
    if length(action_vec) >= 2
        fprintf('DRL action for initial state: action1=%.3f, action2=%.3f\n', action_vec(1), action_vec(2));
    else
        fprintf('Warning: Action vector has only %d elements\n', length(action_vec));
        disp(action_vec);
    end

    % --- 与VFI结果比较的思路 ---
    % 1. 选择相同的固定宏观参数 (fixedMacro)
    % 2. 使用您的 main_olg_v8_utils.HHSolution_VFI_Huggett 求解VFI下的最优策略
    %    [cPolM_vfi, kPolM_vfi, cPpsPolM_vfi, valM_vfi] =
    %       main_olg_v8_utils.HHSolution_VFI_Huggett(fixedMacro.R_k_net_factor_hh,
    %       fixedMacro.MPL_gross, fixedMacro.TR_total, fixedMacro.bV_payg,
    %       paramS_for_vfi_under_fixed_macro, cS_v8);
    %    (需要构建 paramS_for_vfi_under_fixed_macro，包含固定的tau_l等)
    %
    % 3. 选择相同的状态点 (k, k_pps, eps_idx, age_idx)
    % 4. 从VFI策略矩阵中提取对应的 k_prime_vfi 和 c_pps_vfi
    % 5. 将相同的状态点（归一化后）输入到DRL的行动者网络，得到 k_prime_drl 和 c_pps_drl
    % 6. 比较 (k_prime_vfi, c_pps_vfi) 和 (k_prime_drl, c_pps_drl)
    %
    % 7. 可以绘制策略函数的切片图进行比较，类似于您在V8主脚本中做的那样。
    %    但要注意DRL的策略是基于神经网络的连续函数（或其离散化近似），
    %    而VFI的策略是定义在离散网格点上的。

    fprintf('To compare with VFI:\n');
    fprintf('  1. Run VFI solution (e.g., one iteration of solve_K_tau_l) with the SAME fixedMacro parameters.\n');
    fprintf('  2. Extract policy functions from both DRL agent and VFI solution.\n');
    fprintf('  3. Plot and compare them over a range of states.\n');

    % 绘制DRL策略的一个切片
    fprintf('\n--- Plotting DRL Policy Functions ---\n');
    
    % 设置绘图参数
    plot_a_idx_rl = min(round(cS_v8.aD_new / 2), cS_v8.aD_new); % 中年
    if plot_a_idx_rl == 0, plot_a_idx_rl = 1; end
    plot_ie_idx_rl = round(cS_v8.nw / 2); % 中等收入冲击
    plot_k_pps_val_rl = (cS_v8.kppsMax - cS_v8.kppsMin) / 2; % 中等PPS资产
    
    % 创建k网格
    k_grid_plot = linspace(cS_v8.kMin, cS_v8.kMax, 20);
    action1_drl_plot = zeros(size(k_grid_plot)); % k_prime相关的动作
    action2_drl_plot = zeros(size(k_grid_plot)); % c_pps相关的动作
    
    fprintf('  Analyzing policy for age_idx=%d, eps_idx=%d, k_pps=%.2f\n', ...
        plot_a_idx_rl, plot_ie_idx_rl, plot_k_pps_val_rl);
    
    for i_k = 1:length(k_grid_plot)
        current_k = k_grid_plot(i_k);
        
        % 构建状态向量
        state_unnorm = [current_k; plot_k_pps_val_rl; plot_ie_idx_rl; plot_a_idx_rl];
        
                 % 归一化状态（使用环境的方法）
         try
             obs_norm = env.normalizeObservation(state_unnorm);
         catch
             % 如果环境方法不可用，使用本地实现
             obs_norm = normalizeState(state_unnorm, cS_v8);
         end
        
                 % 获取动作
         action_norm = getAction(agent, obs_norm);
         
         % 检查返回的动作格式并正确提取
         if iscell(action_norm)
             action_vec = action_norm{1};
         else
             action_vec = action_norm;
         end
         
         action1_drl_plot(i_k) = action_vec(1);
         action2_drl_plot(i_k) = action_vec(2);
    end
    
    % 绘制策略函数
    figure('Name', 'DRL Policy Functions', 'Position', [100, 100, 1000, 400]);
    
    subplot(1,2,1); 
    plot(k_grid_plot, action1_drl_plot, 'b-', 'LineWidth', 2); 
    title(sprintf('DRL Action 1 (k'' related) at age=%d, eps=%d', plot_a_idx_rl, plot_ie_idx_rl));
    xlabel('k (current assets)'); 
    ylabel('Normalized Action 1');
    grid on;
    
    subplot(1,2,2); 
    plot(k_grid_plot, action2_drl_plot, 'r-', 'LineWidth', 2); 
    title(sprintf('DRL Action 2 (c_{pps} related) at age=%d, eps=%d', plot_a_idx_rl, plot_ie_idx_rl));
    xlabel('k (current assets)'); 
    ylabel('Normalized Action 2');
    grid on;
    
    fprintf('  Policy plots generated successfully!\n');
    
    %% 8. 保存RL结果
    fprintf('\n--- 8. Saving RL Results ---\n');
    
    % 创建结果保存文件夹
    result_folder = 'rl_result';
    if ~exist(result_folder, 'dir')
        mkdir(result_folder);
        fprintf('  Created directory: %s\n', result_folder);
    end
    
    % 保存RL智能体
    agent_file = fullfile(result_folder, 'trained_agent.mat');
    save(agent_file, 'agent');
    fprintf('  Saved trained agent to: %s\n', agent_file);
    
    % 保存环境参数
    env_params_file = fullfile(result_folder, 'environment_params.mat');
    save(env_params_file, 'cS_v8', 'fixedMacro', 'envConstantParams');
    fprintf('  Saved environment parameters to: %s\n', env_params_file);
    
    % 保存训练统计（如果存在）
    if exist('trainingStats', 'var')
        training_file = fullfile(result_folder, 'training_stats.mat');
        save(training_file, 'trainingStats');
        fprintf('  Saved training statistics to: %s\n', training_file);
    end
    
    % 保存策略函数数据
    policy_data = struct();
    policy_data.k_grid_plot = k_grid_plot;
    policy_data.action1_drl_plot = action1_drl_plot;
    policy_data.action2_drl_plot = action2_drl_plot;
    policy_data.plot_a_idx_rl = plot_a_idx_rl;
    policy_data.plot_ie_idx_rl = plot_ie_idx_rl;
    policy_data.plot_k_pps_val_rl = plot_k_pps_val_rl;
    
    policy_file = fullfile(result_folder, 'policy_data.mat');
    save(policy_file, 'policy_data');
    fprintf('  Saved policy function data to: %s\n', policy_file);
    
    % 保存当前时间戳
    timestamp = struct();
    timestamp.datetime = datetime('now');
    timestamp.datestr = datestr(now);
    
    timestamp_file = fullfile(result_folder, 'timestamp.mat');
    save(timestamp_file, 'timestamp');
    fprintf('  Saved timestamp to: %s\n', timestamp_file);
    
    fprintf('  ✓ All RL results saved successfully!\n');


else
    fprintf('No agent trained or loaded. Skipping evaluation.\n');
end

fprintf('\n--- DRL Demo Finished ---\n');

%% 本地辅助函数
function obs_norm = normalizeState(state_unnorm, cS)
    % 本地归一化函数
    obs_norm = zeros(4,1);
    
    % 安全地获取边界值
    if isfield(cS, 'kMin')
        kMin_val = cS.kMin;
    else
        kMin_val = 0;
    end
    
    if isfield(cS, 'kMax')
        kMax_val = cS.kMax;
    else
        kMax_val = 10;
    end
    
    if isfield(cS, 'kppsMin')
        kppsMin_val = cS.kppsMin;
    else
        kppsMin_val = 0;
    end
    
    if isfield(cS, 'kppsMax')
        kppsMax_val = cS.kppsMax;
    else
        kppsMax_val = 5;
    end
    
    if isfield(cS, 'nw')
        nw_val = cS.nw;
    else
        nw_val = 3;
    end
    
    if isfield(cS, 'aD_new')
        aD_new_val = cS.aD_new;
    else
        aD_new_val = 10;
    end
    
    obs_norm(1) = (state_unnorm(1) - kMin_val) / (kMax_val - kMin_val + 1e-6);
    obs_norm(2) = (state_unnorm(2) - kppsMin_val) / (kppsMax_val - kppsMin_val + 1e-6);
    obs_norm(3) = (state_unnorm(3)-1) / (nw_val -1 + 1e-6);
    if nw_val == 1, obs_norm(3) = 0.5; end
    obs_norm(4) = (state_unnorm(4)-1) / (aD_new_val -1 + 1e-6);
    obs_norm = max(0, min(1, obs_norm));
end

