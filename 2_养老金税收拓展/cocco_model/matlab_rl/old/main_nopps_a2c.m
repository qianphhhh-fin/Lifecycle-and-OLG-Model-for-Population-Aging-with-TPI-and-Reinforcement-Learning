
% =========================================================================
%               生命周期模型求解与模拟 (A2C 深度强化学习方法)
% =========================================================================
% 描述:
% 本脚本使用 Advantage Actor-Critic (A2C) 算法求解并模拟一个生命周期模型。
% A2C 是一种深度强化学习方法，其中:
% - Actor (策略网络): 直接学习从状态到行动的映射 (决策函数)。
% - Critic (值网络): 学习评估状态的优劣 (值函数)。
% 这种方法不依赖于状态空间的离散网格，而是使用神经网络作为函数逼近器。
%
% 模型特点与原 VFI 版本一致:
% - CRRA 效用函数
% - 持久性和暂时性收入冲击
% - 工作期和退休期
% - 无风险和风险资产投资
% - 工资税 (tau_y)
%
% 结构:
% 1. setup_parameters: 初始化所有模型和 DRL 训练参数。
% 2. solve_model_drl: 搭建并训练 Actor 和 Critic 神经网络。
% 3. simulate_model: 使用训练好的 Actor 网络进行前向蒙特卡洛模拟。
%
% 作者: Gemini
% 日期: 2025-09-15
% =========================================================================


% --- 初始化 ---
clear;
close all;
clc;

% --- 主流程控制器 ---
fprintf('开始执行生命周期模型 (A2C 深度强化学习)...\n');

% 1. 设置所有模型参数
cS = setup_parameters();
save_path = cS.save_path;
file_path = fullfile(save_path, 'drl_results_crra_nopps_solo.mat');

% 【新增】加载VFI结果作为对比基准
vfi_file_path = fullfile(save_path, 'vfi_results_crra_nopps_solo.mat');
if exist(vfi_file_path, 'file')
    fprintf('加载VFI策略用于对比...\n');
    vfi_cS = load(vfi_file_path,'cS');
    load(vfi_file_path,'vfi_results');
    vfi_results.cS = vfi_cS.cS;
else
    fprintf('警告: 未找到VFI策略文件 (%s)，将不进行奖励对比。\n', vfi_file_path);
    vfi_results = []; % 传入空值
end

% 2. 求解模型 (训练 DRL 网络)
fprintf('\n===== 步骤 2: 求解模型 (DRL 训练) =====\n');
if ~exist(file_path,'file')
    % 【修改】将vfi_results传入求解器
    drl_results = solve_model_drl(cS, vfi_results);
    if ~exist(save_path, 'dir'), mkdir(save_path); end
    save(file_path, 'drl_results', 'cS');
    fprintf('DRL 模型训练结果已保存到 %s\n', file_path);
else
    load(file_path, 'drl_results', 'cS');
    fprintf('已从 %s 加载预训练的 DRL 模型。\n', file_path);
end

% 3. 模拟模型
fprintf('\n===== 步骤 3: 模拟模型 =====\n');
simulate_model(drl_results, cS);

fprintf('\n模型执行完毕。\n');

%% =====================================================================
%                       1. 参数设置函数
%  =====================================================================
function cS = setup_parameters()
fprintf('===== 步骤 1: 设置模型参数 =====\n');

cS = struct();

% --- A. 生命周期与人口结构 ---
cS.tb = 22;       % 初始工作年龄
cS.tr = 25;       % 退休年龄
cS.td = 30;      % 最高寿命
cS.tn = cS.td - cS.tb + 1; % 总期数

% --- B. 经济主体偏好 (CRRA) ---
cS.beta = 0.95;     % 折现因子
cS.gamma = 3.84;    % 相对风险规避系数

% --- C. 收入过程 ---
cS.aa = -2.170042 + 2.700381;
cS.b1 = 0.16818;
cS.b2 = -0.0323371 / 10;
cS.b3 = 0.0019704 / 100;
cS.smay =  sqrt(0.169993); % 暂时冲击标准差 (u_t)
cS.smav =  sqrt(0.112572); % 持久冲击标准差 (z_t)

% --- D. 资产与回报率 ---
cS.rf = 1.02;     % 无风险总回报率
cS.mu = 0.04;     % 风险资产超额回报率
cS.sigr = 0.27;   % 风险资产回报率标准差

% --- E. 养老金与税收 ---
cS.ret_fac = 0.6827;      % 退休后基础养老金 (替代率)
cS.tau_y = 0.06;          % 工资税率

% --- F. DRL 和神经网络参数 ---
cS.drl.num_episodes = 5000;         % 训练的总模拟生命周期次数
cS.drl.actor_lr = 1e-4;             % Actor学习率
cS.drl.critic_lr = 1e-3;            % Critic学习率
cS.drl.state_dim = 3;               % 状态维度 [normalized_wealth, normalized_income, normalized_age]
cS.drl.action_dim = 2;              % 动作维度 [consumption_ratio, alpha]
cS.drl.actor_hidden_units = [64, 64]; % Actor网络隐藏层结构
cS.drl.critic_hidden_units = [64, 64];% Critic网络隐藏层结构
cS.drl.exploration_std = [0.3; 0.5];% 消费使用较小噪声(0.2)，alpha使用较大噪声(0.6)
cS.drl.wealth_scale = 10;           % 用于归一化财富/现金状态的缩放因子
% 【！！！修改！！！】为收入设置一个独立的、更合适的缩放因子
cS.drl.income_scale = 2.0;            % 收入的尺度通常在1附近，所以用一个较小的因子

% --- G. 其他数值参数 ---
cS.save_path = 'result/'; % 结果保存路径
cS.n_shocks = 5;      % 【新增】离散化的随机冲击节点数 (与VFI保持一致)
cS.corr_z_epsilon = 0; % 【新增】与VFI一致，即使为0也要定义
cS.corr_u_epsilon = 0; % 【新增】与VFI一致，即使为0也要定义

% --- H. 派生参数与预计算 ---
% 生存概率
cS.survprob = ones(cS.tn - 1, 1); % 简化生存概率为1

% 确定性收入剖面
cS.f_y = zeros(cS.tr - cS.tb + 1, 1);
for i1 = cS.tb:cS.tr
    age = i1;
    cS.f_y(age - cS.tb + 1) = exp(cS.aa + cS.b1*age + cS.b2*age^2 + cS.b3*age^3);
end

% 【新增】冲击离散化 (与VFI版本完全一致)
tauchenoptions.parallel=0;
[grid, weig_matrix] = discretizeAR1_Tauchen(0, 0, 1, cS.n_shocks, 2, tauchenoptions);
cS.shock_grid = grid;
cS.shock_weig = diag(weig_matrix);

% 风险资产回报率网格
cS.gret = cS.rf + cS.mu + cS.shock_grid * cS.sigr;

% 状态转移概率 (三维冲击)
cS.nweig1 = zeros(cS.n_shocks, cS.n_shocks, cS.n_shocks);
for i6 = 1:cS.n_shocks
    for i7 = 1:cS.n_shocks
        for i8 = 1:cS.n_shocks
            cS.nweig1(i6, i7, i8) = cS.shock_weig(i6) * cS.shock_weig(i7) * cS.shock_weig(i8);
        end
    end
end
% 收入过程网格 (与VFI版本完全一致)
cS.yh = zeros(cS.n_shocks, cS.n_shocks); % 暂时冲击 exp(u_t)
cS.yp = zeros(cS.n_shocks, cS.n_shocks); % 持久冲击增长 exp(z_t)
for i1 = 1:cS.n_shocks
    grid2_u = cS.shock_grid(i1) .* cS.corr_u_epsilon + cS.shock_grid .* sqrt(1 - cS.corr_u_epsilon^2);
    cS.yh(:, i1) = exp(grid2_u * cS.smay);

    grid2_z = cS.shock_grid(i1) .* cS.corr_z_epsilon + cS.shock_grid .* sqrt(1 - cS.corr_z_epsilon^2);
    cS.yp(:, i1) = exp(grid2_z * cS.smav);
end


cS.gyp = zeros(cS.n_shocks, cS.n_shocks, cS.tn - 1); % 总收入增长

% 工作期收入增长
work_periods = cS.tr - cS.tb;
for t = 1:work_periods
    G_t = cS.f_y(t+1) / cS.f_y(t); % 确定性增长
    % 【修正】使用与VFI版本完全相同的构造方式
    cS.gyp(:,:,t) = repmat(G_t, cS.n_shocks, cS.n_shocks) .* cS.yp;
end

% 退休期收入增长 (无增长)
cS.gyp(:,:,(work_periods+1):(cS.tn-1)) = 1.0;


fprintf('参数设置完成。\n');
end

%% =====================================================================
%                       2. DRL 求解器函数
%  =====================================================================

%% =====================================================================
%                       2. DRL 求解器函数
%  =====================================================================

function results = solve_model_drl(cS, vfi_results)

% --- 1. 创建 Actor 和 Critic 网络 ---
[actor, critic] = createActorCriticNetworks(cS);

% 为 Adam 优化器初始化一阶和二阶矩估计
actor_m = []; actor_v = [];
critic_m = []; critic_v = [];

% --- 2. 准备VFI策略插值 (用于对比) ---
do_comparison = ~isempty(vfi_results);
if do_comparison
    fprintf('正在准备VFI策略插值函数用于对比...\n');
    C_policy_vfi = vfi_results.C_policy;
    A_policy_vfi = vfi_results.A_policy;
    C_interp_vfi = cell(1, cS.tn);
    A_interp_vfi = cell(1, cS.tn);
    for t_idx = 1:cS.tn
        C_interp_vfi{t_idx} = griddedInterpolant(vfi_results.cS.gcash, C_policy_vfi(:,t_idx), 'spline', 'linear');
        A_interp_vfi{t_idx} = griddedInterpolant(vfi_results.cS.gcash, A_policy_vfi(:,t_idx), 'spline', 'linear');
    end
end

% --- 3. 初始化训练历史和实时监控图 ---
fprintf('开始 A2C 训练...\n');
history.episode_rewards = zeros(1, cS.drl.num_episodes);
history.vfi_episode_rewards = zeros(1, cS.drl.num_episodes);
history.mean_actor_loss = zeros(1, cS.drl.num_episodes);
history.mean_critic_loss = zeros(1, cS.drl.num_episodes);
history.mean_advantage = zeros(1, cS.drl.num_episodes);
history.mean_c_ratio = zeros(1, cS.drl.num_episodes);
history.mean_alpha = zeros(1, cS.drl.num_episodes);

% 创建图形窗口和子图
fig = figure('Name', 'A2C Training Monitor', 'Position', [100, 100, 1400, 800]);
subplot_handles.rewards = subplot(2, 3, 1);
title('Episode Rewards'); xlabel('Episode'); ylabel('Total Discounted Reward'); grid on; hold on;
subplot_handles.losses = subplot(2, 3, 2);
title('Mean Actor & Critic Loss'); xlabel('Episode'); ylabel('Loss'); grid on; hold on;
subplot_handles.advantage = subplot(2, 3, 3);
title('Mean Advantage (Normalized)'); xlabel('Episode'); ylabel('Advantage'); grid on; hold on;
subplot_handles.actions = subplot(2, 3, 4);
title('Mean Actions'); xlabel('Episode'); ylabel('Ratio/Share'); grid on; hold on;

% 初始化绘图句柄
plot_handles.reward_ma = plot(subplot_handles.rewards, NaN, NaN, 'b-', 'LineWidth', 2, 'DisplayName', 'A2C Reward (MA)');
if do_comparison
    plot_handles.vfi_reward_ma = plot(subplot_handles.rewards, NaN, NaN, 'r--', 'LineWidth', 1.5, 'DisplayName', 'VFI Reward (MA)');
end
legend(subplot_handles.rewards, 'show', 'Location', 'southeast');

plot_handles.actor_loss = plot(subplot_handles.losses, NaN, NaN, 'g-', 'DisplayName', 'Actor Loss (MA)');
plot_handles.critic_loss = plot(subplot_handles.losses, NaN, NaN, 'm-', 'DisplayName', 'Critic Loss (MA)');
legend(subplot_handles.losses, 'show', 'Location', 'northeast');

plot_handles.advantage = plot(subplot_handles.advantage, NaN, NaN, 'k-', 'DisplayName', 'Advantage (MA)');
legend(subplot_handles.advantage, 'show', 'Location', 'northeast');

plot_handles.c_ratio = plot(subplot_handles.actions, NaN, NaN, 'c-', 'LineWidth', 2, 'DisplayName', 'Consumption Ratio (MA)');
plot_handles.alpha = plot(subplot_handles.actions, NaN, NaN, 'y-', 'LineWidth', 2, 'DisplayName', 'Risky Share \alpha (MA)');
legend(subplot_handles.actions, 'show', 'Location', 'northeast');

global_step_counter = 0;
update_freq = 20;
reward_window = 100;

% --- 4. 训练循环 ---
for episode = 1:cS.drl.num_episodes
    
    % -- 初始化此episode的存储 --
    episode_buffer.states = cell(1, cS.tn);
    episode_buffer.rewards = zeros(1, cS.tn);
    episode_buffer.action_means = cell(1, cS.tn);
    episode_buffer.action_samples = cell(1, cS.tn);
    
    current_w = 0.1;
    cumulative_g = 1.0; 
    
    if do_comparison, current_vfi_w = current_w; vfi_episode_reward = 0; end
    
    work_periods = cS.tr - cS.tb;
    u_shock_indices = randi(cS.n_shocks, cS.tn, 1);
    z_shock_indices = randi(cS.n_shocks, cS.tn, 1);
    r_shock_indices = randi(cS.n_shocks, cS.tn, 1);
    
    % -- A. 前向模拟并收集经验 --
    actual_episode_length = 0;
    for t = 1:cS.tn
        actual_episode_length = t;
        
        u_shock_val = cS.yh(u_shock_indices(t), 1);
        if t <= work_periods, income = (1 - cS.tau_y) * u_shock_val;
        else, income = cS.ret_fac; end
        cash_on_hand = current_w + income;
        
        % 【！！！修改！！！】使用独立的 income_scale 来归一化收入
        state = dlarray([current_w / cS.drl.wealth_scale; ...
                         income / cS.drl.income_scale; ...
                         (t - 1) / (cS.tn - 1)], 'CB');

        action_mean = forward(actor, state);
        
        action_sample = action_mean + dlarray(randn(cS.drl.action_dim, 1) .* cS.drl.exploration_std, 'CB');
        
        action_sample_data = extractdata(action_sample);
        c_ratio_logit = action_sample_data(1);
        alpha_logit = action_sample_data(2);
        
        c_ratio = 1 / (1 + exp(-c_ratio_logit)); 
        alpha = 1 / (1 + exp(-alpha_logit));
        c_ratio = max(1e-6, min(0.999, c_ratio));
        
        consumption = c_ratio * cash_on_hand;
        consumption = max(consumption, 1e-6);
        sav = cash_on_hand - consumption;
        
        reward_scale = abs(1 - cS.gamma);
        if cS.gamma == 1
            instant_utility = log(consumption);
            reward = instant_utility + log(cumulative_g); 
        else
            instant_utility = (consumption.^(1 - cS.gamma) - 1) / (1 - cS.gamma);
            reward = (cumulative_g^(1 - cS.gamma)) * instant_utility;
        end
        
        episode_buffer.states{t} = state;
        episode_buffer.rewards(t) = reward / reward_scale;
        episode_buffer.action_means{t} = action_mean;
        episode_buffer.action_samples{t} = action_sample;
        
        if do_comparison
             vfi_cash_on_hand = current_vfi_w + income;
             vfi_c = C_interp_vfi{t}(vfi_cash_on_hand);
             vfi_a = A_interp_vfi{t}(vfi_cash_on_hand);
             vfi_c = max(1e-6, min(vfi_c, vfi_cash_on_hand*0.9999)); vfi_a = max(0, min(vfi_a, 1.0));
             vfi_sav = vfi_cash_on_hand - vfi_c;
             if cS.gamma == 1, vfi_reward = log(vfi_c); else, vfi_reward = (vfi_c.^(1 - cS.gamma)) / (1 - cS.gamma); end
             vfi_episode_reward = vfi_episode_reward + (cS.beta^(t-1)) * vfi_reward;
        end

        next_w = 0; next_vfi_w = 0;
        gpy = 1.0; 
        if t < cS.tn
            yp_val = cS.yp(z_shock_indices(t), 1); gret_val = cS.gret(r_shock_indices(t));
            if t < work_periods
                G_t = cS.f_y(t+1) / cS.f_y(t); 
                gpy = G_t * yp_val;
            end
            portfolio_return = (1 - alpha) * cS.rf + alpha * gret_val;
            next_w = portfolio_return * sav / gpy;
            if do_comparison, vfi_portfolio_return = (1 - vfi_a) * cS.rf + vfi_a * gret_val; next_vfi_w = vfi_portfolio_return * vfi_sav / gpy; end
        end
        current_w = next_w;
        cumulative_g = cumulative_g * gpy;

        if do_comparison, current_vfi_w = next_vfi_w; end
        if current_w <= 1e-6, break; end
    end
    
    % -- B. 计算回报并进行网络更新 --
    rewards_to_process = episode_buffer.rewards(1:actual_episode_length);

    returns = zeros(1, actual_episode_length);
    returns(actual_episode_length) = rewards_to_process(actual_episode_length);
    for t = (actual_episode_length-1):-1:1
        returns(t) = rewards_to_process(t) + cS.beta * cS.survprob(t) * returns(t+1);
    end
    
    batch_states = cat(2, episode_buffer.states{1:actual_episode_length});
    batch_action_means = cat(2, episode_buffer.action_means{1:actual_episode_length});
    batch_action_samples = cat(2, episode_buffer.action_samples{1:actual_episode_length});
    
    global_step_counter = global_step_counter + 1;
    [actor_grads, critic_grads, actor_loss, critic_loss, advantages] = ...
        dlfeval(@computeBatchGradients, actor, critic, batch_states, dlarray(returns, 'CB'), ...
                batch_action_means, batch_action_samples, cS.drl.exploration_std);

    [actor.Learnables, actor_m, actor_v] = adamupdate(actor.Learnables, actor_grads, actor_m, actor_v, global_step_counter, cS.drl.actor_lr);
    [critic.Learnables, critic_m, critic_v] = adamupdate(critic.Learnables, critic_grads, critic_m, critic_v, global_step_counter, cS.drl.critic_lr);
    
    % -- C. 记录和更新监控图 --
    original_rewards = zeros(1, actual_episode_length);
    for t_hist = 1:actual_episode_length
         action_sample_data = extractdata(episode_buffer.action_samples{t_hist});
         state_data = extractdata(episode_buffer.states{t_hist});
         % 在历史记录中也要正确反推
         norm_w = state_data(1) * cS.drl.wealth_scale;
         norm_y = state_data(2) * cS.drl.income_scale;
         cash_on_hand_hist = norm_w + norm_y;
         c_ratio_hist = 1 / (1 + exp(-action_sample_data(1)));
         c_hist = c_ratio_hist * cash_on_hand_hist;
         c_hist = max(1e-6, c_hist);
         if cS.gamma == 1
            original_rewards(t_hist) = log(c_hist);
         else
            original_rewards(t_hist) = (c_hist.^(1 - cS.gamma)) / (1 - cS.gamma);
         end
    end
    history.episode_rewards(episode) = sum(original_rewards .* (cS.beta.^(0:actual_episode_length-1)));

    if do_comparison, history.vfi_episode_rewards(episode) = vfi_episode_reward; end
    history.mean_actor_loss(episode) = extractdata(mean(actor_loss));
    history.mean_critic_loss(episode) = extractdata(mean(critic_loss));
    history.mean_advantage(episode) = extractdata(mean(advantages));
    
    temp_c = zeros(1, actual_episode_length); temp_a = zeros(1, actual_episode_length);
    for t=1:actual_episode_length
        action_sample_data = extractdata(episode_buffer.action_samples{t});
        temp_c(t) = 1 / (1 + exp(-action_sample_data(1)));
        temp_a(t) = 1 / (1 + exp(-action_sample_data(2)));
    end
    history.mean_c_ratio(episode) = mean(temp_c);
    history.mean_alpha(episode) = mean(temp_a);

    if mod(episode, update_freq) == 0 || episode == cS.drl.num_episodes
        start_idx = max(1, episode - reward_window + 1);
        
        ma_rewards = movmean(history.episode_rewards(1:episode), [reward_window-1, 0]);
        ma_actor_loss = movmean(history.mean_actor_loss(1:episode), [reward_window-1, 0]);
        ma_critic_loss = movmean(history.mean_critic_loss(1:episode), [reward_window-1, 0]);
        ma_advantage = movmean(history.mean_advantage(1:episode), [reward_window-1, 0]);
        ma_c_ratio = movmean(history.mean_c_ratio(1:episode), [reward_window-1, 0]);
        ma_alpha = movmean(history.mean_alpha(1:episode), [reward_window-1, 0]);

        set(plot_handles.reward_ma, 'XData', 1:episode, 'YData', ma_rewards);
        if do_comparison
            ma_vfi_rewards = movmean(history.vfi_episode_rewards(1:episode), [reward_window-1, 0]);
            set(plot_handles.vfi_reward_ma, 'XData', 1:episode, 'YData', ma_vfi_rewards);
        end
        set(plot_handles.actor_loss, 'XData', 1:episode, 'YData', ma_actor_loss);
        set(plot_handles.critic_loss, 'XData', 1:episode, 'YData', ma_critic_loss);
        set(plot_handles.advantage, 'XData', 1:episode, 'YData', ma_advantage);
        set(plot_handles.c_ratio, 'XData', 1:episode, 'YData', ma_c_ratio);
        set(plot_handles.alpha, 'XData', 1:episode, 'YData', ma_alpha);
        
        drawnow;
        
        fprintf('Episode: %d | Avg Reward (MA): %.4f', episode, ma_rewards(end));
        if do_comparison, fprintf(' (VFI: %.4f)', ma_vfi_rewards(end)); end
        fprintf(' | Avg Advantage (MA): %.4f\n', ma_advantage(end));
    end
end

fprintf('DRL 训练完成。\n');
results.actor = actor;
results.critic = critic;
end


%% =====================================================================
%                  DRL - 网络创建与梯度计算
% =====================================================================

function [actor, critic] = createActorCriticNetworks(cS)
% 创建 Actor 和 Critic 网络的辅助函数

% --- Actor 网络 ---
% 输入: 2x1 (财富, 年龄) -> 输出: 2x1 (消费比例的logit, alpha的logit)
actor_layers = [
    featureInputLayer(cS.drl.state_dim, 'Name', 'state', 'Normalization', 'none')
    fullyConnectedLayer(cS.drl.actor_hidden_units(1), 'Name', 'afc1')
    tanhLayer('Name', 'atanh1')
    fullyConnectedLayer(cS.drl.actor_hidden_units(2), 'Name', 'afc2')
    tanhLayer('Name', 'atanh2')
    fullyConnectedLayer(cS.drl.action_dim, 'Name', 'action_mean')
];
actor = dlnetwork(layerGraph(actor_layers));

% --- Critic 网络 ---
% 输入: 2x1 (财富, 年龄) -> 输出: 1x1 (状态价值)
critic_layers = [
    featureInputLayer(cS.drl.state_dim, 'Name', 'state', 'Normalization', 'none')
    fullyConnectedLayer(cS.drl.critic_hidden_units(1), 'Name', 'cfc1')
    reluLayer('Name', 'crelu1')
    fullyConnectedLayer(cS.drl.critic_hidden_units(2), 'Name', 'cfc2')
    reluLayer('Name', 'crelu2')
    fullyConnectedLayer(1, 'Name', 'value')
];
critic = dlnetwork(layerGraph(critic_layers));
end



% 【新增】用于批处理梯度计算的辅助函数
function [actor_grads, critic_grads, actor_loss, critic_loss, advantages_for_logging] = computeBatchGradients(actor, critic, states, returns, action_means, action_samples, exploration_std)

    % Critic 更新
    state_values = forward(critic, states);
    
    advantages = returns - state_values;
    
    critic_loss = huber(advantages);

    % Actor 更新
    log_probs = -0.5 * sum(((action_samples - action_means) ./ exploration_std).^2, 1);
    
    % 【！！！修改！！！】 移除了对advantage的标准化步骤
    % normalized_advantages = (advantages - mean(advantages)) / (std(advantages) + 1e-8);
    % 我们直接使用原始的advantages来计算actor_loss
    
    actor_loss = -mean(log_probs .* extractdata(advantages));
    
    % 计算梯度
    critic_grads = dlgradient(mean(critic_loss), critic.Learnables);
    actor_grads = dlgradient(actor_loss, actor.Learnables);

    % 为了日志记录，我们仍然可以记录标准化后的advantage，以便在图中保持一致的尺度
    advantages_for_logging = (advantages - mean(advantages)) / (std(advantages) + 1e-8);
    % advantages_for_logging = advantages;
end


%% =====================================================================
%                       3. 模拟器函数
% =====================================================================
function simulate_model(results, cS)

fprintf('模拟开始...\n');
rng(42); % 可复现性
nsim = 10000;

% 解包训练好的Actor网络
actor_network = results.actor;

% 初始化模拟数组
simGPY = ones(cS.tn, nsim);
simY_norm_shock = zeros(cS.tn, nsim);
simR = zeros(cS.tn, nsim);
simW_norm = zeros(cS.tn + 1, nsim); % 状态变量: 归一化财富 (期初)
simC = zeros(cS.tn, nsim);          % 决策变量: 绝对消费
simA = zeros(cS.tn, nsim);          % 决策变量: 风险资产比例
simCash = zeros(cS.tn, nsim);      % For debugging or other plots

% --- 1. 模拟外生冲击 ---
fprintf('  生成收入和回报率路径...\n');
work_periods = cS.tr - cS.tb;

for i1 = 1:nsim
    z_shocks = randn(work_periods, 1) * cS.smav;
    u_shocks_log = randn(cS.tn, 1) * cS.smay;
    r_shocks = randn(cS.tn, 1) * cS.sigr;
    
    if work_periods > 0
        simGPY(2:work_periods+1, i1) = cS.f_y(2:end) ./ cS.f_y(1:end-1) .* exp(z_shocks);
    end
    simY_norm_shock(:, i1) = exp(u_shocks_log);
    simR(:, i1) = cS.rf + cS.mu + r_shocks;
end

% --- 2. 迭代模拟生命周期决策 (使用Actor网络) ---
fprintf('  模拟生命周期决策 (使用Actor网络)...\n');
simW_norm(1, :) = 0.1; % 初始财富

for t = 1:cS.tn
    % 计算当期收入和可支配现金
    if t <= work_periods
        income = (1-cS.tau_y) * simY_norm_shock(t, :);
    else
        income = repmat(cS.ret_fac, 1, nsim);
    end
    cash_on_hand = simW_norm(t, :) + income;
    
    % 【！！！修改！！！】准备三维状态输入 (batch) 并使用正确的 income_scale
    state_norm_w = simW_norm(t, :) / cS.drl.wealth_scale;
    state_norm_y = income / cS.drl.income_scale;
    state_norm_age = repmat((t - 1) / (cS.tn - 1), 1, nsim);
    current_states_dl = dlarray([state_norm_w; state_norm_y; state_norm_age], 'CB');

    % 使用Actor网络预测动作 (使用均值进行模拟，不加探索噪声)
    action_means = predict(actor_network, current_states_dl);
    
    % 动作后处理
    c_ratios = sigmoid(action_means(1,:));
    simA(t, :) = sigmoid(action_means(2,:));
    
    % 计算消费和储蓄
    simC(t, :) = c_ratios .* cash_on_hand;
    simC(t, :) = min(simC(t, :), cash_on_hand * 0.9999); % 确保消费可行
    liquid_sav_norm = cash_on_hand - simC(t, :);
    simCash(t, :) = cash_on_hand;
    % 更新下一期期初财富
    if t < cS.tn
        portfolio_return_next = (1-simA(t, :)) .* cS.rf + simA(t, :) .* simR(t, :);
        simW_norm(t+1, :) = portfolio_return_next .* liquid_sav_norm ./ simGPY(t+1, :);
    end
end
simW_norm = simW_norm(1:cS.tn, :);

% --- 3. 结果汇总与绘图 ---
fprintf('  汇总归一化结果并绘图...\n');
meanC = mean(simC, 2);
meanW = mean(simW_norm, 2);
meanA = mean(simA, 2);
meanY = zeros(cS.tn, 1);
work_periods = cS.tr - cS.tb; % Redefine for clarity
meanY(1:work_periods) = mean((1-cS.tau_y) * simY_norm_shock(1:work_periods,:), 2);
meanY(work_periods+1:cS.tn) = cS.ret_fac;

ages = cS.tb:cS.td;

figure('Position', [100, 100, 1200, 800], 'Name', 'A2C Simulation Results');

subplot(2,2,1);
plot(ages, meanW, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Liquid Wealth (w)');
hold on;
plot(ages, meanC, 'r-', 'LineWidth', 2, 'DisplayName', 'Consumption (c)');
plot(ages, meanY, 'k:', 'LineWidth', 2, 'DisplayName', 'Income (y)');
xline(cS.tr, '--', 'Retirement', 'LineWidth', 1.5);
hold off;
title('Mean Lifecycle Profiles (Normalized by P_t)'); xlabel('Age');
legend('show', 'Location', 'northwest'); grid on; xlim([cS.tb, cS.td]);

subplot(2,2,2);
plot(ages, meanA, 'm-', 'LineWidth', 2, 'DisplayName', 'Risky Share (\alpha)');
hold on;
xline(cS.tr, '--', 'Retirement', 'LineWidth', 1.5);
hold off;
title('Mean Portfolio Choices'); xlabel('Age'); ylabel('Share');
legend('show', 'Location', 'best'); grid on; xlim([cS.tb, cS.td]); ylim([-0.05, 1.05]);

subplot(2,2,4);
simY_for_ratio = zeros(size(simY_norm_shock));
simY_for_ratio(1:work_periods,:) = (1-cS.tau_y) * simY_norm_shock(1:work_periods,:);
simY_for_ratio(work_periods+1:end,:) = cS.ret_fac;
simWY_norm = simW_norm ./ simY_for_ratio;
simWY_norm(isinf(simWY_norm) | isnan(simWY_norm) | simY_for_ratio < 1e-6) = NaN;
meanWY_norm = nanmean(simWY_norm, 2);
plot(ages, meanWY_norm, 'LineWidth', 2);
xline(cS.tr, '--', 'Retirement', 'LineWidth', 1.5);
title('Mean Wealth-to-Income Ratio (w/y)'); xlabel('Age');
grid on; xlim([cS.tb, cS.td]);

sgtitle('Lifecycle Simulation Results (A2C Deep Reinforcement Learning)', 'FontSize', 16, 'FontWeight', 'bold');

output_filename = fullfile(cS.save_path, 'drl_simulation_results_a2c.png');
print(gcf, output_filename, '-dpng', '-r300');
fprintf('模拟图形已保存到 %s\n', output_filename);
end


% =========================================================================
%                辅助函数
% =========================================================================
function loss = huber(x, delta)
% Huber loss function.
if nargin < 2, delta = 1; end
abs_x = abs(x);
loss = sum(0.5 * (abs_x.^2) .* (abs_x <= delta) + delta * (abs_x - 0.5 * delta) .* (abs_x > delta));
end


function [z_grid, P] = discretizeAR1_Tauchen(mew, rho, sigma, znum, Tauchen_q, ~)
% Tauchen方法离散化AR(1)过程
if znum == 1
    z_grid = mew / (1 - rho);
    P = 1;
    return;
end

zstar = mew / (1 - rho);
sigmaz = sigma / sqrt(1 - rho^2);

z_grid = zstar + linspace(-Tauchen_q * sigmaz, Tauchen_q * sigmaz, znum)';
omega = z_grid(2) - z_grid(1);

P = zeros(znum, znum);
for i = 1:znum
    for j = 1:znum
        if j == 1
            P(i, j) = normcdf((z_grid(j) + omega/2 - rho * z_grid(i) - mew) / sigma);
        elseif j == znum
            P(i, j) = 1 - normcdf((z_grid(j) - omega/2 - rho * z_grid(i) - mew) / sigma);
        else
            P(i, j) = normcdf((z_grid(j) + omega/2 - rho * z_grid(i) - mew) / sigma) - ...
                normcdf((z_grid(j) - omega/2 - rho * z_grid(i) - mew) / sigma);
        end
    end
end
end