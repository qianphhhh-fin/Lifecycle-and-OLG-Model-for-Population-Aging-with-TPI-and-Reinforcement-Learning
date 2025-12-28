%% test_env_rewards.m
% 描述:
% 这是一个用于测试 CoccoEnv 环境奖励分布的独立脚本。
% 它通过执行一个完全随机的策略来模拟一个未经训练的智能体，
% 并分析其获得的单步奖励（效用）的分布。
%
% 目的:
% 1. 检查奖励是否存在极端的负值，这可能阻碍学习。
% 2. 了解奖励的尺度和方差，为RL训练的稳定性提供参考。
% 3. 验证环境在随机输入下的行为是否符合经济学直觉。

% --- 初始化 ---
clear;
close all;
clc;
rng('default'); % 使用默认随机种子

fprintf('开始测试 CoccoEnv 的奖励分布...\n');

% --- 1. 参数设置 ---
N_EPISODES = 2000; % 要模拟的生命周期数量
DISABLE_DEATH = true; % 在有随机死亡的真实环境中测试

% --- 2. 创建环境实例 ---
env = CoccoEnv('disable_random_death', DISABLE_DEATH);
tn = env.tn;
% env = NormalizedCoccoEnv(env, 'norm_obs_idx', [1]); % 只归一化财富和永久收入
actInfo = getActionInfo(env);
low = actInfo.LowerLimit;
low(1)=0.4;
high = actInfo.UpperLimit;

% --- 3. 运行模拟并收集奖励 ---
all_rewards = []; % 用于存储所有时间步的奖励
fprintf('正在使用随机策略模拟 %d 个生命周期...\n', N_EPISODES);

h = waitbar(0, '正在模拟...');
Y = ones(N_EPISODES,tn-1);
for i_episode = 1:N_EPISODES
    % 重置环境
    [obs,info] = env.reset();
    is_done = false;
    ii=1;
    Y(i_episode,ii) = info.AbsoluteIncome;
    W(i_episode,ii) = obs(1);
    P(i_episode,ii) = obs(2);
    while ~is_done
        % a. 生成一个完全随机的动作
        % 动作是 [消费比例, 风险资产比例]
        random_action = low + rand(size(low)) .* (high - low);
        
        % b. 在环境中执行动作
        [obs, reward, is_done, info] = env.step(random_action);
        
        % c. 存储奖励值
        all_rewards(end+1) = reward;
        ii = ii + 1;
        Y(i_episode,ii) = info.AbsoluteIncome_next;
        W(i_episode,ii) = obs(1);
        P(i_episode,ii) = obs(2);
    end
    ii=0;
    waitbar(i_episode / N_EPISODES, h);
end
close(h);

fprintf('模拟完成，共收集到 %d 个单步奖励。\n', length(all_rewards));

% --- 4. 分析与可视化奖励分布 ---

% a. 计算描述性统计量
mean_reward = mean(all_rewards);
std_reward = std(all_rewards);
min_reward = min(all_rewards);
max_reward = max(all_rewards);
median_reward = median(all_rewards);
percentile_1 = prctile(all_rewards, 1);  % 1%分位数
percentile_99 = prctile(all_rewards, 99);% 99%分位数

fprintf('\n--- 奖励分布统计 ---\n');
fprintf('均值 (Mean):      %.4f\n', mean_reward);
fprintf('标准差 (Std Dev): %.4f\n', std_reward);
fprintf('最小值 (Min):       %.4f\n', min_reward);
fprintf('最大值 (Max):       %.4f\n', max_reward);
fprintf('中位数 (Median):    %.4f\n', median_reward);
fprintf('1%% 分位数:        %.4f\n', percentile_1);
fprintf('99%% 分位数:       %.4f\n', percentile_99);
fprintf('-----------------------\n');

% b. 绘制直方图
figure('Position', [200, 200, 1000, 600]);
histogram(all_rewards, 100, 'Normalization', 'pdf');
hold on;

% 绘制均值和中位数的垂直线
xline(mean_reward, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Mean: %.2f', mean_reward));
xline(median_reward, 'g-', 'LineWidth', 2, 'DisplayName', sprintf('Median: %.2f', median_reward));

title('CoccoEnv: 单步奖励 (效用) 分布 (随机策略)');
xlabel('奖励值 (Utility)');
ylabel('概率密度');
legend('show');
grid on;
set(gca, 'YScale', 'log'); % 使用对数刻度以便更好地观察尾部
xlim([percentile_1 * 1.1, percentile_99 * 0.9]); % 限制X轴范围以关注核心分布

fprintf('\n测试完成。请查看生成的图表。\n');