%% 养老金模型敏感性分析脚本
% 此脚本自动运行不同参数组合的模拟，用于敏感性分析
% 使用方法: matlab -nosplash -nodesktop -r "sensitivity_analysis;exit"

clear;
close all;

% 设置所有参数值
retirement_ages = [61, 62, 63, 64, 65];
ret_fac_values = [0, 0.2, 0.4, 0.6, 0.8];
rp_values = [1.00, 1.02, 1.04, 1.06, 1.08];
gamma_values = [2.0, 3.0, 4.0, 5.0, 7.0];
beta_values = [0.92, 0.94, 0.95, 0.96, 0.98];
sigr_values = [0.15, 0.20, 0.27, 0.35, 0.45];
smav_values = [sqrt(0), sqrt(0.2), sqrt(0.4), sqrt(0.5), sqrt(0.8)];
smay_values = [sqrt(0), sqrt(0.2), sqrt(0.4), sqrt(0.5), sqrt(0.8)];
corr_values = [0.1, 0.2, 0.5, 0.8];

% 计算总模拟次数
total_sims = length(retirement_ages) + length(ret_fac_values) + length(rp_values) + ...
             length(gamma_values) + length(beta_values) + length(sigr_values) + ...
             length(smav_values) + length(smay_values) + length(corr_values);

% 添加日志和结果保存目录
root_dir = 'result/sensitivity';
if ~exist(root_dir, 'dir')
    mkdir(root_dir);
end

% 创建日志文件
% logFile = fullfile(root_dir, sprintf('sensitivity_log_%s.txt', datestr(now, 'yyyymmdd_HHMMSS')));
% diary(logFile);

% 记录分析开始时间
fprintf('敏感性分析开始时间: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

% 创建参数记录文件
% params_file = fullfile(root_dir, 'parameter_summary.csv');
% fid = fopen(params_file, 'w');
% fprintf(fid, '分析ID,参数类别,参数名称,参数值,保存路径,运行时间(秒),完成状态\n');
% fclose(fid);

save_path = root_dir;

% 设置默认状态为成功
status = '成功';

% 初始化进度计数器和总时间
sim_count = 0;
total_start_time = tic;

% 更新进度显示函数
function update_progress(sim_count, total_sims, total_start_time)
    elapsed_time = toc(total_start_time);
    fprintf('\r总进度: %d/%d (%.1f%%) - 已用时间: %.1f分钟', ...
        sim_count, total_sims, sim_count/total_sims*100, elapsed_time/60);
end

%% 1. 模拟不同退休年龄 (tr)
fprintf('\n========== 退休年龄敏感性分析 ==========\n');
param_category = '退休年龄';

for i = 1:length(retirement_ages)
    tr = retirement_ages(i);
    analysis_id = sprintf('TR_%d', tr);
    param_name = 'tr';
    param_value = tr;
    
    % 设置保存路径
    filename = sprintf('retirement_age_tr_%d.mat', tr);
    
    fprintf('运行模拟: 退休年龄 = %d\n', tr);
    fprintf('结果将保存到: %s/%s\n', save_path, filename);
    
    % 记录开始时间
    start_time = tic;
    
    % try
        % 运行模型
        batch_baseline_VFI('tr', tr, 'save_path', save_path, 'filename', filename);

        % 计算运行时间
        run_time = toc(start_time);
        % status = '成功';
    % catch ME
    %     run_time = toc(start_time);
    %     status = sprintf('失败: %s', ME.message);
    %     fprintf('错误: %s\n', ME.message);
    % end
    
    % 记录参数到CSV
    % fid = fopen(params_file, 'a');
    % fprintf(fid, '%s,%s,%s,%.4f,%s,%.2f,%s\n', ...
    %     analysis_id, param_category, param_name, param_value, save_path, run_time, status);
    % fclose(fid);
    
    fprintf('完成模拟: 退休年龄 = %d, 耗时 %.2f 秒 (%.2f 分钟)\n\n', ...
        tr, run_time, run_time/60);
    
    % 更新进度
    sim_count = sim_count + 1;
    update_progress(sim_count, total_sims, total_start_time);
end

%% 2. 模拟不同基本养老金水平 (ret_fac)
fprintf('\n========== 基本养老金水平敏感性分析 ==========\n');
param_category = '基本养老金';

for i = 1:length(ret_fac_values)
    ret_fac = ret_fac_values(i);
    analysis_id = sprintf('RF_%.1f', ret_fac);
    param_name = 'ret_fac';
    param_value = ret_fac;
    
    % 设置保存路径
    filename = sprintf('basic_pension_ret_fac_%.1f.mat', ret_fac);
    
    fprintf('运行模拟: 基本养老金 = %.1f\n', ret_fac);
    fprintf('结果将保存到: %s/%s\n', save_path, filename);
    
    % 记录开始时间
    start_time = tic;
    
    % try
        % 运行模型
        batch_baseline_VFI('ret_fac', ret_fac, 'save_path', save_path, 'filename', filename);

        % 计算运行时间
        run_time = toc(start_time);
        % status = '成功';
    % catch ME
    %     run_time = toc(start_time);
    %     status = sprintf('失败: %s', ME.message);
    %     fprintf('错误: %s\n', ME.message);
    % end
    
    % 记录参数到CSV
    % fid = fopen(params_file, 'a');
    % fprintf(fid, '%s,%s,%s,%.4f,%s,%.2f,%s\n', ...
    %     analysis_id, param_category, param_name, param_value, save_path, run_time, status);
    % fclose(fid);
    
    fprintf('完成模拟: 基本养老金 = %.1f, 耗时 %.2f 秒 (%.2f 分钟)\n\n', ...
        ret_fac, run_time, run_time/60);
    
    % 更新进度
    sim_count = sim_count + 1;
    update_progress(sim_count, total_sims, total_start_time);
end

%% 3. 模拟不同个人养老金收益率 (rp)
fprintf('\n========== 个人养老金收益率敏感性分析 ==========\n');
param_category = '养老金收益率';

for i = 1:length(rp_values)
    rp = rp_values(i);
    analysis_id = sprintf('RP_%.2f', rp);
    param_name = 'rp';
    param_value = rp;
    
    % 设置保存路径
    filename = sprintf('pension_return_rp_%.2f.mat', rp);
    
    fprintf('运行模拟: 养老金收益率 = %.2f\n', rp);
    fprintf('结果将保存到: %s/%s\n', save_path, filename);
    
    % 记录开始时间
    start_time = tic;
    
    % try
        % 运行模型
        batch_baseline_VFI('rp', rp, 'save_path', save_path, 'filename', filename);

        % 计算运行时间
        run_time = toc(start_time);
        % status = '成功';
    % catch ME
    %     run_time = toc(start_time);
    %     status = sprintf('失败: %s', ME.message);
    %     fprintf('错误: %s\n', ME.message);
    % end
    
    % 记录参数到CSV
    % fid = fopen(params_file, 'a');
    % fprintf(fid, '%s,%s,%s,%.4f,%s,%.2f,%s\n', ...
    %     analysis_id, param_category, param_name, param_value, save_path, run_time, status);
    % fclose(fid);
    
    fprintf('完成模拟: 养老金收益率 = %.2f, 耗时 %.2f 秒 (%.2f 分钟)\n\n', ...
        rp, run_time, run_time/60);
    
    % 更新进度
    sim_count = sim_count + 1;
    update_progress(sim_count, total_sims, total_start_time);
end

%% 4. 模拟不同风险厌恶系数 (gamma)
fprintf('\n========== 风险厌恶系数敏感性分析 ==========\n');
param_category = '风险厌恶';

for i = 1:length(gamma_values)
    gamma = gamma_values(i);
    analysis_id = sprintf('GAMMA_%.1f', gamma);
    param_name = 'gamma';
    param_value = gamma;
    
    % 设置保存路径
    filename = sprintf('risk_aversion_gamma_%.1f.mat', gamma);
    
    fprintf('运行模拟: 风险厌恶系数 = %.1f\n', gamma);
    fprintf('结果将保存到: %s/%s\n', save_path, filename);
    
    % 记录开始时间
    start_time = tic;
    
    % try
        % 运行模型
        batch_baseline_VFI('gamma', gamma, 'save_path', save_path, 'filename', filename);

        % 计算运行时间
        run_time = toc(start_time);
        % status = '成功';
    % catch ME
    %     run_time = toc(start_time);
    %     status = sprintf('失败: %s', ME.message);
    %     fprintf('错误: %s\n', ME.message);
    % end
    
    % 记录参数到CSV
    % fid = fopen(params_file, 'a');
    % fprintf(fid, '%s,%s,%s,%.4f,%s,%.2f,%s\n', ...
    %     analysis_id, param_category, param_name, param_value, save_path, run_time, status);
    % fclose(fid);
    
    fprintf('完成模拟: 风险厌恶系数 = %.1f, 耗时 %.2f 秒 (%.2f 分钟)\n\n', ...
        gamma, run_time, run_time/60);
    
    % 更新进度
    sim_count = sim_count + 1;
    update_progress(sim_count, total_sims, total_start_time);
end

%% 5. 模拟不同时间偏好系数 (beta)
fprintf('\n========== 时间偏好系数敏感性分析 ==========\n');
param_category = '时间偏好';

for i = 1:length(beta_values)
    beta = beta_values(i);
    analysis_id = sprintf('BETA_%.2f', beta);
    param_name = 'beta';
    param_value = beta;
    
    % 设置保存路径
    filename = sprintf('time_preference_beta_%.2f.mat', beta);
    
    fprintf('运行模拟: 时间偏好系数 = %.2f\n', beta);
    fprintf('结果将保存到: %s/%s\n', save_path, filename);
    
    % 记录开始时间
    start_time = tic;
    
    % try
        % 运行模型
        batch_baseline_VFI('beta', beta, 'save_path', save_path, 'filename', filename);

        % 计算运行时间
        run_time = toc(start_time);
        % status = '成功';
    % catch ME
    %     run_time = toc(start_time);
    %     status = sprintf('失败: %s', ME.message);
    %     fprintf('错误: %s\n', ME.message);
    % end
    
    % 记录参数到CSV
    % fid = fopen(params_file, 'a');
    % fprintf(fid, '%s,%s,%s,%.4f,%s,%.2f,%s\n', ...
    %     analysis_id, param_category, param_name, param_value, save_path, run_time, status);
    % fclose(fid);
    
    fprintf('完成模拟: 时间偏好系数 = %.2f, 耗时 %.2f 秒 (%.2f 分钟)\n\n', ...
        beta, run_time, run_time/60);
    
    % 更新进度
    sim_count = sim_count + 1;
    update_progress(sim_count, total_sims, total_start_time);
end

%% 7. 模拟风险资产波动率变化 (sigr)
fprintf('\n========== 风险资产波动率敏感性分析 ==========\n');
param_category = '风险波动率';

for i = 1:length(sigr_values)
    sigr = sigr_values(i);
    analysis_id = sprintf('SIGR_%.2f', sigr);
    param_name = 'sigr';
    param_value = sigr;
    
    % 设置保存路径
    filename = sprintf('risk_volatility_sigr_%.2f.mat', sigr);
    
    fprintf('运行模拟: 风险资产波动率 = %.2f\n', sigr);
    fprintf('结果将保存到: %s/%s\n', save_path, filename);
    
    % 记录开始时间
    start_time = tic;
    
    % try
        % 运行模型
        batch_baseline_VFI('sigr', sigr, 'save_path', save_path, 'filename', filename);

        % 计算运行时间
        run_time = toc(start_time);
        % status = '成功';
    % catch ME
    %     run_time = toc(start_time);
    %     status = sprintf('失败: %s', ME.message);
    %     fprintf('错误: %s\n', ME.message);
    % end
    
    % 记录参数到CSV
    % fid = fopen(params_file, 'a');
    % fprintf(fid, '%s,%s,%s,%.4f,%s,%.2f,%s\n', ...
    %     analysis_id, param_category, param_name, param_value, save_path, run_time, status);
    % fclose(fid);
    
    fprintf('完成模拟: 风险资产波动率 = %.2f, 耗时 %.2f 秒 (%.2f 分钟)\n\n', ...
        sigr, run_time, run_time/60);
    
    % 更新进度
    sim_count = sim_count + 1;
    update_progress(sim_count, total_sims, total_start_time);
end

%% 8. 模拟不同smav值
fprintf('\n========== 持续冲击smav敏感性分析 ==========\n');
param_category = 'smav';

for i = 1:length(smav_values)
    smav = smav_values(i);
    analysis_id = sprintf('SMAV_%.3f', smav);
    param_name = 'smav';
    param_value = smav;
    
    % 设置保存路径
    filename = sprintf('smav_%.3f.mat', smav);
    
    fprintf('运行模拟: smav = %.3f\n', smav);
    fprintf('结果将保存到: %s/%s\n', save_path, filename);
    
    % 记录开始时间
    start_time = tic;
    
    % try
        % 运行模型
        batch_baseline_VFI('smav', smav, 'save_path', save_path, 'filename', filename);

        % 计算运行时间
        run_time = toc(start_time);
        % status = '成功';
    % catch ME
    %     run_time = toc(start_time);
    %     status = sprintf('失败: %s', ME.message);
    %     fprintf('错误: %s\n', ME.message);
    % end
    
    % 记录参数到CSV
    % fid = fopen(params_file, 'a');
    % fprintf(fid, '%s,%s,%s,%.4f,%s,%.2f,%s\n', ...
    %     analysis_id, param_category, param_name, param_value, save_path, run_time, status);
    % fclose(fid);
    
    fprintf('完成模拟: smav = %.3f, 耗时 %.2f 秒 (%.2f 分钟)\n\n', ...
        smav, run_time, run_time/60);
    
    % 更新进度
    sim_count = sim_count + 1;
    update_progress(sim_count, total_sims, total_start_time);
end

%% 9. 模拟不同smay值
fprintf('\n========== 临时性冲击smay敏感性分析 ==========\n');
param_category = 'smay';

for i = 1:length(smay_values)
    smay = smay_values(i);
    analysis_id = sprintf('SMAY_%.3f', smay);
    param_name = 'smay';
    param_value = smay;
    
    % 设置保存路径
    filename = sprintf('smay_%.3f.mat', smay);
    
    fprintf('运行模拟: smay = %.3f\n', smay);
    fprintf('结果将保存到: %s/%s\n', save_path, filename);
    
    % 记录开始时间
    start_time = tic;
    
    % try
        % 运行模型
        batch_baseline_VFI('smay', smay, 'save_path', save_path, 'filename', filename);

        % 计算运行时间
        run_time = toc(start_time);
        % status = '成功';
    % catch ME
    %     run_time = toc(start_time);
    %     status = sprintf('失败: %s', ME.message);
    %     fprintf('错误: %s\n', ME.message);
    % end
    
    % 记录参数到CSV
    % fid = fopen(params_file, 'a');
    % fprintf(fid, '%s,%s,%s,%.4f,%s,%.2f,%s\n', ...
    %     analysis_id, param_category, param_name, param_value, save_path, run_time, status);
    % fclose(fid);
    
    fprintf('完成模拟: smay = %.3f, 耗时 %.2f 秒 (%.2f 分钟)\n\n', ...
        smay, run_time, run_time/60);
    
    % 更新进度
    sim_count = sim_count + 1;
    update_progress(sim_count, total_sims, total_start_time);
end

%% 10. 模拟不同corr_z_epsilon值
fprintf('\n========== corr_z_epsilon敏感性分析 ==========\n');
param_category = 'corr_z_epsilon';

for i = 1:length(corr_values)
    corr = corr_values(i);
    analysis_id = sprintf('CORR_%.2f', corr);
    param_name = 'corr_z_epsilon';
    param_value = corr;
    
    % 设置保存路径
    filename = sprintf('corr_z_epsilon_%.2f.mat', corr);
    
    fprintf('运行模拟: corr_z_epsilon = %.2f\n', corr);
    fprintf('结果将保存到: %s/%s\n', save_path, filename);
    
    % 记录开始时间
    start_time = tic;
    
    % try
        % 运行模型
        batch_baseline_VFI('corr_z_epsilon', corr, 'save_path', save_path, 'filename', filename);

        % 计算运行时间
        run_time = toc(start_time);
        % status = '成功';
    % catch ME
    %     run_time = toc(start_time);
    %     status = sprintf('失败: %s', ME.message);
    %     fprintf('错误: %s\n', ME.message);
    % end
    
    % 记录参数到CSV
    % fid = fopen(params_file, 'a');
    % fprintf(fid, '%s,%s,%s,%.4f,%s,%.2f,%s\n', ...
    %     analysis_id, param_category, param_name, param_value, save_path, run_time, status);
    % fclose(fid);
    
    fprintf('完成模拟: corr_z_epsilon = %.2f, 耗时 %.2f 秒 (%.2f 分钟)\n\n', ...
        corr, run_time, run_time/60);
    
    % 更新进度
    sim_count = sim_count + 1;
    update_progress(sim_count, total_sims, total_start_time);
end

% 显示总运行时间
total_time = toc(total_start_time);
fprintf('\n\n敏感性分析完成！\n');
fprintf('总运行时间: %.2f 秒 (%.2f 分钟, %.2f 小时)\n', ...
    total_time, total_time/60, total_time/3600);

