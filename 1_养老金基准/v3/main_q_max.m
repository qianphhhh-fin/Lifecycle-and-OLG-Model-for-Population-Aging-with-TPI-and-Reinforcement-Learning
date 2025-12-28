%% 批量求解不同 Q_max 设定下的生命周期模型
%
% 描述:
%   本脚本旨在生成一系列在不同个人养老金年度缴费上限(Q_max)设定下的
%   模型求解结果(VFI results)。脚本会循环遍历一组预定义的、以人民币计价的
%   缴费上限，为每个值调用VFI求解器，并将结果分别保存到不同的 .mat 文件中。
%
%   这些生成的文件将作为后续进行效用成本分析(utility cost analysis)的
%   基础数据。
%
%   本脚本假设 Q_max=0 (对应 vfi_results_nopps.mat) 和
%   Q_max=99999 (对应 vfi_results.mat) 的结果已经存在，因此不会重复计算。
%
% 输出:
%   - result/vfi_results_qmax_6000.mat
%   - result/vfi_results_qmax_12000.mat
%   - result/vfi_results_qmax_24000.mat
%   - result/vfi_results_qmax_36000.mat
%   (以及对应的 cS 结构体)

% --- 初始化 ---
clear;
close all;
clc;
warning off;

% --- 1. 定义要进行敏感性分析的 Q_max 值 (人民币绝对值) ---
% 我们选择当前政策(12000元)及其倍数作为分析点，以评估政策变化的潜在影响。
absolute_q_max_yuan_vec = [3000, 6000, 12000, 24000, 36000];  % 

% --- 2. 循环求解每个 Q_max 设定 ---
for i = 1:length(absolute_q_max_yuan_vec)

    current_q_max_abs = absolute_q_max_yuan_vec(i);
    fprintf('\n\n======================================================\n');
    fprintf('===== 开始处理 Q_max = %d 元的设定 =====\n', current_q_max_abs);
    fprintf('======================================================\n');

    % --- a. 获取基础参数 ---
    % 每次循环都重新获取，以确保 cS 是一个干净的基准。
    cS = utils_pps.setup_parameters();

    % --- b. 计算并覆盖当前循环的归一化 Q_max ---
    % 注意: 这里假设 calculate_q_max 也是 utils_pps 类的一个静态方法。
    % 如果不是，请确保该函数在MATLAB路径上。
    % 锚定收入需要与 setup_parameters 中的设定保持一致。
    anchor_income_yuan = 60000;
    try
        cS.Q_max = utils_pps.calculate_q_max(current_q_max_abs, anchor_income_yuan, cS.f_y);
    catch ME
        error('无法调用 utils_pps.calculate_q_max。请确保该辅助函数已移入 utils_pps.m 静态类中。错误信息: %s', ME.message);
    end

    % --- c. 定义当前设定的文件保存路径 ---
    save_path = cS.save_path;

    file_name = sprintf('vfi_results_qmax_%d.mat', current_q_max_abs);
    if cS.Q_max==0
        file_name = 'vfi_results_nopps.mat';
    end
    file_path = fullfile(save_path, file_name);

    % --- d. 求解或加载模型 ---
    fprintf('\n--- 步骤 2: 求解模型 (VFI) for Q_max = %d ---\n', current_q_max_abs);

    if ~exist(file_path, 'file')
        fprintf('未找到结果文件，开始进行值函数迭代...\n');
        vfi_results = utils_pps.solve_model_vfi(cS);

        if ~exist(save_path, 'dir')
            mkdir(save_path);
        end

        save(file_path, 'vfi_results', 'cS');
        fprintf('模型求解结果已成功保存到: %s\n', file_path);
    else
        fprintf('已在路径 %s 找到对应的结果文件，跳过计算。\n', file_path);
    end
end

fprintf('\n\n所有指定的 Q_max 值的模型求解已全部完成。\n');