%% 批量求解不同参数设定下的生命周期模型 (敏感性分析)
%
% 描述:
%   本脚本旨在为一系列关键模型参数生成VFI求解结果，用于后续的
%   敏感性分析。脚本会遍历每个指定参数的一组典型值，为每个值调用
%   VFI求解器，并将结果分别保存到 'result/riskyrp/' 文件夹下。
%
%   脚本会自动跳过基准参数值的计算，以避免重复工作。
%

%
% 输出:
%   - result/riskyrp/vfi_results_riskyrp_sigrrp_0p1000.mat
%   - result/riskyrp/vfi_results_riskyrp_sigrrp_0p2000.mat


% --- 初始化 ---
clear;
close all;
clc;
warning off;

% --- 1. 定义敏感性分析的配置 ---
% 每个cell元素代表一个要分析的参数
% 'name': 用于文件命名的简称
% 'field': 在cS结构体中对应的字段名(字符串形式)
% 'values': 需要进行分析的一系列数值
sensitivity_params = {
    % 参数 1: 
    struct('name', 'sigrrp', ...
           'field', 'cS.sigrrp', ...
           'values', [0, 0.05, 0.1, 0.2,0.27, 0.3]), ...
};

% --- 2. 循环遍历每个参数及其设定值 ---
for p_idx = 1:length(sensitivity_params)
    param_info = sensitivity_params{p_idx};
    param_name = param_info.name;
    param_field = param_info.field;
    param_values = param_info.values;

    for v_idx = 1:length(param_values)
        current_value = param_values(v_idx);

        fprintf('\n\n======================================================\n');
        fprintf('===== 开始敏感性分析 =====\n');
        fprintf('===== 参数: %s, 当前值: %f =====\n', upper(param_name), current_value);
        fprintf('======================================================\n');

        % --- a. 获取基准参数 ---
        % 每次都重新加载，以确保从一个干净的基准开始
        cS = utils_pps.setup_parameters('sensitivity');
        cS.ncash = 51;
        cS.nfund = 51;

        % --- c. 覆盖当前分析所需的参数值 ---
        % 使用 eval 动态修改 cS 结构体中对应的参数
        eval(sprintf('%s = %f;', param_field, current_value));


        % --- d. 定义文件保存路径和名称 ---
        save_dir = fullfile('result', 'riskyrp');
        if ~exist(save_dir, 'dir')
            mkdir(save_dir);
        end

        value_str = strrep(sprintf('%.4f', current_value), '.', 'p');
        file_name = sprintf('vfi_results_riskyrp_%s_%s.mat', param_name, value_str);
        file_path = fullfile(save_dir, file_name);

        % --- e. 求解或加载模型 ---
        fprintf('\n--- 正在为 %s = %f 求解模型 ---\n', param_name, current_value);
        cS = utils_pps.setup_process(cS); % 派生参数设置
        if ~exist(file_path, 'file')
            fprintf('未找到结果文件，开始进行值函数迭代...\n');
            vfi_results = utils_pps.solve_model_vfi_riskyrp(cS);

            % utils_pps.simulate_model(vfi_results, cS);

            save(file_path, 'vfi_results', 'cS');
            fprintf('模型求解结果已成功保存到: %s\n', file_path);
        else
            fprintf('已在路径 %s 找到对应的结果文件，跳过计算。\n', file_path);
        end
    end
end

fprintf('\n\n所有指定的敏感性分析已全部完成。\n');