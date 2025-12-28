%% 批量求解不同参数设定下的生命周期模型 (敏感性分析)
%
% 描述:
%   本脚本旨在为一系列关键模型参数生成VFI求解结果，用于后续的
%   敏感性分析。脚本会遍历每个指定参数的一组典型值，为每个值调用
%   VFI求解器，并将结果分别保存到 'result/sensitivity/' 文件夹下。
%
%   脚本会自动跳过基准参数值的计算，以避免重复工作。
%
% 分析的参数:
%   1. 永久性收入冲击标准差 (smav)
%   2. 相对风险规避系数 (gamma)
%   3. 风险资产波动率 (sigr)
%   4. 个人养老金账户回报率 (rp)
%
% 输出:
%   - result/sensitivity/vfi_results_smav_0p1000.mat
%   - result/sensitivity/vfi_results_gamma_2p0000.mat
%   - ... (以此类推)

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
    % 参数 1: 永久性收入冲击标准差 (smav)，基于方差设定
    struct('name', 'smav', ...
           'field', 'cS.smav', ...
           'values', [sqrt(0.02), sqrt(0.2), sqrt(0.4), sqrt(0.6)]), ...

    % 参数 2: 相对风险规避系数 (gamma)
    struct('name', 'gamma', ...
           'field', 'cS.gamma', ...
           'values', [2.0, 5.0, 7.5, 10.0]), ...

    % 参数 3: 风险资产波动率 (sigr)
    struct('name', 'sigr', ...
           'field', 'cS.sigr', ...
           'values', [0.01, 0.15, 0.35, 0.45]), ...

    % 参数 4: 个人养老金账户回报率 (rp)
    struct('name', 'rp', ...
           'field', 'cS.rp', ...
           'values', [1.01, 1.02,1.06, 1.08]),...
               
    % 参数 5: 退休年龄 (tr)           
    struct('name', 'tr', ...
           'field', 'cS.tr', ...
           'values', [58, 62, 63, 64, 65]),...
              
           % 参数 6: 个人所得税率           
    struct('name', 'tau_y', ...
           'field', 'cS.tau_y', ...
           'values', [0.03, 0.05, 0.07, 0.08, 0.1, 0.15])
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
        cS.ncash = 31;
        cS.nfund = 31;


        % --- b. 检查当前值是否为基准值，若是则跳过 ---
        % temp = strsplit(param_field, '.');
        % baseline_value = getfield(cS, temp{2} ); % 动态获取基准值
        % if abs(current_value - baseline_value) < 1e-9
        %     fprintf('当前值为此参数的基准值 (%f)，跳过计算。\n', baseline_value);
        %     continue;
        % end

        % --- c. 覆盖当前分析所需的参数值 ---
        % 使用 eval 动态修改 cS 结构体中对应的参数
        eval(sprintf('%s = %f;', param_field, current_value));


        % --- d. 定义文件保存路径和名称 ---
        save_dir = fullfile('result', 'sensitivity');
        if ~exist(save_dir, 'dir')
            mkdir(save_dir);
        end

        value_str = strrep(sprintf('%.4f', current_value), '.', 'p');
        file_name = sprintf('vfi_results_%s_%s.mat', param_name, value_str);
        file_path = fullfile(save_dir, file_name);

        % --- e. 求解或加载模型 ---
        fprintf('\n--- 正在为 %s = %f 求解模型 ---\n', param_name, current_value);
        cS = utils_pps.setup_process(cS); % 派生参数设置
        if ~exist(file_path, 'file')
            fprintf('未找到结果文件，开始进行值函数迭代...\n');
            vfi_results = utils_pps.solve_model_vfi(cS);

            % utils_pps.simulate_model(vfi_results, cS);

            save(file_path, 'vfi_results', 'cS');
            fprintf('模型求解结果已成功保存到: %s\n', file_path);
        else
            fprintf('已在路径 %s 找到对应的结果文件，跳过计算。\n', file_path);
        end
    end
end

fprintf('\n\n所有指定的敏感性分析已全部完成。\n');