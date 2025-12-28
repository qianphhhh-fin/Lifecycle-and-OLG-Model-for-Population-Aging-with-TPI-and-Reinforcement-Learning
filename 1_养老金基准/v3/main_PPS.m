

% --- 初始化 ---
clear;
close all;
clc;
warning off

% --- 主流程控制器 ---
fprintf('开始执行生命周期模型...\n');

% 1. 设置所有模型参数
cS = utils_pps.setup_parameters();
save_path = cS.save_path;
file_path = fullfile(save_path, 'vfi_results_rikyrp.mat');

% 2. 求解模型 (值函数迭代)
fprintf('\n===== 步骤 2: 求解模型 (VFI) =====\n');
if ~exist(file_path,'file')
    vfi_results = utils_pps.solve_model_vfi(cS);
    if ~exist(save_path, 'dir'), mkdir(save_path); end
    save(file_path, 'vfi_results', 'cS');
    fprintf('模型求解结果已保存到 %s\n', file_path);
else
    load(file_path)
end

% 保存结果


% 3. 模拟模型
fprintf('\n===== 步骤 3: 模拟模型 =====\n');
utils_pps.simulate_model(vfi_results, cS);

fprintf('\n模型执行完毕。\n');

