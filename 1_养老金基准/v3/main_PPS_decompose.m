

% --- 初始化 ---
clear;
close all;
clc;
warning off

% --- 主流程控制器 ---
fprintf('开始执行生命周期模型...\n');

%%  无税收激励
% 1. 设置所有模型参数
cS = utils_pps_notax_incentive.setup_parameters();
save_path = cS.save_path;
file_path = fullfile(save_path, 'breakdown','vfi_results_notaxincentive.mat');

% 2. 求解模型 (值函数迭代)
fprintf('\n===== 步骤 2: 求解模型 (VFI) =====\n');
if ~exist(file_path,'file')
    vfi_results = utils_pps_notax_incentive.solve_model_vfi(cS);
    if ~exist(save_path, 'dir'), mkdir(save_path); end
    save(file_path, 'vfi_results', 'cS');
    fprintf('模型求解结果已保存到 %s\n', file_path);
else
    load(file_path)
end

%%  完全流动性个人养老金
% 1. 设置所有模型参数
cS = utils_pps_fullliquidity.setup_parameters();
save_path = cS.save_path;
file_path = fullfile(save_path, 'breakdown','vfi_results_fullliquidity.mat');

% 2. 求解模型 (值函数迭代)
fprintf('\n===== 步骤 2: 求解模型 (VFI) =====\n');
if ~exist(file_path,'file')
    vfi_results = utils_pps_fullliquidity.solve_model_vfi(cS);
    if ~exist(save_path, 'dir'), mkdir(save_path); end
    save(file_path, 'vfi_results', 'cS');
    fprintf('模型求解结果已保存到 %s\n', file_path);
else
    load(file_path)
end


%% 无收入不确定性

% 1. 设置所有模型参数
cS = utils_pps.setup_parameters('test');
cS.smav = 0;
cS.smay = 0;
cS = utils_pps.setup_process(cS);
save_path = cS.save_path;
file_path = fullfile(save_path, 'breakdown','vfi_results_noincomerisk.mat');

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



