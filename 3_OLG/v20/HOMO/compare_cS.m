% =========================================================================
% == SCRIPT: compare_cS.m
% == 版本: [v1.1 - 全字段自动扫描版]
% ==
% == 目的:
% ==   - [核心修正] 不再使用预定义的字段列表，而是自动获取一个cS
% ==     结构体的所有字段名，并以此为基础进行全面比较。
% ==   - 能够检测并报告在一个情景中存在但在另一个中不存在的参数。
% ==   - 使诊断过程更完整、更具鲁棒性。
% =========================================================================

clear; close all; clc;
fprintf('=== cS 参数对比诊断工具 (v1.1 - 全字段扫描版) ===\n\n');

%% --- 1. 加载数据 ---
fprintf('--- 1. 加载数据 ---\n');
try
    data_nopps = load('SS/data_for_transition_nopps.mat', 'data_for_transition');
    cS_nopps = data_nopps.data_for_transition.cS;
    fprintf('   ✅ 已加载 "无PPS" 情景参数。\n');

    data_pps = load('SS/data_for_transition_pps.mat', 'data_for_transition');
    cS_pps = data_pps.data_for_transition.cS;
    fprintf('   ✅ 已加载 "有PPS" 情景参数。\n\n');
catch ME
    error('加载数据失败，请确保已成功运行两个情景的 main_run_ss.m。错误: %s', ME.message);
end

%% --- 2. [核心修正] 获取全字段列表并进行比较 ---
% 以 cS_nopps 的字段为基准，并合并 cS_pps 中独有的字段
fields_nopps = fieldnames(cS_nopps);
fields_pps = fieldnames(cS_pps);
all_fields = union(fields_nopps, fields_pps, 'stable'); % 'stable' 保持原有顺序

fprintf('--- 2. 参数比较结果 ---\n');
fprintf('%-30s | %-30s | %-30s | %-12s\n', '参数字段', '无PPS (nopps)', '有PPS (pps)', '是否相同?');
fprintf('%s\n', repmat('-', 1, 110));

for i = 1:numel(all_fields)
    field = all_fields{i};
    
    val_nopps_str = '(字段不存在)';
    val_pps_str = '(字段不存在)';
    is_same_str = 'N/A';
    
    has_nopps = isfield(cS_nopps, field);
    has_pps = isfield(cS_pps, field);
    
    if has_nopps
        val_nopps = cS_nopps.(field);
        % 对大型数组或复杂类型进行特殊处理，避免命令行刷屏
        if numel(val_nopps) > 10 && ~ischar(val_nopps)
             val_nopps_str = sprintf('[%dx%d %s]', size(val_nopps,1), size(val_nopps,2), class(val_nopps));
        else
             val_nopps_str = mat2str(val_nopps, 4);
        end
    end
    
    if has_pps
        val_pps = cS_pps.(field);
        if numel(val_pps) > 10 && ~ischar(val_pps)
            val_pps_str = sprintf('[%dx%d %s]', size(val_pps,1), size(val_pps,2), class(val_pps));
        else
            val_pps_str = mat2str(val_pps, 4);
        end
    end
    
    if has_nopps && has_pps
        % 使用 isequaln 来比较，可以处理 NaN
        if isequaln(cS_nopps.(field), cS_pps.(field))
            is_same_str = '是';
        else
            is_same_str = '*** 否 ***';
        end
    else
        % 如果一个存在一个不存在，那肯定不同
        is_same_str = '*** 否 ***';
    end
    
    fprintf('%-30s | %-30s | %-30s | %-12s\n', field, val_nopps_str, val_pps_str, is_same_str);
end
fprintf('%s\n', repmat('-', 1, 110));

fprintf('\n--- 诊断脚本执行完毕 ---\n');