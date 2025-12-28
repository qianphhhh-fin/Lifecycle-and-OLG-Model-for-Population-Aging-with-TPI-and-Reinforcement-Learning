% *************************************************************************
% *                    在这里插入下面的调试代码                           *
% *************************************************************************

% --- [开始] 策略路径一致性检验代码 ---
fprintf('\n--- [调试检查] 正在验证所有时期的策略路径(Pol_path)是否恒定... ---\n');

% 设置一个非常小的容忍度，以应对浮点数计算的微小差异
tolerance = 1e-8; 
all_policies_are_constant = true;

% 遍历每个家庭类型 h
for h = 1:nH
    fprintf('   --- 检查类型 h = %d ---\n', h);
    
    % 获取该类型的策略路径和终期稳态策略
    policy_path_h = Pol_path_h{h};
    policy_final_h = polF_h{h};
    
    % 遍历转轨路径上的每一个时期 t (从后向前检查更有启发性)
    for t = (T-1):-1:1
        policy_t_h = policy_path_h{t};
        
        % 获取策略结构体的所有字段名
        fields = fieldnames(policy_t_h);
        
        % 遍历每一个策略字段 (如 'c', 'k_prime' 等)
        for i = 1:length(fields)
            field = fields{i};
            
            % 遍历每一个年龄 a
            for a = 1:cS.aD_new
                % 提取 t 时期的策略矩阵和终期稳态的策略矩阵
                matrix_t = policy_t_h(a).(field);
                matrix_F = policy_final_h(a).(field);
                
                % 计算两个矩阵之间的最大绝对差值
                max_abs_diff = max(abs(matrix_t - matrix_F), [], 'all');
                
                % 检查差值是否超过容忍度
                if max_abs_diff > tolerance
                    fprintf('      [!!! 失败 !!!] 在 t=%3d, 类型 h=%d, 年龄 a=%2d, 策略 ''%s'' 偏离终态！\n', t, h, a, field);
                    fprintf('      最大绝对差异: %.4e\n', max_abs_diff);
                    
                    % (可选) 如果你想在发现第一个差异时就暂停，取消下面的注释
                    % keyboard; 
                    
                    all_policies_are_constant = false;
                    % 找到一个错误就足够了，可以跳出内层循环以提高效率
                    break; 
                end
            end
            if ~all_policies_are_constant, break; end
        end
        if ~all_policies_are_constant, break; end
    end
    if ~all_policies_are_constant, break; end
end

% 打印最终的检验结论
if all_policies_are_constant
    fprintf('\n   [✅ 成功] 所有家庭类型在所有时期的策略路径均与终期稳态策略保持一致。\n');
else
    fprintf('\n   [❌ 失败] 策略路径存在时变性！请检查上方报告的第一个偏离点。\n');
    % 强制程序在此处停止，以便你可以在命令窗口中检查变量
    error('策略路径一致性检验失败，终止执行。');
end
fprintf('--- [调试检查] 验证结束 ---\n\n');

% --- [结束] 策略路径一致性检验代码 ---