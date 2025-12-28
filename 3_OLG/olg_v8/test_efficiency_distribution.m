% 测试修改后的劳动效率分布
clc; clear;

fprintf('=== 测试修改后的劳动效率分布 ===\n');

% 加载参数
cS = main_olg_v8_utils.ParameterValues_HuggettStyle();

% 生成劳动效率网格
[leLogGridV, leTrProbM, leProb1V] = main_olg_v8_utils.EarningProcess_olgm(cS);
leGridV = exp(leLogGridV(:));

fprintf('\n劳动效率状态网格值：\n');
for i = 1:length(leGridV)
    fprintf('状态 %d: %.4f\n', i, leGridV(i));
end

fprintf('\n初始分布概率：\n');
for i = 1:length(leProb1V)
    fprintf('状态 %d: %.4f\n', i, leProb1V(i));
end

fprintf('\n效率比较：\n');
fprintf('最高/最低效率比值: %.2f\n', leGridV(end)/leGridV(1));
fprintf('平均效率: %.4f\n', sum(leGridV .* leProb1V));

% 检查分布是否更合理
if leGridV(end)/leGridV(1) < 10
    fprintf('\n✓ 效率分布合理（最高最低比值 < 10）\n');
else
    fprintf('\n✗ 效率分布仍然过于极端（最高最低比值 >= 10）\n');
end 