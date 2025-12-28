function val = fun_valuefunc_retired(policy, W, F, gamma)
% 退休期值函数计算（策略为比例形式）
% 输入:
%   policy: 策略向量 [C, A]
%     C - 消费比例（占总可用财富）
%     A - 风险资产投资比例（在第3期无意义）
%   W: 当前现金财富
%   F: 养老金收入（在退休期直接作为收入）
%   gamma: 相对风险厌恶系数
%
% 输出:
%   val: 当前状态下的值函数值

% 从策略向量中提取决策变量
c_ratio = policy(1);  % 消费比例

% 快速约束检查
if c_ratio > 1 || c_ratio < 0
    val = -1e10;
    return;
end

% 计算总可用财富和消费量
total_wealth = W + F;
C = c_ratio * total_wealth;

% 计算效用值
if C <= 0
    val = -1e10;
else
    val = (C^(1-gamma)) / (1-gamma);
end

% 注意：在退休期第三期（模型最后一期），没有未来期望效用
% 因此值函数仅由当期效用决定

end 