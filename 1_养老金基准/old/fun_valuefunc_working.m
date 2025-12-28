function val = fun_valuefunc_working(policy, W, F, next_V, Y_grid, R_grid, joint_prob, Rf, Rp, gamma, beta, gW, gF, t)
% 工作期值函数计算（策略为比例形式）
% 输入:
%   policy: 策略向量 [C, A, Q]
%       C - 消费比例（占当前财富）
%       A - 风险资产投资比例（占总投资）
%       Q - 养老金购买比例（占当前财富）
%   W: 当前现金财富
%   F: 当前养老金余额
%   next_V: 下一期值函数矩阵
%   Y_grid: 劳动收入网格
%   R_grid: 风险资产收益率网格
%   joint_prob: 联合概率
%   Rf: 无风险收益率
%   Rp: 养老金账户收益率
%   gamma: 相对风险厌恶系数
%   beta: 折现因子
%   gW: 现金财富网格
%   gF: 养老金余额网格
%   t: 当前时期
%
% 输出:
%   val: 当前状态下的值函数值

% 从策略向量中提取决策变量
c_ratio = policy(1);  % 消费比例
a_ratio = policy(2);  % 风险资产投资比例
q_ratio = policy(3);  % 养老金购买比例

% 检查比例约束
if c_ratio + a_ratio + q_ratio > 1 || c_ratio < 0 || a_ratio < 0 || q_ratio < 0
    val = -1e10;  % 如果约束不满足，返回很大的负值
    return;
end

% 计算绝对值
C = c_ratio * W;  % 消费绝对值
Q = q_ratio * W;  % 养老金购买绝对值

% 当期效用
if C <= 0
    val = -1e10;  % 避免非正消费
    return;
end
current_utility = (C^(1-gamma)) / (1-gamma);

% 计算投资金额和下一期养老金余额
S = W - C - Q;  % 总投资金额
risk_invest = a_ratio * W;  % 风险资产投资金额
safe_invest = W - C - Q - risk_invest;  % 无风险资产投资金额
F_next = (F + Q) * Rp;  % 下一期养老金余额

% 检查F_next是否在网格范围内
F_next = min(max(F_next, gF(1)), gF(end));

% 获取网格维度
[nY, nR] = size(joint_prob);

% 创建网格矩阵用于vectorized计算
[Y_mat, R_mat] = meshgrid(Y_grid, R_grid);
prob_mat = joint_prob';  % 转置以匹配meshgrid的维度结构

% 计算所有可能的未来状态的现金财富（矩阵运算）
W_next_mat = risk_invest * R_mat(:) + safe_invest * Rf + Y_mat(:);

% 限制W_next在网格范围内
W_next_mat = min(max(W_next_mat, gW(1)), gW(end));

% 使用interp2进行网格插值
% 创建用于插值的坐标和值
[gW_mat, gF_mat] = meshgrid(gW, gF);
V_mat = next_V';  % 转置以匹配meshgrid的维度结构

% 为所有可能的未来状态计算值函数
V_next_all = interp2(gW_mat, gF_mat, V_mat, W_next_mat, F_next * ones(size(W_next_mat)), 'linear', 0);

% 重塑为与joint_prob相匹配的尺寸
V_next_all = reshape(V_next_all, size(prob_mat));

% 计算期望值函数（矩阵计算）
expected_utility = sum(prob_mat(:) .* V_next_all(:));

% 计算总值函数值
val = current_utility + beta * expected_utility;

end 