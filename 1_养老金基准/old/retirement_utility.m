function u = retirement_utility(x, w, y3, gamma)
% RETIREMENT_UTILITY 计算退休期效用
% 输入:
%   x - 决策变量 [消费比例, 风险资产配置比例]
%   w - 现金财富
%   y3 - 退休期收入（养老金账户余额）
%   gamma - 风险规避系数
% 输出:
%   u - 效用值

c_ratio = x(1);  % 消费比例
consumption = c_ratio * w;  % 实际消费金额

% 计算消费效用
u = period_utility(consumption, gamma);

end 