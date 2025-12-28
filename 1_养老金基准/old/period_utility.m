function u = period_utility(consumption, gamma)
% PERIOD_UTILITY 计算每期消费的效用
% 输入:
%   consumption - 消费水平
%   gamma - 风险规避系数
% 输出:
%   u - 效用值

if consumption <= 0
    u = -1e10;  % 对于零或负消费，返回一个很大的负值
    return;
end

if gamma == 1
    u = log(consumption);  % 对数效用
else
    u = (consumption^(1-gamma)) / (1-gamma);  % CRRA效用函数
end

end 