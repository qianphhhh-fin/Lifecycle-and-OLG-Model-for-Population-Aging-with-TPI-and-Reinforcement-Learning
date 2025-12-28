

function interp_value = interp_bilinear(value_grid, gcash, gfund, cash_value, fund_value)
% 双线性插值函数
% 参数:
% value_grid - 待插值的二维值函数网格
% gcash - 现金网格
% gfund - 养老基金网格
% cash_value - 要插值的现金值
% fund_value - 要插值的养老基金值
%
% 返回:
% interp_value - 插值结果

% 限制输入值在网格范围内
cash_value = max(min(cash_value, gcash(end)), gcash(1));
fund_value = max(min(fund_value, gfund(end)), gfund(1));

% 找到在gcash网格中的位置
cash_idx = find(gcash <= cash_value, 1, 'last');
if isempty(cash_idx)
    cash_idx = 1;
elseif cash_idx == length(gcash)
    cash_idx = length(gcash) - 1;
end
cash_weight = (cash_value - gcash(cash_idx)) / (gcash(cash_idx + 1) - gcash(cash_idx));

% 找到在gfund网格中的位置
fund_idx = find(gfund <= fund_value, 1, 'last');
if isempty(fund_idx)
    fund_idx = 1;
elseif fund_idx == length(gfund)
    fund_idx = length(gfund) - 1;
end
fund_weight = (fund_value - gfund(fund_idx)) / (gfund(fund_idx + 1) - gfund(fund_idx));

% 双线性插值
v00 = value_grid(cash_idx, fund_idx);
v01 = value_grid(cash_idx, fund_idx + 1);
v10 = value_grid(cash_idx + 1, fund_idx);
v11 = value_grid(cash_idx + 1, fund_idx + 1);

interp_value = (1 - cash_weight) * (1 - fund_weight) * v00 + ...
              (1 - cash_weight) * fund_weight * v01 + ...
              cash_weight * (1 - fund_weight) * v10 + ...
              cash_weight * fund_weight * v11;
end 