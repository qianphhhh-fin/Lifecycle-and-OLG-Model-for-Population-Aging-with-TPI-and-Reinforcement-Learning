function idx = find_nearest_index(grid, value)
% FIND_NEAREST_INDEX 在网格上查找最接近指定值的索引
% 输入:
%   grid - 单调递增的网格点
%   value - 要查找的值
% 输出:
%   idx - 最接近值的索引

% 处理边界情况
if value <= grid(1)
    idx = 1;
    return;
elseif value >= grid(end)
    idx = length(grid);
    return;
end

% 二分查找
left = 1;
right = length(grid);

while left <= right
    mid = floor((left + right) / 2);
    
    if grid(mid) == value
        idx = mid;
        return;
    elseif grid(mid) < value
        left = mid + 1;
    else
        right = mid - 1;
    end
end

% 确定最接近的索引
if left > length(grid)
    idx = length(grid);
elseif left == 1
    idx = 1;
else
    % 计算与左右两个点的距离，取较近的点
    if abs(grid(left) - value) < abs(grid(left-1) - value)
        idx = left;
    else
        idx = left - 1;
    end
end

end 