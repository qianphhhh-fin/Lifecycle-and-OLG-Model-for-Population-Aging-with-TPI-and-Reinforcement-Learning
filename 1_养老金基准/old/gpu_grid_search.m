function [x_opt, f_opt] = gpu_grid_search(fun, lb, ub, grid_size, varargin)
% GPU网格搜索优化函数，用于替代fmincon
%
% 参数:
% fun - 目标函数句柄
% lb - 下界 [lb1, lb2, ...] 
% ub - 上界 [ub1, ub2, ...]
% grid_size - 每个维度的网格点数量 [size1, size2, ...]
% varargin - 传递给目标函数的额外参数
%
% 返回:
% x_opt - 最优参数
% f_opt - 最优函数值

% 获取变量维度
dim = length(lb);

% 如果grid_size是标量，则所有维度使用相同的网格大小
if length(grid_size) == 1
    grid_size = repmat(grid_size, 1, dim);
end

% 创建网格点
grid_cell = cell(dim, 1);
for i = 1:dim
    grid_cell{i} = gpuArray(linspace(lb(i), ub(i), grid_size(i)));
end

% 创建完整网格
[grid_arrays{1:dim}] = ndgrid(grid_cell{:});

% 预分配结果数组
values = gpuArray(inf(grid_size));

% 使用循环代替arrayfun，避免嵌套函数调用
switch dim
    case 1
        % 1D情况
        for i = 1:grid_size(1)
            x = grid_arrays{1}(i);
            values(i) = fun([x], varargin{:});
        end
    case 2
        % 2D情况
        for i = 1:grid_size(1)
            for j = 1:grid_size(2)
                x = grid_arrays{1}(i, j);
                y = grid_arrays{2}(i, j);
                values(i, j) = fun([x, y], varargin{:});
            end
        end
    case 3
        % 3D情况
        for i = 1:grid_size(1)
            for j = 1:grid_size(2)
                for k = 1:grid_size(3)
                    x = grid_arrays{1}(i, j, k);
                    y = grid_arrays{2}(i, j, k);
                    z = grid_arrays{3}(i, j, k);
                    values(i, j, k) = fun([x, y, z], varargin{:});
                end
            end
        end
    otherwise
        error('目前仅支持1-3维网格搜索');
end

% 找到最优点
[f_opt, idx] = min(values(:));
f_opt = gather(f_opt);

% 将线性索引转换为多维索引
idx_sub = cell(1, dim);
[idx_sub{:}] = ind2sub(size(values), idx);

% 构建最优参数
x_opt = zeros(1, dim);
for i = 1:dim
    x_opt(i) = gather(grid_cell{i}(idx_sub{i}));
end

end 