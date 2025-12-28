% 辅助函数：检查值函数范围并在必要时进行缩放
function [scaled_V, scaling_factor, scaled] = scale_value_function(V_matrix, min_range)
% 检查值函数的范围，如果范围小于阈值则进行缩放
%
% 参数:
% V_matrix - 要检查的值函数矩阵
% min_range - 最小可接受的值函数范围（可选，默认为0.01）
%
% 返回:
% scaled_V - 缩放后的值函数
% scaling_factor - 使用的缩放因子
% scaled - 是否进行了缩放

if nargin < 2
    min_range = 0.01;
end

v_min = min(V_matrix(:));
v_max = max(V_matrix(:));
v_range = v_max - v_min;

% 如果值函数范围小于阈值，则进行缩放
if v_range < min_range
    % 计算需要的缩放因子，使得缩放后的范围至少为min_range
    scaling_factor = min_range / v_range;
    % 缩放值函数，保持相对关系不变
    scaled_V = v_min + (V_matrix - v_min) * scaling_factor;
    scaled = true;
else
    % 如果范围足够大，则不需要缩放
    scaled_V = V_matrix;
    scaling_factor = 1.0;
    scaled = false;
end
end