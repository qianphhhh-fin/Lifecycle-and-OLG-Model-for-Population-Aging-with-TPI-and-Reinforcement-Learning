function [Z_grid, P] = tauchen(mu, rho, sigma, N, m)
% TAUCHEN 使用Tauchen方法离散化AR(1)过程
% 输入:
%   mu    - AR(1)过程的常数项
%   rho   - AR(1)过程的自回归系数
%   sigma - 随机扰动项的标准差
%   N     - 离散的状态数
%   m     - 状态空间范围的倍数（通常取3）
%
% 输出:
%   Z_grid - 状态变量的网格点
%   P      - 状态转移概率矩阵

% 计算无条件方差
if abs(rho) < 1
    var_z = sigma^2 / (1 - rho^2);
else
    var_z = sigma^2; % 对于单位根过程
end

% 计算状态空间的范围
std_z = sqrt(var_z);
z_max = mu / (1 - rho) + m * std_z;
z_min = mu / (1 - rho) - m * std_z;

% 生成状态空间网格
Z_grid = linspace(z_min, z_max, N);
step = (z_max - z_min) / (N - 1);

% 计算转移概率矩阵
P = zeros(N, N);
for i = 1:N
    for j = 1:N
        % 计算上下边界
        if j == 1
            % 最左边的区间
            upper_bound = Z_grid(j) + step/2;
            P(i, j) = normcdf((upper_bound - mu - rho * Z_grid(i)) / sigma);
        elseif j == N
            % 最右边的区间
            lower_bound = Z_grid(j) - step/2;
            P(i, j) = 1 - normcdf((lower_bound - mu - rho * Z_grid(i)) / sigma);
        else
            % 中间的区间
            upper_bound = Z_grid(j) + step/2;
            lower_bound = Z_grid(j) - step/2;
            P(i, j) = normcdf((upper_bound - mu - rho * Z_grid(i)) / sigma) - ...
                      normcdf((lower_bound - mu - rho * Z_grid(i)) / sigma);
        end
    end
    
    % 归一化，确保每行概率和为1
    P(i, :) = P(i, :) / sum(P(i, :));
end
end 