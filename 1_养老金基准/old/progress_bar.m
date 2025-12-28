function progress_bar(current, total, bar_width, message)
% PROGRESS_BAR 显示命令行进度条
% 输入:
%   current - 当前进度
%   total - 总任务数
%   bar_width - 进度条宽度（字符数）
%   message - 显示的消息前缀

if nargin < 3
    bar_width = 30;
end

if nargin < 4
    message = '';
end

% 计算完成百分比
percentage = current/total;
filled_width = round(bar_width * percentage);
empty_width = bar_width - filled_width;

% 构建进度条字符串
bar = ['[', repmat('=', 1, filled_width), repmat(' ', 1, empty_width), ']'];
percentage_str = sprintf('%3d%%', round(percentage * 100));

% 构建完整的消息
full_message = sprintf('%s %s %s (%d/%d)', message, bar, percentage_str, current, total);

% 回车不换行，更新同一行
fprintf('%s\r', full_message);

% 如果是最后一个，添加换行符
if current == total
    fprintf('\n');
end

end 