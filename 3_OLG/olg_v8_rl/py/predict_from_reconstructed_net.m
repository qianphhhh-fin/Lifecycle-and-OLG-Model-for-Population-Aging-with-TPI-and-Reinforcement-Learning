function action = predict_from_reconstructed_net(observation)
% 该函数用于从Python调用，以测试重建的RL网络。
% 它加载网络，进行一次预测，并返回结果。
% 使用 persistent 变量来缓存网络，避免重复加载。

persistent net params;

% --- 1. 仅在第一次调用时加载和初始化网络 ---
if isempty(net)
    disp('MATLAB: 第一次调用，正在加载和重建网络...');
    
    % 加载模型参数 (用于动作缩放)
    paramsFile = './py/best_model_sb3/best_model_params.json';
    if ~isfile(paramsFile)
        error('MATLAB: RL参数文件未找到: %s', paramsFile);
    end
    params = jsondecode(fileread(paramsFile));
    
    % 从ONNX文件加载网络
    onnxFile = './py/best_model_sb3/best_model.onnx';
    if ~isfile(onnxFile)
        error('MATLAB: ONNX模型文件未找到: %s', onnxFile);
    end
    
    try
        disp(['MATLAB: 正在从 ' onnxFile ' 导入ONNX网络...']);
        % 更新: 使用推荐的 importNetworkFromONNX 函数, 并移除无效的 'TargetNetwork' 参数
        net = importNetworkFromONNX(onnxFile, "InputDataFormats", "BC", "OutputDataFormats", "BC");
        disp('MATLAB: ONNX网络导入成功!');
    catch ME
        disp('MATLAB: ONNX网络导入失败。');
        disp('这可能是由于MATLAB的ONNX解析器版本不兼容。');
        disp('错误信息:');
        disp(ME.message);
        rethrow(ME);
    end
    
    disp('MATLAB: 网络加载并缓存成功!');
else
    disp('MATLAB: 使用已缓存的网络。');
end

% --- 2. 进行预测 ---
% 将输入转换为 single 类型
observation_vector = single(observation(:)); 

% 调用封装的预测函数，该函数处理dlarray转换和动作缩放
action = predict_for_olg(net, observation_vector, params.action_space_low, params.action_space_high);

% --- 3. 将结果转换为double类型的行向量，以便与Python兼容 ---
action = double(action(:)');

end


%% HELPER FUNCTION
% =========================================================================

function action = predict_for_olg(net, observation, low, high)
    % 封装的预测函数，用于OLG模拟
    
    % 为ONNX导入的网络，输入需要是 (Batch, Features) 格式
    observation_row = single(observation(:)'); % 确保是行向量 (1, num_features)
    dlX = dlarray(observation_row, 'BC');      % 'BC' for (Batch, Channel/Features)
    
    % 进行预测
    dlY = predict(net, dlX);
    raw_action = extractdata(dlY);
    
    % 将动作从 [-1, 1] 缩放到 [low, high]
    % 注意：要确保 low 和 high 是列向量以便正确广播
    action = low(:) + (high(:) - low(:)) .* (raw_action(:) + 1) / 2;
end 