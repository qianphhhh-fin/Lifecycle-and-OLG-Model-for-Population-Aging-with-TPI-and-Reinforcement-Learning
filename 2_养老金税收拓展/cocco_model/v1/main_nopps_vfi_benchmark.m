%% 生命周期模型最优基准求解 (VFI in Absolute Units)
% 描述:
% 本脚本使用值函数迭代(VFI)求解Cocco(2005)生命周期模型的最优策略。
% **核心特点**: 本脚本的经济环境被设计为与 `CoccoEnv.m` RL环境完全一致，
% 以便为RL智能体的性能评估提供一个严谨、可比的“最优”基准。
%
% - 使用绝对单位的状态变量 (Cash-on-hand)，不进行收入归一化。
% - 参数与 CoccoEnv.m 完全匹配。
% - 模型为确定性版本 (无收入或回报冲击)，与RL的简化环境对应。
%
% 流程:
% 1. setup_parameters_vfi: 初始化与 CoccoEnv 匹配的参数和网格。
% 2. solve_model_vfi_abs: 使用VFI在绝对单位状态空间上逆向求解模型。
% 3. 保存结果供 AgentEvaluator.m 使用。

% --- 初始化 ---
clear;
close all;
clc;

% --- 主流程 ---
fprintf('开始执行 VFI 基准求解...\n');

% 1. 设置模型参数
cS = setup_parameters_vfi();
if ~exist(cS.save_path, 'dir'), mkdir(cS.save_path); end
file_path = fullfile(cS.save_path, cS.save_filename);

% 2. 求解模型 (如果结果文件不存在)
if ~exist(file_path, 'file')
    fprintf('\n===== 步骤 2: 求解模型 (VFI) =====\n');
    vfi_results = solve_model_vfi_abs(cS);
    save(file_path, 'vfi_results', 'cS');
    fprintf('VFI 基准求解结果已保存到 %s\n', file_path);
else
    fprintf('\n已找到现有的 VFI 结果文件，跳过求解: %s\n', file_path);
end

fprintf('\nVFI 基准求解脚本执行完毕。\n');
fprintf('现在你可以在 AgentEvaluator.m 中使用 "%s" 作为 VFI 路径进行对比评估。\n', file_path);


%% =====================================================================
%                       1. 参数设置函数
%  =====================================================================
% main_nopps_vfi_benchmark.m

%% =====================================================================
%                       1. 参数设置函数
%  =====================================================================
function cS = setup_parameters_vfi()
    fprintf('===== 步骤 1: 设置与 CoccoEnv 匹配的参数 =====\n');
    cS = struct();

    % --- A. 从 CoccoEnv.m 复制的核心参数 ---
    cS.tb = 22; cS.tr = 26; cS.td = 30;
    cS.tn = cS.td - cS.tb + 1;
    cS.beta = 0.95; cS.gamma = 3.84;
    cS.rf = 1.02; cS.mu = 0.04;
    cS.ret_fac = 0.6827; cS.tau_y = 0;

    cS.sigr = 0; cS.smay = 0; cS.smav = 0;
    
    % --- B. 收入过程 ---
    cS.f_y = zeros(cS.tn, 1);
    aa = -2.170042 + 2.700381; b1 = 0.16818;
    b2 = -0.0323371 / 10; b3 = 0.0019704 / 100;
    for i = 1:cS.tn
        age_loop = cS.tb + i - 1;
        if age_loop <= cS.tr
            cS.f_y(i) = exp(aa + b1*age_loop + b2*age_loop^2 + b3*age_loop^3 - 3);
        end
    end
    
    % --- C. 数值求解参数 ---
    cS.nX = 100;         
    cS.maxX = 50;        
    % **核心修改**: 将网格的下限从 0 调整为一个很小的正数，以避免 lb > ub 的问题
    cS.minX = 0.25;      
    
    cS.save_path = 'result/';
    cS.save_filename = 'vfi_benchmark_results.mat';

    % --- D. 派生参数与预计算 ---
    % 创建非线性状态网格，现在从 minX (一个小的正数) 开始
    % cS.x_grid = linspace(cS.minX^0.5, cS.maxX^0.5, cS.nX)'.^2;
    
    lgX = linspace(log(cS.minX), log(cS.maxX), cS.nX)';
    cS.x_grid = exp(lgX);

    cS.survprob = ones(cS.tn - 1, 1);
    
    cS.n_shocks = 1;
    cS.shock_grid = 0; cS.shock_weig = 1;
    cS.gret = cS.rf + cS.mu + cS.shock_grid * cS.sigr;
    cS.yh = exp(cS.shock_grid * cS.smay);
    
    fprintf('参数设置完成。\n');
end
%% =====================================================================
%               2. VFI 求解器函数 (绝对单位)
%  =====================================================================
% main_nopps_vfi_benchmark.m

%% =====================================================================
%               2. VFI 求解器函数 (绝对单位, fmincon)
%  =====================================================================
function results = solve_model_vfi_abs(cS)
    % 初始化
    V = zeros(cS.nX, cS.tn);
    C_policy = zeros(cS.nX, cS.tn);
    A_policy = zeros(cS.nX, cS.tn);

    % --- 终端期 (t=T) ---
    C_policy(:, cS.tn) = cS.x_grid;
    A_policy(:, cS.tn) = 0.0; % 最后一期无储蓄
    if cS.gamma == 1
        V(:, cS.tn) = log(C_policy(:, cS.tn));
    else
        V(:, cS.tn) = (C_policy(:, cS.tn).^(1 - cS.gamma)) / (1 - cS.gamma);
    end
    V(C_policy(:, cS.tn) <= 0, cS.tn) = -inf; % 防止log(0)

    % --- 逆向归纳求解 ---
    % **核心修改**: 使用与 main_noPPS.m 一致的 fmincon 设置
    options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point', ...
        'MaxFunctionEvaluations', 1000, 'MaxIterations', 400, ...
        'OptimalityTolerance', 1e-10, 'StepTolerance', 1e-10, ...
        'ConstraintTolerance', 1e-10);
    
    % -- 退休期 (t = T-1, ..., K) --
    fprintf('求解退休期 (绝对单位, 使用 fmincon)...\n');
    for t = (cS.tn-1):-1:(cS.tr - cS.tb + 1)
        fprintf('  t = %d (年龄 %d)\n', t, t + cS.tb - 1);
        V_next_interp = griddedInterpolant(cS.x_grid, V(:, t+1), 'spline', 'linear');
        
        for i_x = 1:cS.nX
            X_t = cS.x_grid(i_x);
            
            % 在确定性且mu>rf-1时，最优alpha恒为1.0
            alpha_opt = 1.0; 
            
            % 我们只优化消费 C
            % 定义目标函数，只接受一个输入 C
            obj_fun = @(C) -objective_ret(C, alpha_opt, X_t, V_next_interp, cS, t);
            
            % 设置优化变量 C 的边界和初始值
            lb = 1e-6; % C > 0
            ub = X_t;  % C <= X_t
            if i_x > 1
                x0 = C_policy(i_x-1, t); % 使用上一个网格点的解作为热启动
            else
                x0 = X_t * 0.5;
            end
            x0 = max(lb, min(ub, x0)); % 确保初始值在边界内

            % 使用 fmincon 进行一维优化
            [C_opt, V_opt_neg] = fmincon(obj_fun, x0, [], [], [], [], lb, ub, [], options);
            
            C_policy(i_x, t) = C_opt;
            A_policy(i_x, t) = alpha_opt;
            V(i_x, t) = -V_opt_neg;
        end
    end

    % -- 工作期 (t = K-1, ..., 1) --
    fprintf('求解工作期 (绝对单位, 使用 fmincon)...\n');
    P_retirement_final = 0; 
    for t = (cS.tr - cS.tb):-1:1
        fprintf('  t = %d (年龄 %d)\n', t, t + cS.tb - 1);
        V_next_interp = griddedInterpolant(cS.x_grid, V(:, t+1), 'spline', 'linear');
        
        for i_x = 1:cS.nX
            X_t = cS.x_grid(i_x);
            
            alpha_opt = 1.0;
            
            obj_fun = @(C) -objective_work(C, alpha_opt, X_t, V_next_interp, cS, t);
            
            lb = 1e-6;
            ub = X_t;
            if i_x > 1
                x0 = C_policy(i_x-1, t);
            else
                x0 = X_t * 0.5;
            end
            x0 = max(lb, min(ub, x0));

            [C_opt, V_opt_neg] = fmincon(obj_fun, x0, [], [], [], [], lb, ub, [], options);
            
            C_policy(i_x, t) = C_opt;
            A_policy(i_x, t) = alpha_opt;
            V(i_x, t) = -V_opt_neg;
        end
        if (t + cS.tb - 1) == (cS.tr - 1), P_retirement_final = 0; end
    end

    results.V = V;
    results.C_policy = C_policy;
    results.A_policy = A_policy;
    results.P_retirement_final = P_retirement_final;
end
function value = objective_ret(C, alpha, X_t, V_next_interp, cS, t)
    if C <= 0 || C > X_t, value = -1e20; return; end
    
    % 即期效用
    if cS.gamma == 1, u = log(C); else, u = (C.^(1 - cS.gamma)) / (1 - cS.gamma); end

    % 下一期状态
    W_next = (X_t - C) * (alpha * cS.gret + (1 - alpha) * cS.rf);
    Y_next = cS.ret_fac * cS.f_y(cS.tr - cS.tb + 1); % 退休收入是固定的
    X_next = W_next + Y_next;

    % 期望价值 (这里是确定性的)
    EV_next = V_next_interp(X_next);
    
    value = u + cS.beta * cS.survprob(t) * EV_next;
end

function value = objective_work(C, alpha, X_t, V_next_interp, cS, t)
    if C <= 0 || C > X_t, value = -1e20; return; end
    
    % 即期效用
    if cS.gamma == 1, u = log(C); else, u = (C.^(1 - cS.gamma)) / (1 - cS.gamma); end

    % 下一期状态
    W_next = (X_t - C) * (alpha * cS.gret + (1 - alpha) * cS.rf);
    Y_next = cS.f_y(t + 1) * cS.yh; % P=1, 只有确定性增长和暂时冲击
    X_next = W_next + Y_next;

    % 期望价值 (这里是确定性的)
    EV_next = V_next_interp(X_next);
    
    value = u + cS.beta * cS.survprob(t) * EV_next;
end