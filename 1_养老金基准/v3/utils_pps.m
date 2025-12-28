classdef utils_pps
    methods (Static)


        %% =====================================================================
        %                       1. 参数设置函数
        %  =====================================================================
        function cS = setup_parameters(varargin)
            fprintf('===== 步骤 1: 设置模型参数 =====\n');

            cS = struct();

            % --- A. 生命周期与人口结构 ---
            cS.tb = 22;       % 初始工作年龄
            cS.tr = 61;       % 退休年龄
            cS.td = 100;      % 最高寿命
            cS.tn = cS.td - cS.tb + 1; % 总期数

            % --- B. 经济主体偏好 (Epstein-Zin) ---
            cS.beta = 0.95;     % 折现因子
            cS.gamma = 3.84;    % 相对风险规避系数
            cS.psi = 0.15;      % 跨期替代弹性

            % --- C. 收入过程 ---
            % 确定性部分 g(t) = exp(aa + b1*t + b2*t^2 + b3*t^3)
            cS.aa = -2.170042 + 2.700381;
            cS.b1 = 0.16818;
            cS.b2 = -0.0323371 / 10;
            cS.b3 = 0.0019704 / 100;
            % 随机部分
            cS.smay = sqrt(0.169993); % 暂时冲击标准差 (u_t)
            cS.smav = sqrt(0.112572); % 持久冲击标准差 (z_t)
            cS.corr_z_epsilon = 0; % 持久冲击与风险收益的相关性 (未使用)
            cS.corr_u_epsilon = 0; % 暂时冲击与风险收益的相关性 (未使用)

            % --- D. 资产与回报率 ---
            cS.rf = 1.02;     % 无风险总回报率
            cS.mu = 0.04;     % 风险资产超额回报率
            cS.sigr = 0.27;   % 风险资产回报率标准差
            cS.rp = 1.04;     % 个人养老金账户(PPS)回报率
            cS.sigrrp = 0.27; % 个人养老金账户(PPS)风险资产回报率标准差

            % --- E. 养老金与税收 ---
            cS.ret_fac = 0.6827;      % 退休后基础养老金 (替代率)
            cS.tau_y = 0.06;          % 工资税率
            cS.tau_q = 0.03;          % 个人养老金提取税率
            % cS.Q_max 将在下方通过新函数计算


            % --- F. 数值求解参数 ---
            cS.n_shocks = 5;      % 离散化的随机冲击节点数
            cS.ncash = 51;        % 现金持有量网格数
            cS.nfund = 51;        % 养老金账户网格数
            cS.maxcash = 200;     % 现金网格上限 (归一化)
            cS.mincash = 0.25;    % 现金网格下限 (归一化)
            cS.maxfund = 200;     % 养老金网格上限 (归一化)
            cS.minfund = 0.01;    % 养老金网格下限 (归一化)
            cS.save_path = 'result/'; % 结果保存路径

                        % 生存概率
            cS.survprob = zeros(cS.tn - 1, 1);
            cS.survprob_data = [0.99975, 0.99974, 0.99973, 0.99972, 0.99971, 0.99969, 0.99968, 0.99966, 0.99963, 0.99961, ...
                0.99958, 0.99954, 0.99951, 0.99947, 0.99943, 0.99938, 0.99934, 0.99928, 0.99922, 0.99916, ...
                0.99908, 0.999, 0.99891, 0.99881, 0.9987, 0.99857, 0.99843, 0.99828, 0.99812, 0.99794, ...
                0.99776, 0.99756, 0.99735, 0.99713, 0.9969, 0.99666, 0.9964, 0.99613, 0.99583, 0.99551, ...
                0.99515, 0.99476, 0.99432, 0.99383, 0.9933, 0.9927, 0.99205, 0.99133, 0.99053, 0.98961, ...
                0.98852, 0.98718, 0.98553, 0.98346, 0.98089, 0.97772, 0.97391, 0.96943, 0.96429, 0.95854, ...
                0.95221, 0.94537, 0.93805, 0.93027, 0.92202, 0.91327, 0.90393, 0.89389, 0.88304, 0.87126, ...
                0.85846, 0.84452, 0.82935, 0.81282, 0.79485, 0.77543, 0.75458, 0.7324, 0.70893, 0.68424, ...
                0.68424, 0.68424];

                        % *** 新增: 调用函数计算归一化的 Q_max ***
            % cS.Q_max = utils_pps.calculate_q_max(absolute_q_max_yuan, anchor_income_yuan, cS.f_y);
            cS.Q_max =99999;

            if isempty(varargin) % 便于敏感性分析时改动参数而不进一步计算派生参数
                cS = utils_pps.setup_process(cS);
            end
           
        end

        function cS = setup_process(cS)
            % --- G. 派生参数与预计算 ---
            cS.pension_rate = 1 / (cS.td - cS.tr + 1); % 退休后年金化给付比率
            % 效用函数参数
            cS.psi_1 = 1.0 - 1.0 / cS.psi;
            cS.psi_2 = 1.0 / cS.psi_1;
            cS.theta = (1.0 - cS.gamma) / (1.0 - cS.psi_1);


            for i = 1:min(length(cS.survprob_data), length(cS.survprob))
                cS.survprob(i) = cS.survprob_data(i);
            end

            % 离散化正态分布 (Tauchen方法)
            tauchenoptions.parallel=0;
            [grid, weig_matrix] = utils_pps.discretizeAR1_Tauchen(0, 0, 1, cS.n_shocks, 2, tauchenoptions);
            cS.shock_grid = grid;
            cS.shock_weig = diag(weig_matrix);

            % 风险资产回报率网格
            cS.gret = cS.rf + cS.mu + cS.shock_grid * cS.sigr;
            cS.gret_rp = cS.rp + cS.shock_grid * cS.sigrrp;

            % 状态转移概率 (三维冲击)
            cS.nweig1 = zeros(cS.n_shocks, cS.n_shocks, cS.n_shocks);
            for i6 = 1:cS.n_shocks
                for i7 = 1:cS.n_shocks
                    for i8 = 1:cS.n_shocks
                        cS.nweig1(i6, i7, i8) = cS.shock_weig(i6) * cS.shock_weig(i7) * cS.shock_weig(i8);
                    end
                end
            end

            % 状态变量网格
            lgcash = linspace(log(cS.mincash), log(cS.maxcash), cS.ncash)';
            cS.gcash = exp(lgcash);
            lgfund = linspace(log(cS.minfund), log(cS.maxfund), cS.nfund)';
            cS.gfund = exp(lgfund);
            cS.gfund(1) = 0; % 确保可以从0开始

            % --- 收入过程计算 ---
            % 确定性收入剖面
            cS.f_y = zeros(cS.tr - cS.tb + 1, 1);
            for i1 = cS.tb:cS.tr
                age = i1;
                cS.f_y(age - cS.tb + 1) = exp(cS.aa + cS.b1*age + cS.b2*age^2 + cS.b3*age^3);
            end


            % 收入过程网格
            cS.yh = zeros(cS.n_shocks, cS.n_shocks); % 暂时冲击 exp(u_t)
            cS.yp = zeros(cS.n_shocks, cS.n_shocks); % 持久冲击增长 exp(z_t)
            cS.gyp = zeros(cS.n_shocks, cS.n_shocks, cS.tn - 1); % 总收入增长

            for i1 = 1:cS.n_shocks
                grid2_u = cS.shock_grid(i1) .* cS.corr_u_epsilon + cS.shock_grid .* sqrt(1 - cS.corr_u_epsilon^2);
                cS.yh(:, i1) = exp(grid2_u * cS.smay);

                grid2_z = cS.shock_grid(i1) .* cS.corr_z_epsilon + cS.shock_grid .* sqrt(1 - cS.corr_z_epsilon^2);
                cS.yp(:, i1) = exp(grid2_z * cS.smav);
            end

            % 工作期收入增长
            work_periods = cS.tr - cS.tb;
            for t = 1:work_periods
                G_t = cS.f_y(t+1) / cS.f_y(t); % 确定性增长
                cS.gyp(:,:,t) = repmat(G_t, cS.n_shocks, cS.n_shocks) .* cS.yp;
            end

            % 退休期收入增长 (无增长)
            cS.gyp(:,:,(work_periods+1):(cS.tn-1)) = 1.0;

            fprintf('参数设置完成。\n');
        end

        %
        %% =====================================================================
        %                       2. VFI 求解器函数 (修正版)
        %  =====================================================================
        function results = solve_model_vfi(cS)
            % --- 初始化 ---
            % 启动并行池
            if isempty(gcp('nocreate'))
                parpool('local');
            end

            % 初始化策略和价值函数数组
            V = zeros(cS.ncash, cS.nfund, cS.tn);
            C_policy = zeros(cS.ncash, cS.nfund, cS.tn);
            A_policy = zeros(cS.ncash, cS.nfund, cS.tn);
            Q_policy = zeros(cS.ncash, cS.nfund, cS.tn);

            % --- 终端期 (t=T) ---
            % Epstein-Zin效用下的终端价值 V_T = ((1-beta)*C_T^(1-1/psi))^(1/(1-1/psi))
            % 其中 C_T = W_T (即 cS.gcash), V_T = (1-beta)^psi_2 * W_T
            for i_fund = 1:cS.nfund
                V(:, i_fund, cS.tn) = cS.gcash * (1-cS.beta)^cS.psi_2;
                C_policy(:, i_fund, cS.tn) = cS.gcash;
            end
            A_policy(:, :, cS.tn) = 0.0;
            Q_policy(:, :, cS.tn) = 0.0;

            % --- 逆向归纳求解 ---
            fprintf('开始逆向归纳求解...\n');
            options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point', ...
                'MaxFunctionEvaluations', 1000, 'MaxIterations', 400, ...
                'OptimalityTolerance', 1e-10, 'StepTolerance', 1e-10, ...
                'ConstraintTolerance', 1e-10);

            % -- 退休期 (t = T-1, ..., K) --
            fprintf('求解退休期...\n');
            retire_start_time = tic;
            for t = (cS.tn-1):-1:(cS.tr-cS.tb+1)
                fprintf('  t = %d (年龄 %d)\n', t, t + cS.tb - 1);
                V_next = squeeze(V(:,:,t+1)); % 下一期值函数

                % *** 修改: 创建 griddedInterpolant 对象 ***
                % V_next 维度 (ncash, nfund)。行对应gcash, 列对应gfund
                % griddedInterpolant({Y,X}, V)


                temp_C = zeros(cS.ncash, cS.nfund);
                temp_A = zeros(cS.ncash, cS.nfund);
                temp_V = zeros(cS.ncash, cS.nfund);

                parfor i_fund = 1:cS.nfund

                    local_C = zeros(cS.ncash, 1);
                    local_A = zeros(cS.ncash, 1);
                    local_V = zeros(cS.ncash, 1);

                    gfund_at_retirement = cS.gfund(i_fund);
                    constant_pension_payout = (1 - cS.tau_q) * gfund_at_retirement * cS.pension_rate;

                    V_next_interp_obj = griddedInterpolant(cS.gcash, V_next(:,i_fund),  'spline', 'linear');

                    for i_cash = 1:cS.ncash
                        gcash_val = cS.gcash(i_cash);

                        lb = [1e-6, 0];
                        x0 = [0.5 * gcash_val, 1-1e-06];
                        if i_cash > 1 && local_C(i_cash-1) > 0
                            x0 = [local_C(i_cash-1), local_A(i_cash-1)];
                            % lb = [local_C(i_cash-1)-1e-06, 0];
                        end

                        ub = [gcash_val, 1];

                        % *** 修改: 将插值对象传递给目标函数 ***
                        obj_fun = @(x) utils_pps.fun_valuefunc_retired(x, gcash_val, constant_pension_payout, i_fund, V_next_interp_obj, cS, t);

                        [policy, fval] = fmincon(obj_fun, x0, [], [], [], [], lb, ub, [], options);

                        local_C(i_cash) = policy(1);
                        % if abs(policy(1)-gcash_val)<1e-05
                        %     local_A(i_cash) = nan; % C = Cash
                        % else
                        local_A(i_cash) = policy(2);
                        % end
                        local_V(i_cash) = -fval;
                    end
                    temp_C(:, i_fund) = local_C;
                    temp_A(:, i_fund) = local_A;
                    temp_V(:, i_fund) = local_V;
                end
                C_policy(:, :, t) = temp_C;
                A_policy(:, :, t) = temp_A;
                V(:, :, t) = temp_V;
            end
            fprintf('退休期求解完成，耗时 %.2f 分钟\n', toc(retire_start_time)/60);

            % -- 工作期 (t = K-1, ..., 1) --
            fprintf('求解工作期...\n');
            work_start_time = tic;
            for t = (cS.tr - cS.tb):-1:1
                fprintf('  t = %d (年龄 %d)\n', t, t + cS.tb - 1);
                V_next = squeeze(V(:,:,t+1));

                % *** 修改: 创建 griddedInterpolant 对象 ***
                [X,Y] = ndgrid(cS.gcash,cS.gfund);
                V_next_interp_obj = griddedInterpolant(X, Y, V_next,  'spline','linear');

                temp_C = zeros(cS.ncash, cS.nfund);
                temp_A = zeros(cS.ncash, cS.nfund);
                temp_Q = zeros(cS.ncash, cS.nfund);
                temp_V = zeros(cS.ncash, cS.nfund);

                parfor i_fund = 1:cS.nfund
                    local_C = zeros(cS.ncash, 1);
                    local_A = zeros(cS.ncash, 1);
                    local_Q = zeros(cS.ncash, 1);
                    local_V = zeros(cS.ncash, 1);

                    gfund_val = cS.gfund(i_fund);
                    for i_cash = 1:cS.ncash
                        gcash_val = cS.gcash(i_cash);

                        % 2. 放弃在 lb 中施加单调性硬约束，Q 的下界永远是 0
                        lb = [1e-3, 0, 0];
                        % *** 核心修改 ***
                        if i_cash > 1
                            x0 = [local_C(i_cash-1), local_A(i_cash-1), local_Q(i_cash-1)];
                            % lb = [local_C(i_cash-1)-1e-6, 0, local_Q(i_cash-1)-1e-6];
                        else
                            x0 = [gcash_val * 0.5, 1-1e-06, 1e-6]; % C, alpha, Q
                        end

                        % 3. ub 只包含物理和政策约束
                        ub = [gcash_val, 1, min(cS.Q_max,gcash_val/(1 - cS.tau_y))];


                        % *** 修改: 将插值对象传递给目标函数 ***
                        obj_fun = @(x) utils_pps.fun_valuefunc_work(x, gcash_val, gfund_val, V_next_interp_obj, cS, t);

                        A_ineq = [1, 0, 1 - cS.tau_y];
                        b_ineq = gcash_val;

                        [policy, fval] = fmincon(obj_fun, x0, A_ineq, b_ineq, [], [], lb, ub, [], options);

                        local_C(i_cash) = policy(1);
                        % if abs(A_ineq*policy'-b_ineq)<1e-05
                        %     local_A(i_cash) = nan; % C+（1-tau_y）Q = Cash
                        % else
                        local_A(i_cash) = policy(2);
                        % end
                        local_Q(i_cash) = policy(3);
                        % if abs(policy(3))<0.0001
                        %     policy(3)
                        % end
                        local_V(i_cash) = -fval;
                    end
                    temp_C(:, i_fund) = local_C;
                    temp_A(:, i_fund) = local_A;
                    temp_Q(:, i_fund) = local_Q;
                    temp_V(:, i_fund) = local_V;
                end
                C_policy(:, :, t) = temp_C;
                A_policy(:, :, t) = temp_A;
                Q_policy(:, :, t) = temp_Q;
                V(:, :, t) = temp_V;
            end
            fprintf('工作期求解完成，耗时 %.2f 分钟\n', toc(work_start_time)/60);

            % --- 封装结果 ---
            results.V = V;
            results.C_policy = C_policy;
            results.A_policy = A_policy;
            results.Q_policy = Q_policy;
            % 平滑处理
            % results.A_policy(find(abs(cS.gcash-results.C_policy-(1-cS.tau_y)*results.Q_policy)<1e-05))=nan;
        end
        
        
        function results = solve_model_vfi_riskyrp(cS)
            % --- 初始化 ---
            % 启动并行池
            if isempty(gcp('nocreate'))
                parpool('local');
            end

            % 初始化策略和价值函数数组
            V = zeros(cS.ncash, cS.nfund, cS.tn);
            C_policy = zeros(cS.ncash, cS.nfund, cS.tn);
            A_policy = zeros(cS.ncash, cS.nfund, cS.tn);
            Q_policy = zeros(cS.ncash, cS.nfund, cS.tn);

            % --- 终端期 (t=T) ---
            % Epstein-Zin效用下的终端价值 V_T = ((1-beta)*C_T^(1-1/psi))^(1/(1-1/psi))
            % 其中 C_T = W_T (即 cS.gcash), V_T = (1-beta)^psi_2 * W_T
            for i_fund = 1:cS.nfund
                V(:, i_fund, cS.tn) = cS.gcash * (1-cS.beta)^cS.psi_2;
                C_policy(:, i_fund, cS.tn) = cS.gcash;
            end
            A_policy(:, :, cS.tn) = 0.0;
            Q_policy(:, :, cS.tn) = 0.0;

            % --- 逆向归纳求解 ---
            fprintf('开始逆向归纳求解...\n');
            options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point', ...
                'MaxFunctionEvaluations', 1000, 'MaxIterations', 400, ...
                'OptimalityTolerance', 1e-10, 'StepTolerance', 1e-10, ...
                'ConstraintTolerance', 1e-10);

            % -- 退休期 (t = T-1, ..., K) --
            fprintf('求解退休期...\n');
            retire_start_time = tic;
            for t = (cS.tn-1):-1:(cS.tr-cS.tb+1)
                fprintf('  t = %d (年龄 %d)\n', t, t + cS.tb - 1);
                V_next = squeeze(V(:,:,t+1)); % 下一期值函数

                % *** 修改: 创建 griddedInterpolant 对象 ***
                % V_next 维度 (ncash, nfund)。行对应gcash, 列对应gfund
                % griddedInterpolant({Y,X}, V)


                temp_C = zeros(cS.ncash, cS.nfund);
                temp_A = zeros(cS.ncash, cS.nfund);
                temp_V = zeros(cS.ncash, cS.nfund);

                parfor i_fund = 1:cS.nfund

                    local_C = zeros(cS.ncash, 1);
                    local_A = zeros(cS.ncash, 1);
                    local_V = zeros(cS.ncash, 1);

                    gfund_at_retirement = cS.gfund(i_fund);
                    constant_pension_payout = (1 - cS.tau_q) * gfund_at_retirement * cS.pension_rate;

                    V_next_interp_obj = griddedInterpolant(cS.gcash, V_next(:,i_fund),  'spline', 'linear');

                    for i_cash = 1:cS.ncash
                        gcash_val = cS.gcash(i_cash);

                        lb = [1e-6, 0];
                        x0 = [0.5 * gcash_val, 1-1e-06];
                        if i_cash > 1 && local_C(i_cash-1) > 0
                            x0 = [local_C(i_cash-1), local_A(i_cash-1)];
                            % lb = [local_C(i_cash-1)-1e-06, 0];
                        end

                        ub = [gcash_val, 1];

                        % *** 修改: 将插值对象传递给目标函数 ***
                        obj_fun = @(x) utils_pps.fun_valuefunc_retired(x, gcash_val, constant_pension_payout, i_fund, V_next_interp_obj, cS, t);

                        [policy, fval] = fmincon(obj_fun, x0, [], [], [], [], lb, ub, [], options);

                        local_C(i_cash) = policy(1);
                        % if abs(policy(1)-gcash_val)<1e-05
                        %     local_A(i_cash) = nan; % C = Cash
                        % else
                        local_A(i_cash) = policy(2);
                        % end
                        local_V(i_cash) = -fval;
                    end
                    temp_C(:, i_fund) = local_C;
                    temp_A(:, i_fund) = local_A;
                    temp_V(:, i_fund) = local_V;
                end
                C_policy(:, :, t) = temp_C;
                A_policy(:, :, t) = temp_A;
                V(:, :, t) = temp_V;
            end
            fprintf('退休期求解完成，耗时 %.2f 分钟\n', toc(retire_start_time)/60);

            % -- 工作期 (t = K-1, ..., 1) --
            fprintf('求解工作期...\n');
            work_start_time = tic;
            for t = (cS.tr - cS.tb):-1:1
                fprintf('  t = %d (年龄 %d)\n', t, t + cS.tb - 1);
                V_next = squeeze(V(:,:,t+1));

                % *** 修改: 创建 griddedInterpolant 对象 ***
                [X,Y] = ndgrid(cS.gcash,cS.gfund);
                V_next_interp_obj = griddedInterpolant(X, Y, V_next,  'spline','linear');

                temp_C = zeros(cS.ncash, cS.nfund);
                temp_A = zeros(cS.ncash, cS.nfund);
                temp_Q = zeros(cS.ncash, cS.nfund);
                temp_V = zeros(cS.ncash, cS.nfund);

                parfor i_fund = 1:cS.nfund
                    local_C = zeros(cS.ncash, 1);
                    local_A = zeros(cS.ncash, 1);
                    local_Q = zeros(cS.ncash, 1);
                    local_V = zeros(cS.ncash, 1);

                    gfund_val = cS.gfund(i_fund);
                    for i_cash = 1:cS.ncash
                        gcash_val = cS.gcash(i_cash);

                        % 2. 放弃在 lb 中施加单调性硬约束，Q 的下界永远是 0
                        lb = [1e-3, 0, 0];
                        % *** 核心修改 ***
                        if i_cash > 1
                            x0 = [local_C(i_cash-1), local_A(i_cash-1), local_Q(i_cash-1)];
                            % lb = [local_C(i_cash-1)-1e-6, 0, local_Q(i_cash-1)-1e-6];
                        else
                            x0 = [gcash_val * 0.5, 1-1e-06, 1e-6]; % C, alpha, Q
                        end

                        % 3. ub 只包含物理和政策约束
                        ub = [gcash_val, 1, min(cS.Q_max,gcash_val/(1 - cS.tau_y))];


                        % *** 修改: 将插值对象传递给目标函数 ***
                        obj_fun = @(x) utils_pps.fun_valuefunc_work_riskrp(x, gcash_val, gfund_val, V_next_interp_obj, cS, t);

                        A_ineq = [1, 0, 1 - cS.tau_y];
                        b_ineq = gcash_val;

                        [policy, fval] = fmincon(obj_fun, x0, A_ineq, b_ineq, [], [], lb, ub, [], options);

                        local_C(i_cash) = policy(1);
                        % if abs(A_ineq*policy'-b_ineq)<1e-05
                        %     local_A(i_cash) = nan; % C+（1-tau_y）Q = Cash
                        % else
                        local_A(i_cash) = policy(2);
                        % end
                        local_Q(i_cash) = policy(3);
                        % if abs(policy(3))<0.0001
                        %     policy(3)
                        % end
                        local_V(i_cash) = -fval;
                    end
                    temp_C(:, i_fund) = local_C;
                    temp_A(:, i_fund) = local_A;
                    temp_Q(:, i_fund) = local_Q;
                    temp_V(:, i_fund) = local_V;
                end
                C_policy(:, :, t) = temp_C;
                A_policy(:, :, t) = temp_A;
                Q_policy(:, :, t) = temp_Q;
                V(:, :, t) = temp_V;
            end
            fprintf('工作期求解完成，耗时 %.2f 分钟\n', toc(work_start_time)/60);

            % --- 封装结果 ---
            results.V = V;
            results.C_policy = C_policy;
            results.A_policy = A_policy;
            results.Q_policy = Q_policy;
            % 平滑处理
            % results.A_policy(find(abs(cS.gcash-results.C_policy-(1-cS.tau_y)*results.Q_policy)<1e-05))=nan;
        end

        %% =====================================================================
        %                VFI - 退休期值函数
        %  =====================================================================
        function v = fun_valuefunc_retired(x, cash, pension_payout, i_fund_K, V_next_interp_obj, cS, t)
            % x = [消费比例 c_rate, 风险资产比例 alpha]
            % pension_payout: 这是一个根据F_K计算出的固定值
            % i_fund_K: 这是F_K在网格上的索引，用于插值V_next，因为它在退休期间不变

            alpha = x(2);

            consumption = x(1);
            if consumption <= 1e-8, v = 1e10; return; end

            sav = cash - consumption;
            if sav < 0, v = 1e10; return; end

            % Epstein-Zin 瞬时效用项
            u_term = (1 - cS.beta) * (consumption^cS.psi_1);

            % 下一期现金持有量 (归一化)
            ret_fac_norm = cS.ret_fac;
            cash_1 = (cS.rf * (1 - alpha) + cS.gret * alpha) * sav + ret_fac_norm + pension_payout;
            % cash_1 = max(min(cash_1, cS.gcash(end)), cS.gcash(1));

            % 下一期养老金状态 F_{t+1} 就是 F_K
            fund_1 = repmat(cS.gfund(i_fund_K), cS.n_shocks, 1);

            % *** 修改: 使用 griddedInterpolant 对象进行插值 ***
            % 调用方式 F({Yq, Xq})
            int_V = V_next_interp_obj(cash_1);
            int_V(int_V <= 0) = 1e-9;

            % 计算期望
            expected_term = cS.shock_weig' * (int_V.^(1 - cS.gamma));

            % 完整的 Epstein-Zin 值函数 (退休后 gyp=1)
            value = (u_term + cS.beta * cS.survprob(t) * (expected_term)^(cS.psi_1 / (1 - cS.gamma)))^cS.psi_2;

            v = -value; % fmincon 最小化
        end

        %% =====================================================================
        %                VFI - 工作期值函数
        %  =====================================================================

        function v = fun_valuefunc_work(x, cash, fund, V_next_interp_obj, cS, t)
            % x = [绝对消费 C, 风险资产比例 alpha, 绝对养老金缴费 Q]
            consumption = x(1);
            alpha       = x(2);
            pension_contrib = x(3);

            if consumption <= 1e-8
                v = 1e10;
                return;
            end

            liquid_sav = cash - consumption - (1 - cS.tau_y) * pension_contrib;

            if liquid_sav < -1e-6
                v = 1e10;
                return;
            end

            % Epstein-Zin 瞬时效用项
            u_term = (1 - cS.beta) * (consumption^cS.psi_1);

            % 构造下一期状态的网格 (三维冲击)
            n = cS.n_shocks;
            gyp_3d = repmat(reshape(cS.gyp(:,:,t), [n, 1, n]), [1, n, 1]);
            gret_3d = repmat(reshape(cS.gret, [1, 1, n]), [n, n, 1]);
            yh_3d = repmat(reshape(cS.yh, [1, n, n]), [n, 1, 1]);

            % 下一期归一化状态变量
            portfolio_return = cS.rf * (1 - alpha) + gret_3d * alpha;
            cash_1 = portfolio_return .* (liquid_sav ./ gyp_3d) + (1-cS.tau_y) .* yh_3d;
            fund_1 = (fund + pension_contrib) * cS.rp ./ gyp_3d;

            % cash_1 = max(min(cash_1, cS.gcash(end)), cS.gcash(1));
            % fund_1 = max(min(fund_1, cS.gfund(end)), cS.gfund(1));

            % *** 修改: 使用 griddedInterpolant 对象进行插值 ***
            % 调用方式 F({Yq, Xq})
            int_V = V_next_interp_obj(cash_1, fund_1);
            int_V(int_V <= 0) = 1e-9;

            % 构造期望项 (包含归一化因子的增长)
            term_inside_exp = (gyp_3d .* int_V).^(1 - cS.gamma);
            expected_term = sum(cS.nweig1 .* term_inside_exp, 'all');

            % 完整的 Epstein-Zin 值函数
            value = (u_term + cS.beta * cS.survprob(t) * (expected_term)^(cS.psi_1 / (1 - cS.gamma)))^cS.psi_2;

            v = -value; % fmincon 最小化
        end


        function v = fun_valuefunc_work_riskrp(x, cash, fund, V_next_interp_obj, cS, t)
            % x = [绝对消费 C, 风险资产比例 alpha, 绝对养老金缴费 Q]
            consumption = x(1);
            alpha       = x(2);
            pension_contrib = x(3);

            if consumption <= 1e-8
                v = 1e10;
                return;
            end

            liquid_sav = cash - consumption - (1 - cS.tau_y) * pension_contrib;

            if liquid_sav < -1e-6
                v = 1e10;
                return;
            end

            % Epstein-Zin 瞬时效用项
            u_term = (1 - cS.beta) * (consumption^cS.psi_1);

            % 构造下一期状态的网格 (三维冲击)
            n = cS.n_shocks;
            gyp_3d = repmat(reshape(cS.gyp(:,:,t), [n, 1, n]), [1, n, 1]);
            gret_3d = repmat(reshape(cS.gret, [1, 1, n]), [n, n, 1]);
            gret_rp_3d = repmat(reshape(cS.gret_rp, [1, 1, n]), [n, n, 1]);
            yh_3d = repmat(reshape(cS.yh, [1, n, n]), [n, 1, 1]);

            % 下一期归一化状态变量
            portfolio_return = cS.rf * (1 - alpha) + gret_3d * alpha;
            cash_1 = portfolio_return .* (liquid_sav ./ gyp_3d) + (1-cS.tau_y) .* yh_3d;
            fund_1 = (fund + pension_contrib) * gret_rp_3d ./ gyp_3d;

            % cash_1 = max(min(cash_1, cS.gcash(end)), cS.gcash(1));
            % fund_1 = max(min(fund_1, cS.gfund(end)), cS.gfund(1));

            % *** 修改: 使用 griddedInterpolant 对象进行插值 ***
            % 调用方式 F({Yq, Xq})
            int_V = V_next_interp_obj(cash_1, fund_1);
            int_V(int_V <= 0) = 1e-9;

            % 构造期望项 (包含归一化因子的增长)
            term_inside_exp = (gyp_3d .* int_V).^(1 - cS.gamma);
            expected_term = sum(cS.nweig1 .* term_inside_exp, 'all');

            % 完整的 Epstein-Zin 值函数
            value = (u_term + cS.beta * cS.survprob(t) * (expected_term)^(cS.psi_1 / (1 - cS.gamma)))^cS.psi_2;

            v = -value; % fmincon 最小化
        end

        %% =====================================================================
        %                       3. 模拟器函数
        %  =====================================================================
        function simulate_model(results, cS)
            % --- 初始化 ---
            fprintf('模拟开始...\n');
            rng(42); % 可复现性
            nsim = 10000; % 模拟个体数量

            % 解包策略函数
            C_policy = results.C_policy;
            A_policy = results.A_policy;
            Q_policy = results.Q_policy;

            % 初始化模拟数组 (所有变量首先在归一化单位下计算)
            % 冲击和外生过程
            simGPY = ones(cS.tn, nsim);       % 持久收入增长因子 G_t * exp(z_t)
            simY_norm = zeros(cS.tn, nsim);    % 归一化暂时/基础收入 exp(u_t) 或 ret_fac
            simR = zeros(cS.tn, nsim);         % 风险回报率

            % 归一化单位下的决策和状态变量
            simW_norm = zeros(cS.tn, nsim);    % 归一化流动性财富 w_t = W_t/P_t
            simF_norm = zeros(cS.tn, nsim);    % 归一化养老金 f_t = F_t/P_t
            simC_norm = zeros(cS.tn, nsim);    % 归一化消费 c_t = C_t/P_t
            simQ_norm = zeros(cS.tn, nsim);    % 归一化缴费 q_t = Q_t/P_t
            simA = zeros(cS.tn, nsim);         % 风险资产配置 alpha_t
            simS_norm = zeros(cS.tn, nsim);    % 归一化风险资产 s_t
            simB_norm = zeros(cS.tn, nsim);    % 归一化无风险资产 b_t
            simQ_Cash_ratio = zeros(cS.tn, nsim); % 缴费与现金流比率 q/cash

            % --- 1. 模拟外生冲击和收入过程 ---
            fprintf('  生成收入和回报率路径...\n');
            work_periods = cS.tr - cS.tb;

            % 使用 Antithetic Variates 方法减少模拟误差
            for i1 = 1:floor(nsim/2)
                z_shocks = randn(work_periods, 1) * cS.smav;
                u_shocks_log = randn(work_periods, 1) * cS.smay;
                r_shocks = randn(cS.tn, 1) * cS.sigr;

                % 正向冲击路径
                simGPY(2:work_periods+1, i1) = cS.f_y(2:work_periods+1) ./ cS.f_y(1:work_periods) .* exp(z_shocks);
                simY_norm(1:work_periods, i1) = exp(u_shocks_log);
                simR(:, i1) = cS.rf + cS.mu + r_shocks;

                % 反向冲击路径 (镜像)
                i2 = nsim/2 + i1;
                simGPY(2:work_periods+1, i2) = cS.f_y(2:work_periods+1) ./ cS.f_y(1:work_periods) .* exp(-z_shocks);
                simY_norm(1:work_periods, i2) = exp(-u_shocks_log);
                simR(:, i2) = cS.rf + cS.mu - r_shocks;
            end

            % 退休后基础养老金 (替代率)
            norm_base_pension = cS.ret_fac ;
            simY_norm(work_periods+1:cS.tn, :) = norm_base_pension;

            % --- 2. 迭代模拟生命周期决策 (在归一化单位下) ---
            fprintf('  模拟生命周期决策 (归一化单位)...\n');
            for t = 1:cS.tn
                if t <= work_periods % 工作期
                    norm_cash = simW_norm(t, :) + (1-cS.tau_y) * simY_norm(t, :);
                else % 退休期
                    norm_pension_payout = (1-cS.tau_q) * simF_norm(t, :) * cS.pension_rate;
                    norm_cash = simW_norm(t, :) + simY_norm(t, :) + norm_pension_payout;
                end

                % 为当前期 t 的所有策略函数创建插值器
                Fc_interpolant = griddedInterpolant({cS.gcash, cS.gfund}, C_policy(:,:,t), 'spline', 'linear');
                Fa_interpolant = griddedInterpolant({cS.gcash, cS.gfund}, A_policy(:,:,t), 'spline', 'linear');
                Fq_interpolant = griddedInterpolant({cS.gcash, cS.gfund}, Q_policy(:,:,t), 'spline', 'linear');

                % 对所有 nsim 个体进行插值，得到决策
                simC_norm(t, :) = Fc_interpolant(norm_cash, simF_norm(t, :));
                simA(t, :) = Fa_interpolant(norm_cash, simF_norm(t, :));
                simQ_norm(t, :) = Fq_interpolant(norm_cash, simF_norm(t, :));

                % 应用约束和边界条件
                simA(t, :) = max(min(simA(t, :), 1), 0);
                simC_norm(t, :) = max(min(simC_norm(t, :), norm_cash), 0);
                simQ_norm(t, :) = max(min(simQ_norm(t, :), cS.Q_max), 0);

                if t <= work_periods % 工作期
                    total_outflow = simC_norm(t, :) + (1-cS.tau_y)*simQ_norm(t, :);
                    idx_exceed = total_outflow > norm_cash;
                    if any(idx_exceed)
                        scale_factor = norm_cash(idx_exceed) ./ total_outflow(idx_exceed) * 0.9999;
                        simC_norm(t, idx_exceed) = simC_norm(t, idx_exceed) .* scale_factor;
                        simQ_norm(t, idx_exceed) = simQ_norm(t, idx_exceed) .* scale_factor;
                    end
                    simQ_Cash_ratio(t, :) = simQ_norm(t, :) ./ norm_cash;
                    simQ_Cash_ratio(t, isnan(simQ_Cash_ratio(t,:)) | isinf(simQ_Cash_ratio(t,:))) = 0;
                    liquid_sav_norm = norm_cash - simC_norm(t, :) - (1-cS.tau_y)*simQ_norm(t, :);
                else % 退休期
                    simQ_norm(t,:) = 0;
                    simC_norm(t,:) = min(simC_norm(t,:), norm_cash * 0.9999);
                    liquid_sav_norm = norm_cash - simC_norm(t, :);
                end

                simS_norm(t, :) = simA(t, :) .* liquid_sav_norm;
                simB_norm(t, :) = liquid_sav_norm - simS_norm(t, :);

                % 更新下一期状态 (归一化单位)
                if t < cS.tn
                    portfolio_return_next = simB_norm(t, :) * cS.rf + simS_norm(t, :) .* simR(t, :);
                    simW_norm(t+1, :) = portfolio_return_next ./ simGPY(t+1, :);

                    if t < work_periods
                        simF_norm(t+1, :) = ((simF_norm(t, :) + simQ_norm(t, :)) * cS.rp) ./ simGPY(t+1, :);
                    elseif t == work_periods % 退休前最后一年
                        simF_norm(t+1, :) = ((simF_norm(t, :) + simQ_norm(t, :)) * cS.rp) ./ simGPY(t+1, :);
                    else % 退休后
                        simF_norm(t+1, :) = simF_norm(t, :);
                    end
                end
            end
            fprintf('  归一化单位模拟完成。\n');

            % --- 3. 将结果转换为以人民币计价的真实数值 ---
            fprintf('  将模拟结果转换为真实人民币单位...\n');

            % **锚定参数**
            initial_income_yuan = 60000; % 假设22岁时平均年持久性收入为 60,000 元

            % 初始化真实单位数组
            simP_real = zeros(cS.tn, nsim); % 真实持久性收入 (CNY)
            simY_real = zeros(cS.tn, nsim); % 真实总收入 (CNY)
            simW_real = zeros(cS.tn, nsim); % 真实流动性财富 (CNY)
            simF_real = zeros(cS.tn, nsim); % 真实养老金财富 (CNY)
            simC_real = zeros(cS.tn, nsim); % 真实消费 (CNY)
            simQ_real = zeros(cS.tn, nsim); % 真实养老金缴费 (CNY)

            % 模拟真实持久性收入路径 P_t
            simP_real(1, :) = initial_income_yuan; % 锚定初始值
            for t = 1:(cS.tn - 1)
                simP_real(t+1, :) = simP_real(t, :) .* simGPY(t+1, :);
            end

            % "去归一化" 所有相关变量
            for t = 1:cS.tn
                simW_real(t, :) = simW_norm(t, :) .* simP_real(t, :);
                simF_real(t, :) = simF_norm(t, :) .* simP_real(t, :);
                simC_real(t, :) = simC_norm(t, :) .* simP_real(t, :);
                simQ_real(t, :) = simQ_norm(t, :) .* simP_real(t, :);
                % 计算真实总收入
                if t <= work_periods
                    simY_real(t, :) = simP_real(t, :) .* simY_norm(t, :);
                else
                    simY_real(t, :) = simP_real(t, :);
                end
            end

            % --- 4. 汇总真实结果并绘图 (单位: 万元) ---
            fprintf('  汇总真实结果并绘图 (单位: 万元)...\n');
            yuan_to_wanyuan = 1 / 10000;
            ages = (cS.tb:cS.td)';
            x_fill = [ages; flipud(ages)]; % 用于填充置信区间的x坐标

            % --- 数据处理 ---
            % 均值
            meanC_plot = mean(simC_real, 2) * yuan_to_wanyuan;
            meanY_plot = mean(simY_real, 2) * yuan_to_wanyuan;
            meanTotalW_plot = mean(simW_real + simF_real, 2) * yuan_to_wanyuan;
            meanLiquidW_plot = mean(simW_real, 2) * yuan_to_wanyuan;
            meanF_plot = mean(simF_real, 2) * yuan_to_wanyuan;
            meanA_plot = nanmean(simA, 2);
            meanQ_plot = mean(simQ_real, 2) * yuan_to_wanyuan;

            % 标准差
            stdC_plot = std(simC_real, 0, 2) * yuan_to_wanyuan;
            stdY_plot = std(simY_real, 0, 2) * yuan_to_wanyuan;
            stdTotalW_plot = std(simW_real + simF_real, 0, 2) * yuan_to_wanyuan;
            stdLiquidW_plot = std(simW_real, 0, 2) * yuan_to_wanyuan;
            stdF_plot = std(simF_real, 0, 2) * yuan_to_wanyuan;
            stdA_plot = nanstd(simA, 0, 2);
            stdQ_plot = std(simQ_real, 0, 2) * yuan_to_wanyuan;

            % --- 绘图 ---
            figure('Position', [100, 100, 1400, 900]);

            % 子图 1: 生命周期剖面均值
            subplot(2,2,1);
            hold on;
            % 颜色定义
            color_total_w = [0, 0.4470, 0.7410]; % Blue
            color_liquid_w = [0.3010, 0.7450, 0.9330]; % Cyan
            color_c = [0.8500, 0.3250, 0.0980]; % Red
            color_y = [0, 0, 0]; % Black

            % 绘制置信区间 (fill)
            fill(x_fill, [meanTotalW_plot - stdTotalW_plot; flipud(meanTotalW_plot + stdTotalW_plot)], color_total_w, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
            fill(x_fill, [meanLiquidW_plot - stdLiquidW_plot; flipud(meanLiquidW_plot + stdLiquidW_plot)], color_liquid_w, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
            fill(x_fill, [meanC_plot - stdC_plot; flipud(meanC_plot + stdC_plot)], color_c, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
            fill(x_fill, [meanY_plot - stdY_plot; flipud(meanY_plot + stdY_plot)], color_y, 'FaceAlpha', 0.1, 'EdgeColor', 'none');

            % 绘制均值线
            plot(ages, meanTotalW_plot, 'Color', color_total_w, 'LineWidth', 2, 'DisplayName', '总财富 (流动+养老金)');
            plot(ages, meanLiquidW_plot, 'Color', color_liquid_w, 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', '流动性财富');
            plot(ages, meanC_plot, 'Color', color_c, 'LineWidth', 2, 'DisplayName', '消费');
            plot(ages, meanY_plot, 'Color', color_y, 'LineStyle', ':', 'LineWidth', 2, 'DisplayName', '收入');
            xline(cS.tr, '--', '退休', 'LineWidth', 1.5);
            hold off;
            title('生命周期剖面 (均值 \pm 1\sigma)'); xlabel('年龄'); ylabel('人民币 (万元)');
            legend('show', 'Location', 'northwest'); grid on; xlim([cS.tb, cS.td]);

            % 子图 2: 投资与养老金缴费决策
            subplot(2,2,2);
            hold on;
            % 颜色定义
            color_a = [0.4940, 0.1840, 0.5560]; % Purple
            color_q = [0.4660, 0.6740, 0.1880]; % Green

            % 左轴
            yyaxis left;
            % 置信区间
            fill(x_fill, [meanA_plot - stdA_plot; flipud(meanA_plot + stdA_plot)], color_a, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
            fill(x_fill, [meanQ_plot - stdQ_plot; flipud(meanQ_plot + stdQ_plot)], color_q, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
            % 均值线
            plot(ages, meanA_plot, 'Color', color_a, 'LineWidth', 2, 'DisplayName', '风险资产配置比例 (\alpha)');
            plot(ages(1:work_periods), meanQ_plot(1:work_periods), 'Color', color_q, 'LineWidth', 2, 'DisplayName', '个人养老金年缴费 (万元)');
            ylabel('比例 / 人民币 (万元)');
            ylim([-0.05, 1.05]);

            % 右轴
            yyaxis right;
            simQY_ratio_real = simQ_real ./ simY_real;
            simQY_ratio_real(isinf(simQY_ratio_real) | isnan(simQY_ratio_real)) = NaN;
            meanQY_ratio_plot = nanmean(simQY_ratio_real, 2);
            stdQY_ratio_plot = nanstd(simQY_ratio_real, 0, 2);
            % 置信区间
            fill(x_fill(1:2*work_periods), [meanQY_ratio_plot(1:work_periods) - stdQY_ratio_plot(1:work_periods); flipud(meanQY_ratio_plot(1:work_periods) + stdQY_ratio_plot(1:work_periods))], color_q, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
            % 均值线
            plot(ages(1:work_periods), meanQY_ratio_plot(1:work_periods), 'Color', color_q, 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', '缴费/收入比');
            ylabel('缴费收入比');

            hold off;
            xline(cS.tr, '--', '退休', 'LineWidth', 1.5);
            title('投资与养老金缴费决策 (均值 \pm 1\sigma)'); xlabel('年龄');
            grid on; xlim([cS.tb, cS.td]);
            legend('show', 'Location', 'best');

            % 子图 3: 个人养老金账户财富动态
            subplot(2,2,3);
            hold on;
            color_f = [0.9290, 0.6940, 0.1250]; % Orange
            % 置信区间
            fill(x_fill, [meanF_plot - stdF_plot; flipud(meanF_plot + stdF_plot)], color_f, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
            % 均值线
            plot(ages, meanF_plot, 'Color', color_f, 'LineWidth', 2, 'DisplayName', '养老金账户资产');
            xline(cS.tr, '--', '退休', 'LineWidth', 1.5);
            hold off;
            title('个人养老金账户财富动态 (均值 \pm 1\sigma)'); xlabel('年龄'); ylabel('人民币 (万元)');
            legend('show', 'Location', 'northwest'); grid on; xlim([cS.tb, cS.td]);

            % 子图 4: 总财富-收入比均值
            subplot(2,2,4);
            hold on;
            color_ratio = [0.6350, 0.0780, 0.1840]; % Maroon
            simWY_ratio_real = (simW_real + simF_real) ./ simY_real;
            simWY_ratio_real(isinf(simWY_ratio_real) | isnan(simWY_ratio_real)) = NaN;
            meanWY_ratio_plot = nanmean(simWY_ratio_real, 2);
            stdWY_ratio_plot = nanstd(simWY_ratio_real, 0, 2);
            % 置信区间
            fill(x_fill, [meanWY_ratio_plot - stdWY_ratio_plot; flipud(meanWY_ratio_plot + stdWY_ratio_plot)], color_ratio, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
            % 均值线
            plot(ages, meanWY_ratio_plot, 'Color', color_ratio, 'LineWidth', 2);
            xline(cS.tr, '--', '退休', 'LineWidth', 1.5);
            hold off;
            title('总财富-收入比 (均值 \pm 1\sigma)'); xlabel('年龄'); ylabel('比率');
            grid on; xlim([cS.tb, cS.td]);

            sgtitle('生命周期模拟结果 (以2024年人民币计价)', 'FontSize', 16, 'FontWeight', 'bold');

            output_filename = fullfile(cS.save_path, 'vfi_simulation_results_real_cny_with_ci.png');
            print(gcf, output_filename, '-dpng', '-r300');
            fprintf('真实单位模拟图形已保存到 %s\n', output_filename);
        end
        % =====================================================================
        %                辅助函数: 计算归一化的养老金缴费上限
        % =====================================================================
        function q_max_normalized = calculate_q_max(absolute_q_max_yuan, anchor_income_yuan, f_y_profile)
            % 描述:
            % 本函数根据现实世界的绝对缴费上限(如12000元)和代表性收入水平，
            % 计算出适用于模型归一化单位的缴费上限比率 cS.Q_max。
            %
            % 输入:
            %   absolute_q_max_yuan - 政策规定的年度个人养老金缴费绝对上限 (人民币)
            %   anchor_income_yuan  - 用于锚定模型收入水平的初始工作年收入 (人民币)
            %   f_y_profile         - 模型计算出的整个工作生涯的确定性收入剖面 (无单位)
            %
            % 输出:
            %   q_max_normalized    - 归一化后的缴费上限 (相对于平均持久性收入的比率)

            fprintf('  正在计算归一化的养老金缴费上限 (Q_max)...\n');

            % 1. 检查输入是否有效
            if isempty(f_y_profile) || f_y_profile(1) <= 0
                error('输入的收入剖面 f_y_profile 无效。');
            end

            % 2. 计算转换因子，将模型的抽象收入单位与人民币挂钩
            % 因子 = 人民币 / 模型单位
            yuan_per_model_unit = anchor_income_yuan / f_y_profile(1);
            fprintf('    锚定收入: 初始年收入 %.0f 元 对应模型单位 %.4f\n', anchor_income_yuan, f_y_profile(1));

            % 3. 计算整个工作生涯的确定性收入路径 (人民币计价)
            work_income_yuan = f_y_profile * yuan_per_model_unit;

            % 4. 计算工作生涯的平均年收入 (代表性的持久性收入水平)
            average_work_income_yuan = mean(work_income_yuan);
            fprintf('    模型隐含的工作生涯平均年收入为: %.0f 元\n', average_work_income_yuan);

            % 5. 计算归一化的 Q_max
            q_max_normalized = absolute_q_max_yuan / average_work_income_yuan;
            fprintf('    政策上限 %.0f 元 对应的归一化 Q_max 为: %.4f\n', absolute_q_max_yuan, q_max_normalized);

        end
        %  =====================================================================
        %  =====================================================================
        function [z_grid, P] = discretizeAR1_Tauchen(mew, rho, sigma, znum, Tauchen_q, ~)
            % Tauchen方法离散化AR(1)过程
            if znum == 1
                z_grid = mew / (1 - rho);
                P = 1;
                return;
            end

            zstar = mew / (1 - rho);
            sigmaz = sigma / sqrt(1 - rho^2);

            z_grid = zstar + linspace(-Tauchen_q * sigmaz, Tauchen_q * sigmaz, znum)';
            omega = z_grid(2) - z_grid(1);

            P = zeros(znum, znum);
            for i = 1:znum
                for j = 1:znum
                    if j == 1
                        P(i, j) = normcdf((z_grid(j) + omega/2 - rho * z_grid(i) - mew) / sigma);
                    elseif j == znum
                        P(i, j) = 1 - normcdf((z_grid(j) - omega/2 - rho * z_grid(i) - mew) / sigma);
                    else
                        P(i, j) = normcdf((z_grid(j) + omega/2 - rho * z_grid(i) - mew) / sigma) - ...
                            normcdf((z_grid(j) - omega/2 - rho * z_grid(i) - mew) / sigma);
                    end
                end
            end
        end


    end

end