function household()
%% 生命周期消费/投资/养老金缴费模型
% 基于 Cocco, Gomes, and Maenhout (2005) 的框架
% 扩展模型包含一个个人养老金账户 (IRA)

clear all;
close all;

% --- 模型参数设定 ---
disp('1. 设置模型参数...');
p = set_parameters();

% --- 状态变量网格和冲击过程离散化 ---
disp('2. 构建网格和离散化冲击过程...');
g = create_grids_and_shocks(p);

% --- 价值函数迭代 (VFI) ---
disp('3. 开始价值函数迭代 (逆向归纳法)...');
vfi_results = value_function_iteration(p, g);

% --- 生命周期模拟 ---
disp('4. 进行生命周期模拟...');
simulation_results = life_cycle_simulation(p, g, vfi_results);

% --- 可视化结果 (可选) ---
disp('5. 绘制结果...');
plot_results(p, simulation_results);

disp('计算完成。');

end

%% 1. 参数设置函数
function p = set_parameters()
    % 年龄和生命周期
    p.tb = 20;         % 工作起始年龄
    p.tr = 66;         % 退休年龄
    p.td = 100;         % 最高年龄
    p.tn = p.td - p.tb + 1; % 总期数
    p.work_life = p.tr - p.tb; % 工作期长度

    % 偏好
    p.rho = 10.0;      % 相对风险厌恶系数 (gamma in paper)
    p.psi = 0.5;       % 跨期替代弹性 (EIS)
    p.delta = 0.97;    % 时间折扣因子 (beta in paper)

    % 资产收益
    p.r = 1.015;       % 无风险利率 (Rf in paper)
    p.mu = 0.04;       % 风险资产平均超额收益
    p.sigr = 0.20;     % 风险资产收益标准差
    p.Rp = 1.015;      % 养老金账户投资回报率 (设为无风险利率)

    % 收入过程 (log-income: ln(Y_t) = f(t) + p_t + u_t, p_t = p_{t-1} + z_t)
    p.smay = 0.1;      % 暂时性收入冲击的标准差 (sigma_u)
    p.smav = 0.1;      % 持久性收入冲击的标准差 (sigma_z)
    p.corr_y = 0.0;    % 收入冲击与股票收益的相关性 (设为0)
    % 确定性收入曲线参数 (来自 life_cycle.m)
    p.aa = -2.170042 + 2.700381;
    p.b1 = 0.16818;
    p.b2 = -0.0323371 / 10;
    p.b3 = 0.0019704 / 100;
    p.ret_fac = 0.68212; % 退休后收入替代率 (相对于退休前最后一个时期的持久性收入)

    % 政策参数 (新加入)
    p.tau_y = 0.20;    % 劳动收入税率
    p.tau_q = 0.15;    % 养老金提取税率
    p.Q_max = 0.10;    % 个人养老金缴费率上限

    % 数值方法参数
    % 状态变量网格点数
    p.nw = 41;         % 标准化流动性财富 w = W/P
    p.nf = 21;         % 标准化养老金财富 f = F/P
    % 决策变量网格点数
    p.nc = 21;         % 消费
    p.nq = 5;          % 养老金缴费率
    p.na = 11;         % 风险资产配置比例
    % 随机冲击离散化节点数
    p.n_perm = 5;      % 持久性收入冲击 z
    p.n_tran = 5;      % 暂时性收入冲击 u
    p.n_ret = 5;       % 资产收益冲击

    % 生存概率 (来自 life_cycle.m)
    survprob_data = [0.99845,0.99839,0.99833,0.9983,0.99827,0.99826,0.99824,0.9982,0.99813,0.99804,0.99795,0.99785,0.99776,0.99766,0.99755,0.99743,0.9973,0.99718,0.99707,0.99696,0.99685,0.99672,0.99656,0.99635,0.9961,0.99579,0.99543,0.99504,0.99463,0.9942,0.9937,0.99311,0.99245,0.99172,0.99091,0.99005,0.98911,0.98803,0.9868,0.98545,0.98409,0.9827,0.98123,0.97961,0.97786,0.97603,0.97414,0.97207,0.9697,0.96699,0.96393,0.96055,0.9569,0.9531,0.94921,0.94508,0.94057,0.9357,0.93031,0.92424,0.91717,0.90922,0.90089,0.89282,0.88503,0.87622,0.86576,0.8544,0.8423,0.82942,0.8154,0.80002,0.78404,0.76842,0.75382,0.73996,0.72464,0.71057,0.6961,0.6809];
    p.survprob = survprob_data(1:p.tn-1)';

    % 效用函数相关计算
    p.psi_1 = 1.0 - 1.0 / p.psi;
    if abs(p.psi_1) < 1e-9, p.psi_1 = 1e-9; end % 避免除以零
    p.psi_2 = 1.0 / p.psi_1;
    p.theta = (1.0 - p.rho) / p.psi_1;

    % 模拟参数
    p.nsim = 5000; % 模拟个体数量
end

%% 2. 网格和冲击过程生成函数
function g = create_grids_and_shocks(p)
    % 状态变量网格 (w, f)
    w_max = 20; w_min = 0;
    g.w_grid = linspace(w_min, w_max, p.nw)';
    f_max = 20; f_min = 0;
    g.f_grid = linspace(f_min, f_max, p.nf)';
    [g.W_mesh, g.F_mesh] = meshgrid(g.w_grid, g.f_grid); % 用于插值

    % 决策变量网格 (alpha, q)
    g.alpha_grid = linspace(0, 1, p.na)';
    g.q_grid = linspace(0, p.Q_max, p.nq)';

    % 冲击离散化 (使用 tauchen 函数)
    % 持久性收入冲击 z_t, N(0, smav^2), iid (rho=0)
    [z_nodes, z_prob_mat] = tauchen(p.n_perm, 0, p.smav, 0, 3);
    z_prob = z_prob_mat(1, :)'; % 对 iid 冲击, 每行概率都一样
    % 暂时性收入冲击 u_t, N(0, smay^2), iid (rho=0)
    [u_nodes, u_prob_mat] = tauchen(p.n_tran, 0, p.smay, 0, 3);
    u_prob = u_prob_mat(1, :)';
    % 风险资产收益冲击 eps_r, N(0, 1), iid (rho=0)
    [ret_shk_nodes, ret_shk_prob_mat] = tauchen(p.n_ret, 0, 1, 0, 3);
    ret_shk_prob = ret_shk_prob_mat(1, :)';

    % 组合所有冲击
    [z_idx, u_idx, r_idx] = ndgrid(1:p.n_perm, 1:p.n_tran, 1:p.n_ret);
    g.shock_nodes = [z_nodes(z_idx(:)), exp(u_nodes(u_idx(:))), ret_shk_nodes(r_idx(:))];
    
    % 构造联合概率分布向量 (125x1)
    g.shock_probs = kron(ret_shk_prob, kron(u_prob, z_prob));
    g.n_shocks = size(g.shock_nodes, 1);
    
    % 确定性收入增长率 gy
    f_y = zeros(p.work_life + 1, 1);
    for age_idx = 1:(p.work_life + 1)
        age = p.tb + age_idx - 1;
        f_y(age_idx) = exp(p.aa + p.b1*age + p.b2*age^2 + p.b3*age^3);
    end
    g.gy = f_y(2:end) ./ f_y(1:end-1);

    % 风险资产收益率 R = Rf + mu + eps_r*sigr
    g.R_tilde = p.r + p.mu + g.shock_nodes(:,3) * p.sigr;
end
%% 3. 价值函数迭代
function vfi_results = value_function_iteration(p, g)
    % --- 初始化 ---
    % 价值函数 V(f_idx, w_idx, t)
    V = zeros(p.nf, p.nw, p.tn + 1);
    % 政策函数 C, Q, alpha
    vfi_results.C_policy = zeros(p.nf, p.nw, p.tn);
    vfi_results.Q_policy = zeros(p.nf, p.nw, p.tn);
    vfi_results.alpha_policy = zeros(p.nf, p.nw, p.tn);

    % --- 退休期 (t = tn, ..., work_life + 1) ---
    % 状态: (w_t = W_t/P_{tr-1}, f_K = F_{tr}/P_{tr-1})
    P_bar_norm_vec = (1 - p.tau_q) * g.f_grid / (p.td - p.tr + 1);
    
    % t = tn (td), 最终期价值函数 (基于Epstein-Zin偏好)
    t = p.tn;
    age = p.td;
    fprintf('求解 VFI: 年龄 %d (t=%d)\n', age, t);
    for f_idx = 1:p.nf
        cash_on_hand_norm = g.w_grid + p.ret_fac + P_bar_norm_vec(f_idx);
        % V_T = ((1-delta)*C_T^psi_1)^psi_2 = C_T * (1-delta)^psi_2
        V(f_idx, :, t) = cash_on_hand_norm * (1 - p.delta)^p.psi_2;
    end

    % t = tn-1, ..., work_life + 1
    for t = p.tn-1:-1:p.work_life + 1
        age = p.tb + t - 1;
        fprintf('求解 VFI: 年龄 %d (t=%d)\n', age, t);
        
        V_next_interp = griddedInterpolant({g.f_grid, g.w_grid}, V(:,:,t+1), 'linear', 'linear');
        
        for f_idx = 1:p.nf
            P_bar_norm = P_bar_norm_vec(f_idx);
            for w_idx = 1:p.nw
                cash_on_hand_norm = g.w_grid(w_idx) + p.ret_fac + P_bar_norm;
                
                max_v = -inf;
                best_c = 0;
                best_alpha = 0;
                
                c_min = 1e-6;
                c_max = cash_on_hand_norm;
                c_grid_ind = linspace(c_min, c_max, p.nc);

                for c_idx = 1:p.nc
                    c_val = c_grid_ind(c_idx);
                    sav_norm = cash_on_hand_norm - c_val;
                    if sav_norm < 0, continue; end
                    
                    for a_idx = 1:p.na
                        alpha_val = g.alpha_grid(a_idx);
                        
                        R_port = (1 - alpha_val) * p.r + alpha_val * g.R_tilde;
                        w_next = sav_norm * R_port;
                        f_next = g.f_grid(f_idx) * ones(g.n_shocks, 1);
                        
                        V_next_points = V_next_interp(f_next, w_next);
                        V_next_points(isnan(V_next_points)) = 1e-9; % 处理外插点
                        V_next_points = max(V_next_points, 1e-9);   % V必须为正

                        % Epstein-Zin 价值函数计算
                        E_V_term = g.shock_probs' * (V_next_points .^ (1 - p.rho));
                        u = (1-p.delta) * (c_val ^ p.psi_1);
                        certainty_equivalent = (p.survprob(t) * E_V_term)^(1/p.theta);
                        current_v = (u + p.delta * certainty_equivalent)^p.psi_2;
                        auxVV(a_idx,c_idx) = current_v;
                        if current_v > max_v
                            max_v = current_v;
                            best_c = c_val;
                            best_alpha = alpha_val;
                        end
                    end
                end
                
                V(f_idx, w_idx, t) = max_v;
                vfi_results.C_policy(f_idx, w_idx, t) = best_c;
                vfi_results.alpha_policy(f_idx, w_idx, t) = best_alpha;
                vfi_results.Q_policy(f_idx, w_idx, t) = 0;
            end
        end
    end
    
    % --- 工作期 (t = work_life, ..., 1) ---
    for t = p.work_life:-1:1
        age = p.tb + t - 1;
        fprintf('求解 VFI: 年龄 %d (t=%d)\n', age, t);

        V_next_interp = griddedInterpolant({g.f_grid, g.w_grid}, V(:,:,t+1), 'linear', 'none');

        for f_idx = 1:p.nf
            for w_idx = 1:p.nw
                max_v = -inf;
                best_c = 0; best_q = 0; best_alpha = 0;
                y_t = 1.0; % 简化：在暂时性收入冲击的均值处求解决策

                for q_idx = 1:p.nq
                    q_val = g.q_grid(q_idx);
                    cash_on_hand_norm = g.w_grid(w_idx) + (1 - p.tau_y) * (y_t - q_val * y_t);
                    if cash_on_hand_norm <= 0, continue; end
                    
                    c_min = 1e-6;
                    c_max = cash_on_hand_norm;
                    c_grid_ind = linspace(c_min, c_max, p.nc);

                    for c_val = c_grid_ind
                        sav_norm = cash_on_hand_norm - c_val;
                        if sav_norm < 0, continue; end
                        
                        for a_idx = 1:p.na
                            alpha_val = g.alpha_grid(a_idx);
                            
                            G_next = g.gy(t) * exp(g.shock_nodes(:,1));
                            R_port = (1 - alpha_val) * p.r + alpha_val * g.R_tilde;
                            w_next = sav_norm * R_port ./ G_next;
                            f_next = (g.f_grid(f_idx) * p.Rp + q_val * y_t) ./ G_next;
                            
                            V_next_points = V_next_interp(f_next, w_next);
                            V_next_points(isnan(V_next_points)) = 1e-9;
                            V_next_points = max(V_next_points, 1e-9);

                            % Epstein-Zin 价值函数计算 (包含收入增长)
                            E_V_term = g.shock_probs' * ((V_next_points .* G_next) .^ (1 - p.rho));
                            u = (1-p.delta) * (c_val ^ p.psi_1);
                            certainty_equivalent = (p.survprob(t) * E_V_term)^(1/p.theta);
                            current_v = (u + p.delta * certainty_equivalent)^p.psi_2;
                            
                            if current_v > max_v
                                max_v = current_v;
                                best_c = c_val;
                                best_q = q_val;
                                best_alpha = alpha_val;
                            end
                        end
                    end
                end 
                
                V(f_idx, w_idx, t) = max_v;
                vfi_results.C_policy(f_idx, w_idx, t) = best_c;
                vfi_results.Q_policy(f_idx, w_idx, t) = best_q;
                vfi_results.alpha_policy(f_idx, w_idx, t) = best_alpha;

            end
        end
    end
end
% 辅助函数: 找到在y_t=1(暂时性冲击为均值)时的最优策略，用于模拟
function [best_c, best_q, best_alpha] = find_policy_for_sim(w_idx, f_idx, V_next_interp, p, g, t)
    max_v = -inf;
    best_c = 0; best_q = 0; best_alpha = 0;
    y_t = 1.0; % 暂时性收入冲击的均值

    for q_idx = 1:p.nq
        q_val = g.q_grid(q_idx);
        cash_on_hand_norm = g.w_grid(w_idx) + (1 - p.tau_y) * (y_t - q_val * y_t);
        if cash_on_hand_norm <= 0, continue; end
        c_min = 1e-6; c_max = cash_on_hand_norm;
        c_grid_ind = linspace(c_min, c_max, p.nc);

        for c_val = c_grid_ind
            sav_norm = cash_on_hand_norm - c_val;
            if sav_norm < 0, continue; end
            for a_idx = 1:p.na
                alpha_val = g.alpha_grid(a_idx);
                G_next = g.gy(t) * exp(g.shock_nodes(:,1));
                R_port = (1 - alpha_val) * p.r + alpha_val * g.R_tilde;
                w_next = sav_norm * R_port ./ G_next;
                f_next = (g.f_grid(f_idx) * p.Rp + q_val * y_t) ./ G_next;
                V_next_points = V_next_interp(f_next, w_next);
                V_next_points(isnan(V_next_points)) = utility(w_next(isnan(V_next_points)), p);
                E_V_next = g.shock_probs' * V_next_points;
                current_v = utility(c_val, p) + p.delta * p.survprob(t) * E_V_next;
                if current_v > max_v
                    max_v = current_v;
                    best_c = c_val;
                    best_q = q_val;
                    best_alpha = alpha_val;
                end
            end
        end
    end
end


%% 4. 生命周期模拟
function sim = life_cycle_simulation(p, g, vfi_results)
    % --- 初始化模拟矩阵 ---
    sim.W = zeros(p.tn, p.nsim); % 流动性财富
    sim.F = zeros(p.tn, p.nsim); % 养老金财富
    sim.P = ones(p.tn, p.nsim);  % 持久性收入
    sim.Y = zeros(p.tn, p.nsim);  % 总收入 Y = P * u
    sim.C = zeros(p.tn, p.nsim);  % 消费
    sim.Q = zeros(p.tn, p.nsim);  % 养老金缴费额
    sim.alpha = zeros(p.tn, p.nsim); % 风险资产比例

    % --- 生成冲击序列 ---
    rng('default'); % for reproducibility
    z_shocks = randn(p.tn, p.nsim) * p.smav;
    u_shocks = exp(randn(p.tn, p.nsim) * p.smay - 0.5 * p.smay^2); % E[u]=1
    r_shocks = randn(p.tn, p.nsim);

    % --- 创建插值对象 ---
    C_interp = griddedInterpolant({g.f_grid, g.w_grid, (1:p.tn)'}, vfi_results.C_policy, 'linear');
    Q_interp = griddedInterpolant({g.f_grid, g.w_grid, (1:p.tn)'}, vfi_results.Q_policy, 'linear');
    alpha_interp = griddedInterpolant({g.f_grid, g.w_grid, (1:p.tn)'}, vfi_results.alpha_policy, 'linear');
    
    % --- 时间向前推进 ---
    for t = 1:p.tn-1
        fprintf('模拟... t=%d\n', t);
        age = p.tb + t - 1;
        
        % 确定持久性收入 P
        if t > 1
            sim.P(t, :) = sim.P(t-1, :) .* g.gy(t-1) .* exp(z_shocks(t,:));
        end
        
        for i = 1:p.nsim
            if age < p.tr % 工作期
                % 状态变量标准化
                w_norm = sim.W(t,i) / sim.P(t,i);
                f_norm = sim.F(t,i) / sim.P(t,i);

                % 插值得到决策
                c_norm = C_interp(f_norm, w_norm, t);
                q_rate = Q_interp(f_norm, w_norm, t);
                alpha = alpha_interp(f_norm, w_norm, t);
                
                % 边界检查
                c_norm = max(0, c_norm);
                q_rate = max(0, min(p.Q_max, q_rate));
                alpha = max(0, min(1.0, alpha));

                % 真实收入
                Y_real = sim.P(t,i) * u_shocks(t,i);
                
                % 真实决策变量
                C_real = c_norm * sim.P(t,i);
                Q_real = q_rate * Y_real;

                % 确保消费可行
                cash_after_tax_pension = sim.W(t,i) + (1-p.tau_y)*(Y_real - Q_real);
                C_real = min(C_real, cash_after_tax_pension * 0.999);
                
                % 储蓄
                saving = cash_after_tax_pension - C_real;
                
                % 更新下一期财富
                R_risky = p.r + p.mu + r_shocks(t,i) * p.sigr;
                R_port = (1-alpha)*p.r + alpha*R_risky;
                
                sim.W(t+1,i) = saving * R_port;
                sim.F(t+1,i) = (sim.F(t,i) + Q_real) * p.Rp;

                % 存储当期信息
                sim.Y(t,i) = Y_real;
                sim.C(t,i) = C_real;
                sim.Q(t,i) = Q_real;
                sim.alpha(t,i) = alpha;

            else % 退休期
                if age == p.tr % 刚退休，计算年金
                   sim.P_bar(i) = (1-p.tau_q) * sim.F(t,i) / (p.td - p.tr + 1);
                   sim.P_retire(i) = sim.P(t,i); % 锁定退休时的持久性收入
                end
                
                P_retire_norm = sim.P_retire(i);
                % 状态变量标准化 (用退休时的P)
                w_norm = sim.W(t,i) / P_retire_norm;
                f_norm = sim.F(t,i) / P_retire_norm;

                % 插值得到决策
                c_norm = C_interp(f_norm, w_norm, t);
                alpha = alpha_interp(f_norm, w_norm, t);

                % 边界检查
                c_norm = max(0, c_norm);
                alpha = max(0, min(1.0, alpha));

                % 真实决策变量
                C_real = c_norm * P_retire_norm;
                
                % 收入
                Y_retire_real = p.ret_fac * P_retire_norm;
                
                cash_on_hand = sim.W(t,i) + Y_retire_real + sim.P_bar(i);
                C_real = min(C_real, cash_on_hand * 0.999);
                
                saving = cash_on_hand - C_real;

                % 更新下一期财富
                R_risky = p.r + p.mu + r_shocks(t,i) * p.sigr;
                R_port = (1-alpha)*p.r + alpha*R_risky;
                
                sim.W(t+1,i) = saving * R_port;
                sim.F(t+1,i) = sim.F(t,i); % 养老金账户在年金化后不再变化

                % 存储当期信息
                sim.Y(t,i) = Y_retire_real + sim.P_bar(i);
                sim.C(t,i) = C_real;
                sim.Q(t,i) = 0;
                sim.alpha(t,i) = alpha;
            end
        end
    end
end


%% 辅助函数: 效用函数
function u = utility(c, p)
    c = max(c, 1e-9); % 避免c为0或负
    if abs(p.psi_1) < 1e-8 % 对数效用
        u = log(c);
    else
        u = (c .^ p.psi_1) / p.psi_1;
    end
end

%% Tauchen 函数 (由用户提供)
function [y_grid_out, trProbM_out] = tauchen(N, rho, sigma, mu, m)
    if N == 1, y_grid_out = 0; trProbM_out = 1; return; end
    std_y = sqrt(sigma^2 / (1-rho^2)); 
    if abs(1-rho^2) < 1e-9, std_y = sigma; end % i.i.d case
    y_max = m*std_y; y_min = -y_max;
    y = linspace(y_min, y_max, N); d = y(2)-y(1);
    trProbM_out = zeros(N,N);
    for j=1:N
        for k=1:N
            m_k = rho*y(j) + mu;
            if k==1
                trProbM_out(j,k) = normcdf((y(1)-m_k+d/2)/sigma);
            elseif k==N
                trProbM_out(j,k) = 1 - normcdf((y(N)-m_k-d/2)/sigma);
            else
                trProbM_out(j,k) = normcdf((y(k)-m_k+d/2)/sigma) - normcdf((y(k)-m_k-d/2)/sigma);
            end
        end
    end
    y_grid_out = y(:);
end

%% 可视化结果
function plot_results(p, sim)
    ages = p.tb : p.td;
    
    % 计算均值
    mean_W = mean(sim.W, 2);
    mean_F = mean(sim.F, 2);
    mean_C = mean(sim.C, 2);
    mean_Y = mean(sim.Y, 2);
    mean_alpha = mean(sim.alpha, 2);
    mean_Q_rate = mean(sim.Q ./ sim.Y, 2, 'omitnan');

    figure('Name', '生命周期财富、收入与消费');
    plot(ages, mean_W, 'b-', 'LineWidth', 2);
    hold on;
    plot(ages, mean_F, 'r-', 'LineWidth', 2);
    plot(ages, mean_C, 'g--', 'LineWidth', 2);
    plot(ages, mean_Y, 'k:', 'LineWidth', 2);
    xline(p.tr, '--', '退休');
    title('生命周期路径 (均值)');
    xlabel('年龄');
    ylabel('金额');
    legend('流动性财富 (W)', '养老金财富 (F)', '消费 (C)', '收入 (Y)');
    grid on;

    figure('Name', '生命周期投资与缴费决策');
    subplot(2,1,1);
    plot(ages, mean_alpha, 'LineWidth', 2);
    xline(p.tr, '--', '退休');
    title('风险资产配置比例 (\alpha)');
    xlabel('年龄');
    ylabel('比例');
    grid on;

    subplot(2,1,2);
    plot(ages(1:p.work_life), mean_Q_rate(1:p.work_life), 'LineWidth', 2);
    title('养老金缴费率 (Q/Y)');
    xlabel('年龄');
    ylabel('比例');
    grid on;

end