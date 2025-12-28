classdef CoccoEnv < rl.env.MATLABEnvironment
    % 描述:
    % Cocco et al. (2005) 生命周期模型的 MATLAB 强化学习环境。
    % **核心重构**:
    % 1. 状态变量明确为 X_cash_on_hand，与论文保持一致。
    % 2. 严格按照论文定义的时序和公式进行状态转移。

    properties
        % --- A. 生命周期与人口结构 ---
        tb = 22; tr = 26; td = 30; tn = 9;

        % --- B. 经济主体偏好 ---
        beta = 0.95; gamma = 3.84;

        % --- C. 收入过程参数 ---
        aa = -2.170042 + 2.700381; b1 = 0.16818;
        b2 = -0.0323371 / 10; b3 = 0.0019704 / 100;
        smay = 0; smav = 0;

        % --- D. 资产与回报率 ---
        rf = 1.02; mu = 0.04;
        sigr = 0;

        % --- E. 养老金与税收 ---
        ret_fac = 0.6827; tau_y = 0;

        % --- F. 生存概率与确定性收入 ---
        survprob; f_y;

        % --- G. 离散化随机过程 ---
        n_shocks = 5;
        shock_grid; shock_weig; gret; yh; yp;

        % --- H. 状态变量 ---
        % **核心修正**: 状态变量从 W (财富) 改为 X (现金)
        X_cash_on_hand (1,1) double = 0.0;
        P (1,1) double = 0.0;
        age (1,1) double = 22;
        P_retirement (1,1) double = 0.0;

        % --- I. 模拟控制 ---
        disable_random_death (1,1) logical = false
        is_done (1,1) logical = false
    end

    methods
        % --- 构造函数 ---
        function this = CoccoEnv(varargin)
            % **核心修正**: 观测的第一个维度现在是 X (Cash-on-hand)
            obsInfo = rlNumericSpec([2 1], 'LowerLimit', [0;  -1], 'UpperLimit', [inf; 1]);
            obsInfo.Name = 'Cocco State (Cash-on-hand, P, Age)';
            actInfo = rlNumericSpec([2 1], 'LowerLimit', [0.1; -0.1], 'UpperLimit', [1.1; 1.1]);
            actInfo.Name = 'Cocco Action (c_prop, alpha)';
            this@rl.env.MATLABEnvironment(obsInfo, actInfo);

            p = inputParser;
            addParameter(p, 'disable_random_death', false);
            parse(p, varargin{:});
            this.disable_random_death = p.Results.disable_random_death;

            this.tn = this.td - this.tb + 1;
            this.precompute_survival_and_income();
            this.precompute_discrete_shocks();
        end

        % --- 重置函数 ---
        function [InitialObservation, info] = reset(this)
            this.age = this.tb;
            this.is_done = false;
            this.P_retirement = 0;

            % 1. 初始化 P_1
            idx_z = randsample(this.n_shocks, 1, true, this.shock_weig);
            this.P = this.yp(idx_z);

            % 2. 初始财富 W_1 = 0
            W1 = 0;

            % 3. 实现初始收入 Y_1
            idx_u = randsample(this.n_shocks, 1, true, this.shock_weig);
            Y1 = exp(this.f_y(1) + this.P + this.yh(idx_u));

            % 4. 初始化状态变量 X_1 = W_1 + Y_1
            this.X_cash_on_hand = W1 + Y1;

            InitialObservation = this.getObs();
            info = struct('AbsoluteIncome', Y1, 'AbsoluteWealth', W1);
        end

        % --- 步进函数 ---
        function [NextObs, Reward, IsDone, LoggedSignal] = step(this, Action)
            c_prop = max(0, min(1, Action(1)));
            alpha  = max(0, min(1, Action(2)));
            
            % 当前状态为 (X_t, P_t, age_t)
            current_X = this.X_cash_on_hand;
            current_age_idx = this.age - this.tb + 1;

            % 1. 智能体决策: 计算消费 C_t 和储蓄 S_t
            C = c_prop * current_X;
            C = min(C, current_X * 0.99999);
            if C < 1e-6, C = 1e-6; end
            S = current_X - C;

            % 2. 计算奖励
            utility = (C^(1 - this.gamma)) / (1 - this.gamma)+1;
            if this.gamma == 1, utility = log(C); end
            Reward = max(utility, -10);

            % 3. 计算下一期期初财富 W_{t+1}
            idx_r = randsample(this.n_shocks, 1, true, this.shock_weig);
            portfolio_return = alpha * this.gret(idx_r) + (1 - alpha) * this.rf;
            W_next = S * portfolio_return;

            % 4. 演进状态 age 和 P 到 t+1
            next_age = this.age + 1;
            idx_z = randsample(this.n_shocks, 1, true, this.shock_weig);
            if next_age <= this.tr
                this.P = this.P + this.yp(idx_z);
                if next_age == this.tr
                    this.P_retirement = this.P;
                end
            end
            this.age = next_age;

            % 5. 实现下一期收入 Y_{t+1}
            Y_next = 0;
            idx_u = randsample(this.n_shocks, 1, true, this.shock_weig);
            next_age_idx = this.age - this.tb + 1;
            if this.age <= this.tr
                Y_next = exp(this.f_y(next_age_idx) + this.P + this.yh(idx_u));
            else
                Y_next = this.ret_fac * exp(this.f_y(this.tr - this.tb + 1) + this.P_retirement);
            end

            % 6. 更新下一期状态：现金 X_{t+1} = W_{t+1} + Y_{t+1}
            this.X_cash_on_hand = W_next + Y_next;

            % 7. 检查终止条件
            died = false;
            if this.age > this.td
                died = true;
            elseif ~this.disable_random_death && rand() > this.survprob(current_age_idx)
                died = true;
            end
            this.is_done = died;

            IsDone = this.is_done;
            NextObs = this.getObs();
            
            % 为日志信号计算当期收入 Y_t
            % Y_t = X_t - W_t = X_t - S_{t-1}*R_t (这很复杂)
            % 为了简单，我们只记录下一期的财富和收入
            LoggedSignal = struct('C_prop', C/current_X, 'AbsoluteConsumption', C, ...
                                'AbsoluteIncome_next', Y_next, 'AbsoluteWealth_next', W_next, ...
                                'RiskyShare', alpha);
        end

        % --- 辅助函数: 获取当前观测 ---
        function Obs = getObs(this)
            normalized_age = 2 * (this.age - this.tb) / (this.td - this.tb) - 1;
            Obs = [this.X_cash_on_hand; normalized_age];
        end
    end

    methods (Access = private)
        % precompute_discrete_shocks, precompute_survival_and_income, discretizeAR1_Tauchen
        % 这三个函数保持不变，无需修改
        function precompute_discrete_shocks(this)
            [grid, weig_matrix] = this.discretizeAR1_Tauchen(0, 0, 1, this.n_shocks, 2);
            this.shock_grid = grid;
            this.shock_weig = diag(weig_matrix);
            this.gret = this.rf + this.mu + this.shock_grid * this.sigr;
            this.yh = this.shock_grid * this.smay;
            this.yp = this.shock_grid * this.smav;
        end
        function precompute_survival_and_income(this)
            survprob_data = [0.99975, 0.99974, 0.99973, 0.99972, 0.99971, 0.99969, 0.99968, 0.99966, 0.99963, 0.99961, 0.99958, 0.99954, 0.99951, 0.99947, 0.99943, 0.99938, 0.99934, 0.99928, 0.99922, 0.99916, 0.99908, 0.999, 0.99891, 0.99881, 0.9987, 0.99857, 0.99843, 0.99828, 0.99812, 0.99794, 0.99776, 0.99756, 0.99735, 0.99713, 0.9969, 0.99666, 0.9964, 0.99613, 0.99583, 0.99551, 0.99515, 0.99476, 0.99432, 0.99383, 0.9933, 0.9927, 0.99205, 0.99133, 0.99053, 0.98961, 0.98852, 0.98718, 0.98553, 0.98346, 0.98089, 0.97772, 0.97391, 0.96943, 0.96429, 0.95854, 0.95221, 0.94537, 0.93805, 0.93027, 0.92202, 0.91327, 0.90393, 0.89389, 0.88304, 0.87126, 0.85846, 0.84452, 0.82935, 0.81282, 0.79485, 0.77543, 0.75458, 0.7324, 0.70893, 0.68424, 0.68424, 0.68424];
            this.survprob = ones(this.tn - 1, 1);
            copy_len = min(length(this.survprob), length(survprob_data));
            this.survprob(1:copy_len) = survprob_data(1:copy_len);
            this.f_y = zeros(this.tn, 1);
            for i = 1:this.tn
                age_loop = this.tb + i - 1;
                if age_loop <= this.tr
                    this.f_y(i) = (this.aa + this.b1*age_loop + this.b2*age_loop^2 + this.b3*age_loop^3)-3;
                end
            end
        end
        function [z_grid, P] = discretizeAR1_Tauchen(~, mew, rho, sigma, znum, Tauchen_q)
            if znum == 1, z_grid = mew / (1 - rho); P = 1; return; end
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