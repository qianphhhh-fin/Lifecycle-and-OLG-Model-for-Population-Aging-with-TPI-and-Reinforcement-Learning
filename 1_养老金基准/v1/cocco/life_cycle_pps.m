%% Life-Cycle Consumption/Portfolio/Pension Choice Problem
% Based on Gomes (2020) and extended with a Private Pension System (PPS)
% Solves the model described in the provided TeX documentation.

clear all;
close all;
warning off
%% ----------------------------------------------------------------
%  1. Model Parameters
%  ----------------------------------------------------------------
output_dir = 'results';
file_path = fullfile(output_dir, 'vfi_results.mat');
% Life-cycle parameters
tb     = 20;         % Beginning of economic life
tr     = 66;         % Retirement age
td     = 100;        % Maximum age
tn     = td-tb+1;    % Total periods
work_life = tr-tb;   % Length of working life

% Preference parameters
delta       = 0.97;     % Subjective discount factor (beta in TeX)
rho         = 10.0;     % Coefficient of relative risk aversion (gamma in TeX)
psi         = 0.5;      % Elasticity of intertemporal substitution (EIS)

% Income process parameters (matches Gomes)
aa          = -2.170042+2.700381; % Constant term for income profile
b1          = 0.16818;             % Linear term for income profile
b2          = -0.0323371/10;       % Quadratic term for income profile
b3          = 0.0019704/100;       % Cubic term for income profile
smay        = 0.1;      % Std. dev. of transitory shock (sigma_u)
smav        = 0.1;      % Std. dev. of permanent shock (sigma_z)
ret_fac     = 0.68212;  % Retirement income replacement rate (rho in TeX)

% Asset return parameters
r           = 1.015;    % Risk-free rate (R_f_bar)
mu          = 0.04;     % Equity premium
sigr        = 0.2;      % Std. dev. of risky asset return (sigma_epsilon)
pps_mu =0.02; % PPS premium

% Tax and Pension parameters
tau_y       = 0.06;     % Labor income tax rate
tau_q       = 0.03;     % Pension withdrawal tax rate
max_Q_hat   = 0.2;      % Max pension contribution (normalized by permanent income)

% Numerical approximation parameters
n           = 5;        % Number of nodes for quadrature
na          = 21;       % Grid size for risky share (alpha)
nc          = 21;       % Grid size for consumption choice
nq          = 21;       % Grid size for pension contribution choice

nw          = 31;       % Grid size for cash-on-hand (W_hat)
max_w       = 200.0;    % Max value for W_hat grid
min_w       = 0.25;     % Min value for W_hat grid

nf          = 21;       % Grid size for pension wealth (F_hat)
max_f       = 20.0;     % Max value for F_hat grid

nsim        = 10000;    % Number of individuals for simulation


%% ----------------------------------------------------------------
%  2. Pre-computation and Grid Setup
%  ----------------------------------------------------------------

% Derived preference parameters
theta = (1.0-rho)/(1.0-1.0/psi);
psi_1 = 1.0-1.0/psi;
psi_2 = 1.0/psi_1;

% Survival probabilities (from Gomes)
survprob = load_survprob(tn);

% Discretize shocks using Gauss-Hermite quadrature
[grid, weig] = gauss_hermite(n);

% Income profile
f_y = zeros(work_life+1,1);
for i1=tb:tr
    f_y(i1-tb+1,1) = exp(aa+b1*i1+b2*i1^2+b3*i1^3);
end

% Deterministic income growth G_t = exp(g(t+1)-g(t))
Gy = zeros(work_life,1);
for i1=1:work_life
    Gy(i1,1) = f_y(i1+1,1)/f_y(i1,1);
end

% Shocks grids
z_shocks = grid * smav;         % Permanent income shocks
u_shocks = grid * smay;         % Transitory income shocks
r_shocks = grid * sigr;         % Risky return shocks

% Returns grid
gret = r + mu + r_shocks; % Gross returns on risky asset

% State and choice grids
gW = logspace(log10(min_w), log10(max_w), nw)'; % Cash-on-hand grid (W_hat)
gF = linspace(0, max_f, nf)';                   % Pension wealth grid (F_hat)
gA = linspace(0, 1, na)';                       % Risky share grid (alpha)
gQ = linspace(0, max_Q_hat, nq)';               % Pension contribution grid (Q_hat)

% Combine weights for 3D integration (z, u, r shocks)
[wZ, wU, wR] = ndgrid(weig, weig, weig);
nweig = wZ .* wU .* wR;


%% ----------------------------------------------------------------
%  3. Value Function Iteration (Backward Induction)
%  ----------------------------------------------------------------
if ~exist(file_path,'file')
% Policy and Value Functions Storage
% During working life (2D state space: W_hat, F_hat)
V_work      = zeros(nw, nf, work_life+1);
C_work_pol  = zeros(nw, nf, work_life+1);
A_work_pol  = zeros(nw, nf, work_life+1);
Q_work_pol  = zeros(nw, nf, work_life+1);

% During retirement (2D state space: W_hat, F_hat_K)
V_ret       = zeros(nw, nf, tn - work_life);
C_ret_pol   = zeros(nw, nf, tn - work_life);
A_ret_pol   = zeros(nw, nf, tn - work_life);

% --- Terminal Period (Age td) ---
t_idx = tn - work_life; % Index for the last period in the V_ret array
num_ret_periods = tn - work_life;
P_annuity_hat = (1 - tau_q) * gF / num_ret_periods; % [nf x 1] vector of annuity flows

% At t=td, all remaining wealth is consumed. V_T = W_T (cash-on-hand).
% Cash-on-hand at this final stage is start-of-period liquid wealth (gW) plus final annuity payment.
% Note: The base retirement income ret_fac is already part of the budget constraint that calculates
% the *next* period's W, so it's implicitly handled. Here we just add the annuity part.
V_ret(:, :, t_idx) = repmat(gW, 1, nf) + repmat(P_annuity_hat', nw, 1);
C_ret_pol(:, :, t_idx) = V_ret(:, :, t_idx);
A_ret_pol(:, :, t_idx) = 0; % No investment in the last period


% --- Retirement Periods (Age td-1 down to tr) ---
fprintf('Solving Retirement Periods (2D VFI)...\n');
n_shocks = length(weig);

for t = (tn-1):-1:(work_life+1)
    t_idx_current = t - work_life;
    V_next = V_ret(:, :, t_idx_current + 1);
    V_next_interp = griddedInterpolant({gW, gF}, V_next, 'linear', 'linear');
    
    % Pre-reshape for broadcasting
    gA_b = reshape(gA, [1, na]);
    gc_b = reshape(linspace(0, 1, nc)', [nc, 1]);
    
    % Loop over the 2D state space (W_hat, F_hat_K)
    parfor i_f = 1:nf
        % Create temporary slices for parfor
        V_slice = zeros(nw,1); 
        C_slice = zeros(nw,1);
        A_slice = zeros(nw,1);
        
        current_F_hat = gF(i_f);
        current_P_annuity_hat = P_annuity_hat(i_f);
        
        for i_w = 1:nw
            % Total cash-on-hand available for consumption and saving
            % This is start-of-period assets (gW) + base retirement income (ret_fac) + annuity
            % The VFI state is gW (assets), but total resources include income.
            cash_on_hand = gW(i_w) + ret_fac + current_P_annuity_hat;
            
            % Consumption grid for this state
            minc = 1e-5;
            maxc = cash_on_hand;
            gc = minc + gc_b * (maxc - minc);
            liquid_sav = cash_on_hand - gc; % [nc x 1] vector of savings

            % --- Vectorized expectation calculation ---
            
            % Portfolio returns, size [1 x na x n_shocks]
            R_portfolio = (1 - gA_b) .* r + gA_b .* reshape(gret, [1, 1, n_shocks]);
            
            % Next period's start-of-period assets (state variable W_hat for V_next)
            W_next = reshape(liquid_sav, [nc, 1, 1]) .* R_portfolio;
            
            % F is fixed during retirement, so F_next is just current_F_hat
            F_next = repmat(current_F_hat, [nc, na, n_shocks]);
            
            % Interpolate to get continuation value
            V_interp = V_next_interp(min(max(W_next, gW(1)), gW(nw)), F_next);
            
            % Integrate over return shocks
            E_V_term = sum(reshape(weig, [1, 1, n_shocks]) .* (V_interp.^(1-rho)), 3); % Result is [nc x na]

            % Current value for all (c, a) choices
            V_candidates = ((1-delta).*(gc.^(psi_1)) + delta.*survprob(t) .* (E_V_term.^((1-1/psi)/(1-rho)))).^psi_2;
            
            % Handle infeasible choices
            V_candidates(liquid_sav < 0, :) = -1e20;

            % Find the optimum
            [best_V, linear_idx] = max(V_candidates(:));
            [idx_c, idx_a] = ind2sub(size(V_candidates), linear_idx);
            
            V_slice(i_w) = best_V;
            C_slice(i_w) = gc(idx_c);
            A_slice(i_w) = gA(idx_a);
        end
        
        % Store results from this parallel worker
        V_ret(:, i_f, t_idx_current) = V_slice;
        C_ret_pol(:, i_f, t_idx_current) = C_slice;
        A_ret_pol(:, i_f, t_idx_current) = A_slice;
    end
    fprintf('Solved for retirement age %d\n', t+tb-1);
end


% --- Store retirement solution for interpolation ---
% This now correctly uses the 2D value function from the first period of retirement
V_work(:,:,work_life+1) = V_ret(:,:,1);
C_work_pol(:,:,work_life+1) = C_ret_pol(:,:,1);
A_work_pol(:,:,work_life+1) = A_ret_pol(:,:,1);
Q_work_pol(:,:,work_life+1) = 0; % No contribution in retirement

% --- Working Periods (Age tb to tr-1) ---
fprintf('Solving Working Periods (Fully Matrixized)...\n');

% Define FIXED PROPORTIONAL grids for choices outside the main loop
gC_prop = linspace(1e-5, 1.0, nc)'; % Proportional consumption grid [nc x 1]
gQ_prop = linspace(0, 1.0, nq);     % Proportional Q contribution grid [1 x nq]


for t = work_life:-1:1
    t1=tic;
    V_next = V_work(:,:,t+1);

    % Create 2D interpolant for next period's value function
    V_next_interp = griddedInterpolant({gW, gF}, V_next, 'linear', 'linear');

    % Pre-reshape choice grids and shock grids for broadcasting
    gA_b = reshape(gA, [1, 1, na]);

    [Z_grid_orig, U_grid_orig, R_grid_orig] = ndgrid(z_shocks, u_shocks, gret);
    [wZ_orig, wU_orig, wR_orig] = ndgrid(weig, weig, weig);
    nweig_mat = reshape(wZ_orig .* wU_orig .* wR_orig, [1, 1, 1, n, n, n]);
    Z_grid = reshape(Z_grid_orig, [1, 1, 1, n, n, n]);
    U_grid = reshape(U_grid_orig, [1, 1, 1, n, n, n]);
    R_grid = reshape(R_grid_orig, [1, 1, 1, n, n, n]);

    % Use parfor for the outer loop over F_hat grid
    parfor i_f = 1:nf

        % Create temporary slices for parfor
        V_slice = zeros(nw,1); C_slice = zeros(nw,1);
        A_slice = zeros(nw,1); Q_slice = zeros(nw,1);

        current_F = gF(i_f);

        for i_w = 1:nw
            current_W = gW(i_w);
            
            % 1. Calculate the maximum possible Q for the current state W
            max_q_for_this_w = min(max_Q_hat, current_W / (1 - tau_y));
            
            % 2. Generate actual Q choices based on the FIXED proportional grid
            gQ_choices = gQ_prop * max_q_for_this_w;
            
            % 3. For each Q choice, determine the max possible consumption
            max_c_after_q = current_W - (1 - tau_y) * gQ_choices;
            
            % --- MODIFICATION START: Smoother Consumption Grid Generation ---
            % 4. Generate the full [nc x nq] matrix of consumption choices.
            % Use max() to ensure the range is always non-negative, preventing
            % discontinuities in the choice grid construction.
            gc = gC_prop * max(max_c_after_q, 1e-5);
            % --- MODIFICATION END ---

            % 5. Calculate savings for all (c, q) combinations
            sav = current_W - (1 - tau_y) * gQ_choices - gc;

            % 6. Reshape for broadcasting
            sav_b = reshape(sav, [nc, nq, 1, 1, 1, 1]);

            % 7. Calculate next-period states
            R_portfolio = (1 - gA_b) * r + gA_b .* R_grid;
            G_t_z = Gy(t) * exp(Z_grid);
            W_next = sav_b .* R_portfolio ./ G_t_z + (1-tau_y) * exp(U_grid);
            
            gQ_b_reshaped = reshape(gQ_choices, [1, nq]);
            F_next_shocks = (current_F + gQ_b_reshaped) * (r + pps_mu)./ G_t_z;

            % 8. Interpolate
            F_next = repmat(reshape(F_next_shocks, [1, nq, 1, n, n, n]), [nc, 1, na, 1, 1, 1]);
            V_interp = V_next_interp(min(max(W_next, gW(1)), gW(nw)), min(max(F_next, gF(1)), gF(nf)));
            V_interp(isnan(V_interp)) = 1e-6;

            % 9. Calculate Expectation
            term = (G_t_z .* V_interp).^(1 - rho);
            E_V_term = squeeze(sum(nweig_mat .* term, [4, 5, 6]));

            % 10. Calculate value for all choices
            gc_b = reshape(gc, [nc, nq, 1]);
            V_candidates = ((1-delta).*(gc_b.^(psi_1)) + delta.*survprob(t).*(E_V_term.^((1-1/psi)/(1-rho)))).^psi_2;

            % Infeasibility is handled cleanly by the savings calculation
            is_infeasible = (sav < 0);
            V_candidates(repmat(is_infeasible, [1, 1, na])) = -1e20;

            % 11. Find the optimum
            [best_V, linear_idx] = max(V_candidates(:));
            [idx_c, idx_q, idx_a] = ind2sub(size(V_candidates), linear_idx);

            % Store results in slices
            V_slice(i_w) = best_V;
            C_slice(i_w) = gc(idx_c, idx_q);
            A_slice(i_w) = gA(idx_a);
            Q_slice(i_w) = gQ_choices(idx_q); 
        end

        % Store results from this parallel worker
        V_work(:, i_f, t) = V_slice;
        C_work_pol(:, i_f, t) = C_slice;
        A_work_pol(:, i_f, t) = A_slice;
        Q_work_pol(:, i_f, t) = Q_slice;
    end
        fprintf('Solved for age %d, %f s\n', t+tb-1, toc(t1));
end
%% ----------------------------------------------------------------
%  3.5. Save VFI Results
%  ----------------------------------------------------------------
fprintf('VFI complete. Saving results...\n');

if ~exist(output_dir, 'dir'), mkdir(output_dir); end
save(file_path, ...
    'V_work', 'C_work_pol', 'A_work_pol', 'Q_work_pol', ...
    'V_ret', 'C_ret_pol', 'A_ret_pol', ...
    'gW', 'gF', 'work_life', 'tn', '-v7.3');
fprintf('Results saved to results/vfi_results.mat\n');

else
vfi_results_file = fullfile('results', 'vfi_results.mat');
if exist(vfi_results_file, 'file')
    fprintf('Loading pre-computed VFI results...\n');
    load(vfi_results_file);
else
    error('VFI results file not found. Please run the VFI section first.');
end
end

%% ----------------------------------------------------------------
%  4. Life-Cycle Simulation
%  ----------------------------------------------------------------

% --- Load VFI results if they exist ---


fprintf('Starting simulation...\n');

% Create interpolants for policy functions from loaded data
C_pol_interp = griddedInterpolant({gW, gF, (1:work_life)'}, C_work_pol(:,:,1:work_life), 'linear', 'linear');
A_pol_interp = griddedInterpolant({gW, gF, (1:work_life)'}, A_work_pol(:,:,1:work_life), 'linear', 'linear');
Q_pol_interp = griddedInterpolant({gW, gF, (1:work_life)'}, Q_work_pol(:,:,1:work_life), 'linear', 'linear');

% For retirement, need a 2D interpolant for C and A
C_ret_pol_interp = griddedInterpolant({gW, gF, (work_life+1:tn)'}, C_ret_pol, 'linear', 'linear');
A_ret_pol_interp = griddedInterpolant({gW, gF, (work_life+1:tn)'}, A_ret_pol, 'linear', 'linear');

% Simulation storage
simAssets = zeros(tn + 1, nsim); % Start-of-period liquid wealth (un-normalized)
simF      = zeros(tn + 1, nsim); % Start-of-period pension wealth (un-normalized)
simM      = zeros(tn, nsim);     % Permanent income component (log)
simY      = zeros(tn, nsim);     % Labor income (un-normalized)
simC      = zeros(tn, nsim);     % Consumption (un-normalized)
simA      = zeros(tn, nsim);     % Risky share
simQ      = zeros(tn, nsim);     % Pension contribution

% Initial conditions (agents start with zero wealth)
simM(1,:) = 0; 

% Generate all shocks at once
z_shocks_sim = randn(tn, nsim) * smav;
u_shocks_sim = randn(tn, nsim) * smay;
r_shocks_sim = randn(tn, nsim) * sigr;

% Loop over time for simulation
for t = 1:tn
    fprintf('Simulating for age %d\n', t+tb-1);
    
    % --- Working Period ---
    if t <= work_life
        % 1. Update permanent income component M_t
        if t > 1
            simM(t,:) = simM(t-1,:) + z_shocks_sim(t,:);
        else % t=1, M_0 = 0, so M_1 = z_1
            simM(t,:) = z_shocks_sim(t,:);
        end
        
        % 2. Calculate normalization factor and income
        P_t = exp(log(f_y(t)) + simM(t,:));
        simY(t,:) = P_t .* exp(u_shocks_sim(t,:));
        
        % 3. Construct the state variables for interpolation
        % The VFI state gW represents start-of-period assets. For policy lookup,
        % we use total cash-on-hand, which is assets + after-tax income.
        % However, the VFI grid gW represents total cash-on-hand. Let's align.
        % The working period VFI uses `current_W` from gW, which is cash-on-hand.
        cash_on_hand_unnorm = simAssets(t,:) + (1-tau_y) * simY(t,:);

        % Normalize the state variables
        W_hat = cash_on_hand_unnorm ./ P_t;
        F_hat = simF(t,:) ./ P_t;
        
        % 4. Interpolate to get optimal policies (normalized)
        W_hat_clipped = min(max(W_hat, gW(1)), gW(nw));
        F_hat_clipped = min(max(F_hat, gF(1)), gF(nf));
        
        C_hat = C_pol_interp(W_hat_clipped, F_hat_clipped, t*ones(1,nsim));
        A_hat = A_pol_interp(W_hat_clipped, F_hat_clipped, t*ones(1,nsim));
        Q_hat = Q_pol_interp(W_hat_clipped, F_hat_clipped, t*ones(1,nsim));
        
        % 5. De-normalize choices
        simC(t,:) = C_hat .* P_t;
        simQ(t,:) = Q_hat .* P_t;
        simA(t,:) = A_hat;
        
        % 6. Apply budget constraints (in un-normalized terms)
        % Ensure Q is feasible
        max_q_unnorm = max_Q_hat * P_t;
        simQ(t,:) = min(simQ(t,:), (1-tau_y)*simY(t,:)); % Cannot contribute more than after-tax income
        simQ(t,:) = min(simQ(t,:), max_q_unnorm);
        simQ(t,:) = max(simQ(t,:), 0);
        
        % Total resources available for C and liquid savings, after Q is chosen
        resources_for_C_and_Sav = cash_on_hand_unnorm - (1-tau_y)*simQ(t,:);
        
        % Ensure C is feasible
        simC(t,:) = min(simC(t,:), resources_for_C_and_Sav - 1e-6);
        simC(t,:) = max(simC(t,:), 1e-6);
        
        liquid_sav = resources_for_C_and_Sav - simC(t,:);
        
        % 7. Evolve wealth to the start of the next period
        if t < tn
            R_portfolio = (1-simA(t,:)).*r + simA(t,:).*(r+mu+r_shocks_sim(t+1,:));
            simAssets(t+1,:) = liquid_sav .* R_portfolio;
            simF(t+1,:) = (simF(t,:) + simQ(t,:)) * (r + pps_mu);
        end
        
    % --- Retirement Period ---
    else
        P_annuity = 0;
        % At the start of retirement, calculate the annuity payment
        if t == work_life + 1
            num_ret_periods = tn - work_life;
            simF_at_retirement = simF(t,:);
            P_annuity_flow = (1-tau_q) * simF_at_retirement / num_ret_periods;
            simM(t,:) = simM(work_life,:); % Permanent income stops evolving
        else
            simM(t,:) = simM(t-1,:);
        end

        % Permanent income level at retirement is the normalization factor
        P_retire = exp(log(f_y(work_life+1)) + simM(work_life,:));
        simY(t,:) = ret_fac * P_retire;
        
        % Total resources (cash-on-hand) for retirees
        total_resources_unnorm = simAssets(t,:) + simY(t,:) + P_annuity_flow;
        
        % Normalize total resources to get the state variable for interpolation
        % The state for retirement VFI is (start-of-period assets, F_at_retirement)
        W_hat_ret = simAssets(t,:) ./ P_retire;
        F_hat_ret = simF(work_life+1,:) ./ P_retire; % F is fixed at the retirement value
        
        W_hat_clipped = min(max(W_hat_ret, gW(1)), gW(nw)); 
        F_hat_clipped = min(max(F_hat_ret, gF(1)), gF(nf));
        
        % Get policies (normalized consumption)
        C_hat = C_ret_pol_interp(W_hat_clipped, F_hat_clipped, t*ones(1,nsim));
        A_hat = A_ret_pol_interp(W_hat_clipped, F_hat_clipped, t*ones(1,nsim));
        
        % De-normalize and apply policies
        simC(t,:) = C_hat .* P_retire;
        simA(t,:) = A_hat;
        simQ(t,:) = 0;
        
        % Apply budget constraints
        simC(t,:) = min(simC(t,:), total_resources_unnorm - 1e-6);
        simC(t,:) = max(simC(t,:), 1e-6);
        liquid_sav = total_resources_unnorm - simC(t,:);
        
        % Evolve wealth to the start of the next period
        if t < tn
            R_portfolio = (1-simA(t,:)).*r + simA(t,:).*(r+mu+r_shocks_sim(t+1,:));
            simAssets(t+1,:) = liquid_sav .* R_portfolio;
            simF(t+1,:) = simF(t,:); % Pension wealth is depleted via annuity, stock for policy lookup remains
        end
    end
end

%% ----------------------------------------------------------------
%  5. Plotting Results
%  ----------------------------------------------------------------

% Use wealth at the beginning of the period for plotting
simAssets_plot = simAssets(1:tn,:);
simF_plot = simF(1:tn,:);

% Calculate mean profiles
meanAssets = mean(simAssets_plot, 2);
meanF = mean(simF_plot, 2);
meanTotalW = meanAssets + meanF;
meanC = mean(simC, 2);
meanY = mean(simY, 2);
meanA = mean(simA, 2);
meanQ = mean(simQ, 2);

% Avoid division by zero for ages with zero income
meanY_for_ratio = meanY;
meanY_for_ratio(meanY_for_ratio == 0) = NaN;
meanW_Y = meanTotalW ./ meanY_for_ratio;

ages = (tb:td)';

fig = figure('Name', 'Life-Cycle Profiles with PPS', 'Position', [100, 100, 1200, 800]);

% Panel 1: Wealth, Consumption, Income
subplot(2,2,1);
plot(ages, meanTotalW, 'b-', 'LineWidth', 2, 'DisplayName', 'Total Wealth ($\bar{W}+\bar{F}$)');
hold on;
plot(ages, meanC, 'r--', 'LineWidth', 2, 'DisplayName', 'Consumption ($\bar{C}$)');
plot(ages, meanY, 'k:', 'LineWidth', 2, 'DisplayName', 'Income ($\bar{Y}$)');
xline(tr, 'k:', 'Retirement', 'LineWidth', 1.5, 'HandleVisibility', 'off');
hold off;
title('Mean Lifecycle Profiles'); xlabel('Age'); ylabel('Value');
grid on; xlim([tb, td]); legend('Location', 'northwest');

% Panel 2: Risky Asset Allocation
subplot(2,2,2);
plot(ages, meanA, 'm-', 'LineWidth', 2);
xline(tr, 'k:', 'Retirement', 'LineWidth', 1.5);
title('Mean Risky Asset Allocation'); xlabel('Age'); ylabel('Risky Share ($\bar{\alpha}$)');
ylim([-0.05, 1.05]); grid on; xlim([tb, td]);

% Panel 3: Pension Wealth and Contribution
subplot(2,2,3);
yyaxis left
plot(ages, meanF, 'g-', 'LineWidth', 2, 'DisplayName', 'Pension Wealth ($\bar{F}$)');
ylabel('Pension Wealth');
ylim_left = get(gca, 'YLim');
if ylim_left(2) == 0, set(gca, 'YLim', [0, 1]); end % Handle case of zero pension wealth
hold on;
yyaxis right
plot(ages(1:work_life), meanQ(1:work_life), 'b--', 'LineWidth', 2, 'DisplayName', 'Pension Contr. ($\bar{Q}$)');
ylabel('Pension Contribution');
ylim_right = get(gca, 'YLim');
if ylim_right(2) == 0, set(gca, 'YLim', [0, 0.1]); end % Handle case of zero contribution
xline(tr, 'k:', 'Retirement', 'LineWidth', 1.5, 'HandleVisibility', 'off');
hold off;
title('Pension Account Dynamics'); xlabel('Age');
grid on; xlim([tb, td]); legend('Location', 'northwest');

% Panel 4: Wealth-to-Income Ratio
subplot(2,2,4);
plot(ages, meanW_Y, 'c-', 'LineWidth', 2);
xline(tr, 'k:', 'Retirement', 'LineWidth', 1.5);
title('Mean Wealth-to-Income Ratio'); xlabel('Age'); ylabel('$(\bar{W}+\bar{F}) / \bar{Y}$');
grid on; xlim([tb, td]);

sgtitle('Life-Cycle Model with Private Pension System', 'FontSize', 16, 'FontWeight', 'bold');

% Save figure
output_dir = 'fig';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
output_filename = fullfile(output_dir, 'life_cycle_pps.png');
print(fig, output_filename, '-dpng', '-r300');


%% Helper Functions

function survprob = load_survprob(tn)
% Conditional Survival Probabilities from Gomes (2020)
survprob_data = [0.99845, 0.99839, 0.99833, 0.9983, 0.99827, 0.99826, 0.99824, 0.9982, 0.99813, 0.99804, 0.99795, 0.99785, 0.99776, 0.99766, 0.99755, 0.99743, 0.9973, 0.99718, 0.99707, 0.99696, 0.99685, 0.99672, 0.99656, 0.99635, 0.9961, 0.99579, 0.99543, 0.99504, 0.99463, 0.9942, 0.9937, 0.99311, 0.99245, 0.99172, 0.99091, 0.99005, 0.98911, 0.98803, 0.9868, 0.98545, 0.98409, 0.9827, 0.98123, 0.97961, 0.97786, 0.97603, 0.97414, 0.97207, 0.9697, 0.96699, 0.96393, 0.96055, 0.9569, 0.9531, 0.94921, 0.94508, 0.94057, 0.9357, 0.93031, 0.92424, 0.91717, 0.90922, 0.90089, 0.89282, 0.88503, 0.87622, 0.86576, 0.8544, 0.8423, 0.82942, 0.8154, 0.80002, 0.78404, 0.76842, 0.75382, 0.73996, 0.72464, 0.71057, 0.6961, 0.6809];
survprob = ones(tn-1, 1);
len = min(tn-1, length(survprob_data));
survprob(1:len) = survprob_data(1:len);
end

function [x, w] = gauss_hermite(n)
% Computes nodes and weights for Gauss-Hermite quadrature
i = 1:n-1;
a = sqrt(i/2);
CM = diag(a,1) + diag(a,-1);
[V, L] = eig(CM);
[x, ind] = sort(diag(L));
V = V(:,ind)';
w = sqrt(pi) * V(:,1).^2;
end