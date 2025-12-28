%% Test Script for Pension Contribution (Q) Policy Function
%  This script isolates the one-period decision problem for a typical
%  working-age individual to analyze the sensitivity of the pension
%  contribution policy (Q_hat) to preference parameters (rho, psi).

clear all;
close all;
clc;

fprintf('Starting policy function sensitivity analysis...\n');

%% ----------------------------------------------------------------
%  1. Setup: Model Environment (Consistent with life_cycle_pps.m)
%  ----------------------------------------------------------------

% Life-cycle and parameter setup
tb = 20; tr = 66;
delta = 0.97;
tau_y = 0.06;
max_Q_hat = 0.2;
r = 1.015; mu = 0.04; sigr = 0.2; pps_mu = 0.02;
smav = 0.1; smay = 0.1;
aa = -2.170042+2.700381; b1 = 0.16818; b2 = -0.0323371/10; b3 = 0.0019704/100;

% Numerical approximation grids
n = 5; na = 21; nc = 21; nq = 21;
nw = 31; max_w = 200.0; min_w = 0.25;
nf = 21; max_f = 20.0;

% Define the age for testing
test_age = 45;
t = test_age - tb + 1; % Model period

% Define parameter scenarios to test
scenarios = {
    struct('name', 'Baseline',            'rho', 10.0, 'psi', 0.5), ...
    struct('name', 'Low Risk Aversion',   'rho', 3.0,  'psi', 0.5), ...
    struct('name', 'High EIS',            'rho', 10.0, 'psi', 1.5), ...
    struct('name', 'CRRA Utility (rho=3)',  'rho', 3.0,  'psi', 1/3)
};

%% ----------------------------------------------------------------
%  2. Pre-computation and Grid Setup
%  ----------------------------------------------------------------

% Survival probabilities
survprob_t = 0.99579; % Approx. for age 45 from Gomes data

% Discretize shocks using Gauss-Hermite quadrature
[grid, weig] = gauss_hermite(n);

% Income profile and growth
f_y_t = exp(aa+b1*(t+tb-1)+b2*(t+tb-1)^2+b3*(t+tb-1)^3);
f_y_t_plus_1 = exp(aa+b1*(t+tb)+b2*(t+tb)^2+b3*(t+tb)^3);
Gy_t = f_y_t_plus_1 / f_y_t;

% Shocks grids
z_shocks = grid * smav;
u_shocks = grid * smay;
r_shocks = grid * sigr;
gret = r + mu + r_shocks;

% State and choice grids
gW = logspace(log10(min_w), log10(max_w), nw)';
gF = linspace(0, max_f, nf)';
gA = linspace(0, 1, na)';

% Combine weights for 3D integration
[wZ_orig, wU_orig, wR_orig] = ndgrid(weig, weig, weig);
nweig_mat = reshape(wZ_orig .* wU_orig .* wR_orig, [1, 1, 1, n, n, n]);

[Z_grid_orig, U_grid_orig, R_grid_orig] = ndgrid(z_shocks, u_shocks, gret);
Z_grid = reshape(Z_grid_orig, [1, 1, 1, n, n, n]);
U_grid = reshape(U_grid_orig, [1, 1, 1, n, n, n]);
R_grid = reshape(R_grid_orig, [1, 1, 1, n, n, n]);

G_t_z = Gy_t * exp(Z_grid);

%% ----------------------------------------------------------------
%  3. Main Loop: Solve One-Period Problem for Each Scenario
%  ----------------------------------------------------------------

Q_policy_results = cell(1, length(scenarios));

for s = 1:length(scenarios)
    % --- Get current scenario parameters ---
    rho = scenarios{s}.rho;
    psi = scenarios{s}.psi;
    fprintf('Testing Scenario: %s (rho=%.2f, psi=%.2f)\n', scenarios{s}.name, rho, psi);
    
    % Derived preference parameters
    theta = (1.0-rho)/(1.0-1.0/psi);
    psi_1 = 1.0-1.0/psi;
    psi_2 = 1.0/psi_1;
    
    % --- Simplification for V_next ---
    % To test the mechanism without solving the entire model backwards, we
    % assume a simple, plausible form for the next-period value function.
    % Here, we assume V_next is a linear function of total wealth, meaning
    % future liquid and illiquid wealth are perfectly substitutable. This
    % is a strong assumption that FAVORS saving in the illiquid account.
    V_next = repmat(gW, 1, nf) + repmat(gF', nw, 1);
    V_next_interp = griddedInterpolant({gW, gF}, V_next, 'linear', 'linear');
    
    % Storage for this scenario's policy function
    Q_policy = zeros(nw, nf);
    
    % --- Solve for optimal Q at each state (W_hat, F_hat) ---
    for i_f = 1:nf
        current_F = gF(i_f);
        
        for i_w = 1:nw
            current_W = gW(i_w);
            
            % Dynamic Q grid
            max_q_for_this_w = min(max_Q_hat, current_W / (1 - tau_y));
            gQ_dynamic = linspace(0, max_q_for_this_w, nq)';
            gQ_b = reshape(gQ_dynamic, [1, nq]);
            
            % Consumption grid
            max_c = current_W - (1 - tau_y) * gQ_b;
            unit_grid = linspace(0, 1, nc)';
            range_c = max_c - 1e-5;
            gc = 1e-5 + unit_grid * range_c;
            gc(:, max_c <= 1e-5) = 1e-5;
            
            % Savings
            sav = current_W - (1 - tau_y) * gQ_b - gc;
            
            % Reshape for broadcasting
            sav_b = reshape(sav, [nc, nq, 1, 1, 1, 1]);
            gA_b = reshape(gA, [1, 1, na]);
            
            % --- STATE TRANSITION EQUATION (Using the THEORETICALLY CORRECT form) ---
            R_portfolio = (1 - gA_b) * r + gA_b .* R_grid;
            W_next = sav_b .* R_portfolio ./ G_t_z + (1-tau_y) * exp(U_grid);
            F_next_shocks = (current_F + gQ_b) * (r + pps_mu)./ G_t_z;
            
            % Interpolate
            F_next = repmat(reshape(F_next_shocks, [1, nq, 1, n, n, n]), [nc, 1, na, 1, 1, 1]);
            V_interp = V_next_interp(min(max(W_next, gW(1)), gW(nw)), min(max(F_next, gF(1)), gF(nf)));
            V_interp(isnan(V_interp)) = 1e-6;
            
            % Calculate Expectation
            term = (G_t_z .* V_interp).^(1 - rho);
            E_V_term = squeeze(sum(nweig_mat .* term, [4, 5, 6]));
            
            % Calculate value for all choices
            gc_b = reshape(gc, [nc, nq, 1]);
            V_candidates = ((1-delta).*(gc_b.^(psi_1)) + delta.*survprob_t.*(E_V_term.^((1-1/psi)/(1-rho)))).^psi_2;
            
            % Handle infeasible choices
            is_infeasible = (sav < 0) | reshape(max_c <= 1e-5, [1, nq]);
            V_candidates(repmat(is_infeasible, [nc, 1, na])) = -1e20;
            
            % Find the optimum
            [~, linear_idx] = max(V_candidates(:));
            [~, idx_q, ~] = ind2sub(size(V_candidates), linear_idx);
            
            % Store result
            Q_policy(i_w, i_f) = gQ_dynamic(idx_q);
        end
    end
    Q_policy_results{s} = Q_policy;
end
fprintf('Analysis complete. Plotting results...\n');

%% ----------------------------------------------------------------
%  4. Plotting Results
%  ----------------------------------------------------------------
fig = figure('Name', 'Pension Contribution Policy Sensitivity', 'Position', [50, 50, 1400, 1000]);
sgtitle(sprintf('Pension Contribution Policy ($\\hat{Q}$) at Age %d', test_age), 'Interpreter', 'latex', 'FontSize', 16);

for s = 1:length(scenarios)
    subplot(2, 2, s);
    imagesc(gF, gW, Q_policy_results{s});
    
    set(gca, 'YDir', 'normal'); % Set origin to bottom-left
    colorbar;
    xlabel('Pension Wealth ($\hat{F}$)', 'Interpreter', 'latex');
    ylabel('Cash-on-Hand ($\hat{W}$)', 'Interpreter', 'latex');
    title(sprintf('%s ($\\rho=%.2f, \\psi=%.2f$)', scenarios{s}.name, scenarios{s}.rho, scenarios{s}.psi), 'Interpreter', 'latex');
    grid on;
end

output_dir = 'fig';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
output_filename = fullfile(output_dir, 'test_Q_policy_sensitivity.png');
print(fig, output_filename, '-dpng', '-r300');


%% Helper Function
function [x, w] = gauss_hermite(n)
    i = 1:n-1;
    a = sqrt(i/2);
    CM = diag(a,1) + diag(a,-1);
    [V, L] = eig(CM);
    [x, ind] = sort(diag(L));
    V = V(:,ind)';
    w = sqrt(pi) * V(:,1).^2;
end