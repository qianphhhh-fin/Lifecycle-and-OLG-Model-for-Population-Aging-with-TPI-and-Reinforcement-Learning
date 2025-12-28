%% Test script for vectorizing the VFI shock integration
% This script compares three versions of the VFI calculation for a single period:
% 1. Loopy Version: Nested loops for everything.
% 2. Matrixized Shocks: Inner loops over shocks (z, u, r) are vectorized.
% 3. Fully Matrixized: Loops over choices (q, a) are also vectorized.
%
% GOAL: Verify correctness and measure progressive speedup.

clear all;
clc;
close all;

%% ----------------------------------------------------------------
%  1. Setup Minimal Environment
%  ----------------------------------------------------------------
fprintf('1. Setting up minimal environment...\n');

% Minimal grid sizes for quick testing
nw = 31;     % Grid size for cash-on-hand (W_hat)
nf = 7;     % Grid size for pension wealth (F_hat)
na = 21;     % Grid size for risky share (alpha)
nc = 21;    % Grid size for consumption choice
nq = 21;     % Grid size for pension contribution choice
n  = 5;     % Number of quadrature nodes for EACH shock (Increased to show speedup)

% Model Parameters (placeholders for calculation)
delta = 0.97;
rho   = 10.0;
psi   = 0.5;
r     = 1.015;
mu    = 0.04;
sigr  = 0.2;
smav  = 0.1;
smay  = 0.1;
tau_y = 0.06;
Gy_t  = 1.02; % Assume a constant deterministic income growth
survprob_t = 0.99;

% Derived preference parameters
theta = (1.0 - rho) / (1.0 - 1.0 / psi);
psi_1 = 1.0 - 1.0 / psi;
psi_2 = 1.0 / psi_1;

% State and Choice Grids
gW = linspace(0.5, 10, nw)';
gF = linspace(0, 20, nf)';
gA = linspace(0, 1, na)';
gQ = linspace(0, 0.2, nq)';

% Shocks
[grid_nodes, weig] = gauss_hermite(n); % Use proper quadrature

z_shocks = grid_nodes * smav;
u_shocks = grid_nodes * smay;
r_shocks = grid_nodes * sigr;
gret = r + mu + r_shocks;

% Create a dummy V_next and its interpolant
[W_grid_3d, F_grid_3d] = ndgrid(gW, gF);
V_next = log(W_grid_3d + F_grid_3d + 1);
V_next_interp = griddedInterpolant({gW, gF}, V_next, 'linear', 'none');

% Storage for results
V_loop = zeros(nw, nf); C_loop = zeros(nw, nf); A_loop = zeros(nw, nf); Q_loop = zeros(nw, nf);
V_mat1 = zeros(nw, nf); C_mat1 = zeros(nw, nf); A_mat1 = zeros(nw, nf); Q_mat1 = zeros(nw, nf);
V_mat2 = zeros(nw, nf); C_mat2 = zeros(nw, nf); A_mat2 = zeros(nw, nf); Q_mat2 = zeros(nw, nf);

%% ----------------------------------------------------------------
%  2. Original Loopy Version (Ground Truth)
%  ----------------------------------------------------------------
fprintf('\n2. Running Original Loopy Version...\n');
tic;

[wZ, wU, wR] = ndgrid(weig, weig, weig);
nweig1 = wZ .* wU .* wR;

for i_f = 1:nf
    for i_w = 1:nw
        current_W = gW(i_w); current_F = gF(i_f);
        best_V = -1e20; best_c = 0; best_a = 0; best_q = 0;

        for i_q = 1:nq
            q_choice = gQ(i_q);
            max_c = current_W - (1 - tau_y) * q_choice;
            if max_c <= 1e-5, continue; end
            gc = linspace(1e-5, max_c, nc)';

            for i_a = 1:na
                alpha = gA(i_a);
                sav = current_W - (1 - tau_y) * q_choice - gc;

                E_V_term = zeros(nc, 1);
                for i_z = 1:n
                    G_t_z = Gy_t * exp(z_shocks(i_z));
                    for i_u = 1:n
                        for i_r = 1:n
                            R_portfolio = (1 - alpha) * r + alpha * gret(i_r);
                            W_next_vec = sav * R_portfolio / G_t_z + (1 - tau_y) * exp(u_shocks(i_u));
                            F_next_scalar = (current_F + q_choice) * r / G_t_z;
                            F_next_vec = repmat(F_next_scalar, size(W_next_vec));

                            V_interp = V_next_interp(min(max(W_next_vec, gW(1)), gW(nw)), min(max(F_next_vec, gF(1)), gF(nf)));
                            V_interp(isnan(V_interp)) = 1e-6;

                            term = (G_t_z * V_interp).^(1 - rho);
                            E_V_term = E_V_term + nweig1(i_z, i_u, i_r) * term;
                        end
                    end
                end

                V_candidates = ((1 - delta) .* (gc.^(psi_1)) + delta .* survprob_t .* (E_V_term.^((1 - 1/psi) / (1 - rho)))).^psi_2;

                [max_val, max_idx] = max(V_candidates);
                if max_val > best_V
                    best_V = max_val; best_c = gc(max_idx); best_a = alpha; best_q = q_choice;
                end
            end
        end
        V_loop(i_w, i_f) = best_V; C_loop(i_w, i_f) = best_c; A_loop(i_w, i_f) = best_a; Q_loop(i_w, i_f) = best_q;
    end
end
time_loop = toc;
fprintf('   Done. Elapsed time: %.4f seconds.\n', time_loop);

%% ----------------------------------------------------------------
%  3. Matrixized Shocks Version
%  ----------------------------------------------------------------
fprintf('\n3. Running Matrixized Shocks Version...\n');
tic;

[Z_grid_orig, U_grid_orig, R_grid_orig] = ndgrid(z_shocks, u_shocks, gret);
[wZ_orig, wU_orig, wR_orig] = ndgrid(weig, weig, weig);
Z_grid = reshape(Z_grid_orig, [1, n, n, n]); U_grid = reshape(U_grid_orig, [1, n, n, n]);
R_grid = reshape(R_grid_orig, [1, n, n, n]); nweig_mat = reshape(wZ_orig .* wU_orig .* wR_orig, [1, n, n, n]);

for i_f = 1:nf
    for i_w = 1:nw
        current_W = gW(i_w); current_F = gF(i_f);
        best_V = -1e20; best_c = 0; best_a = 0; best_q = 0;

        for i_q = 1:nq
            q_choice = gQ(i_q);
            max_c = current_W - (1 - tau_y) * q_choice;
            if max_c <= 1e-5, continue; end
            gc = linspace(1e-5, max_c, nc)';

            for i_a = 1:na
                alpha = gA(i_a);
                sav = current_W - (1 - tau_y) * q_choice - gc;
                sav_b = reshape(sav, [nc, 1, 1, 1]);

                G_t_z = Gy_t * exp(Z_grid);
                R_portfolio = (1 - alpha) * r + alpha * R_grid;
                W_next = sav_b .* R_portfolio ./ G_t_z + (1 - tau_y) * exp(U_grid);
                F_next_shocks = (current_F + q_choice) * r ./ G_t_z;
                F_next = repmat(F_next_shocks, [nc, 1, 1, 1]);

                V_interp = V_next_interp(min(max(W_next, gW(1)), gW(nw)), min(max(F_next, gF(1)), gF(nf)));
                V_interp(isnan(V_interp)) = 1e-6;

                term = (G_t_z .* V_interp).^(1-rho);
                E_V_term_mat = squeeze(sum(nweig_mat .* term, [2, 3, 4]));

                V_candidates = ((1 - delta) .* (gc.^(psi_1)) + delta .* survprob_t .* (E_V_term_mat.^((1 - 1/psi) / (1 - rho)))).^psi_2;

                [max_val, max_idx] = max(V_candidates);
                if max_val > best_V
                    best_V = max_val; best_c = gc(max_idx); best_a = alpha; best_q = q_choice;
                end
            end
        end
        V_mat1(i_w, i_f) = best_V; C_mat1(i_w, i_f) = best_c; A_mat1(i_w, i_f) = best_a; Q_mat1(i_w, i_f) = best_q;
    end
end
time_mat1 = toc;
fprintf('   Done. Elapsed time: %.4f seconds.\n', time_mat1);

%% ----------------------------------------------------------------
%  4. Fully Matrixized Version (Choices and Shocks)
%  ----------------------------------------------------------------
fprintf('\n4. Running Fully Matrixized Version...\n');
tic;

% Reshape choice grids for broadcasting
% gQ_b -> [1 x nq]
% gA_b -> [1 x 1 x na]
gQ_b = reshape(gQ, [1, nq]);
gA_b = reshape(gA, [1, 1, na]);

% Pre-reshape shock grids and weights for broadcasting
% Final size for all will be [1 x 1 x 1 x n x n x n]
[Z_grid_orig, U_grid_orig, R_grid_orig] = ndgrid(z_shocks, u_shocks, gret);
[wZ_orig, wU_orig, wR_orig] = ndgrid(weig, weig, weig);
nweig_mat = reshape(wZ_orig .* wU_orig .* wR_orig, [1, 1, 1, n, n, n]);
Z_grid = reshape(Z_grid_orig, [1, 1, 1, n, n, n]);
U_grid = reshape(U_grid_orig, [1, 1, 1, n, n, n]);
R_grid = reshape(R_grid_orig, [1, 1, 1, n, n, n]);


for i_f = 1:nf
    for i_w = 1:nw
        current_W = gW(i_w);
        current_F = gF(i_f);
        
        % --- Vectorization over Q, A, C starts here ---
        
        % 1. Create consumption grids for all Q choices simultaneously
        % max_c is [1 x nq]
        max_c = current_W - (1 - tau_y) * gQ_b;
        
        % Correctly generate the [nc x nq] consumption grid matrix
        unit_grid = linspace(0, 1, nc)'; % [nc x 1]
        range_c = max_c - 1e-5;          % [1 x nq]
        gc = 1e-5 + unit_grid * range_c; % Broadcasts to [nc x nq]
        
        % Invalidate grids where max_c was non-positive from the start
        gc(:, max_c <= 1e-5) = 1e-5;
        
        % 2. Calculate savings for all (c, q) combinations
        % sav is [nc x nq]
        sav = current_W - (1 - tau_y) * gQ_b - gc;
        
        % 3. Reshape for broadcasting with alpha and shocks
        % sav_b -> [nc x nq x 1 x 1 x 1 x 1]
        sav_b = reshape(sav, [nc, nq, 1, 1, 1, 1]);
        
        % 4. Calculate next-period states for all choices and shocks
        % R_portfolio -> [1 x 1 x na x n x n x n]
        R_portfolio = (1 - gA_b) * r + gA_b .* R_grid;
        G_t_z = Gy_t * exp(Z_grid); % [1 x 1 x 1 x n x n x n]
        
        % W_next broadcasts to [nc x nq x na x n x n x n]
        W_next = sav_b .* R_portfolio ./ G_t_z + (1-tau_y) * exp(U_grid);
        
        % F_next broadcasts to [1 x nq x 1 x n x n x n]
        F_next_shocks = (current_F + gQ_b) * r ./ G_t_z;
        
        % 5. Interpolate
        % To match W_next, F_next needs to be [nc x nq x na x n x n x n]
        F_next = repmat(reshape(F_next_shocks, [1, nq, 1, n, n, n]), [nc, 1, na, 1, 1, 1]);
        V_interp = V_next_interp(min(max(W_next, gW(1)), gW(nw)), min(max(F_next, gF(1)), gF(nf)));
        V_interp(isnan(V_interp)) = 1e-6;
        
        % 6. Calculate Expectation
        % term -> [nc x nq x na x n x n x n]
        term = (G_t_z .* V_interp).^(1 - rho);
        % E_V_term -> [nc x nq x na]
        E_V_term = squeeze(sum(nweig_mat .* term, [4, 5, 6]));
        
        % 7. Calculate value for all choices
        % V_candidates -> [nc x nq x na]
        gc_b = reshape(gc, [nc, nq, 1]);
        V_candidates = ((1-delta).*(gc_b.^(psi_1)) + delta.*survprob_t.*(E_V_term.^((1-1/psi)/(1-rho)))).^psi_2;
        
        % Invalidate choices leading to negative savings or where q was infeasible
        is_infeasible = (sav < 0) | reshape(max_c <= 1e-5, [1, nq]);
        V_candidates(is_infeasible) = -1e20;
        
        % 8. Find the optimum
        [best_V, linear_idx] = max(V_candidates(:));
        [idx_c, idx_q, idx_a] = ind2sub(size(V_candidates), linear_idx);
        
        V_mat2(i_w, i_f) = best_V;
        C_mat2(i_w, i_f) = gc(idx_c, idx_q);
        A_mat2(i_w, i_f) = gA(idx_a);
        Q_mat2(i_w, i_f) = gQ(idx_q);
    end
end
time_mat2 = toc;
fprintf('   Done. Elapsed time: %.4f seconds.\n', time_mat2);
%% ----------------------------------------------------------------
%  5. Super Matrixized Version (Choices, Shocks, and State W)
%  ----------------------------------------------------------------
fprintf('\n5. Running Super Matrixized Version...\n');
tic;

% Pre-allocate storage for this version's results
V_mat3 = zeros(nw, nf); C_mat3 = zeros(nw, nf); A_mat3 = zeros(nw, nf); Q_mat3 = zeros(nw, nf);

% Define the logical order of dimensions for explicit broadcasting:
% 1:W, 2:C, 3:Q, 4:A, 5:Z, 6:U, 7:R

% Reshape choice and state grids into 7D for explicit broadcasting
gW_b = reshape(gW, [nw, 1, 1, 1, 1, 1, 1]);
gQ_b = reshape(gQ, [1, 1, nq, 1, 1, 1, 1]);
gA_b = reshape(gA, [1, 1, 1, na, 1, 1, 1]);

% Reshape shock grids and weights into 7D
[wZ_orig, wU_orig, wR_orig] = ndgrid(weig, weig, weig);
nweig_mat = reshape(wZ_orig .* wU_orig .* wR_orig, [1, 1, 1, 1, n, n, n]);

Z_grid = reshape(z_shocks, [1, 1, 1, 1, n, 1, 1]);
U_grid = reshape(u_shocks, [1, 1, 1, 1, 1, n, 1]);
R_shocks_grid = reshape(r_shocks, [1, 1, 1, 1, 1, 1, n]);

% This version only loops over the F_hat grid
for i_f = 1:nf
    current_F = gF(i_f);
    
    % --- Vectorization over W, Q, A, C starts here ---
    
    % 1. Create consumption grids for all (W, Q) combinations
    % max_c becomes a [nw x 1 x nq] tensor
    max_c = gW_b - (1 - tau_y) * gQ_b;
    
    % Generate a [nw x nc x nq] consumption grid tensor 'gc'
    unit_grid = reshape(linspace(0, 1, nc), [1, nc, 1]);
    range_c = max_c - 1e-5;
    gc = 1e-5 + unit_grid .* range_c;
    gc(gc <= 1e-5) = 1e-5;
    
    % 2. Calculate savings for all (W, C, Q) combinations
    % sav becomes a [nw x nc x nq] tensor
    sav = max_c - gc;
    
    % 3. Reshape sav into the 7D space
    sav_b = reshape(sav, [nw, nc, nq, 1, 1, 1, 1]);
    
    % 4. Calculate next-period states for all choices and shocks
    R_portfolio = (1 - gA_b) .* r + gA_b .* (r + mu + R_shocks_grid);
    G_t_z = Gy_t * exp(Z_grid);
    
    % All components will now broadcast correctly to the final size [nw, nc, nq, na, n, n, n]
    W_next = sav_b .* R_portfolio ./ G_t_z + (1-tau_y) * exp(U_grid);
    
    % F_next also calculated in the 7D space
    F_next_shocks = (current_F + gQ_b) * r ./ G_t_z;
    
    % 5. Interpolate
    % Replicate F_next to match W_next's dimensions for interpolation
    F_next = repmat(F_next_shocks, [nw, nc, 1, na, 1, n, n]);
    V_interp = V_next_interp(min(max(W_next, gW(1)), gW(nw)), min(max(F_next, gF(1)), gF(nf)));
    V_interp(isnan(V_interp)) = 1e-6;
    
    % 6. Calculate Expectation
    % term -> [nw x nc x nq x na x n x n x n]
    term = (G_t_z .* V_interp).^(1 - rho);
    % E_V_term -> [nw x nc x nq x na]
    E_V_term = squeeze(sum(nweig_mat .* term, [5, 6, 7]));
    
    % 7. Calculate value for all choices
    % V_candidates -> [nw x nc x nq x na]
    V_candidates = ((1-delta).*(reshape(gc, [nw, nc, nq, 1]).^(psi_1)) + delta.*survprob_t.*(E_V_term.^((1-1/psi)/(1-rho)))).^psi_2;
    
    % Invalidate choices leading to negative savings or where q was infeasible
    is_infeasible = (sav < 0) | (max_c <= 1e-5);
    V_candidates(is_infeasible) = -1e20;
    
    % 8. Find the optimum for each W state simultaneously
    V_candidates_reshaped = reshape(V_candidates, nw, nc * nq * na);
    [best_V_vec, linear_idx_vec] = max(V_candidates_reshaped, [], 2);
    
    [idx_c_vec, idx_q_vec, idx_a_vec] = ind2sub([nc, nq, na], linear_idx_vec);
    
    % Extract optimal policies using the indices
    % For gc, which is [nw x nc x nq], we need to perform indexed subscription
    linear_gc_indices = sub2ind(size(gc), (1:nw)', idx_c_vec, idx_q_vec);
    
    V_mat3(:, i_f) = best_V_vec;
    C_mat3(:, i_f) = gc(linear_gc_indices);
    A_mat3(:, i_f) = gA(idx_a_vec);
    Q_mat3(:, i_f) = gQ(idx_q_vec);
end
time_mat3 = toc;
fprintf('   Done. Elapsed time: %.4f seconds.\n', time_mat3);


%% ----------------------------------------------------------------
%  6. Comparison
%  ----------------------------------------------------------------
fprintf('\n6. Comparing Results...\n');

fprintf('   Comparison between Loop and Matrix-Shocks:\n');
fprintf('   Max absolute difference in V: %e\n', max(abs(V_loop(:) - V_mat1(:))));
fprintf('   Max absolute difference in C: %e\n', max(abs(C_loop(:) - C_mat1(:))));
fprintf('   Max absolute difference in A: %e\n', max(abs(A_loop(:) - A_mat1(:))));
fprintf('   Max absolute difference in Q: %e\n', max(abs(Q_loop(:) - Q_mat1(:))));

fprintf('\n   Comparison between Loop and Full-Matrix (Choices):\n');
fprintf('   Max absolute difference in V: %e\n', max(abs(V_loop(:) - V_mat2(:))));
fprintf('   Max absolute difference in C: %e\n', max(abs(C_loop(:) - C_mat2(:))));
fprintf('   Max absolute difference in A: %e\n', max(abs(A_loop(:) - A_mat2(:))));
fprintf('   Max absolute difference in Q: %e\n', max(abs(Q_loop(:) - Q_mat2(:))));

fprintf('\n   Comparison between Loop and Super-Matrix (Choices+States):\n');
fprintf('   Max absolute difference in V: %e\n', max(abs(V_loop(:) - V_mat3(:))));
fprintf('   Max absolute difference in C: %e\n', max(abs(C_loop(:) - C_mat3(:))));
fprintf('   Max absolute difference in A: %e\n', max(abs(A_loop(:) - A_mat3(:))));
fprintf('   Max absolute difference in Q: %e\n', max(abs(Q_loop(:) - Q_mat3(:))));

fprintf('\nPerformance Summary:\n');
fprintf('   Loopy version time:          %.4f s\n', time_loop);
fprintf('   Matrix-Shocks version time:  %.4f s\n', time_mat1);
fprintf('   Full-Matrix version time:    %.4f s\n', time_mat2);
fprintf('   Super-Matrix version time:   %.4f s\n', time_mat3);

if time_mat1 > 0, fprintf('   Speedup (Matrix-Shocks):     %.2f x\n', time_loop / time_mat1); end
if time_mat2 > 0, fprintf('   Speedup (Full-Matrix):       %.2f x\n', time_loop / time_mat2); end
if time_mat3 > 0, fprintf('   Speedup (Super-Matrix):      %.2f x\n', time_loop / time_mat3); end


%% Helper Function
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