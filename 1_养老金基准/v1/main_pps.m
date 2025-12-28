%% Main script for the Lifecycle Model with Pension Savings
%  This script sets up the model, solves it using value function iteration,
%  simulates lifecycle paths, and plots the results.

clear all;
close all;
clc;

fprintf('Starting the lifecycle model solution and simulation...\n\n');

% -------------------------------------------------------------------------
% 1. Set Model Parameters
% -------------------------------------------------------------------------
fprintf('1. Setting model parameters...\n');
tic;
cS = utils.set_parameters();
fprintf('   Parameters set. Elapsed time: %.2f seconds.\n\n', toc);

output_dir = 'debug';
if ~exist(output_dir, 'dir')
    mkdir(output_dir); 
end
save_dir = fullfile(output_dir, 'vfi_results_extended_model.mat');

% -------------------------------------------------------------------------
% 2. Solve the Model using Value Function Iteration
% -------------------------------------------------------------------------
if ~exist(save_dir, 'file')
    fprintf('2. Solving the household problem via Value Function Iteration...\n');
tic;
[polS, valS] = utils.value_function_iteration_matrix(cS);
fprintf('   VFI solved. Elapsed time: %.2f seconds.\n\n', toc);
fprintf('   VFI results saved to %s.\n\n', save_dir);
    save(save_dir, 'polS', 'valS', 'cS');
else
fprintf('2. Loading the household problem via Value Function Iteration...\n');
load(save_dir)
end


% Save the results to avoid re-computation



% -------------------------------------------------------------------------
% 3. Simulate Lifecycle Paths
% -------------------------------------------------------------------------
fprintf('3. Simulating lifecycle paths for a panel of households...\n');
tic;
simS = utils.simulation(polS, cS);
fprintf('   Simulation complete. Elapsed time: %.2f seconds.\n\n', toc);


% -------------------------------------------------------------------------
% 4. Plot Simulation Results
% -------------------------------------------------------------------------
fprintf('4. Generating and saving result plots...\n');
tic;
utils.plot_results(simS, cS);
fprintf('   Plots generated and saved to the /fig directory. Elapsed time: %.2f seconds.\n\n', toc);

fprintf('Lifecycle model execution finished successfully.\n');