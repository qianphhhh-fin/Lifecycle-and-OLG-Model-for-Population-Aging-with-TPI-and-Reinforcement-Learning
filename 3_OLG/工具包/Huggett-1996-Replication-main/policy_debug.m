% --- START OF FILE policy_debug.m ---

% Script to load and analyze saved policy functions from policy_debug.mat

clear; clc; close all;

fprintf('Loading policy functions for debugging...\n');

% Define the folder and filename
debugFolderName = 'tmp';
saveFilename = fullfile(debugFolderName, 'policy_debug.mat');

% Check if the file exists
if ~exist(saveFilename, 'file')
    error('File not found: %s\nPlease run the main model script first to generate the file.', saveFilename);
end

% Load the data
try
    load(saveFilename, 'kPolM', 'pPolM', 'qPolM', 'cPolM'); % Load only policy functions
    fprintf('Successfully loaded: kPolM, pPolM, qPolM, cPolM\n');
catch ME
    error('Failed to load data from %s: %s', saveFilename, ME.message);
end

% --- Basic Information ---
% Need some parameters to interpret the data. Let's assume standard 16 periods.
% You might need to load cS or define key parameters manually if needed.
aD = size(kPolM, 4);
nK = size(kPolM, 1);
nP = size(kPolM, 2);
nE = size(kPolM, 3);
retirementAgeIdx = 9; % Assuming standard 16-period setup
workingAgeMaxIdx = retirementAgeIdx - 1;
fprintf('Detected dimensions: nK=%d, nP=%d, nE=%d, aD=%d\n', nK, nP, nE, aD);

% --- Check for NaNs or Infs ---
fprintf('\nChecking for NaNs/Infs in policies:\n');
fprintf('  kPolM: %d NaNs, %d Infs\n', sum(isnan(kPolM(:))), sum(isinf(kPolM(:))));
fprintf('  pPolM: %d NaNs, %d Infs\n', sum(isnan(pPolM(:))), sum(isinf(pPolM(:))));
fprintf('  qPolM: %d NaNs, %d Infs\n', sum(isnan(qPolM(:))), sum(isinf(qPolM(:))));
fprintf('  cPolM: %d NaNs, %d Infs\n', sum(isnan(cPolM(:))), sum(isinf(cPolM(:))));
% Note: The code already replaces NaNs during creation/saving, so this should be 0.

% --- Analyze q Policy (Pension Contribution Rate) ---
fprintf('\nAnalyzing qPolM (Pension Contribution Rate):\n');
q_threshold = 1e-4; % Threshold for considering q as positive
total_working_states = nK * nP * nE * workingAgeMaxIdx;
total_positive_q = 0;

avg_q_by_age = zeros(workingAgeMaxIdx, 1);
positive_q_ratio_by_age = zeros(workingAgeMaxIdx, 1);

for a = 1:workingAgeMaxIdx
    q_slice = qPolM(:,:,:,a);
    valid_q_slice = q_slice(isfinite(q_slice)); % Exclude NaNs/Infs if any
    avg_q_by_age(a) = mean(valid_q_slice(:));
    num_positive_q_age = sum(valid_q_slice(:) > q_threshold);
    positive_q_ratio_by_age(a) = num_positive_q_age / numel(valid_q_slice);
    total_positive_q = total_positive_q + num_positive_q_age;
    fprintf('  Age %d: Avg q = %.6f, %% Positive q (>%.1e) = %.2f%%\n', ...
            a, avg_q_by_age(a), q_threshold, positive_q_ratio_by_age(a)*100);
end
fprintf('  Overall %% Positive q states (working age): %.2f%%\n', ...
        100 * total_positive_q / total_working_states);

% Plot average q rate by working age
figure('Name', 'Average q Policy');
bar(1:workingAgeMaxIdx, avg_q_by_age);
title('Average Pension Contribution Rate (q) by Working Age');
xlabel('Working Age Group (1 to 8)'); ylabel('Average q'); grid on;
xticks(1:workingAgeMaxIdx);

% --- Analyze p Policy (Next Period Pension Assets) ---
fprintf('\nAnalyzing pPolM (Next Period Pension Assets):\n');

% Check behavior at high p values (e.g., last p grid point)
p_high_idx = nP;
k_mid_idx = floor(nK/2)+1;
e_mid_idx = floor(nE/2)+1;

fprintf('  Checking p'' policy at high current p (pIdx=%d):\n', p_high_idx);
p_next_high_p_work = zeros(workingAgeMaxIdx, 1);
p_next_high_p_retire = zeros(aD - workingAgeMaxIdx, 1);

for a = 1:workingAgeMaxIdx
    % Get p' for mid k, high p, mid e
    p_next_high_p_work(a) = pPolM(k_mid_idx, p_high_idx, e_mid_idx, a);
    fprintf('    Work Age %d (mid k, high p, mid e): p'' = %.2f\n', a, p_next_high_p_work(a));
end
for a = retirementAgeIdx:aD
    % Get p' for mid k, high p, mid e
    retire_idx = a - workingAgeMaxIdx;
    p_next_high_p_retire(retire_idx) = pPolM(k_mid_idx, p_high_idx, e_mid_idx, a);
     fprintf('    Retire Age %d (mid k, high p, mid e): p'' = %.2f\n', a, p_next_high_p_retire(retire_idx));
end

% Visualize a slice of pPolM (e.g., mid e, fixed age)
target_age_p_plot = 5; % Example working age
target_e_p_plot = e_mid_idx;

if target_age_p_plot <= aD
    figure('Name', sprintf('pPolM Slice (Age %d, e %d)', target_age_p_plot, target_e_p_plot));
    pPolM_slice = pPolM(:,:,target_e_p_plot, target_age_p_plot);
    % Need kGridV and pGridV - load them or define roughly
    % For visualization, let's assume standard grid names might be available
    % If not, you'll need to load cS or define kGridV/pGridV here
    try
        load('tmp/policy_debug.mat', 'kPolM'); % Reload to potentially get grids if saved? No.
        % Define placeholder grids if necessary
        kGridV_viz = linspace(0, 200000, nK); % Use rough bounds if cS not loaded
        pGridV_viz = linspace(0, 200000, nP);
        [Kmesh, Pmesh] = meshgrid(kGridV_viz, pGridV_viz); % Create meshgrid

        surf(Kmesh, Pmesh, pPolM_slice'); % Transpose slice to match meshgrid
        xlabel('Current k'); ylabel('Current p'); zlabel('Next p''');
        title(sprintf('p'' Policy Slice (Age=%d, e=%d)', target_age_p_plot, target_e_p_plot));
        colorbar;
        view(3); % 3D view
        fprintf('  Plotted pPolM slice for Age %d, e %d.\n', target_age_p_plot, target_e_p_plot);
    catch ME_plot
        warning('Could not plot pPolM slice: %s\nNeed kGridV/pGridV definitions or load cS.', ME_plot.message);
    end
else
    fprintf('  Target age %d for plotting pPolM is out of bounds (max %d).\n', target_age_p_plot, aD);
end


% --- Analyze k Policy (Next Period Capital Assets) ---
fprintf('\nAnalyzing kPolM (Next Period Capital Assets):\n');
% Example: Check if k' generally increases with k for a given slice
k_policy_slice = kPolM(:,:,e_mid_idx, target_age_p_plot); % Use same slice as p' plot
diff_k_policy = diff(k_policy_slice, 1, 1); % Difference along k dimension
if all(diff_k_policy(:) >= -1e-6) % Allow for small numerical noise
    fprintf('  k'' policy appears generally non-decreasing with k for slice (Age %d, e %d).\n', target_age_p_plot, target_e_p_plot);
else
    fprintf('  Warning: k'' policy shows decreases with k for slice (Age %d, e %d).\n', target_age_p_plot, target_e_p_plot);
end

% --- Analyze c Policy (Consumption) ---
fprintf('\nAnalyzing cPolM (Consumption):\n');
% Example: Check if consumption is always positive
min_c_policy = min(cPolM(:));
max_c_policy = max(cPolM(:));
fprintf('  Min consumption in cPolM: %.4e\n', min_c_policy);
fprintf('  Max consumption in cPolM: %.4e\n', max_c_policy);
if min_c_policy < 1e-9
    fprintf('  Warning: Very low or non-positive consumption values found in cPolM.\n');
end


fprintf('\nPolicy Debug Analysis Complete.\n');

% --- END OF FILE policy_debug.m ---