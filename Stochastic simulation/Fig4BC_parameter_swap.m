%% Biochemical Reaction Network Parameter Swap Code
% This code is used for parameter swap
% Uses Gillespie algorithm for stochastic simulation


%% ================== Define Reaction Network ==================
% Molecular species indices:
% 1:phiR73, 2:dCpf1, 3:gRNA, 4:dCg-phiR73, 5:mphiR73, 6:mdc, 
% 7:PphiR73, 8:Pdc, 9:Pg, 10:PphiR73-A, 11:Pdc-A, 12:Pg-A, 13:GFP, 14:mCherry

reaction_matrix = defineReactionMatrix();

%% ================== Simulation Parameter Settings ==================
simulation_params = struct();
simulation_params.max_steps = 8e5;           
simulation_params.num_cells = 50;          
simulation_params.equilibrium_time = 4e5;    
simulation_params.stimulation_time = 2e5;    
simulation_params.recording_interval = 0.1 * 3600;  
simulation_params.time_delay = 0.5;          

 %% ================== Parameter Swap Main Loop ==================
% Base parameters
base_params = [8.92078E-08, 1.6858E-5, 0.0717, 0.922, 0.712, 0.0148, 11.716, 0.309];

% Perturbation settings
num_params = 20;           % Number of parameters to perturb
num_perturb = 21;          % Number of perturbation rates
pert_rate_list = logspace(-1, 1, num_perturb); % Perturbation strength from 0.1 to 10

% Initialize results storage
paraswap_recorder = zeros(num_params, num_perturb, 2);

% Progress tracking
fprintf('Starting parameter fitting...\n');
fprintf('Total iterations: %d\n', num_params * num_perturb);
start_time = tic;

% Main fitting loop
for param_index = 1:num_params
    fprintf('\nProcessing parameter %d/%d\n', param_index, num_params);
    
    for pert_rate_index = 1:num_perturb
        pert_rate = pert_rate_list(pert_rate_index);
        
        % Progress indicator
        current_iter = (param_index - 1) * num_perturb + pert_rate_index;
        total_iter = num_params * num_perturb;
        fprintf('  Perturbation %d/%d (rate: %.3f) - Progress: %.1f%%\n', ...
                pert_rate_index, num_perturb, pert_rate, ...
                100 * current_iter / total_iter);
        
        % Run simulation
        simulation_results = runSimulation(base_params, simulation_params, ...
                                         param_index, pert_rate);
        
        % Process simulation data
        paraswap_result = ParaSwapProcess(simulation_results);
        
        % Store results
        paraswap_recorder(param_index, pert_rate_index, :) = paraswap_result;

    end
    
    % Progress update
    elapsed_time = toc(start_time);
    estimated_total = elapsed_time * num_params / param_index;
    remaining_time = estimated_total - elapsed_time;
    fprintf('  Parameter %d completed. Elapsed: %.1fs, Estimated remaining: %.1fs\n', ...
            param_index, elapsed_time, remaining_time);
end

adp_score = paraswap_recorder(:, :, 1);
bif_score = paraswap_recorder(:, :, 2);
% save('../results/adp_score.mat', 'adp_score')
% save('../results/bif_score.mat', 'bif_score')

%% ================== Visualize Results ==================
% load('../results/adp_score.mat')
% load('../results/bif_score.mat')
pert_rate_log_list = linspace(-1, 1, num_perturb);
Rsquared_adp = corr(pert_rate_log_list', adp_score').^2;
Rsquared_bif = corr(pert_rate_log_list', bif_score').^2;

% Fig4B - Response Sensitivity Plot
figure('Position', [100, 100, 800, 600]);
semilogy(pert_rate_log_list, adp_score(1:16, :), 'LineWidth', 1.5);
xlabel('log_{10}(Para_{swap}/Para_{ori})', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Node A - Max Response', 'FontSize', 12, 'FontWeight', 'bold');
title('Response Sensitivity', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 11);
% saveas(gcf, '../results/Response_sensitivity.png');

% Fig4B - Proportion Sensitivity Plot  
figure('Position', [100, 100, 800, 600]);
plot(pert_rate_log_list, bif_score(1:16, :), 'LineWidth', 1.5);
xlabel('log_{10}(Para_{swap}/Para_{ori})', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Node B - High Proportion', 'FontSize', 12, 'FontWeight', 'bold');
title('Proportion Sensitivity', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 11);
% saveas(gcf, '../results/Proportion_sensitivity.png');

% Fig4C - Response Sensitivity R^2 vs Proportion Sensitivity R^2
figure('Position', [100, 100, 600, 600]);
scatter(Rsquared_adp(1:16), Rsquared_bif(1:16), 50, 'LineWidth', 2);
xlabel('Response Sensitivity R^2', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Proportion Sensitivity R^2', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 11);
axis square
% saveas(gcf, '../results/R_squared.png');


%% ================== Helper Functions ==================
function reaction_matrix = defineReactionMatrix()
    % Define reaction matrix
    % Each row represents a reaction, each column represents stoichiometric coefficient of a molecular species
    reaction_matrix = zeros(20, 14);
    
    % Binding reactions
    reaction_matrix(1,:) = [-2, 0, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0];  % phiR73 + Pdc → Pdc-A
    reaction_matrix(2,:) = [-2, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0];  % phiR73 + Pg → Pg-A
    reaction_matrix(3,:) = [-2, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0];  % phiR73 + PphiR73 → PphiR73-A
    
    % Dissociation reactions
    reaction_matrix(4,:) = [2, 0, 0, 0, 0, 1, 0, 1, 0, 0, -1, 0, 0, 0];   % Pdc-A → phiR73 + Pdc + mdc
    reaction_matrix(5,:) = [2, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0];   % Pg-A → phiR73 + Pg + gRNA
    reaction_matrix(6,:) = [2, 0, 0, 0, 1, 0, 1, 0, 0, -1, 0, 0, 0, 0];   % PphiR73-A → phiR73 + PphiR73 + mphiR73
    
    % Transcription reactions
    reaction_matrix(7,:) = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];    % mdc → mdc + dCpf1
    reaction_matrix(8,:) = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];    % mphiR73 → mphiR73 + phiR73
    
    % Degradation reactions
    reaction_matrix(9,:) = [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0];   % mdc → ∅
    reaction_matrix(10,:) = [0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0];  % mphiR73 → ∅
    reaction_matrix(11,:) = [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];  % gRNA → ∅
    reaction_matrix(12,:) = [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];  % phiR73 → ∅
    reaction_matrix(13,:) = [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];  % dCpf1 → ∅
    
    % CRISPR complex formation
    reaction_matrix(14,:) = [0, -1, -1, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0]; % dCpf1 + gRNA + PphiR73 → dCg-phiR73
    reaction_matrix(15,:) = [0, 1, 1, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0];   % dCg-phiR73 → PphiR73 + dCpf1 + gRNA
    
    % Reporter gene expression
    reaction_matrix(16,:) = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0];    % PphiR73 → mphiR73 + PphiR73
    reaction_matrix(17,:) = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0];    % gRNA → gRNA + GFP
    reaction_matrix(18,:) = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1];    % mphiR73 → mphiR73 + mCherry
    reaction_matrix(19,:) = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0];   % GFP → ∅
    reaction_matrix(20,:) = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1];   % mCherry → ∅
end

function results = runSimulation(params, sim_params, pert_index, pert_rate)
    % Run Gillespie simulation
    
    % Setup cumate induction function
    cumate_params = [2.5, 1946.6, 0.013, 1.5];
    cumate_concentrations = [5, 1000];  % Low/High cumate concentration
    
    % Calculate induction strength
    induction_strength = calculateInductionStrength(cumate_params, cumate_concentrations);
    
    % Initialize molecular concentrations
    initial_molecules = [5, 0, 0, 0, 0, 0, 40, 300, 60, 0, 0, 0, 0, 0];
    molecules = repmat(initial_molecules, sim_params.num_cells, 1);
    
    % Run simulation
    recorded_data = gillespieSimulation(molecules, params, sim_params, induction_strength, pert_index, pert_rate);
    results = struct();
    results.recorded_data = recorded_data;
    
end

function induction_strength = calculateInductionStrength(params, concentrations)
    % Calculate Hill function induction strength
    induction_strength = params(1) * (concentrations.^params(4)) ./ ...
                        (params(2) + concentrations.^params(4)) + params(3);
end

function [recorded_data] = gillespieSimulation(molecules, params, sim_params, induction_strength, pert_index, pert_rate)
    % Gillespie stochastic simulation algorithm
    
    num_cells = sim_params.num_cells;
    max_steps = sim_params.max_steps;
    
    % Initialize time and recorders
    time = zeros(num_cells, 1);
    current_induction = induction_strength(1) * ones(num_cells, 1);
    
    % Reaction matrix
    record_space = 4000;
    recorded_data = struct();
    recorded_data.pre_stimulation = zeros(num_cells,record_space, 4);
    recorded_data.during_stimulation = zeros(num_cells,record_space,4);
    recorded_data.post_stimulation = zeros(num_cells,record_space, 4);
    
    % Reaction matrix
    reaction_matrix = defineReactionMatrix();
    
    for step = 1:max_steps
        % Display progress
        if mod(step, 80000) == 0
            fprintf('Step %d/%d\n', step, max_steps);
        end
        
        % Update induction strength
        current_induction = updateInductionStrength(time, sim_params, induction_strength);
        
        % Calculate reaction rates
        reaction_rates = calculateReactionRates(molecules, current_induction, params, pert_index, pert_rate);

        % Execute Gillespie step
        [molecules, time] = gillespieStep(molecules, time, reaction_rates, reaction_matrix);
        
        % Record data
        recorded_data = recordData(recorded_data, molecules, time, sim_params);
    end
end

function current_induction = updateInductionStrength(time, sim_params, induction_strength)
    % Update induction strength
    current_induction = induction_strength(1) * ones(length(time), 1);
    
    stimulation_start = sim_params.equilibrium_time;
    stimulation_end = stimulation_start + sim_params.stimulation_time;
    
    % Use high induction strength during stimulation
    stim_mask = (time >= stimulation_start) & (time < stimulation_end);
    current_induction(stim_mask) = induction_strength(2);
end

function reaction_rates = calculateReactionRates(molecules, induction, params, pert_index, pert_rate)
    % Calculate rates for all reactions
    num_cells = size(molecules, 1);
    num_reactions = 20;
    reaction_rates = zeros(num_cells, num_reactions);
    % Set reaction rate parameters
    rate_params = [params(1)*induction, params(1)*induction, params(2)*ones(num_cells,1), ...
                   params(3)*ones(num_cells,1), params(3)*ones(num_cells,1), params(3)*ones(num_cells,1), ...
                   params(4)*ones(num_cells,1), params(4)*ones(num_cells,1), params(5)*ones(num_cells,1), ...
                   params(5)*ones(num_cells,1), params(5)*ones(num_cells,1), params(6)*ones(num_cells,1), ...
                   params(6)*ones(num_cells,1), params(7)*ones(num_cells,1), params(8)*ones(num_cells,1), ...
                   1e-7*ones(num_cells,1), params(4)*ones(num_cells,1), params(4)*ones(num_cells,1), ...
                   params(6)*ones(num_cells,1), params(6)*ones(num_cells,1)];

    rate_params(:, pert_index) = rate_params(:, pert_index,:)*pert_rate;

    % Calculate rate for each reaction
    for i = 1:num_reactions
        reaction_rates(:, i) = calculateSingleReactionRate(molecules, rate_params(:, i), i);
    end
end

function rate = calculateSingleReactionRate(molecules, k, reaction_index)
    switch reaction_index
        case 1  % phiR73 + Pdc → Pdc-A
            rate = k .* molecules(:,1) .* (molecules(:,1)-1) .* molecules(:,8);
        case 2  % phiR73 + Pg → Pg-A
            rate = k .* molecules(:,1) .* (molecules(:,1)-1) .* molecules(:,9);
        case 3  % phiR73 + PphiR73 → PphiR73-A
            rate = k .* molecules(:,1) .* (molecules(:,1)-1) .* molecules(:,7);
        case 4  % Pdc-A → products
            rate = k .* molecules(:,11);
        case 5  % Pg-A → products
            rate = k .* molecules(:,12);
        case 6  % PphiR73-A → products
            rate = k .* molecules(:,10);
        case 7  % mdc → mdc + dCpf1
            rate = k .* molecules(:,6);
        case 8  % mphiR73 → mphiR73 + phiR73
            rate = k .* molecules(:,5);
        case 9  % mdc degradation
            rate = k .* molecules(:,6);
        case 10 % mphiR73 degradation
            rate = k .* molecules(:,5);
        case 11 % gRNA degradation
            rate = k .* molecules(:,3);
        case 12 % phiR73 degradation
            rate = k .* molecules(:,1);
        case 13 % dCpf1 degradation
            rate = k .* molecules(:,2);
        case 14 % CRISPR complex formation
            rate = k .* molecules(:,2) .* molecules(:,7) .* molecules(:,3);
        case 15 % CRISPR complex dissociation
            rate = k .* molecules(:,4);
        case 16 % PphiR73 transcription
            rate = k .* molecules(:,7);
        case 17 % GFP expression
            rate = k .* molecules(:,3);
        case 18 % mCherry expression
            rate = k .* molecules(:,5);
        case 19 % GFP degradation
            rate = k .* molecules(:,13);
        case 20 % mCherry degradation
            rate = k .* molecules(:,14);
        otherwise
            rate = zeros(size(molecules, 1), 1);
    end
end

function [molecules, time] = gillespieStep(molecules, time, reaction_rates, reaction_matrix)
    % Execute one step of Gillespie algorithm
    num_cells = size(molecules, 1);
    
    % Calculate total reaction rates
    total_rates = sum(reaction_rates, 2);
    
    % Calculate time step
    dt = exp(min(total_rates.^(-1), 1));
    time = time + dt;
    
    % Select reaction
    random_vals = rand(num_cells, 1);
    cumulative_probs = cumsum(reaction_rates, 2) ./ total_rates;
    
    for cell = 1:num_cells
        if total_rates(cell) > 0
            reaction_idx = find(random_vals(cell) <= cumulative_probs(cell, :), 1);
            if ~isempty(reaction_idx)
                molecules(cell, :) = molecules(cell, :) + reaction_matrix(reaction_idx, :);
            end
        end
    end
end

function recorded_data = recordData(recorded_data, molecules, time, sim_params)
    % Record simulation data
    % Data recording logic needs to be implemented here
    % Record data into different arrays based on time periods
    % Simplified version, actual implementation requires more complex logic
    recorded_data.molecules = molecules;
    recorded_data.time = time;

    for times =1: length(time)
        t_cell = time(times);
            if t_cell < sim_params.equilibrium_time && floor(t_cell/sim_params.recording_interval)>recorded_data.pre_stimulation(times,1,4)
                recorded_data.pre_stimulation(times,1,4) = floor(t_cell/sim_params.recording_interval);
                recorded_data.pre_stimulation(times,floor(t_cell/sim_params.recording_interval),1:3) = [molecules(times,14),molecules(times,13),time(times)];
                continue
            end
        t_cell = time(times)-sim_params.equilibrium_time;
            if t_cell>= 0 && t_cell < sim_params.stimulation_time && floor(t_cell/sim_params.recording_interval)>recorded_data.during_stimulation(times,1,4)
                recorded_data.during_stimulation(times,1,4) = floor(t_cell/sim_params.recording_interval);
                recorded_data.during_stimulation(times,floor(t_cell/sim_params.recording_interval),1:3) = [molecules(times,14),molecules(times,13),time(times)];
                continue
            end
        t_cell = time(times)-sim_params.equilibrium_time-sim_params.stimulation_time;
            if t_cell >= 0 && floor(t_cell/sim_params.recording_interval)>recorded_data.post_stimulation(times,1,4)
                recorded_data.during_stimulation(times,1,4) = floor(t_cell/sim_params.recording_interval);
                recorded_data.post_stimulation(times,floor(t_cell/sim_params.recording_interval),1:3) = [molecules(times,14),molecules(times,13),time(times)];
                continue
            end

    end

end

function paraswap_result = ParaSwapProcess(sim_results)

    % This function analyzes simulation data to compute two key metrics:
    % 1. Adaptation score: Maximum response during stimulation
    % 2. Bifurcation score: Fraction of values above threshold in pre-stimulation data
    %
    % Input:
    %   sim_results - Structure containing simulation data with fields:
    %                 - recorded_data.during_stimulation: Data during stimulation phase
    %                 - recorded_data.pre_stimulation: Data before stimulation phase
    %
    % Output:
    %   paraswap_result - 1x2 array containing [adaptation_score, bifurcation_score]

    % ========== ADAPTATION SCORE CALCULATION ==========
    % Extract data from node A during stimulation (3rd dimension, index 2)
    stimulation_data = sim_results.recorded_data.during_stimulation(:, :, 2);
    
    % Find maximum value across columns 1-555 for each row
    max_per_row = max(stimulation_data(:, 1:555), [], 2);
    
    % Calculate adaptation score as the overall maximum
    adaptation_score = max(max_per_row);
    
    % ========== BIFURCATION SCORE CALCULATION ==========
    % Extract data from node B before stimulation (3rd dimension, index 1)
    pre_stimulation_data = sim_results.recorded_data.pre_stimulation(:, :, 1);
    
    % Focus on time points 500-1000 for analysis
    analysis_window = pre_stimulation_data(:, 500:1000);
    
    % Calculate bifurcation score using fixed threshold of 50
    % This represents the fraction of values above threshold
    values_above_threshold = sum(analysis_window(:) > 50);
    total_valid_values = sum(analysis_window(:) >= 0);
    bifurcation_score = values_above_threshold / total_valid_values;
    
    % ========== RETURN RESULTS ==========
    % Combine both scores into output array
    paraswap_result = [adaptation_score, bifurcation_score];
    
end
