%% Biochemical Reaction Network Parameter Fitting Code
% This code is used for fitting parameters of CRISPR-Cas system based biological switches
% Uses Gillespie algorithm for stochastic simulation

clear all;
close all;
clc;
warning('off');

%% ================== Load Experimental Data ==================
% Dataset A: A-node expression data during stimulation
experimental_data_A = [33.21034333, 159.5112333, 233.9074333, 276.2746, 285.7165333, ...
                       277.0111667, 252.9988333, 213.2478, 185.0001, 145.9864667, ...
                       107.04506, 90.17657667, 71.67451, 59.09033333, 52.08817, ...
                       47.95737333, 46.54682333, 46.51112667];

% Dataset B: B-node expression data during stimulation
experimental_data_B = [942.7497333, 962.2061, 923.3546667, 822.5150667, 680.4816333, ...
                       516.8535333, 394.9654, 273.5470667, 208.0181333, 143.5598, ...
                       117.1417, 100.0699133, 83.43736333, 72.91440333, 73.49182, ...
                       70.41011333, 62.96470667, 56.75531667];

% Dataset C: A-node expression data after stimulation
experimental_data_C = [126.0033667, 140.2573, 129.5871667, 125.7795, 141.622, ...
                       179.8475333, 204.1667, 202.3627667, 191.3148333, 162.2389667, ...
                       190.7813333, 203.8730667, 227.8121667, 250.6040667, 252.4869333, ...
                       311.4014667, 321.5470667, 317.6719];

% Time points (hours)
time_points = [0, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];

% Calculate weights (based on experimental data variance)
weights_A = calculateWeights(experimental_data_A);
weights_B = calculateWeights(experimental_data_B);
weights_C = ones(1, length(experimental_data_C)) / length(experimental_data_C);

%% ================== Define Reaction Network ==================
% Molecular species indices:
% 1:phiR73, 2:dCpf1, 3:gRNA, 4:dCg-phiR73, 5:mphiR73, 6:mdc, 
% 7:PphiR73, 8:Pdc, 9:Pg, 10:PphiR73-A, 11:Pdc-A, 12:Pg-A, 13:GFP, 14:mCherry

reaction_matrix = defineReactionMatrix();

%% ================== Simulation Parameter Settings ==================
simulation_params = struct();
simulation_params.max_steps = 8e5;           
simulation_params.num_cells = 1000;          
simulation_params.equilibrium_time = 4e5;    
simulation_params.stimulation_time = 2e5;    
simulation_params.recording_interval = 0.05 * 3600;  
simulation_params.time_delay = 0.5;          

%% ================== Parameter Fitting Main Loop ==================
best_loss = inf;
best_params = [];
num_iterations = 100;

% 基准参数
base_params = [8.92078E-08, 1.6858E-5, 0.0717, 0.922, 0.712, 0.0148, 11.716, 0.309];

fprintf('Starting parameter fitting...\n');

for iteration = 1:num_iterations
    fprintf('Iteration %d/%d\n', iteration, num_iterations);
    
    % Randomly perturb parameters
    current_params = perturbParameters(base_params, 0.1);
    
    % Run simulation
    [loss_total, simulation_results] = runSimulation(current_params, simulation_params, ...
                                                    experimental_data_A, experimental_data_B, ...
                                                    experimental_data_C, weights_A, weights_B, ...
                                                    weights_C, time_points);
    
    % Update best parameters
    if loss_total < best_loss
        best_loss = loss_total;
        best_params = current_params;
        fprintf('Found better parameters, loss function: %.6f\n', best_loss);
    end
    
    % Visualize results
    if mod(iteration, 1) == 0
        visualizeResults(simulation_results, iteration);
    end
end

fprintf('Parameter fitting completed!\n');
fprintf('Best parameters: [%.6e, %.6e, %.6e, %.6e, %.6e, %.6e, %.6e, %.6e]\n', best_params);
fprintf('Minimum loss function: %.6f\n', best_loss);

%% ================== Helper Functions ==================

function weights = calculateWeights(data)
    % Calculate variance-based weights
    weights = 1 ./ (data.^2);
    weights = weights / sum(weights);
end

function perturbed_params = perturbParameters(base_params, perturbation_factor)
    % Randomly perturb parameters
    perturbation = (rand(1, length(base_params)) * 2 * perturbation_factor - perturbation_factor + 1);
    perturbed_params = perturbation .* base_params;
end

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

function [total_loss, results] = runSimulation(params, sim_params, exp_data_A, exp_data_B, exp_data_C, ...
                                              weights_A, weights_B, weights_C, time_points)
    % Run Gillespie simulation
    
    % Setup cumate induction function
    cumate_params = [2.5, 1946.6, 0.013, 1.5];
    cumate_concentrations = [5, 1000];  % Low/High cumate
    
    % Calculate induction strength
    induction_strength = calculateInductionStrength(cumate_params, cumate_concentrations);
    
    % Initialize molecular concentrations
    initial_molecules = [5, 0, 0, 0, 0, 0, 40, 300, 60, 0, 0, 0, 0, 0];
    molecules = repmat(initial_molecules, sim_params.num_cells, 1);
    
    % Run simulation
    [recorded_data] = gillespieSimulation(molecules, params, sim_params, induction_strength);
    
    % Process data and calculate loss function
    [model_A, model_B, model_C] = processSimulationData(recorded_data, sim_params, time_points);
    
    % Calculate fitting loss
    loss_A = calculateFitScore(exp_data_A, weights_A, model_A);
    loss_B = calculateFitScore(exp_data_B, weights_B, model_B);
    loss_C = calculateFitScore(exp_data_C, weights_C, model_C);
    
    total_loss = loss_A + loss_B + loss_C;
    
    % Save results for visualization
    results = struct();
    results.model_A = model_A;
    results.model_B = model_B;
    results.model_C = model_C;
    results.recorded_data = recorded_data;
end

function induction_strength = calculateInductionStrength(params, concentrations)
    % Calculate Hill function induction strength
    induction_strength = params(1) * (concentrations.^params(4)) ./ ...
                        (params(2) + concentrations.^params(4)) + params(3);
end

function [recorded_data] = gillespieSimulation(molecules, params, sim_params, induction_strength)
    % Gillespie stochastic simulation algorithm
    
    num_cells = sim_params.num_cells;
    max_steps = sim_params.max_steps;
    
    % Initialize time and recorders
    time = zeros(num_cells, 1);
    current_induction = induction_strength(1) * ones(num_cells, 1);
    
    % Reaction matrix
    recorded_data = struct();
    recorded_data.pre_stimulation = [];
    recorded_data.during_stimulation = [];
    recorded_data.post_stimulation = [];
    
    % Reaction matrix
    reaction_matrix = defineReactionMatrix();
    
    for step = 1:max_steps
        % Display progress
        if mod(step, 80000) == 0
            fprintf('Step %d/%d\n', step, max_steps);
        end
        
        % Check if all cells have completed simulation
        total_time = sim_params.equilibrium_time + sim_params.stimulation_time + 72*3600;
        if sum(time > total_time) == length(time)
            break;
        end
        
        % Update induction strength
        current_induction = updateInductionStrength(time, sim_params, induction_strength);
        
        % Calculate reaction rates
        reaction_rates = calculateReactionRates(molecules, current_induction, params);
        
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

function reaction_rates = calculateReactionRates(molecules, induction, params)
    % Calculate rates for all reactions
    num_cells = size(molecules, 1);
    num_reactions = 20;
    reaction_rates = zeros(num_cells, num_reactions);
    
    % Set reaction rate parameters
    rate_params = [params(1)*induction, params(1)*induction, params(2)*ones(num_cells,1), ...
                   params(3)*ones(num_cells,1), params(4)*ones(num_cells,1), params(5)*ones(num_cells,1), ...
                   params(6)*ones(num_cells,1), params(7)*ones(num_cells,1), params(8)*ones(num_cells,1), ...
                   params(8)*ones(num_cells,1), params(8)*ones(num_cells,1), params(1)*induction, ...
                   params(1)*induction, params(2)*ones(num_cells,1), params(3)*ones(num_cells,1), ...
                   1e-7*ones(num_cells,1), params(6)*ones(num_cells,1), params(7)*ones(num_cells,1), ...
                   params(8)*ones(num_cells,1), params(8)*ones(num_cells,1)];
    
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
end

function [model_A, model_B, model_C] = processSimulationData(recorded_data, sim_params, time_points)
    % Process simulation data to match experimental time points
    % Data processing logic needs to be implemented here
    % Extract data at corresponding time points from recorded_data

    model_A = log(rand(1, length(time_points)) * 100 + 1);
    model_B = log(rand(1, length(time_points)) * 100 + 1);
    model_C = log(rand(1, length(time_points)) * 100 + 1);
end

function loss = calculateFitScore(experimental_data, weights, model_data)
    % Calculate fit score
    % Linear fitting

    X = [ones(length(model_data), 1), model_data'];
    beta = (X' * diag(weights) * X) \ (X' * diag(weights) * experimental_data');
    
    % Calculate predicted values
    predicted = beta(1) + beta(2) * model_data;
    
    % Calculate mean squared error
    loss = sum((log(predicted) - log(experimental_data)).^2);
end

function visualizeResults(results, iteration)
    % Visualize results
    
    figure('Name', sprintf('Iteration %d Results', iteration));
    
    subplot(2, 2, 1);
    plot(results.model_A, 'b-', 'LineWidth', 2);
    title('Node A Expression Model');
    xlabel('Time Point');
    ylabel('Expression Level');
    
    subplot(2, 2, 2);
    plot(results.model_B, 'r-', 'LineWidth', 2);
    title('Node B Expression Model');
    xlabel('Time Point');
    ylabel('Expression Level');
    
    subplot(2, 2, 3);
    plot(results.model_C, 'g-', 'LineWidth', 2);
    title('Node C Expression Model');
    xlabel('Time Point');
    ylabel('Expression Level');

    drawnow;
end