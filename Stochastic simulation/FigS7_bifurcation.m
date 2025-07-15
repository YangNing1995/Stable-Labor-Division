%% Biochemical Reaction Network Simulation Code - Fixed Cumate Concentration
% This code shows the simulation results of the bifurcation diagram of
% Fig.S7 Panel 1
% Uses Gillespie algorithm for stochastic simulation

%% ================== Simulation Parameter Settings ==================
simulation_params = struct();
simulation_params.max_steps = 8e5;           
simulation_params.num_cells = 100;          
simulation_params.equilibrium_time = 4e5;    
simulation_params.stimulation_time = 5e5;    
simulation_params.recording_interval = 0.05 * 3600;  
simulation_params.time_delay = 0.5;

%% ================== Set Current Parameters ==================
current_params = [8.92078E-08, 1.6858E-5, 0.0717, 0.922, 0.712, 0.0148, 11.716, 0.309];

%% ================== Fixed Cumate Concentrations ==================
cumate_concentrations = logspace(0, 3, 10);
num_concentrations = length(cumate_concentrations);

% Initialize results storage
steady_state_B = zeros(simulation_params.num_cells, num_concentrations);
all_simulation_results = cell(1, num_concentrations);

% Setup cumate induction parameters
cumate_params = [2.5, 1946.6, 0.013, 1.5];

fprintf('Starting simulations with fixed cumate concentrations...\n');

%% ================== Run Simulations for Each Concentration ==================
for i = 1:num_concentrations
    current_concentration = cumate_concentrations(i);
    fprintf('\n--- Running simulation for cumate concentration: %d ---\n', current_concentration);
    
    % Calculate induction strength for current concentration
    induction_strength = calculateInductionStrength(cumate_params, current_concentration);
    
    % Run simulation with fixed concentration
    [simulation_results, steady_state_value] = runFixedCumateSimulation(current_params, simulation_params, induction_strength);
    
    % Store results
    load("beta.mat")   % load the fitted parameters
    beta = [5.25473709749359, 0.0136075199728784];
    steady_state_B(:, i) = beta(1) + beta(2)*steady_state_value; 
    all_simulation_results{i} = simulation_results;
    
    fprintf('Steady state B node concentration: %.6f\n', steady_state_value);
end

%% ================== Save Results ==================
results_data = struct();
results_data.cumate_concentrations = cumate_concentrations;
results_data.steady_state_B = steady_state_B;
results_data.all_simulation_results = all_simulation_results;
results_data.simulation_params = simulation_params;
results_data.current_params = current_params;

% Save to file
save('cumate_concentration_results.mat', 'results_data');
fprintf('\nResults saved to cumate_concentration_results.mat\n');

%% ================== Visualize Results ==================
base_line = 50;
threshold = 6.5;

low_states  = zeros(1, length(cumate_concentrations));
high_states = zeros(1, length(cumate_concentrations));

for j = 1: length(cumate_concentrations)
    column_data = steady_state_B(:, j);

    low_state  = column_data(column_data <= threshold);
    high_state = column_data(column_data > threshold);

    if ~isempty(low_state)
        low_states(j) = mean(low_state);
    else
        low_states(j) = NaN;
    end

    if ~isempty(high_state)
        high_states(j) = mean(high_state);
    else
        high_states(j) = NaN;
    end
end

% Plot cumate concentration vs steady state B node concentration
figure('Name', 'Cumate Concentration vs Steady State B Node', 'Position', [100 100 600 450]);
loglog(cumate_concentrations, exp(low_states) - base_line , ...
    'ro', 'LineWidth', 1, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
hold on;
loglog(cumate_concentrations, exp(high_states) - base_line , ...
    'bs', 'LineWidth', 1, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
grid on;
xlabel('Cumate Concentration');
ylabel('Steady State B Node Expression');
ylim([5e1, 1e4])
set(gca, 'FontSize', 12);
legend({'Low state', 'High state'}, 'Location', 'best');

%% ================== Plot Time Series for Each Concentration ==================
figure('Name', 'Time Series for Different Cumate Concentrations', 'Position', [200 200 1200 800]);

for i = 1:num_concentrations
    subplot(ceil(num_concentrations / 3), 3, i);

    % Get time series data
    time_series = all_simulation_results{i}.time_series;
    time_points = all_simulation_results{i}.time_points;

    plot(time_points, time_series, 'LineWidth', 2);
    grid on;
    xlabel('Time (hours)');
    ylabel('B Node Expression');
    title(sprintf('Cumate Concentration: %d', cumate_concentrations(i)));
    set(gca, 'FontSize', 10);
end

%% ================== Helper Functions ==================

function induction_strength = calculateInductionStrength(params, concentration)
    % Calculate Hill function induction strength for a single concentration
    induction_strength = params(1) * (concentration^params(4)) / ...
                        (params(2) + concentration^params(4)) + params(3);
end

function [simulation_results, steady_state_value] = runFixedCumateSimulation(params, sim_params, induction_strength)
    % Run simulation with fixed cumate concentration
    
    % Initialize molecular concentrations
    initial_molecules = [5, 0, 0, 0, 0, 0, 40, 300, 60, 0, 0, 0, 0, 0];
    molecules = repmat(initial_molecules, sim_params.num_cells, 1);
    
    % Run simulation
    [recorded_data, time_series_data] = gillespieSimulationFixed(molecules, params, sim_params, induction_strength);
    
    % Process results
    simulation_results = struct();
    simulation_results.recorded_data = recorded_data;
    simulation_results.time_series = time_series_data.B_node;
    simulation_results.time_points = time_series_data.time_points;
    
    % Calculate steady state (average of last 50% of simulation)
    steady_state_start = round(0.5 * length(time_series_data.time_points));
    rec = time_series_data.B_node(:, steady_state_start:end);
    steady_state_value =  mean(rec, 2);

end

function [recorded_data, time_series_data] = gillespieSimulationFixed(molecules, params, sim_params, induction_strength)
    % Gillespie stochastic simulation algorithm with fixed induction
    
    num_cells = sim_params.num_cells;
    max_steps = sim_params.max_steps;
    
    % Initialize time and recorders
    time = zeros(num_cells, 1);
    
    % Time series recording
    max_time_points = 1000;
    time_series_data = struct();
    time_series_data.B_node = zeros(num_cells, max_time_points);
    time_series_data.time_points = zeros(1, max_time_points);
    time_series_counter = 1;
    next_recording_time = 0;
    recording_interval = 600; % Record every 10 minute
    
    % Reaction matrix
    record_space = 4000;
    recorded_data = struct();
    recorded_data.molecules = zeros(num_cells, record_space, 14);
    recorded_data.time_points = zeros(1, record_space);
    record_counter = 1;
    
    % Reaction matrix
    reaction_matrix = defineReactionMatrix();
    
    % Total simulation time (longer for steady state)
    total_simulation_time = 96 * 3600; 
    
    for step = 1:max_steps
        % Display progress
        if mod(step, 80000) == 0
            fprintf('Step %d/%d, Time: %.2f hours\n', step, max_steps, mean(time)/3600);
        end
        
        % Check if all cells have completed simulation
        if sum(time > total_simulation_time) == length(time)
            break;
        end
        
        % Calculate reaction rates (fixed induction strength)
        reaction_rates = calculateReactionRates(molecules, induction_strength * ones(num_cells, 1), params);
        
        % Execute Gillespie step
        [molecules, time] = gillespieStep(molecules, time, reaction_rates, reaction_matrix);
        
        % Record time series data
        if mean(time) >= next_recording_time && time_series_counter <= max_time_points
            % Record B node concentration (mCherry, index 14)
            time_series_data.B_node(:, time_series_counter) = molecules(:, 14);
            time_series_data.time_points(time_series_counter) = mean(time) / 3600; % Convert to hours
            time_series_counter = time_series_counter + 1;
            next_recording_time = next_recording_time + recording_interval;
        end
        
        % Record detailed data occasionally
        if mod(step, 10000) == 0 && record_counter <= record_space
            recorded_data.molecules(:, record_counter, :) = molecules;
            recorded_data.time_points(record_counter) = mean(time);
            record_counter = record_counter + 1;
        end
    end
    
    % Trim unused space
    time_series_data.B_node = time_series_data.B_node(:, 1:time_series_counter-1);
    time_series_data.time_points = time_series_data.time_points(1:time_series_counter-1);
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

function reaction_rates = calculateReactionRates(molecules, induction, params)
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