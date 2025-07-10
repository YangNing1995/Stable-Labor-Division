% This program exhaustively enumerates two-node topologies achieving both bistability and adaptation by Ning Yang, 2023/1/10
% Node A is the output node receiving external input, Node B is the control node
% The two-node network has 4 edges: A->A, A->B, B->A, B->B, totaling 3^4 possible topologies
% Each topology is represented by a 4x3 matrix where rows correspond to edges in order: 
% [A->A; A->B; B->A; B->B], and columns use one-hot encoding for regulation types: [activation, repression, no regulation]
% The system uses AND logic transcriptional regulation with Hill equation
% parameters (K, n) and degradation coefficients Gamma
% Added screening for high-to-low state transitions from original program

%% Generate topology matrix
TopologyTensor = zeros(3^4, 4, 3); % Tensor storing all possible topologies (81x4x3)
TopologyIndex = double(dec2base(0:3^4-1, 3, 4)) - 48; % Convert to base-3 encoded indices

% Convert indices to one-hot encoded regulation types
for i = 0:2
    switch i
        case 0
            OneHotEncode = [1,0,0];    % Activation encoding
        case 1
            OneHotEncode = [0,1,0];    % Repression encoding
        case 2
            OneHotEncode = [0,0,1];    % No regulation encoding
    end
    [Row, Col] = find(TopologyIndex==i);
    for j = 1:length(Row)
        TopologyTensor(Row(j), Col(j), :) = OneHotEncode;
    end
end

%% Parameter sampling using Latin Hypercube Sampling 
N_sample = 1e4;       % Number of parameter sets to sample
K_lower = 1e-2;
K_upper = 1e2;
n_lower = 1;
n_upper = 4;
Gamma_lower = 1e-2;
Gamma_upper = 1;

X = lhsdesign(N_sample, 10); 
K_list = K_lower*(K_upper/K_lower).^X(:, 1:4);
n_list = n_lower + (n_upper - n_lower).*X(:, 5:8);
Gamma_list = Gamma_lower*(Gamma_upper/Gamma_lower).^X(:, 9:10);
% save('ParameterSet_1e4.mat', 'Gamma_list',"n_list","K_list");

%% Exhaustive enumeration of topologies
N_sample = 1e4;       % Number of parameter sets to sample
load('ParameterSet_1e4.mat'); % Pre-loaded parameters: K_list, n_list, Lambda_list

Result_final = cell(3^4, 1); % Cell array to store successful parameters for each topology
% Simulation parameters
I0 = 1;                % Base input stimulus
limit = 1e-2;          % Convergence threshold
Bistability_threshold = 5e-2;     % Bistability threshold
Error_threshold = 0.1;        % Error threshold 
Response_threshold = 0.2;       % Response threshold
tspan = 0:100;         % Simulation time points
x0_range = logspace(-2,1,10); % Initial condition range (10^−2 to 10^1)

% Main enumeration loop
for n = 1:3^4
    tic
    topology = squeeze(TopologyTensor(n, :, :)); % Current topology
    Result_temp = {}; % Temporary storage for successful parameters
    
    for m = 1:N_sample
        % Set ODE parameters from sampled set
        p.K_AA = K_list(m, 1); p.K_AB = K_list(m, 2);
        p.K_BA = K_list(m, 3); p.K_BB = K_list(m, 4);
        p.n_AA = n_list(m, 1); p.n_AB = n_list(m, 2);
        p.n_BA = n_list(m, 3); p.n_BB = n_list(m, 4);
        p.Gamma_A = Gamma_list(m, 1); p.Gamma_B = Gamma_list(m, 2);
        
        % Test multiple initial conditions
        x_final = zeros(length(x0_range)^2, 2); % Store final states
        convergence_failures = 0; 
        
        for i = 1:length(x0_range)
            for j = 1:length(x0_range)
                % Simulate system from current initial condition
                [~, x] = ode23s(@(t,x) ODE_circuit(t, x, I0, p, topology), tspan, [x0_range(i), x0_range(j)]);
                
                % Check convergence
                if all(abs(x(end,:) - x(end-1,:)) < limit)
                    x_final(i*(length(x0_range)-1)+j, :) = x(end, :);
                else % Extend simulation time if not converged
                    [~, x] = ode23s(@(t,x) ODE_circuit(t, x, I0, p, topology), 10*tspan, x(end,:));  
                    if all(abs(x(end,:) - x(end-1,:)) < limit)
                        x_final(i*(length(x0_range)-1)+j, :) = x(end, :); 
                    else
                        convergence_failures = convergence_failures + 1;
                    end
                end
            end
        end
        
        % Skip if majority of initial conditions fail to converge
        if convergence_failures/(length(x0_range)^2) > 0.9
            continue
        end
        
        % Bistability screening
        x_min = min(x_final); % Minimum steady states
        x_max = max(x_final); % Maximum steady states
        [~, dominant_node] = max(x_max - x_min); % Node with largest state difference
        
        % Check bistability criteria
        if (x_max(dominant_node) - x_min(dominant_node)) > Bistability_threshold
            % Count steady state populations
            num_ss = sum(abs(x_final(:,dominant_node) - x_min(dominant_node)) < limit) + ...
                     sum(abs(x_final(:,dominant_node) - x_max(dominant_node)) < limit);
            
            if num_ss/(length(x0_range)^2) > 0.9
                % Adaptation screening using high-input stimulus
                high_input = 1e3*I0;
                ss_index = find(x_final(:, dominant_node) == x_max(dominant_node), 1);
                x0_adapt = x_final(ss_index, :);
                
                % Simulate adaptation response
                [~, x] = ode23s(@(t,x) ODE_circuit(t, x, high_input, p, topology), tspan, x0_adapt);
                
                % Calculate adaptation metrics
                x_peak = max(x(:,1));
                Deviation = abs(x(end, 1) - x0_adapt(1));    % Final deviation from baseline
                response = abs(x_peak - x0_adapt(1));     % Peak response magnitude
                error = x(end, dominant_node) - x_min(dominant_node); % Return to low state
                
                % Check adaptation criteria
                if Deviation < limit && response > Response_threshold && error < Error_threshold*(x_max(dominant_node)-x_min(dominant_node))
                    Result_temp{end+1} = p;
                    disp(['Valid parameter set found for topology ', num2str(n)])
                end
            end
        end
    end
    
    Result_final{n} = Result_temp;
    save('Enumeration_test.mat', "Result_final"); % Periodic saving of results
    toc
end

%% ODE System Definition
function dydt = ODE_circuit(t, x, I0, p, topology)
% ODE system for two-node network with AND logic regulation
% Inputs:
%   t: time
%   x: state vector [A; B]
%   I0: input stimulus
%   p: parameter struct (K, n, Lambda)
%   topology: 4x3 matrix defining regulation types

    A = x(1);
    B = x(2);
    
    % Calculate Hill functions for each interaction
    H_AA = A^p.n_AA / (A^p.n_AA + p.K_AA^p.n_AA);
    H_AB = A^p.n_AB / (A^p.n_AB + p.K_AB^p.n_AB);
    H_BA = B^p.n_BA / (B^p.n_BA + p.K_BA^p.n_BA);
    H_BB = B^p.n_BB / (B^p.n_BB + p.K_BB^p.n_BB);
    
    % Calculate regulation terms using topology matrix
    v_AA = topology(1,1)*H_AA + topology(1,2)*(1 - H_AA) + topology(1,3);
    v_AB = topology(2,1)*H_AB + topology(2,2)*(1 - H_AB) + topology(2,3);
    v_BA = topology(3,1)*H_BA + topology(3,2)*(1 - H_BA) + topology(3,3);
    v_BB = topology(4,1)*H_BB + topology(4,2)*(1 - H_BB) + topology(4,3);

    % ODE equations
    dydt = zeros(2,1);
    dydt(1) = p.Gamma_A * (I0 * v_AA * v_BA - A);  % Node A dynamics
    dydt(2) = p.Gamma_B * (v_AB * v_BB - B);       % Node B dynamics
end