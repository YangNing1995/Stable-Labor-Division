%% Bifurcation diagram
% B activates A, A inhibits B, B self-activates
% Plot nullcline for node B
% Hyperparameters
Gamma_A = 1; Gamma_B = 1;        % Degradation rate
Beta_A = 1;  Beta_B = 2;           % Maximal production rate
n_BA = 2; n_AB = 1; n_BB = 4;      % Hill coefficients
K_BA = 50; K_AB = 1; K_BB = 0.9;   % Binding affinities constant

% Find steady states using fsolve
I_seq = logspace(0, 6, 100);
x_high_seq = zeros(2, length(I_seq));      % Record high-state steady values
flag_high_seq = zeros(length(I_seq), 1);   % Track high-state exit flags
x_low_seq = zeros(2, length(I_seq));       % Record low-state steady values
flag_low_seq = zeros(length(I_seq), 1);    % Track low-state exit flags
                                                   
x0_high = [1, 2];                        % Initial guess for high state
x0_low = [0.1, 1];                       % Initial guess for low state
opt = optimset('TolFun',1e-10,'TolX', 1e-10, 'Algorithm' , 'trust-region');  % fsolve options

for i = 1: length(I_seq)
    I = I_seq(i);
    [x_high,~,exitflag_high,~] = fsolve(@(x)circuit(x, I, Gamma_A, Gamma_B, Beta_A, Beta_B, n_BA, n_AB, n_BB, K_BA, K_AB, K_BB), x0_high,opt);
    [x_low,~,exitflag_low,~] = fsolve(@(x)circuit(x, I, Gamma_A, Gamma_B, Beta_A, Beta_B, n_BA, n_AB, n_BB, K_BA, K_AB, K_BB), x0_low,opt);
    
    x_high_seq(:, i) = x_high;
    x_low_seq(:, i) = x_low;
    flag_high_seq(i) = exitflag_high;
    flag_low_seq(i) = exitflag_low; 
end

Index_high = find(flag_high_seq==-2, 1);  % flag=-2 indicates non-convergence; discard subsequent results
Index_low = find(flag_low_seq==-2, 1);   

%---------------------------------- Node B ---------------------------------------%
color_seq = {'#0A5D89','#7F76A2'}; % Blue-green, wisteria-gray

figure('unit','points','PaperUnits','points', 'position',[500 300 400 150]);
p1 = semilogx(I_seq(1:Index_high), x_high_seq(2, 1:Index_high),'-','lineWidth',3, 'Color', color_seq{1});
hold on 
p2 = semilogx(I_seq(1:Index_low), x_low_seq(2, 1:Index_high),'-.','lineWidth',3, 'Color', color_seq{2});
p3 = semilogx(I_seq, zeros(length(I_seq), 1),'-','lineWidth',3, 'Color', color_seq{1});
p4 = plot([I_seq(Index_low), I_seq(Index_low)], [-5 5], '--', 'lineWidth',1.5, 'Color', 'k');
p5 = scatter(I_seq(Index_high), x_high_seq(2, Index_high), 500, '.', 'k');

xticks([1e0, 1e3, 1e6])
yticks([0, 1, 2])
xlim([1 1e6])
ylim([-0.2 2.2])
xlabel('$[I]$','Interpreter','LaTex', 'Units', 'normalized')
ylabel('$[B]$','Interpreter','LaTex', 'Units', 'normalized')
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 16);  
grid on

%---------------------------------- Node A ---------------------------------------%
figure('unit','points','PaperUnits','points', 'position',[500 300 400 150]);
p1 = semilogx(I_seq(1:Index_high), x_high_seq(1, 1:Index_high),'-','lineWidth',3, 'Color',  color_seq{1});
hold on 
p2 = semilogx(I_seq(1:Index_low), x_low_seq(1, 1:Index_high),'-.','lineWidth',3, 'Color',  color_seq{2});
p3 = semilogx(I_seq, zeros(length(I_seq), 1),'-','lineWidth',3, 'Color',  color_seq{1});
p4 = plot([I_seq(Index_low), I_seq(Index_low)], [-5 5], '--', 'lineWidth',1.5, 'Color', 'k');
p5 = scatter(I_seq(Index_high), x_high_seq(1, Index_high), 500, '.', 'k');

xticks([1e0, 1e3, 1e6])
yticks([0, 0.3, 0.6])
xlim([1 1e6])
ylim([-0.05 0.6])
xlabel('$[I]$','Interpreter','LaTex', 'Units', 'normalized')
ylabel('$[A]$','Interpreter','LaTex', 'Units', 'normalized')
leg = legend([p1, p2, p5], 'Stable fixed point', 'Unstable fixed point','SN bifurcation point','FontSize',12,'FontName','Times New Roman','Location','northeast');
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 16);  
grid on

%% 
function F = circuit(x, I, Gamma_A, Gamma_B, Beta_A, Beta_B, n_BA, n_AB, n_BB, K_BA, K_AB, K_BB)
    A = x(1);
    B = x(2);
    F(1) = -Gamma_A*A + Beta_A*I*(B^n_BA/(B^n_BA+K_BA^n_BA));
    F(2) = -Gamma_B*B + Beta_B*(K_AB^n_AB/(A^n_AB+K_AB^n_AB))*(B^n_BB/(B^n_BB+K_BB^n_BB));
end