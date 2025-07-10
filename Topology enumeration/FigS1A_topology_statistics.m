%% FigS1 parameters statistics of the three topologies
load('ParameterSet_1e4.mat');
load('Enumeration.mat');

% Reference values for comparison
n_ref = 2.5;
K_ref = 0;      % log10
Gamma_ref = -1;  % log10

% Define topology configurations
topologies = [10, 37, 64];

% Create figure window for each topology
for topo_idx = 1:length(topologies)
    i = topologies(topo_idx);
    Parameters_final = [];
    
    % Load parameters for current topology
    Temp = Result_final{i};
    for j = 1:length(Temp)
        Parameters_final = [Parameters_final; struct2array(Temp{j})];
    end
    
    % Calculate statistics
    Gamma_mean = mean(log(Parameters_final(:, 1:2))/log(10));
    Gamma_std = std(log(Parameters_final(:, 1:2))/log(10));
    n_mean = mean(Parameters_final(:, 7:10));
    n_std = std(Parameters_final(:, 7:10));
    K_mean = mean(log(Parameters_final(:, 3:6))/log(10));
    K_std = std(log(Parameters_final(:, 3:6))/log(10));

    % Create figure with three subplots
    figure('unit','points','position',[100 100 1200 350], 'Color','w')
    
    % ==========================
    % Hill Coefficients Subplot
    % ==========================
    subplot(1,3,1)
    if ismember(i, [10, 37])
        % Plot all 4 parameters for topologies 10 and 37
        bar_handle = bar(n_mean, 'FaceColor', [0.2 0.6 0.8]);
        hold on
        errorbar(1:4, n_mean, n_std, 'k', 'linestyle', 'none');
        label = {'$n_\mathrm{AA}$', '$n_\mathrm{AB}$', '$n_\mathrm{BA}$', '$n_\mathrm{BB}$'};
        xlim([0.5 4.5])
        xticks(1:4)
    else
        % Plot 3 parameters for topology 64
        bar_handle = bar(n_mean(2:4), 'FaceColor', [0.2 0.6 0.8]);
        hold on
        errorbar(1:3, n_mean(2:4), n_std(2:4), 'k', 'linestyle', 'none');
        label = {'$n_\mathrm{AB}$', '$n_\mathrm{BA}$', '$n_\mathrm{BB}$'};
        xlim([0.25 3.75])
        xticks(1:3)
    end
    
    line(xlim(), [n_ref, n_ref], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--');
    hold off
    
    yticks([0 2 4])
    xticklabels(label);
    set(gca, 'FontSize',12, 'TickLabelInterpreter','latex', 'Box','off', 'FontName', 'Arial');
    ylabel('Hill Coefficients')

    % ==========================
    % Affinity Constants Subplot
    % ==========================
    subplot(1,3,2)
    if ismember(i, [10, 37])
        % Plot all 4 parameters for topologies 10 and 37
        bar_handle = bar(K_mean, 'FaceColor', [0.8 0.4 0.2]);
        hold on
        errorbar(1:4, K_mean, K_std, 'k', 'linestyle', 'none');
        label = {'$\log_{10}K_\mathrm{AA}$', '$\log_{10}K_\mathrm{AB}$',...
                 '$\log_{10}K_\mathrm{BA}$', '$\log_{10}K_\mathrm{BB}$'};
        xlim([0.5 4.5])
        xticks(1:4)
    else
        % Plot 3 parameters for topology 64
        bar_handle = bar(K_mean(2:4), 'FaceColor', [0.8 0.4 0.2]);
        hold on
        errorbar(1:3, K_mean(2:4), K_std(2:4), 'k', 'linestyle', 'none');
        label = {'$\log_{10}K_\mathrm{AB}$', '$\log_{10}K_\mathrm{BA}$',...
                 '$\log_{10}K_\mathrm{BB}$'};
        xlim([0.25 3.75])
        xticks(1:3)
    end
    
    line(xlim(), [K_ref, K_ref], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--');
    hold off
    
    yticks([-2 0 2])
    ylim([-2.5 2.5])
    xticklabels(label);
    set(gca, 'FontSize',12, 'TickLabelInterpreter','latex', 'Box','off', 'FontName', 'Arial');
    ylabel('Affinity Constants')

    % ==========================
    % Degradation Rates Subplot
    % ==========================
    subplot(1,3,3)
    bar_handle = bar(Gamma_mean,  0.45, 'FaceColor', [0.4 0.7 0.4]);
    hold on
    line(xlim(), [Gamma_ref, Gamma_ref], 'Color', 'k', 'LineWidth',1,'LineStyle','--');
    errorbar(1:2, Gamma_mean, Gamma_std, 'k', 'linestyle', 'none');
    hold off
    
    xlim([0.5 2.5])
    yticks([-1 -.5 0 .5])
    ylim([-1.2 0.5])
    xticks(1:2)
    label = {'$\log_{10}\gamma_\mathrm{A}$', '$\log_{10}\gamma_\mathrm{B}$'};
    xticklabels(label);
    ylabel('Degradation Rates')
    set(gca, 'FontSize',12, 'TickLabelInterpreter','latex', 'Box','off', 'FontName', 'Arial');

    
    % Adjust subplot spacing
    set(gcf, 'Position', [100 100 1200 350]);
    print(strcat('FigS1_', int2str(i),'.pdf'), '-dpdf', '-bestfit');
end


