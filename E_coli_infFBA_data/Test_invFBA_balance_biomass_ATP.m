%%
initCobraToolbox
%% Model definition
clear;
fileName = 'ec_iJO1366.mat';
load(fileName);
modelOri = readCbModel(fileName);
[m, n] = size(modelOri.S);

objective1 = 'ATPS4rpp';
objective2 = 'Ec_biomass_iJO1366_WT_53p95M';
obj_tag1 = 'ATP Synthesis';
obj_tag2 = 'Biomass';



%% Forward FBA predictions
N = 3;
forward_model = modelOri;
forward_model.ub(forward_model.ub == Inf | forward_model.ub > 10000) = 1000;
forward_model.lb(forward_model.lb == -Inf | forward_model.lb < -10000) = -1000;
forward_model.lb(9:332) = 0;
forward_model.lb([
    76
    85
    86
   127
   185
   187
   237
   239
   245
   252
   284
   288
   313]) = -.001;
forward_model.lb(164) = 0;
forward_model.lb(82) = -8.5;
forward_model.lb(242) = -8.5;
forward_model.lb(244) = -19;
forward_model.lb(206) = -184;
forward_model.lb(291) = -9.5;
forward_model.lb(263) = -85;
forward_model.lb(233) = -2;
forward_model.lb(128) = -.01;
forward_model.lb(95) = -.001;
forward_model.lb(332) = -.015;
forward_model.lb(74) = -0.03;
forward_model.lb(252) = -1000;

forward_model.lb(7) = 0.1;

reac_idx = [164,174,293];
reac_lim = [7,12,25];
carbon_sources = {'Glucose', 'Glycerol', 'Succinate'};

alphas = 0:0.05:1;
c_L1_atp = zeros(length(alphas), N);
c_L1_bio = zeros(length(alphas), N);
c_L0_atp = zeros(length(alphas), N);
c_L0_bio = zeros(length(alphas), N);

S = forward_model.S;
idx_atp = findRxnIDs(forward_model, objective1);
idx_bio = findRxnIDs(forward_model, objective2);

idx_other = setdiff(1:n, [idx_atp, idx_bio]);
zeros_vec = zeros(size(idx_other));

n_rep = 1;

for i = 1:N
    forward_model.lb(reac_idx(i)) = -reac_lim(i);
    curr_lb = forward_model.lb;
    curr_ub = forward_model.ub;
    
    for k = 1:length(alphas)
        disp(k)
        temp_model = changeObjective(forward_model, {objective1, objective2}, [alphas(k), 1-alphas(k)]);
        sol_a = optimizeCbModel(temp_model, 'max', 'one');
        % sol_a = enumerateOptimalSolutions(temp_model).fluxes(:,1:n_rep);
        sol_all{k, i} = sol_a;
        disp(k)

        if sol_a.stat == ones(n_rep)
            % [cL1, ~, cL0, ~, cLB] = invFBA(sol_a.v, curr_lb, curr_ub, S, true, false,zeros_vec,idx_other);
            [cL1, ~, cL0, ~, cLB] = invFBA(sol_a.v, curr_lb, curr_ub, S, true, false);
            
            c_L1_atp(k, i) = cL1(idx_atp);
            c_L1_bio(k, i) = cL1(idx_bio);
            c_L0_atp(k, i) = cL0(idx_atp);
            c_L0_bio(k, i) = cL0(idx_bio);
            c_LB_atp(k, i) = cLB(idx_atp);
            c_LB_bio(k, i) = cLB(idx_bio);
        end
        
    end
    forward_model.lb(reac_idx(i)) = 0;
end

flux_matrix = cell2mat(cellfun(@(s) s.x, sol_all(:, 1)', 'UniformOutput', false));
%% All figures 
figure;
for i = 1:N
    subplot(N, 3, 3*i - 2);
    plot(alphas, c_LB_atp(:, i), 'b.', alphas, c_LB_bio(:, i), 'r.');
    title(sprintf('Base invFBA: %s (Rxn %d)', carbon_sources{i}, reac_idx(i)));
    xlabel('Alpha');
    ylabel('invFBA Coefficient');
    legend(obj_tag1, obj_tag2);

    subplot(N, 3, 3*i - 1);
    plot(alphas, c_L1_atp(:, i), 'b.', alphas, c_L1_bio(:, i), 'r.');
    title(sprintf('L1 Norm: %s (Rxn %d)', carbon_sources{i}, reac_idx(i)));
    xlabel('Alpha');
    ylabel('invFBA Coefficient');
    legend(obj_tag1, obj_tag2);
    
    subplot(N, 3, 3*i);
    plot(alphas, c_L0_atp(:, i), 'b.', alphas, c_L0_bio(:, i), 'r.');
    title(sprintf('L0 Norm: %s (Rxn %d)', carbon_sources{i}, reac_idx(i)));
    xlabel('Alpha');
    ylabel('invFBA Coefficient');
    legend(obj_tag1, obj_tag2);
end

%% Figures for poster
figure;
NN = 1 %Glucose
for i = 1:NN
    subplot(NN, 3, 3*i - 2);
    plot(alphas, c_LB_atp(:, i), 'b.', alphas, c_LB_bio(:, i), 'r.');
    title(sprintf('invFBA Base'));
    xlabel('α');
    ylabel('invFBA Coefficient');
    legend(obj_tag1, obj_tag2);

    subplot(NN, 3, 3*i - 1);
    plot(alphas, c_L1_atp(:, i), 'b.', alphas, c_L1_bio(:, i), 'r.');
    title(sprintf('invFBA + L1 Norm'));
    xlabel('α');
    ylabel('invFBA Coefficient');
    legend(obj_tag1, obj_tag2);
    
    subplot(NN, 3, 3*i);
    plot(alphas, c_L0_atp(:, i), 'b.', alphas, c_L0_bio(:, i), 'r.');
    title(sprintf('invFBA + L0 Norm'));
    xlabel('α');
    ylabel('invFBA Coefficient');
    legend(obj_tag1, obj_tag2);
end

%% Figures for poster
figure('Color', 'w', 'Position', [100, 100, 1730, 300]);
NN = 1 %Glucose
for i = 1:NN
    subplot(NN, 3, 3*i - 2);
    plot(alphas, c_LB_atp(:, i), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5); hold on;
    plot(alphas, c_LB_bio(:, i), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 4); hold off;
    title(sprintf('invFBA Base'), 'Color', 'k');
    xlabel('Original α');
    ylabel('invFBA α');
    set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'FontSize', 20, 'LineWidth', 1.5);
    set(gca, 'OuterPosition', [0.01, 0, 0.273, 1]);

    subplot(NN, 3, 3*i - 1);
    plot(alphas, c_L1_atp(:, i), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5); hold on;
    plot(alphas, c_L1_bio(:, i), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 4); hold off;
    title(sprintf('invFBA + L1 Norm'), 'Color', 'k');
    xlabel('Original α');
    % ylabel('invFBA Coefficient');
    set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'FontSize', 20, 'LineWidth', 1.5);
    set(gca, 'OuterPosition', [0.29, 0, 0.25, 1]);

    subplot(NN, 3, 3*i);
    plot(alphas, c_L0_atp(:, i), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5); hold on;
    plot(alphas, c_L0_bio(:, i), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 4); hold off;
    title(sprintf('invFBA + L0 Norm'), 'Color', 'k');
    xlabel('Original α');
    % ylabel('invFBA Coefficient');
    set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'FontSize', 20, 'LineWidth', 1.5);
    set(gca, 'OuterPosition', [0.55, 0, 0.25, 1]);

    lgd = legend(obj_tag1, obj_tag2, 'Color', 'w','FontSize', 20, 'TextColor', 'k');
    title(lgd, 'Reaction');
    set(lgd, 'Units', 'normalized');
    set(lgd, 'Position', [0.83, 0.4, 0.08, 0.2]);
    
    
end
