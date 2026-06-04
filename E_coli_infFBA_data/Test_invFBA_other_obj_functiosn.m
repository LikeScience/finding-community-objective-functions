%%
initCobraToolbox

%%
clear;

% load EcoliSimulatedData.mat

fileName = 'ec_iJO1366.mat';

load(fileName);


modelOri = readCbModel(fileName);

[m, n] = size(modelOri.S);


%% Forward model definition

prob = optimproblem('ObjectiveSense', 'minimize');

N = 3;

forward_model = modelOri;
forward_model = changeObjective(forward_model,'Ec_biomass_iJO1366_WT_53p95M');
% forward_model = changeObjective(forward_model,'ATPS4rpp'); forward_model.lb(7)=0.6;
% forward_model = changeObjective(forward_model, {'ATPS4rpp','Ec_biomass_iJO1366_WT_53p95M'},[0.005,0.995]);
% forward_model = changeObjective(forward_model,'EX_glc(e)'); forward_model.lb(7)=0.68; 
% [ecoli.rxns(contains(ecoli.rxnNames, 'atp', 'IgnoreCase', true)), ecoli.rxnNames(contains(ecoli.rxnNames, 'atp', 'IgnoreCase', true))]
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

forward_model.lb(164) = 0; % Glucose
forward_model.lb(82) = -8.5; % Chloride
forward_model.lb(242) = -8.5; % Na+
forward_model.lb(244) = -19; % Ammonium
forward_model.lb(206) = -184; % K+
forward_model.lb(291) = -9.5; % Sulfate
forward_model.lb(263) = -85; % Phosphate
forward_model.lb(233) = -2; % Magnesium
forward_model.lb(128) = -.01; % Fe3+
forward_model.lb(95) = -.001; % Cu2+
forward_model.lb(332) = -.015; % Zinc
forward_model.lb(74) = -0.03; % Calcium
%Aerobic or anaerobic
forward_model.lb(252) = -1000; % O2

reac_idx = [164,174,293];
reac_lim = [7,12,25];
% reac_lim = [14,24,500];

for i = 1:N
% for i = 1:length(reac_idx)
    forward_model.lb(reac_idx(i)) = -reac_lim(i);
    sol = optimizeCbModel(forward_model, 'max', 'one');
    flux(i) = sol.f;
    v2(:, i) = sol.x;
    ub(:, i) = forward_model.ub;
    lb(:, i) = forward_model.lb;
    forward_model.lb(reac_idx(i)) = 0;
    
end

% stdev = 1;
% v2 = (optimizeCbModel(forward_model,'max','one').x).* (1 + stdev * randn(n, N));

S  = forward_model.S;

% clear reac_idx reac_lim m n fileName prob forward_model modelOri;

model = forward_model;

%% Finding Objective Function

for i = 1:N
    [objL1{i}, f1{i}, objL0{i}, f0{i}, objBase{i}] = invFBA(v2(:,i), lb(:,i), ub(:,i), S, true);
end
%% Display objective functions
for i = 1:N
    sprintf("Limiting reaction: %s", model.rxns{reac_idx(i)});
    disp("After L1:");
    idx = objL1{i} > 1e-4;
    disp([model.rxns(idx), num2cell(objL1{i}(idx))])
    disp("After L0:");
    idx = objL0{i} > 1e-4;
    disp([model.rxns(idx), num2cell(objL0{i}(idx))])
end


%% Forward plots

substrate_names = {'glucose', 'glycerol', 'succinate'}

for j = 1:N
    t_model = forward_model
    t_model.lb(reac_idx(j)) = -reac_lim(j);
    label = 'E. coli, objective biomass'
    short_label = 'E.coli_biomass'
    
    objs = [objL1{j}, objL0{j}, objBase{j}];
    labels = {'L1','L0','Base'};
    for i = 1:3
        fModel = changeObjective(t_model, t_model.rxns, objs(:,i));
        sol_f = optimizeCbModel(fModel,'max','one');
        sol_f.f
        fig = figure('Visible', 'off'); 
        scatter(sol.x, sol_f.x,'blue');
        xlabel('Fluxes with original objective');
        ylabel('Fluxes with recovered objective');
        current_label = labels{i};
        title(sprintf('%s (%s)', label, current_label));
        filename = sprintf('flux_plots/fluxes_original_objective_vs_invFBA_%s_%s_%s.png', current_label, short_label, substrate_names{j});
        saveas(fig, filename);
        close(fig);
    end
    t_model.lb(reac_idx(j)) = 0;

end



%% Try with constraint to make L1 norm 1

for i = 1:N
    [objBase{i}, sol_c_pos{i},sol_c_neg{i}] = invFBAmod(v2(:,i), lb(:,i), ub(:,i), S, true);
end
for i = 1:N
    sprintf("Limiting reaction: %s", model.rxns{reac_idx(i)});
    disp("Base:");
    idx = abs(objBase{i}) > 1e-4;
    disp([model.rxns(idx), num2cell(objBase{i}(idx))])
end