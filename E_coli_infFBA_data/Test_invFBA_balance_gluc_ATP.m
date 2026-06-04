%%
initCobraToolbox
%%
clear;
fileName = 'ec_iJO1366.mat';
load(fileName);
modelOri = readCbModel(fileName);
[m, n] = size(modelOri.S);

% objective1 = 'ATPS4rpp';
objective1 = 'Ec_biomass_iJO1366_WT_53p95M';
objective2 = 'EX_glc(e)';
% obj_tag1 = ' ATP Synthase';
obj_tag2 = 'Glucose exchange';
obj_tag1 = 'Biomass';

%% Forward model definition
N = 1;
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
reac_idx = [164];
reac_lim = [7];
carbon_sources = {'Glucose'};

alphas = 0.9:0.001:0.95;
c_L1_atp = zeros(length(alphas), N);
c_L1_bio = zeros(length(alphas), N);
c_L0_atp = zeros(length(alphas), N);
c_L0_bio = zeros(length(alphas), N);

S = forward_model.S;
idx_atp = findRxnIDs(forward_model, objective1);
idx_bio = findRxnIDs(forward_model, objective2);

for i = 1:N
    forward_model.lb(reac_idx(i)) = -reac_lim(i);
    curr_lb = forward_model.lb;
    curr_ub = forward_model.ub;
    
    for k = 1:length(alphas)
        disp(k)
        temp_model = changeObjective(forward_model, {objective1, objective2}, [alphas(k), 1-alphas(k)]);
        sol_a = optimizeCbModel(temp_model, 'max', 'one');
        
        if sol_a.stat == 1
            [cL1, ~, cL0, ~] = invFBA(sol_a.v, curr_lb, curr_ub, S, true);
            
            c_L1_atp(k, i) = cL1(idx_atp);
            c_L1_bio(k, i) = cL1(idx_bio);
            c_L0_atp(k, i) = cL0(idx_atp);
            c_L0_bio(k, i) = cL0(idx_bio);
        end
    end
    forward_model.lb(reac_idx(i)) = 0;
end

figure;
for i = 1:N
    subplot(3, 2, 2*i - 1);
    plot(alphas, c_L1_atp(:, i), 'b.', alphas, c_L1_bio(:, i), 'r.');
    title(sprintf('L1 Norm: %s (Rxn %d)', carbon_sources{i}, reac_idx(i)));
    xlabel('Alpha');
    ylabel('invFBA Coefficient');
    legend(obj_tag1, obj_tag2);
    
    subplot(3, 2, 2*i);
    plot(alphas, c_L0_atp(:, i), 'b.', alphas, c_L0_bio(:, i), 'r.');
    title(sprintf('L0 Norm: %s (Rxn %d)', carbon_sources{i}, reac_idx(i)));
    xlabel('Alpha');
    ylabel('invFBA Coefficient');
    legend(obj_tag1, obj_tag2);
end