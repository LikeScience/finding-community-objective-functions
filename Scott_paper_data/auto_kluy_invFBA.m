%clear all
%Cobra Toolbox v3.33
%MMT 2.0
initCobraToolbox

%% Create community model
modelFolder=[pwd filesep 'model_folder']; 
modelList = {'auto';'kluy'};
joinModelsPairwiseFromList({'auto';'kluy'},modelFolder,'pairwiseModelFolder',modelFolder,'biomasses',{'BIOt','BIOt'}); %This function generates the community model and the pairedModelinfo file

%% Load model
clear
load([pwd filesep 'model_folder' filesep 'pairedModel_auto_kluy.mat']); %pairedModel_auto_kluy.mat is the name of the community model
% changeCobraSolver('glpk', 'LP');
changeCobraSolver('gurobi', 'LP');

%% Assign bounds to the community model
Data= readtable('Bounds_MMT.xlsx'); %Assign bounds to the reactions defined in [u] compartment in the community model 
lb_com=Data{1:2419,2};
ub_com=Data{1:2419,3};
for w = 1:2419
    pairedModel = changeRxnBounds(pairedModel,pairedModel.rxns(w), lb_com(w,1), 'l');
    pairedModel = changeRxnBounds(pairedModel,pairedModel.rxns(w), ub_com(w,1), 'u'); 
end

%% Input data and constraints

pairedModel=changeRxnBounds(pairedModel, 'EX_CO[u]', -4.8, 'l'); %Define CO uptake rate (mmol/h)
% pairedModel=changeRxnBounds(pairedModel, 'EX_CO[u]', -0.5, 'l'); %Define CO uptake rate (mmol/h)
% pairedModel=changeRxnBounds(pairedModel, 'EX_CO[u]', -0.073, 'l'); %Define CO uptake rate (mmol/h)
pairedModel=changeRxnBounds(pairedModel, 'auto_BIOt', 0.021*0.22*0.4, 'u'); %Fixed flux through biomass reactions (growth rate* biomass of each species)
pairedModel=changeRxnBounds(pairedModel, 'kluy_BIOt',0.021*0.22*0.6, 'u');
% pairedModel=changeRxnBounds(pairedModel, 'auto_BIOt',0.02,'u')
% pairedModel=changeRxnBounds(pairedModel, 'kluy_BIOt',0.5,'u')
%pairedModel=changeRxnBounds(pairedModel, 'auto_atpm', 8.4*0.22*0.4, 'l');%ATP maintenance flux scaled to mmol/h. Infeasible
%pairedModel=changeRxnBounds(pairedModel, 'kluy_Rckl725',0.45*0.22*0.6, 'l');
save([pwd filesep 'model_folder' filesep  'pairedModel_auto_kluy.mat'],'pairedModel') %This saves the community model with the added constraints

%% Run simulatePairwiseInteractions

modPath=[pwd filesep 'model_folder'];
[pairwiseInteractions, pairwiseSolutions] = simulatePairwiseInteractions(modPath,'saveSolutionsFlag',true); %Outputs the fluxes and growth rates


%% Export results to excel header
exportData = cell(numel(pairedModel.rxns) + 2, 5);
exportData(2, 1:5) = {'Reaction_ID', 'Reaction_Name', 'Equation', 'Lower_Bound', 'Upper_Bound'};
exportData(3:end, 1:5) = [pairedModel.rxns, pairedModel.rxnNames, printRxnFormula(pairedModel, pairedModel.rxns, false), num2cell(pairedModel.lb), num2cell(pairedModel.ub)];
simIdx = 1; 

%% Run invFBA and save/output results

function [simIdx, exportData, objL1, objL0, objBase] = runinvFBA(fba_sol, pairedModel, zeroArr, targetIdx, titleStr, simIdx, exportData, terminal_output)
    if isempty(zeroArr)
        [objL1, ~, objL0, ~, objBase, ~, epsilon] = invFBA(fba_sol.v, pairedModel.lb, pairedModel.ub, pairedModel.S, true, false);
    else
        [objL1, ~, objL0, ~, objBase, ~, epsilon] = invFBA(fba_sol.v, pairedModel.lb, pairedModel.ub, pairedModel.S, true, false, zeroArr, targetIdx);
    end

    if terminal_output
        fprintf("Epsilon: %g\n", epsilon);
        disp("After L1:");
        idx = abs(objL1) > 0;
        disp([pairedModel.rxns(idx), num2cell(objL1(idx)), num2cell(fba_sol.x(idx))]);
        disp("After L0:");
        idx = abs(objL0) > 0;
        disp([pairedModel.rxns(idx), num2cell(objL0(idx)), num2cell(fba_sol.x(idx))]);
    end

    simIdx = simIdx + 1;
    c = 7 + (simIdx - 1) * 8;
    exportData(1, c) = {titleStr};
    exportData(1, c+5) = {'Epsilon:'};
    exportData(1, c+6) = {epsilon};
    exportData(2, c:c+6) = {'Forward_FBA_fluxes', 'invFBA_obj_before_L1', 'invFBA_obj_after_L1', 'invFBA_obj_after_L1_L0', 'Fwd_FBA_Base', 'Fwd_FBA_L1', 'Fwd_FBA_L0'};
    exportData(3:end, c:c+3) = num2cell([fba_sol.v(:), objBase(:), objL1(:), objL0(:)]);
end

%% Compare obtained fluxes/objs with fluxes/objs from single reaction/metabolite objective 
function [allFluxes, A] = computeFBAs(model, loadFile, filename)
    if nargin > 1 && loadFile && isfile(filename)
        load(filename, 'allFluxes', 'A');
        return;
    end
    A = [speye(size(model.S, 2)), -speye(size(model.S, 2)), (model.S ./ sum(abs(model.S), 2))', -(model.S ./ sum(abs(model.S), 2))'];
    allFluxes = NaN(size(model.S, 2), size(A, 2));
    
    for i = 1:size(A, 2)
        fprintf('Iteration %d / %d\n', i, size(A, 2));
        
        model.c = A(:, i);
        sol = optimizeCbModel(model, 'max', 'one');
        
        if sol.stat == 1 && ~isempty(sol.x)
            allFluxes(:, i) = sol.x;
        end
    end
    save(filename, 'allFluxes', 'A');
end

function [maxDist, matches] = compareFluxes(model, refObjs, allFluxes, A, tol)
    numRefs = size(refObjs, 2);
    numFBAs = size(allFluxes, 2);
    refFluxes = NaN(size(model.S, 2), numRefs);
    
    for j = 1:numRefs
        model.c = refObjs(:, j);
        sol = optimizeCbModel(model, 'max', 'one');
        if sol.stat == 1 && ~isempty(sol.x)
            refFluxes(:, j) = sol.x;
        end
    end
    
    maxDist = NaN(numFBAs, numRefs);
    matches = cell(1, numRefs);
    
    for j = 1:numRefs
        maxDist(:, j) = max(abs(allFluxes - refFluxes(:, j)), [], 1)';
        match_idx = find(maxDist(:, j) < tol);
        
        matches{j} = cell(length(match_idx), 1);
        for k = 1:length(match_idx)
            matches{j}{k} = model.rxns(A(:, match_idx(k)) ~= 0);
        end
    end
end

function [diffVals, matches] = compareObjectives(model, refObjs, allFluxes, A, tol, ignoreZeroObj)
    if nargin < 6, ignoreZeroObj = false; end
    numRefs = size(refObjs, 2);
    numObjs = size(A, 2);
    numRxns = size(model.S, 2);
    numMets = size(model.S, 1);
    
    refFluxes = NaN(numRxns, numRefs);
    for j = 1:numRefs
        model.c = refObjs(:, j);
        sol = optimizeCbModel(model, 'max', 'one');
        if sol.stat == 1 && ~isempty(sol.x)
            refFluxes(:, j) = sol.x;
        end
    end
    
    diffVals = NaN(numRefs, numObjs);
    matches = cell(1, numRefs);
    allObjVals = sum(A .* allFluxes, 1);
    
    for j = 1:numRefs
        diffVals(j, :) = allObjVals - (A' * refFluxes(:, j))';
        match_idx = find(abs(diffVals(j, :)) < tol & (~ignoreZeroObj | allObjVals ~= 0));
        
        matches{j} = cell(length(match_idx), 1);
        for k = 1:length(match_idx)
            idx = match_idx(k);
            if idx <= 2 * numRxns
                matches{j}{k} = model.rxns{mod(idx - 1, numRxns) + 1};
            else
                matches{j}{k} = model.mets{mod(idx - 2 * numRxns - 1, numMets) + 1};
            end
        end
    end
end

%Function to generate presentable flux plots
function generateFluxPlots(pairedModel, objs, sol, labels, label, short_label)
    figWidth = 12; 
    figHeight = 12;
    for i = 1:3
        fModel = changeObjective(pairedModel, pairedModel.rxns, objs(:,i));
        sol_f = optimizeCbModel(fModel, 'max', 'one');
        
        fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', 'Position', [1, 1, figWidth, figHeight]);
        scatter(sol.x, sol_f.x, 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'none');
        hold on;
        
        ax = gca;
        min_v = min([ax.XLim(1), ax.YLim(1)]);
        max_v = max([ax.XLim(2), ax.YLim(2)]);
        plot([min_v, max_v], [min_v, max_v], 'k--', 'LineWidth', 1);
        axis(ax, 'square');
        xlim(ax, [min_v max_v]);
        ylim(ax, [min_v max_v]);
        hold off;
        
        set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'FontSize', 18, 'Box', 'on');
        xlabel('Fluxes with original objective', 'Color', 'k', 'FontSize', 18, 'FontWeight', 'bold');
        ylabel('Fluxes with recovered objective', 'Color', 'k', 'FontSize', 18, 'FontWeight', 'bold');
        % title(sprintf('%s (%s)', label, labels{i}), 'Color', 'k', 'FontSize', 16, 'FontWeight', 'bold');
        
        rmsd = sqrt(mean((sol.x - sol_f.x).^2));
        text(0.05, 0.9, sprintf('RMSD = %.4g', rmsd), 'Units', 'normalized', 'Color', 'k', 'BackgroundColor', 'w', 'EdgeColor', 'k', 'FontSize', 18); % Background box prevents point overlap
        
        exportgraphics(fig, sprintf('plots/fluxes_original_objective_vs_invFBA_%s_%s.png', labels{i}, short_label), 'Resolution', 300);
        close(fig);
    end
end

%% Run invFBA - joint, biomass (0.5, 0.5)
changeCobraSolverParams('LP', 'feasTol', 1e-5);
changeCobraSolverParams('LP', 'optTol', 1e-5);
% Objective

pairedModel = changeObjective(pairedModel, {'auto_BIOt', 'kluy_BIOt'}, [0.5, 0.5]);
disp([num2cell(find(pairedModel.c(:))), pairedModel.rxns(find(pairedModel.c(:))), num2cell(pairedModel.c(find(pairedModel.c(:))))])

sol = optimizeCbModel(pairedModel, 'max', 'one');

fprintf('Flux: %g Biomass Auto: %g Biomass kluy: %g \n', sol.f, sol.x(strcmp(pairedModel.rxns, 'auto_BIOt')), sol.x(strcmp(pairedModel.rxns, 'kluy_BIOt')));

if numel(sol.x) == numel(pairedModel.rxns)
    targetIdx = find(~((startsWith(pairedModel.rxns, 'EX') | startsWith(pairedModel.rxns, 'auto_IEX') | startsWith(pairedModel.rxns, 'kluy_IEX') | startsWith(pairedModel.rxns, 'auto_BIO') | startsWith(pairedModel.rxns, 'kluy_BIO')) & (abs(sol.x) >= 1e-4)));
else
    targetIdx = [];
end
zeroArr = zeros(size(targetIdx));

label = 'Joint model, biomass objective (0.5,0.5) only exchanges/biomass, zeoring';
short_label = 'Joint_biom_0.5_0.5_ez';

[simIdx, exportData, objL1, objL0, objBase] = runinvFBA(sol, pairedModel, zeroArr, targetIdx, label, simIdx, exportData, true);

% [objBM] = invFBAmod(sol.v, pairedModel.lb, pairedModel.ub, pairedModel.S, true, false, zeroArr, targetIdx);

objs = [objL1, objL0, objBase];
labels = {'L1','L0','Base'};

generateFluxPlots(pairedModel, objs, sol, labels, label, short_label)

% for i = 1:3
%     fModel = changeObjective(pairedModel, pairedModel.rxns, objs(:,i));
%     sol_f = optimizeCbModel(fModel,'max','one');
%     fig = figure('Visible', 'off'); 
%     scatter(sol.x, sol_f.x,'blue');
%     xlabel('Fluxes with original objective');
%     ylabel('Fluxes with recovered objective');
%     current_label = labels{i};
%     title(sprintf('%s (%s)', label, current_label));
%     filename = sprintf('plots/fluxes_original_objective_vs_invFBA_%s_%s.png', current_label, short_label);
%     saveas(fig, filename);
%     close(fig);
% end
%% OVA
r_idx = [961,2071];
[C_lb, C_ub, objBase, flagBase, sum_eps,C_lb_full, C_ub_full] = OVA(sol.v, pairedModel.lb, pairedModel.ub, pairedModel.S, r_idx,zeroArr, targetIdx);
%%
[allFluxes, A] = computeFBAs(pairedModel, true, 'precomputed_fluxes.mat');
[diffVals, objMatches] = compareObjectives(pairedModel, objs, allFluxes, A, 1e-5,true);
% [maxDist, matches] = compareFBAs(pairedModel, objs, allFluxes, A, 1e-5)

%% Run invFBA - joint, biomass (0.5, 0.5), no zeroing of reactions with small flux

% Objective

pairedModel = changeObjective(pairedModel, {'auto_BIOt', 'kluy_BIOt'}, [0.5, 0.5]);
disp([num2cell(find(pairedModel.c(:))), pairedModel.rxns(find(pairedModel.c(:))), num2cell(pairedModel.c(find(pairedModel.c(:))))])

sol = optimizeCbModel(pairedModel, 'max', 'one');

fprintf('Flux: %g Biomass Auto: %g Biomass kluy: %g \n', sol.f, sol.x(strcmp(pairedModel.rxns, 'auto_BIOt')), sol.x(strcmp(pairedModel.rxns, 'kluy_BIOt')));

if numel(sol.x) == numel(pairedModel.rxns)
    targetIdx = find(~(startsWith(pairedModel.rxns, 'EX') | startsWith(pairedModel.rxns, 'auto_IEX') | startsWith(pairedModel.rxns, 'kluy_IEX') | startsWith(pairedModel.rxns, 'auto_BIO') | startsWith(pairedModel.rxns, 'kluy_BIO') ));
else
    targetIdx = [];
end
zeroArr = zeros(size(targetIdx));

label = 'Joint model, biomass objective (0.5,0.5) only exchanges/biomass'
short_label = 'Joint_biom_0.5_0.5_e'

[simIdx, exportData,objL1, objL0, objBase] = runinvFBA(sol, pairedModel, zeroArr, targetIdx, label, simIdx, exportData, true);
objs = [objL1, objL0, objBase];
labels = {'L1','L0','Base'};

generateFluxPlots(pairedModel, objs, sol, labels, label, short_label)

%% Run invFBA - joint, biomass (0.5, 0.5)

% Objective

pairedModel = changeObjective(pairedModel, {'auto_BIOt', 'kluy_BIOt'}, [0.5, 0.5]);
disp([num2cell(find(pairedModel.c(:))), pairedModel.rxns(find(pairedModel.c(:))), num2cell(pairedModel.c(find(pairedModel.c(:))))])

sol = optimizeCbModel(pairedModel, 'max', 'one');

fprintf('Flux: %g Biomass Auto: %g Biomass kluy: %g \n', sol.f, sol.x(strcmp(pairedModel.rxns, 'auto_BIOt')), sol.x(strcmp(pairedModel.rxns, 'kluy_BIOt')));

if numel(sol.x) == numel(pairedModel.rxns)
    targetIdx = find(~(abs(sol.x) >= 1e-4));
else
    targetIdx = [];
end
zeroArr = zeros(size(targetIdx));

label = 'Joint model, biomass objective (0.5,0.5), zeroing'
short_label = 'Joint_biom_0.5_0.5_z'

[simIdx, exportData,objL1, objL0, objBase] = runinvFBA(sol, pairedModel, zeroArr, targetIdx, label, simIdx, exportData, true);

objs = [objL1, objL0, objBase];
labels = {'L1','L0','Base'};
generateFluxPlots(pairedModel, objs, sol, labels, label, short_label)

%% Run invFBA - joint, biomass (0.5, 0.5)

% Objective

pairedModel = changeObjective(pairedModel, {'auto_BIOt', 'kluy_BIOt'}, [0.5, 0.5]);
disp([num2cell(find(pairedModel.c(:))), pairedModel.rxns(find(pairedModel.c(:))), num2cell(pairedModel.c(find(pairedModel.c(:))))])

sol = optimizeCbModel(pairedModel, 'max', 'one');

fprintf('Flux: %g Biomass Auto: %g Biomass kluy: %g \n', sol.f, sol.x(strcmp(pairedModel.rxns, 'auto_BIOt')), sol.x(strcmp(pairedModel.rxns, 'kluy_BIOt')));

label = 'Joint model, biomass objective (0.5,0.5)'
short_label = 'Joint_biom_0.5_0.5'

[simIdx, exportData, objL1, objL0, objBase] = runinvFBA(sol, pairedModel, [], [], label, simIdx, exportData, true);

objs = [objL1, objL0, objBase];
labels = {'L1','L0','Base'};

generateFluxPlots(pairedModel, objs, sol, labels, label, short_label)

%% Run invFBA - joint, biomass (1, 0)

% Objective

pairedModel = changeObjective(pairedModel, {'auto_BIOt', 'kluy_BIOt'}, [1, 0]);
disp([num2cell(find(pairedModel.c(:))), pairedModel.rxns(find(pairedModel.c(:))), num2cell(pairedModel.c(find(pairedModel.c(:))))])

sol = optimizeCbModel(pairedModel, 'max', 'one');

fprintf('Flux: %g Biomass Auto: %g Biomass kluy: %g \n', sol.f, sol.x(strcmp(pairedModel.rxns, 'auto_BIOt')), sol.x(strcmp(pairedModel.rxns, 'kluy_BIOt')));

if numel(sol.x) == numel(pairedModel.rxns)
    targetIdx = find(~((startsWith(pairedModel.rxns, 'EX') | startsWith(pairedModel.rxns, 'auto_IEX') | startsWith(pairedModel.rxns, 'kluy_IEX') | startsWith(pairedModel.rxns, 'auto_BIO') | startsWith(pairedModel.rxns, 'kluy_BIO')) & (abs(sol.x) >= 1e-8)));
else
    targetIdx = [];
end
zeroArr = zeros(size(targetIdx));


[simIdx, exportData, objL1, objL0, objBase] = runinvFBA(sol, pairedModel, zeroArr, targetIdx, 'Joint, biomass, (1,0) only exchanges/biomass, zeoring', simIdx, exportData, true);

objs = [objL1, objL0, objBase];

for i = 1:3
    fModel = changeObjective(pairedModel, pairedModel.rxns, objs(i));
    sol_f = optimizeCbModel(fModel,'max','one');
    plot = scatter(sol.x, sol_f.x)
end

%% Run invFBA - joint, biomass (0.7, 0.3)

% Objective

pairedModel = changeObjective(pairedModel, {'auto_BIOt', 'kluy_BIOt'}, [0.7, 0.3]);
disp([num2cell(find(pairedModel.c(:))), pairedModel.rxns(find(pairedModel.c(:))), num2cell(pairedModel.c(find(pairedModel.c(:))))])

sol = optimizeCbModel(pairedModel, 'max', 'one');

fprintf('Flux: %g Biomass Auto: %g Biomass kluy: %g \n', sol.f, sol.x(strcmp(pairedModel.rxns, 'auto_BIOt')), sol.x(strcmp(pairedModel.rxns, 'kluy_BIOt')));

if numel(sol.x) == numel(pairedModel.rxns)
    targetIdx = find(~((startsWith(pairedModel.rxns, 'EX') | startsWith(pairedModel.rxns, 'auto_IEX') | startsWith(pairedModel.rxns, 'kluy_IEX') | startsWith(pairedModel.rxns, 'auto_BIO') | startsWith(pairedModel.rxns, 'kluy_BIO')) & (abs(sol.x) >= 1e-8)));
else
    targetIdx = [];
end
zeroArr = zeros(size(targetIdx));

[simIdx, exportData] = runinvFBA(sol, pairedModel, zeroArr, targetIdx, 'Joint, biomass, (.7,.3) only exchanges/biomass, zeoring', simIdx, exportData, true);



%% Run invFBA - joint, biomass (0.3, 0.7)

% Objective

pairedModel = changeObjective(pairedModel, {'auto_BIOt', 'kluy_BIOt'}, [.3, .7]);
disp([num2cell(find(pairedModel.c(:))), pairedModel.rxns(find(pairedModel.c(:))), num2cell(pairedModel.c(find(pairedModel.c(:))))])

sol = optimizeCbModel(pairedModel, 'max', 'one');

fprintf('Flux: %g Biomass Auto: %g Biomass kluy: %g \n', sol.f, sol.x(strcmp(pairedModel.rxns, 'auto_BIOt')), sol.x(strcmp(pairedModel.rxns, 'kluy_BIOt')));

if numel(sol.x) == numel(pairedModel.rxns)
    targetIdx = find(~((startsWith(pairedModel.rxns, 'EX') | startsWith(pairedModel.rxns, 'auto_IEX') | startsWith(pairedModel.rxns, 'kluy_IEX') | startsWith(pairedModel.rxns, 'auto_BIO') | startsWith(pairedModel.rxns, 'kluy_BIO')) & (abs(sol.x) >= 1e-8)));
else
    targetIdx = [];
end
zeroArr = zeros(size(targetIdx));

[simIdx, exportData] = runinvFBA(sol, pairedModel, zeroArr, targetIdx, 'Joint, biomass, (0.3,0.7) only exchanges/biomass, zeoring', simIdx, exportData, true);
%% Run invFBA - joint, biomass (0, 1)

% Objective

pairedModel = changeObjective(pairedModel, {'auto_BIOt', 'kluy_BIOt'}, [0,1]);
disp([num2cell(find(pairedModel.c(:))), pairedModel.rxns(find(pairedModel.c(:))), num2cell(pairedModel.c(find(pairedModel.c(:))))])

sol = optimizeCbModel(pairedModel, 'max', 'one');

fprintf('Flux: %g Biomass Auto: %g Biomass kluy: %g \n', sol.f, sol.x(strcmp(pairedModel.rxns, 'auto_BIOt')), sol.x(strcmp(pairedModel.rxns, 'kluy_BIOt')));

if numel(sol.x) == numel(pairedModel.rxns)
    targetIdx = find(~((startsWith(pairedModel.rxns, 'EX') | startsWith(pairedModel.rxns, 'auto_IEX') | startsWith(pairedModel.rxns, 'kluy_IEX') | startsWith(pairedModel.rxns, 'auto_BIO') | startsWith(pairedModel.rxns, 'kluy_BIO')) & (abs(sol.x) >= 1e-8)));
else
    targetIdx = [];
end
zeroArr = zeros(size(targetIdx));

[simIdx, exportData] = runinvFBA(sol, pairedModel, zeroArr, targetIdx, 'Joint, biomass, (0,1) only exchanges/biomass, zeoring', simIdx, exportData, true);


%% Run invFBA - kluy alone (auto blocked to 0)
% Set auto rxn bounds to zero
kluyModel = pairedModel
autoIdx = startsWith(pairedModel.rxns, 'auto');
kluyModel.lb(autoIdx) = 0;
kluyModel.ub(autoIdx) = 0;

kluyModel = changeObjective(kluyModel, {'auto_BIOt', 'kluy_BIOt'}, [0,1]);
disp([num2cell(find(kluyModel.c(:))), kluyModel.rxns(find(kluyModel.c(:))), num2cell(kluyModel.c(find(kluyModel.c(:))))])
sol = optimizeCbModel(kluyModel, 'max', 'one');
fprintf('Flux: %g Biomass Auto: %g Biomass kluy: %g \n', sol.f, sol.x(strcmp(kluyModel.rxns, 'auto_BIOt')), sol.x(strcmp(kluyModel.rxns, 'kluy_BIOt')));
if numel(sol.x) == numel(kluyModel.rxns)
    targetIdx = find((~((startsWith(kluyModel.rxns, 'EX') | startsWith(kluyModel.rxns, 'kluy_IEX') | startsWith(kluyModel.rxns, 'kluy_BIO')) & (abs(sol.x) >= 1e-8))) | startsWith(kluyModel.rxns, 'auto'));
else
    targetIdx = [];
end
zeroArr = zeros(size(targetIdx));

[simIdx, exportData] = runinvFBA(sol, kluyModel, zeroArr, targetIdx, 'Kluy, biomass, only exchanges/biomass, zeoring', simIdx, exportData, true);

%% Run invFBA - Auto alone (kluy blocked to 0)
% Set auto rxn bounds to zero
autoModel = pairedModel
kluyIdx = startsWith(autoModel.rxns, 'kluy');

autoModel.lb(kluyIdx) = 0;
autoModel.ub(kluyIdx) = 0;

autoModel = changeObjective(autoModel, {'auto_BIOt', 'kluy_BIOt'}, [1,0]);
disp([num2cell(find(autoModel.c(:))), autoModel.rxns(find(autoModel.c(:))), num2cell(autoModel.c(find(autoModel.c(:))))])
sol = optimizeCbModel(autoModel, 'max', 'one');
fprintf('Flux: %g Biomass Auto: %g Biomass kluy: %g \n', sol.f, sol.x(strcmp(autoModel.rxns, 'auto_BIOt')), sol.x(strcmp(autoModel.rxns, 'kluy_BIOt')));
targetIdx = cell(1, numel(autoModel));
zeroArr = cell(1, numel(autoModel));
targetIdx = find(abs(sol.x) <= 0);
if numel(sol.x) == numel(autoModel.rxns)
    targetIdx = find((~((startsWith(autoModel.rxns, 'EX') | startsWith(autoModel.rxns, 'auto_IEX') | startsWith(autoModel.rxns, 'auto_BIO')) & (abs(sol.x) >= 1e-8))) | startsWith(autoModel.rxns, 'kluy'));
else
    targetIdx = [];
end
zeroArr = zeros(size(targetIdx));
[simIdx, exportData] = runinvFBA(sol, pairedModel, zeroArr, targetIdx, 'Auto, biomass, only exchanges/biomass, zeoring', simIdx, exportData, true);

%% Export excel
% Note for people cloning this repo from Github: this code block frequently leads to errors if Excel is open. Don't
% forget to save all you open excel files before running this block.
system('taskkill /F /IM excel.exe');

excelFileName = 'auto_kluy_invFBA_results_auto_poster.xlsx';

writecell(exportData, excelFileName);
pause(1); % Sync I/O buffer
x = actxserver('Excel.Application');
w = x.Workbooks.Open(fullfile(pwd, excelFileName));

w.Sheets.Item(1).ListObjects.Add(1, w.Sheets.Item(1).Range(['A2:', w.Sheets.Item(1).UsedRange.SpecialCells(11).Address]), [], 1);
% w.Sheets.Item(1).ListObjects.Add(1, w.Sheets.Item(1).UsedRange, [], 1);
w.Sheets.Item(1).Columns.ColumnWidth = 10;

w.Save;
w.Close;
x.Quit;
x.delete;


%% Run TIOBjFInd


