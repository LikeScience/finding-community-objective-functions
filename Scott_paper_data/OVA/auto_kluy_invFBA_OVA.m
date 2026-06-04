%clear all
%Cobra Toolbox v3.33
%MMT 2.0
initCobraToolbox

%% Load model
clear
load(['pairedModel_auto_kluy.mat']); %pairedModel_auto_kluy.mat is the name of the community model
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
pairedModel=changeRxnBounds(pairedModel, 'auto_BIOt', 0.021*0.22*0.4, 'u'); %Fixed flux through biomass reactions (growth rate* biomass of each species)
pairedModel=changeRxnBounds(pairedModel, 'kluy_BIOt',0.021*0.22*0.6, 'u');
save(['pairedModel_auto_kluy.mat'],'pairedModel') %This saves the community model with the added constraints


%% Prepare excel results exporting
exportData = cell(numel(pairedModel.rxns) + 2, 5);
exportData(2, 1:5) = {'Reaction_ID', 'Reaction_Name', 'Equation', 'Lower_Bound', 'Upper_Bound'};
exportData(3:end, 1:5) = [pairedModel.rxns, pairedModel.rxnNames, printRxnFormula(pairedModel, pairedModel.rxns, false), num2cell(pairedModel.lb), num2cell(pairedModel.ub)];
simIdx = 1; 

%% Run invFBA and save/output results

function [simIdx, exportData, objL1, objL0, objBase] = runinvFBA(fba_sol, pairedModel, zeroArr, targetIdx, titleStr, simIdx, exportData, terminal_output,min_objval,min_coef_val)
    if isempty(zeroArr)
        [objL1, ~, objL0, ~, objBase, ~, epsilon] = invFBA(fba_sol, pairedModel.lb, pairedModel.ub, pairedModel.S, true, false,[],[],min_objval,min_coef_val);
    else
        [objL1, ~, objL0, ~, objBase, ~, epsilon] = invFBA(fba_sol, pairedModel.lb, pairedModel.ub, pairedModel.S, true, false, zeroArr, targetIdx,min_objval,min_coef_val);
    end

    if terminal_output
        fprintf("Epsilon: %g\n", epsilon);
        disp("After L1:");
        idx = abs(objL1) > 0;
        disp([pairedModel.rxns(idx), num2cell(objL1(idx)), num2cell(fba_sol(idx))]);
        disp("After L0:");
        idx = abs(objL0) > 0;
        disp([pairedModel.rxns(idx), num2cell(objL0(idx)), num2cell(fba_sol(idx))]);
    end

    simIdx = simIdx + 1;
    c = 7 + (simIdx - 1) * 8;
    exportData(1, c) = {titleStr};
    exportData(1, c+5) = {'Epsilon:'};
    exportData(1, c+6) = {epsilon};
    exportData(2, c:c+6) = {'Forward_FBA_fluxes', 'invFBA_obj_before_L1', 'invFBA_obj_after_L1', 'invFBA_obj_after_L1_L0', 'Fwd_FBA_Base', 'Fwd_FBA_L1', 'Fwd_FBA_L0'};
    exportData(3:end, c:c+3) = num2cell([fba_sol(:), objBase(:), objL1(:), objL0(:)]);
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
pairedModel = changeObjective(pairedModel, {'auto_BIOt', 'kluy_BIOt'}, [0.5, 0.5]);
sol = optimizeCbModel(pairedModel,"max","one")
selected_samples=sol.x;

if size(selected_samples,1) == numel(pairedModel.rxns)
    targetIdx = find(~(abs(selected_samples) >= 1e-4));
else
    targetIdx = [];
end
zeroArr = zeros(size(targetIdx));

label = 'Joint model, biomass objective (0.5,0.5)';
short_label = 'Joint_biom_0.5_0.5';

min_coef_val = -1000; % (No restriction)
% min_coef_val = 0.01;
[simIdx, exportData, objL1, objL0, objBase] = runinvFBA(selected_samples, pairedModel, zeroArr, targetIdx, label, simIdx, exportData, true,0,min_coef_val);

objs = [objL1, objL0, objBase];
labels = {'L1','L0','Base'};
% generateFluxPlots(pairedModel, objs, sol, labels, label, short_label)


%% Select reactions with nonzero flux
targetIdx = find(~(abs(selected_samples) >= 1e-4));
acceptableIdx = find((abs(selected_samples) >= 1e-4));
min_objval = 1e-2;
%% OVA
[C_lb, C_ub, objBase, flagBase, sum_eps,C_lb_full, C_ub_full] = OVA(sol.x, pairedModel.lb, pairedModel.ub, pairedModel.S, acceptableIdx,zeroArr, targetIdx,min_objval);
%% Save
save('OVA_AK_target1e-4_minobj1e-2_results.mat', 'C_lb', 'C_ub', 'objBase', 'flagBase', 'sum_eps', 'C_lb_full', 'C_ub_full');
%% Load OVA results
load('OVA_AK_target1e-4_minobj1e-2_results.mat')
%%
threshold = 0.5; 
meaningful_rxns = pairedModel.rxns(acceptableIdx(C_ub > threshold));
non_meaningful_idx = setdiff(1:numel(pairedModel.rxns), acceptableIdx(C_ub > threshold));
zeroArr_2 = zeros(size(non_meaningful_idx));

label = 'Joint model, biomass objective (0.5,0.5) after OVA'
short_label = 'Joint_biom_0.5_0.5_OVA'
% label = 'Joint model, biomass objective (0.5,0.5) after OVA min coefs 0.01'
% short_label = 'Joint_biom_0.5_0.5_OVA_coef0.01'

% min_coef_val = -1000; % (No restriction)
min_coef_val = 0.01;
[simIdx, exportData, objL1, objL0, objBase] = runinvFBA(selected_samples, pairedModel, zeroArr_2, non_meaningful_idx, label, simIdx, exportData, true,min_objval,min_coef_val);


%% Export excel
% Note for people cloning this repo from Github: this code block frequently leads to errors if Excel is open. Don't
% forget to save all you open excel files before running this block.
system('taskkill /F /IM excel.exe');

excelFileName = 'auto_kluy_invFBA_OVA_results_poster.xlsx';

writecell(exportData, excelFileName);
pause(1); % Sync I/O buffer
x = actxserver('Excel.Application');
x.DisplayAlerts = false;
w = x.Workbooks.Open(fullfile(pwd, excelFileName));

w.Sheets.Item(1).ListObjects.Add(1, w.Sheets.Item(1).Range(['A2:', w.Sheets.Item(1).UsedRange.SpecialCells(11).Address]), [], 1);
% w.Sheets.Item(1).ListObjects.Add(1, w.Sheets.Item(1).UsedRange, [], 1);
w.Sheets.Item(1).Columns.ColumnWidth = 10;

w.Save;
w.Close;
x.Quit;
x.delete;


