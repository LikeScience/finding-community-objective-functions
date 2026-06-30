%clear all
initCobraToolbox

%% Load model created in MIMECO
clear
pairedModel = readCbModel('pairedModel_Lacto_Akker_Bio_ATPM.mat')  %pairedModel_auto_kluy.mat is the name of the community model
% changeCobraSolver('glpk', 'LP');
changeCobraSolver('gurobi', 'LP');

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

%% Function to generate presentable flux plots
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


%% Run invFBA - joint, biomass/ATPM (0.5, 0.5)
changeCobraSolverParams('LP', 'feasTol', 1e-5);
changeCobraSolverParams('LP', 'optTol', 1e-5);
pairedModel = changeObjective(pairedModel, {'Growth:L_plantarum', 'ATPM:A_muciniphila'}, [0.5, 0.5]);
pairedModel.lb(2526)=0.1;
sol = optimizeCbModel(pairedModel,"max","one")
selected_samples=sol.x;

if size(selected_samples,1) == numel(pairedModel.rxns)
    targetIdx = find(~(abs(selected_samples) >= 1e-4));
else
    targetIdx = [];
end
zeroArr = zeros(size(targetIdx));

label = 'Joint model, ATPM_biomass objective (0.5,0.5)';
short_label = 'Joint_ATPM_biom_0.5_0.5';

% label = 'Joint model, ATPM_biomass objective (0.5,0.5)_c0.01';
% short_label = 'Joint_ATPM_biom_0.5_0.5_c0.01';

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
save('OVA_BioATPM_target1e-4_minobj1e-2_results.mat', 'C_lb', 'C_ub', 'objBase', 'flagBase', 'sum_eps', 'C_lb_full', 'C_ub_full');
%% Load OVA results
load('OVA_BioATPM_target1e-4_minobj1e-2_results.mat')
%%
threshold = 0.5; 
meaningful_rxns = pairedModel.rxns(acceptableIdx(C_ub > threshold));
non_meaningful_idx = setdiff(1:numel(pairedModel.rxns), acceptableIdx(C_ub > threshold));
zeroArr_2 = zeros(size(non_meaningful_idx));

% label = 'Joint model, ATP_biomass objective (0.5,0.5) after OVA '
% short_label = 'Joint_ATP_biom_0.5_0.5_OVA'

label = 'Joint model, biomass objective (0.5,0.5) after OVA min coefs 0.01'
short_label = 'Joint_biom_0.5_0.5_OVA_coef0.01'

% min_coef_val = -1000; % (No restriction)
min_coef_val = 0.01;
[simIdx, exportData, objL1, objL0, objBase] = runinvFBA(selected_samples, pairedModel, zeroArr_2, non_meaningful_idx, label, simIdx, exportData, true,min_objval,min_coef_val);

%% Try with other ratios 
for i = 0:10
    ratio = i/10
    pairedModel = changeObjective(pairedModel, {'Growth:L_plantarum', 'Growth:A_muciniphila'}, [ratio, 1-ratio]);
    sol(i+1) = optimizeCbModel(pairedModel,"max","one");
    disp(sol(i+1).f)
    disp(sol(i+1).x(1234))
    disp(sol(i+1).x(2526))
end

% At this point we see that above 0.5 for L_plantarum, A_muciniphila doesnt
% grow, while with 0.5 or less, their growth doesn't vary a lot 

%% Try with sampled point from Pareto Front:

selected_samples = cell2mat(pairedModel.fluxes(:, 2));

if size(selected_samples,1) == numel(pairedModel.rxns)
    targetIdx = find(~(abs(selected_samples) >= 1e-4));
else
    targetIdx = [];
end
zeroArr = zeros(size(targetIdx));

label = 'Joint model, biomass objective (0.6,0.4)';
short_label = 'Joint_biom_0.6_0.4';

% min_coef_val = -1000; % (No restriction)
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
[C_lb, C_ub, objBase, flagBase, sum_eps,C_lb_full, C_ub_full] = OVA(selected_samples, pairedModel.lb, pairedModel.ub, pairedModel.S, acceptableIdx,zeroArr, targetIdx,min_objval);
%% Save
save('OVA_pareto_target1e-4_minobj1e-2_results.mat', 'C_lb', 'C_ub', 'objBase', 'flagBase', 'sum_eps', 'C_lb_full', 'C_ub_full');
%% Load OVA results
load('OVA_pareto_target1e-4_minobj1e-2_results.mat')

%%
threshold = 0.2; 
meaningful_rxns = pairedModel.rxns(acceptableIdx(C_ub > threshold));
non_meaningful_idx = setdiff(1:numel(pairedModel.rxns), acceptableIdx(C_ub > threshold));
zeroArr_2 = zeros(size(non_meaningful_idx));

label = 'Joint model, objective from Pareto front after OVA min coefs 0.01'
short_label = 'Joint_biom_Pareto_OVA_coef0.01'

min_coef_val = -1000; % (No restriction)
% min_coef_val = 0.01;
[simIdx, exportData, objL1, objL0, objBase] = runinvFBA(selected_samples, pairedModel, zeroArr_2, non_meaningful_idx, label, simIdx, exportData, true,min_objval,min_coef_val);




%% Export excel
% Note for people cloning this repo from Github: this code block frequently leads to errors if Excel is open. Don't
% forget to save all you open excel files before running this block.
system('taskkill /F /IM excel.exe');

excelFileName = 'lacto_akker_invFBA_OVA_results_BioATPM.xlsx';

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


