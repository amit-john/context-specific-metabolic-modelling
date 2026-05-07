% Load cobratoolbox
clc;
addpath(genpath('C:/Users/USER/Documents/MATLAB/cobratoolbox'));
initCobraToolbox(0)
changeCobraSolver('gurobi','all')
changeCobraSolverParams('LP', 'feasTol', 1e-8);
addpath(genpath('C:/Users/USER/Desktop/Desktop/personal_docs/IITM_academics/5th_year/Final year Project/Making the MKN GSM model/SprintGapFiller-main/SprintGapFiller-main'));

%% Parent Model
clc;
model=load('C:/Users/USER/Downloads/Recon3D_301 (1)/Recon3D_301/Recon3DModel_301.mat').Recon3DModel;
%
%%
[grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution] = singleRxnDeletion(model, 'FBA');
essentialRxnIndices = find(grRateKO < 1e-6);
%
clc;
biomassTargetID = 'biomass_reaction'; 
bioIdx = find(strcmp(model.rxns, biomassTargetID));
essentialRxnIndices(end+1)=bioIdx;
ATP_id = find(ismember(model.rxns,'DM_atp_c_'));
essentialRxnIndices(end+1)=ATP_id;
%%
clc;
% New exchanges
exchanges = {
'EX_gly[e]',
'EX_ala_L[e]',
'EX_arg_L[e]',
'EX_asn_L[e]',
'EX_asp_L[e]',
'EX_cys_L[e]',
'EX_glu_L[e]',
'EX_gln_L[e]',
'EX_his_L[e]',
'EX_4hpro[e]',
'EX_ile_L[e]',
'EX_leu_L[e]',
'EX_lys_L[e]',
'EX_met_L[e]',
'EX_phe_L[e]',
'EX_pro_L[e]',
'EX_ser_L[e]',
'EX_thr_L[e]',
'EX_trp_L[e]',
'EX_tyr_L[e]',
'EX_val_L[e]',
'EX_ascb_L[e]',
'EX_btn[e]',
'EX_chol[e]',
'EX_pnto_R[e]',
'EX_fol[e]',
'EX_pydxn[e]',
'EX_ribflv[e]',
'EX_thm[e]',
'EX_inost[e]',
'EX_glc_D[e]',
'EX_etha[e]',
'EX_gthrd[e]',
'EX_na1[e]',
'EX_k[e]',
'EX_so4[e]',
'EX_hco3[e]',
'EX_o2[e]'
};

normallimits_0h = [
-0.04220478
-0.031653581
-0.363834293
-0.119899937
-0.04759937
-0.065734281
-0.043066102
-1.916559414
-0.030632499
-0.048326079
-0.120815206
-0.120815206
-0.069188158
-0.031866023
-0.028775984
-0.055049712
-0.090438814
-0.0531993
-0.007758231
-0.035170648
-0.054108692
-0.002733093
-0.000259456
-0.006782911
-0.000165899
-0.000717768
-0.001536582
-0.00016837
-0.000939275
-0.061548636
-3.517064763
-0.006165861
-0.001031061
-42.768369
-1.68891
-0.12888
-7.53658
-1
];


%%
[tf2, loc2] = ismember(model.rxns,exchanges);
idx_model_ids_new = find(tf2);    % Just for checking      
idx_exchanges_new = loc2(tf2); 

ex_index=find(findExcRxns(model));
model = changeRxnBounds(model, model.rxns(ex_index),-0.01, 'l');
model.lb(idx_model_ids_new) = normallimits_0h(idx_exchanges_new); 



%%
%Transcriptomic data
dat=readtable('C:/Users/USER/Downloads/context_data_clean_all_samples.csv',"VariableNamingRule","preserve");
dat1=readtable('C:/Users/USER/Downloads/GES1_control_ENSEMBL_log2.csv',"VariableNamingRule","preserve");
dat2=readtable('C:/Users/USER/Downloads/GES1_HP26695_ENSEMBL_log2.csv',"VariableNamingRule","preserve");
dat3=readtable('C:/Users/USER/Downloads/ges_count_data.csv',"VariableNamingRule","preserve");
%gene mapping
gene_names=readtable("C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Datasets/gene_identifiers.csv","VariableNamingRule","preserve");

% Preprocessing

dat = renamevars(dat, "GeneID", "GENE");
dat  = removevars(dat, {'GSM6102478','GSM6102479','GSM6102480','GSM6102481','GSM6102486','GSM6102487','GSM6102488','GSM6102489'});
dat1 = removevars(dat1, 'COND_MEDIAN');
dat2 = removevars(dat2, 'COND_MEDIAN');
%% Quantile normalisation

%GENE=dat(:,1);
tpm=dat(:,2:end);
%tpml=log2(tpm + 1); 
%tpml_withGene = [GENE, tpml];  % RNA-seq table with gene IDs

micro = outerjoin(dat1, dat2, 'Keys', 'GENE', 'MergeKeys', true);
GENE_micro   = micro(:,1);                        % gene IDs
micro_dat    = table2array(micro(:, 2:end));       % microarray values (log2)
micro_linear = 2.^micro_dat - 1;                  % inverse log2 → linear scale

% Convert back to table with GENE
micro_linear_table = [GENE_micro, array2table(micro_linear, ...
                      'VariableNames', micro.Properties.VariableNames(2:end))];

% RNA-seq - already in TPM (linear scale)
% dat has GENE in col 1, TPM in remaining cols
% NO transformation needed since microarray is now also in linear scale

% Innerjoin on GENE to get common genes
merged = innerjoin(micro_linear_table, dat, 'Keys', 'GENE');

g = merged(:,1);
fprintf('Common genes: %d\n', height(merged));

% Extract numeric data from merged table
n_micro_cols  = width(micro) - 1;         % microarray columns (excluding GENE)
n_rnaseq_cols = width(tpm);              % rnaseq columns

microarray_data = table2array(merged(:, 2 : n_micro_cols + 1));
rnaseq_log      = table2array(merged(:, n_micro_cols + 2 : end));

fprintf('Microarray: %d genes x %d samples\n', size(microarray_data,1), size(microarray_data,2));
fprintf('RNAseq:     %d genes x %d samples\n', size(rnaseq_log,1),      size(rnaseq_log,2));


% Quantile Normalization
combined = [microarray_data, rnaseq_log];
n_micro  = size(microarray_data, 2);
n_total  = size(combined, 2);

sorted   = zeros(size(combined));
sort_idx = zeros(size(combined));
for j = 1:n_total
    [sorted(:,j), sort_idx(:,j)] = sort(combined(:,j), 'ascend');
end

row_means = mean(sorted, 2);

normalized = zeros(size(combined));
for j = 1:n_total
    normalized(sort_idx(:,j), j) = row_means;
end

microarray_norm = normalized(:, 1:n_micro);
rnaseq_norm     = normalized(:, n_micro+1:end);

microarray_norm_table = [g, array2table(microarray_norm, ...
    'VariableNames', merged.Properties.VariableNames(2 : n_micro_cols + 1))];

% Reconstruct rnaseq normalized table
rnaseq_norm_table = [g, array2table(rnaseq_norm, ...
    'VariableNames', merged.Properties.VariableNames(n_micro_cols + 2 : end))];
fprintf('Microarray normalized: %d genes x %d samples\n', size(microarray_norm,1), size(microarray_norm,2));
fprintf('RNAseq normalized:     %d genes x %d samples\n', size(rnaseq_norm,1),     size(rnaseq_norm,2));


%%

m = outerjoin(g, microarray_norm_table , 'Keys', 'GENE', 'MergeKeys', true);
%mergedTable = outerjoin(dat1,dat2 , 'Keys', 'GENE', 'MergeKeys', true);
%mergedTable = outerjoin(dat3,mergedTable , 'Keys', 'GENE', 'MergeKeys', true);
%mergedTable = outerjoin(dat,mergedTable,'Keys','GENE','MergeKeys',true);
%mergedTable = rmmissing(mergedTable);
%mergedTable = dat;
mergedTable = m;
%%
% 1. Prepare the Mapping Table
% Extract only the columns we need: Ensembl (Key) and Entrez (Value)
mapping_subset = gene_names(:, {'ensembl_gene', 'entrez_id'});

% Remove duplicates from the mapping file to prevent data explosion
[~, uniqueIdx] = unique(mapping_subset.ensembl_gene, 'stable');
mapping_subset = mapping_subset(uniqueIdx, :);

% Rename key column to match mergedTable exactly
mapping_subset.Properties.VariableNames{'ensembl_gene'} = 'GENE';

% 2. Perform INNER JOIN (This is the key change)
% 'innerjoin' automatically discards rows that do not have a match in both tables.
% This effectively removes any row where entrez_id would have been NaN.
finalTable = innerjoin(mergedTable, mapping_subset, 'Keys', 'GENE');

% 3. (Optional) Move entrez_id to the second column
if ismember('entrez_id', finalTable.Properties.VariableNames)
    finalTable = movevars(finalTable, 'entrez_id', 'After', 'GENE');
end



%finalTable.entrez_id = string(finalTable.entrez_id) + ".1";

% Convert entrez to string (base IDs, no suffix yet)
base_entrez = string(finalTable.entrez_id);
expression_matrix = table2array(finalTable(:, 3:end));

% 3. Model isoform base IDs
model_gene_strings = string(model.genes);
model_gene_base = extractBefore(model_gene_strings, '.');

% 4. Expand expression matrix to cover ALL isoforms
expanded_genes = {};
expanded_values = [];

for i = 1:numel(base_entrez)
    % Find ALL isoforms in model matching this base Entrez ID
    isoform_idx = find(strcmp(model_gene_base, base_entrez(i)));
    
    if ~isempty(isoform_idx)
        for j = 1:numel(isoform_idx)
            expanded_genes{end+1} = model_gene_strings(isoform_idx(j));
            expanded_values(end+1, :) = expression_matrix(i, :);
        end
    end
end

% 4. View result
disp(['Rows remaining: ', num2str(height(finalTable))]);
disp(head(finalTable));

%% sprintcore
% 5. Build geneExpression struct
geneExpression = struct();
geneExpression.value = expanded_values;
geneExpression.genes = string(expanded_genes)';
geneExpression.context = string(finalTable.Properties.VariableNames(3:end));


fprintf('Original genes in expression data:  %d\n', numel(base_entrez));
fprintf('Expanded genes (all isoforms):       %d\n', numel(geneExpression.genes));
fprintf('Model genes total:                   %d\n', numel(model.genes));
fprintf('Expression matrix size:              %dx%d\n', size(geneExpression.value));
fprintf('Context (samples):                   %d\n', numel(geneExpression.context));

% Verify no isoform is duplicated
assert(numel(unique(geneExpression.genes)) == numel(geneExpression.genes), ...
    'Duplicate isoforms found — check mapping table for duplicate Ensembl IDs');

%%
%ccledata2023_new=readmatrix('C:/Users/USER/Downloads/corrected.csv');

%% sprintcore

%ccledata2023=table2array(finalTable(1:end,3:end));

%geneExpression=struct();
%geneExpression.value = ccledata2023;
%geneExpression.genes = finalTable.entrez_id;
%geneExpression.context = string(finalTable.Properties.VariableNames(3:end));



%%
%[RxnImp, Contexts] = GiniReactionImportance(geneExpression, model,90,10,coreRxn);
[RxnImp, Contexts] = GiniReactionImportance(geneExpression, model,90,10,essentialRxnIndices);
RxnImp(isnan(RxnImp))=0;
%%
%for mkn28_0
%con_id=find(strcmp(Contexts, ''));
curr_core = mean(RxnImp(:, 1:4), 2);
% weights are for those non-core reactions more the weights, less the
% probabilty of getting added
coreRxn=find(curr_core >= 0.5);
weights = zeros(numel(curr_core),1);
weights(curr_core<0.5)=1-curr_core(curr_core<0.5);


lb_before = model.lb;
ub_before = model.ub;
rxns_before = model.rxns;

[model_new1, LPS] = sprintcore(model, coreRxn, 1e-4, [], weights);

% Check if model was mutated
if ~isequal(model.lb, lb_before) || ~isequal(model.ub, ub_before)
    warning('sprintcore mutated the input model bounds!');
elseif ~isequal(model.rxns, rxns_before)
    warning('sprintcore mutated model.rxns!');
else
    disp('Model is safe — sprintcore did not mutate input.');
end

%%
%for mkn28_2
%con_id=find(strcmp(Contexts, ''));
curr_core = mean(RxnImp(:, 5:8), 2);
% weights are for those non-core reactions more the weights, less the
% probabilty of getting added
coreRxn=find(curr_core >= 0.5);
weights = zeros(numel(curr_core),1);
weights(curr_core<0.5)=1-curr_core(curr_core<0.5);
[model_new2,LPS] = sprintcore(model,coreRxn,1e-4,[],weights);

%%
%for mkn28_6
%con_id=find(strcmp(Contexts, ''));
curr_core = mean(RxnImp(:, 9:12), 2);
% weights are for those non-core reactions more the weights, less the
% probabilty of getting added
coreRxn=find(curr_core >= 0.5);
weights = zeros(numel(curr_core),1);
weights(curr_core<0.5)=1-curr_core(curr_core<0.5);
[model_new3,LPS] = sprintcore(model,coreRxn,1e-4,[],weights);

%%

%for ges1_0
curr_core = mean(RxnImp(:, 1:3), 2);
% weights are for those non-core reactions more the weights, less the
% probabilty of getting added
coreRxn=find(curr_core >= 0.5);
weights = zeros(numel(curr_core),1);
weights(curr_core<0.5)=1-curr_core(curr_core<0.5);


lb_before = model.lb;
ub_before = model.ub;
rxns_before = model.rxns;

[model_newg1, LPS] = sprintcore(model, coreRxn, 1e-4, [], weights);

% Check if model was mutated
if ~isequal(model.lb, lb_before) || ~isequal(model.ub, ub_before)
    warning('sprintcore mutated the input model bounds!');
elseif ~isequal(model.rxns, rxns_before)
    warning('sprintcore mutated model.rxns!');
else
    disp('Model is safe — sprintcore did not mutate input.');
end
%%
%for ges1_24
con_id=find(strcmp(Contexts, ''));
curr_core = mean(RxnImp(:, 4:6), 2);
% weights are for those non-core reactions more the weights, less the
% probabilty of getting added
coreRxn=find(curr_core >= 0.5);
weights = zeros(numel(curr_core),1);
weights(curr_core<0.5)=1-curr_core(curr_core<0.5);
[model_newg2,LPS] = sprintcore(model,coreRxn,1e-4,[],weights);
%%

save('C:/Users/USER/Downloads/NEW_model/ges1_0_updated.mat','model_newg1')
save('C:/Users/USER/Downloads/NEW_model/ges1_24_updated.mat','model_newg2')


%%
save('C:/Users/USER/Downloads/NEW_model/mkn28_0.mat','model_new1')
save('C:/Users/USER/Downloads/NEW_model/mkn28_2.mat','model_new2')
save('C:/Users/USER/Downloads/NEW_model/mkn28_6.mat','model_new3')

%[model_new1,LPS] = sprintcore(model,coreRxn,1e-4,[],weights);
%[model_new2,LPS] = sprintcore(model,coreRxn,1e-4,[],weights);
%[model_new3,LPS] = sprintcore(model,coreRxn,1e-4,[],weights);
%[model_new4,LPS] = sprintcore(model,coreRxn,1e-4,[],weights);

%%
%for GES1_non infected
con_id=find(strcmp(Contexts, ''));
curr_core = mean(RxnImp(:, 13:15), 2);
% weights are for those non-core reactions more the weights, less the
% probabilty of getting added
coreRxn=find(curr_core >= 0.5);
weights = zeros(numel(curr_core),1);
weights(curr_core<0.5)=1-curr_core(curr_core<0.5);
[model_new2,LPS] = sprintcore(model,coreRxn,1e-4,[],weights);


%%
biomassTargetID ='biomass_reaction'; 
bioIdx = find(strcmp(model.rxns, biomassTargetID)); % Use strcmp for exact match
model=changeObjective(model,biomassTargetID);
FBA=optimizeCbModel(model,'max');
FBA.stat 
FBA.f
fprintf('biomass flux after constraining %f',FBA.x(bioIdx))

%%

model_new1=load('C:/Users/USER/Downloads/NEW_model/mkn28_0.mat').model_new1;
model_new2=load('C:/Users/USER/Downloads/NEW_model/mkn28_2.mat').model_new2;
model_new3=load('C:/Users/USER/Downloads/NEW_model/mkn28_6.mat').model_new3;

%%
model_newg1=load('C:/Users/USER/Downloads/NEW_model/ges1_0_updated.mat').model_newg1;
model_newg2=load('C:/Users/USER/Downloads/NEW_model/ges1_24_updated.mat').model_newg2;

%%

% Adding Calcium
model_new1 = addReaction(model_new1, 'EX_ca2[e]', ...
    'reactionFormula', 'ca2[e] <=>', ...
    'reactionName',   'Exchange of Calcium', ...
    'subSystem',      'Exchange/demand reaction', ...
    'lowerBound',     -0.134125, ...
    'upperBound',      1000);


% Adding Chlorine
model_new1 = addReaction(model_new1, 'EX_cl[e]', ...
    'reactionFormula', 'cl[e] <=>', ...
    'reactionName',   'Exchange of Chloride', ...
    'subSystem',      'Exchange/demand reaction', ...
    'lowerBound',     -34.571522, ...
    'upperBound',      1000);

% Adding Calcium
model_new2 = addReaction(model_new2, 'EX_ca2[e]', ...
    'reactionFormula', 'ca2[e] <=>', ...
    'reactionName',   'Exchange of Calcium', ...
    'subSystem',      'Exchange/demand reaction', ...
    'lowerBound',     -0.134125, ...
    'upperBound',      1000);


% Adding Chlorine
model_new2 = addReaction(model_new2, 'EX_cl[e]', ...
    'reactionFormula', 'cl[e] <=>', ...
    'reactionName',   'Exchange of Chloride', ...
    'subSystem',      'Exchange/demand reaction', ...
    'lowerBound',     -34.571522, ...
    'upperBound',      1000);

% Adding Calcium
model_new3 = addReaction(model_new3, 'EX_ca2[e]', ...
    'reactionFormula', 'ca2[e] <=>', ...
    'reactionName',   'Exchange of Calcium', ...
    'subSystem',      'Exchange/demand reaction', ...
    'lowerBound',     -0.134125, ...
    'upperBound',      1000);


% Adding Chlorine
model_new3 = addReaction(model_new3, 'EX_cl[e]', ...
    'reactionFormula', 'cl[e] <=>', ...
    'reactionName',   'Exchange of Chloride', ...
    'subSystem',      'Exchange/demand reaction', ...
    'lowerBound',     -34.571522, ...
    'upperBound',      1000);


%%
model_newg1 = addReaction(model_newg1, 'EX_ca2[e]', ...
    'reactionFormula', 'ca2[e] <=>', ...
    'reactionName',   'Exchange of Calcium', ...
    'subSystem',      'Exchange/demand reaction', ...
    'lowerBound',     -0.134125, ...
    'upperBound',      1000);


% Adding Chlorine
model_newg1 = addReaction(model_newg1, 'EX_cl[e]', ...
    'reactionFormula', 'cl[e] <=>', ...
    'reactionName',   'Exchange of Chloride', ...
    'subSystem',      'Exchange/demand reaction', ...
    'lowerBound',     -34.571522, ...
    'upperBound',      1000);

model_newg2 = addReaction(model_newg2, 'EX_ca2[e]', ...
    'reactionFormula', 'ca2[e] <=>', ...
    'reactionName',   'Exchange of Calcium', ...
    'subSystem',      'Exchange/demand reaction', ...
    'lowerBound',     -0.134125, ...
    'upperBound',      1000);


% Adding Chlorine
model_newg2 = addReaction(model_newg2, 'EX_cl[e]', ...
    'reactionFormula', 'cl[e] <=>', ...
    'reactionName',   'Exchange of Chloride', ...
    'subSystem',      'Exchange/demand reaction', ...
    'lowerBound',     -34.571522, ...
    'upperBound',      1000);

%%
save('C:/Users/USER/Downloads/NEW_model/mkn28_0_updated.mat','model_new1');
save('C:/Users/USER/Downloads/NEW_model/mkn28_2_updated.mat','model_new2');
save('C:/Users/USER/Downloads/NEW_model/mkn28_6_updated.mat','model_new3');

%% Convertion to sbml

%load the models
model_new1=load('C:/Users/USER/Downloads/NEW_model/mkn28_0_updated.mat').model_new1;
model_new2=load('C:/Users/USER/Downloads/NEW_model/mkn28_2_updated.mat').model_new2;
model_new3=load('C:/Users/USER/Downloads/NEW_model/mkn28_6_updated.mat').model_new3;

%%

writeCbModel(model_new1, 'format', 'sbml', 'fileName', 'C:/Users/USER/Downloads/NEW_model/mkn28_0_updated.xml');
disp('0 done')
writeCbModel(model_new2, 'format', 'sbml', 'fileName', 'C:/Users/USER/Downloads/NEW_model/mkn28_2_updated.xml');
disp('2 done')
writeCbModel(model_new3, 'format', 'sbml', 'fileName', 'C:/Users/USER/Downloads/NEW_model/mkn28_6_updated.xml');
disp('6 done')
%%

biomassTargetID ='biomass_reaction'; 
bioIdx1 = find(strcmp(model_new1.rxns, biomassTargetID)); % Use strcmp for exact match
model_new1=changeObjective(model_new1,biomassTargetID);
FBA1=optimizeCbModel(model_new1,'max');
FBA1.stat 
FBA1.f
fprintf('biomass flux after constraining %f',FBA1.x(bioIdx1))

%%

biomassTargetID ='biomass_reaction'; 
bioIdx2 = find(strcmp(model_new2.rxns, biomassTargetID)); % Use strcmp for exact match
model_new2=changeObjective(model_new2,biomassTargetID);
FBA2=optimizeCbModel(model_new2,'max');
FBA2.stat 
FBA2.f
fprintf('biomass flux after constraining %f',FBA2.x(bioIdx2))

%%
biomassTargetID ='biomass_reaction'; 
bioIdx = find(strcmp(model_new3.rxns, biomassTargetID)); % Use strcmp for exact match
model_new3=changeObjective(model_new3,biomassTargetID);
FBA=optimizeCbModel(model_new3,'max');
FBA.stat 
FBA.f
fprintf('biomass flux after constraining %f',FBA.x(bioIdx))

%%

biomassTargetID ='biomass_reaction'; 
bioIdx = find(strcmp(model_newg1.rxns, biomassTargetID)); % Use strcmp for exact match
model_newg1=changeObjective(model_newg1,biomassTargetID);
FBA4=optimizeCbModel(model_newg1,'max');
FBA4.stat 
FBA4.f
fprintf('biomass flux after constraining %f',FBA4.x(bioIdx))

%%
bioIdx = find(strcmp(model_newg2.rxns, biomassTargetID)); % Use strcmp for exact match
model_newg2=changeObjective(model_newg2,biomassTargetID);
FBA5=optimizeCbModel(model_newg2,'max');
FBA5.stat 
FBA5.f
fprintf('biomass flux after constraining %f',FBA5.x(bioIdx))

%% Lets start with the FVAs
% MKNs
% 0
[minFlux, maxFlux] = fluxVariability(model_new1);
fluxTable = table(model_new1.rxns, minFlux, maxFlux,'VariableNames', {'ReactionID', 'minFlux', 'maxFlux'});
% Save the Table to a CSV File
% Specify the desired file name. The file will be saved in your current MATLAB folder.
%csvFileName = 'C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/DDP_models/MKN28_hp_0_non_infected_non_intergrated_flux_variability.csv';
csvFileName = 'C:/Users/USER/Downloads/NEW_model/mkn28_0updated_flux_variability.csv';
writetable(fluxTable, csvFileName);
disp(['Flux variability results saved to: ', csvFileName]);

%%
% 2 model

[minFlux, maxFlux] = fluxVariability(model_new2);
fluxTable = table(model_new2.rxns, minFlux, maxFlux,'VariableNames', {'ReactionID', 'minFlux', 'maxFlux'});
% Save the Table to a CSV File
% Specify the desired file name. The file will be saved in your current MATLAB folder.
%csvFileName = 'C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/DDP_models/MKN28_hp_0_non_infected_non_intergrated_flux_variability.csv';
csvFileName = 'C:/Users/USER/Downloads/NEW_model/mkn28_2_updated_flux_variability.csv';
writetable(fluxTable, csvFileName);
disp(['Flux variability results saved to: ', csvFileName]);

%%
% 6
[minFlux, maxFlux] = fluxVariability(model_new3);
fluxTable = table(model_new3.rxns, minFlux, maxFlux,'VariableNames', {'ReactionID', 'minFlux', 'maxFlux'});
% Save the Table to a CSV File
% Specify the desired file name. The file will be saved in your current MATLAB folder.
%csvFileName = 'C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/DDP_models/MKN28_hp_0_non_infected_non_intergrated_flux_variability.csv';
csvFileName = 'C:/Users/USER/Downloads/NEW_model/mkn28_6_updated_flux_variability.csv';
writetable(fluxTable, csvFileName);
disp(['Flux variability results saved to: ', csvFileName]);


%%
%ges1 0

[minFlux, maxFlux] = fluxVariability(model_newg1);
fluxTable = table(model_newg1.rxns, minFlux, maxFlux,'VariableNames', {'ReactionID', 'minFlux', 'maxFlux'});
% Save the Table to a CSV File
% Specify the desired file name. The file will be saved in your current MATLAB folder.
csvFileName = 'C:/Users/USER/Downloads/NEW_model/ges1_0_updated_flux_variability.csv';
writetable(fluxTable, csvFileName);
disp(['Flux variability results saved to: ', csvFileName]);

%%

%ges24
[minFlux, maxFlux] = fluxVariability(model_newg2);
fluxTable = table(model_newg2.rxns, minFlux, maxFlux,'VariableNames', {'ReactionID', 'minFlux', 'maxFlux'});
% Save the Table to a CSV File
% Specify the desired file name. The file will be saved in your current MATLAB folder.
csvFileName = 'C:/Users/USER/Downloads/NEW_model/ges1_24_updated_flux_variability.csv';
writetable(fluxTable, csvFileName);
disp(['Flux variability results saved to: ', csvFileName]);





%%
[minFlux, maxFlux] = fluxVariability(model_new1);
fluxTable = table(model_new1.rxns, minFlux, maxFlux,'VariableNames', {'ReactionID', 'minFlux', 'maxFlux'});
% Save the Table to a CSV File
% Specify the desired file name. The file will be saved in your current MATLAB folder.
%csvFileName = 'C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/DDP_models/MKN28_hp_0_non_infected_non_intergrated_flux_variability.csv';
csvFileName = 'C:/Users/USER/Downloads/mkn28_0_26_3_26_flux_variability.csv';
writetable(fluxTable, csvFileName);
disp(['Flux variability results saved to: ', csvFileName]);


%% Integration step

%load the models
model_new1=load('C:/Users/USER/Downloads/NEW_model/mkn28_0_updated.mat').model_new1;
model_new2=load('C:/Users/USER/Downloads/NEW_model/mkn28_2_updated.mat').model_new2;
model_new3=load('C:/Users/USER/Downloads/NEW_model/mkn28_6_updated.mat').model_new3;

%%

%filename = 'C:/Users/USER/Downloads/Helicobacter_pylori_26695.xml'; % Replace with your actual filename
%pathogen = readCbModel(filename);
pathogen_biomass_name='biomass525';
pathogen=load('C:/Users/USER/Downloads/Helicobacter_pylori_26695 (1).mat').model;

pathogen_biomass_idx = find(strcmp(pathogen.rxns, pathogen_biomass_name));
temp_model = changeObjective(pathogen, 'biomass525');
FBAsol2 = optimizeCbModel(temp_model);
biomass_flux = FBAsol2.x(pathogen_biomass_idx);
fprintf('The pathogen biomass flux is: %f\n', biomass_flux);

%%
ex_index=find(findExcRxns(pathogen));
pathogen = changeRxnBounds(pathogen, pathogen.rxns(ex_index),-1, 'l');
pathogen = changeRxnBounds(pathogen, pathogen.rxns(ex_index), 10, 'u');
pathogen = changeRxnBounds(pathogen, 'EX_o2(e)', 0, 'l');
pathogen = changeRxnBounds(pathogen,'EX_glc_D(e)',-20,'l');

temp_model = changeObjective(pathogen, 'biomass525');
FBAsol2 = optimizeCbModel(temp_model);
biomass_flux = FBAsol2.x(pathogen_biomass_idx);
fprintf('The pathogen biomass flux is: %f\n', biomass_flux);


%%
clc;
[modelJoint1] = createMultipleSpeciesModel( ...
    {pathogen}, ...                      % The pathogen model in a cell array
    {pathogen_biomass_name}, ...         % Its biomass reaction name
    'modelHost', model_new1, ...          % Your human model as the host
    'nameTagsModels', {'pathogen'}, ...  % A tag for the pathogen's reactions
    'nameTagHost', 'human', ...          % A tag for the human's reactions
    'mergeGenesFlag', true ...           % Explicitly set to true to merge genes and GPRs
);

[modelJoint2] = createMultipleSpeciesModel( ...
    {pathogen}, ...                      % The pathogen model in a cell array
    {pathogen_biomass_name}, ...         % Its biomass reaction name
    'modelHost', model_new2, ...          % Your human model as the host
    'nameTagsModels', {'pathogen'}, ...  % A tag for the pathogen's reactions
    'nameTagHost', 'human', ...          % A tag for the human's reactions
    'mergeGenesFlag', true ...           % Explicitly set to true to merge genes and GPRs
);

[modelJoint3] = createMultipleSpeciesModel( ...
    {pathogen}, ...                      % The pathogen model in a cell array
    {pathogen_biomass_name}, ...         % Its biomass reaction name
    'modelHost', model_new3, ...          % Your human model as the host
    'nameTagsModels', {'pathogen'}, ...  % A tag for the pathogen's reactions
    'nameTagHost', 'human', ...          % A tag for the human's reactions
    'mergeGenesFlag', true ...           % Explicitly set to true to merge genes and GPRs
);



%%  Too restrictive
start_prefix = 'humanEX_';
end_suffix = 'b';

modelJoint4=modelJoint3;

% 1. Block the complete blood compartment
start_match = startsWith(modelJoint4.rxns, start_prefix, 'IgnoreCase', false);
end_match = endsWith(modelJoint4.rxns, end_suffix, 'IgnoreCase', false);
final_matching_indices = start_match & end_match;
reactions_list = modelJoint4.rxns(final_matching_indices);
modelJoint4 = changeRxnBounds(modelJoint4, reactions_list , 0, 'b');
[~, logical_indices] = ismember(reactions_list, modelJoint4.rxns);
lower_bounds_fluxes = modelJoint4.lb(logical_indices);

%%

modelJoint4 = modelJoint1;

% Find ALL reactions involving blood compartment
contains_blood = endsWith(modelJoint4.rxns, '[e]b');
reactions_list = modelJoint4.rxns(contains_blood);
% Block them completely
modelJoint4 = changeRxnBounds(modelJoint4, reactions_list, 0, 'b');
% Optional check
[~, logical_indices] = ismember(reactions_list, modelJoint4.rxns);
lower_bounds_fluxes = modelJoint4.lb(logical_indices);


lumen_exchanges = {
'EX_gly[u]',
'EX_ala_L[u]',
'EX_arg_L[u]',
'EX_asn_L[u]',
'EX_asp_L[u]',
'EX_cys_L[u]',
'EX_glu_L[u]',
'EX_gln_L[u]',
'EX_his_L[u]',
'EX_4hpro[u]',
'EX_ile_L[u]',
'EX_leu_L[u]',
'EX_lys_L[u]',
'EX_met_L[u]',
'EX_phe_L[u]',
'EX_pro_L[u]',
'EX_ser_L[u]',
'EX_thr_L[u]',
'EX_trp_L[u]',
'EX_tyr_L[u]',
'EX_val_L[u]',
'EX_ascb_L[u]',
'EX_btn[u]',
'EX_chol[u]',
'EX_pnto_R[u]',
'EX_fol[u]',
'EX_pydxn[u]',
'EX_ribflv[u]',
'EX_thm[u]',
'EX_inost[u]',
'EX_glc_D[u]',
'EX_etha[u]',
'EX_gthrd[u]',
'EX_na1[u]',
'EX_k[u]',
'EX_so4[u]',
'EX_hco3[u]',
'EX_o2[u]'
};


start_prefix1 = 'EX_';
end_suffix1 = '[u]';
start_match1 = startsWith(modelJoint4.rxns, start_prefix1, 'IgnoreCase', false);
end_match1 = endsWith(modelJoint4.rxns, end_suffix1, 'IgnoreCase', false);
final_matching_indices1 = start_match1 & end_match1;
reactions_list1 = modelJoint4.rxns(final_matching_indices1);
modelJoint4 = changeRxnBounds(modelJoint4, reactions_list1 , -0.01, 'l');
[~, logical_indices] = ismember(reactions_list1, modelJoint4.rxns);
lower_bounds_fluxes1 = modelJoint4.lb(logical_indices);

%
[tf2, loc2] = ismember(modelJoint4.rxns,lumen_exchanges);
idx_model_ids_new1 = find(tf2);    % Just for checking      
idx_exchanges_new1 = loc2(tf2); 
modelJoint4.lb(idx_model_ids_new1) = normallimits_0h(idx_exchanges_new1); 
lower_bounds_fluxes1 = modelJoint4.lb(logical_indices);


MOI=100;
% Get pathogen IEX reactions
start_prefix1 = 'pathogenIEX';
end_suffix1 = '[u]tr';
start_match1 = startsWith(modelJoint4.rxns, start_prefix1, 'IgnoreCase', false);
end_match1 = endsWith(modelJoint4.rxns, end_suffix1, 'IgnoreCase', false);
final_matching_indices1 = start_match1 & end_match1;

pathogen_uptake_rxns = modelJoint4.rxns(startsWith(modelJoint4.rxns, 'pathogenIEX','IgnoreCase', false) & endsWith(modelJoint4.rxns, '[u]tr', 'IgnoreCase', false));

modelJoint4 = changeRxnBounds(modelJoint4, modelJoint4.rxns(final_matching_indices1),-0.01*100, 'l');
modelJoint4 = changeRxnBounds(modelJoint4, modelJoint4.rxns(final_matching_indices1), 10*100, 'u');
modelJoint4 = changeRxnBounds(modelJoint4, 'pathogenIEX_o2[u]tr', -5*100, 'l');
modelJoint4 = changeRxnBounds(modelJoint4,'pathogenIEX_glc_D[u]tr',-20*100,'l');

% biomass
%%
modelJoint4=modelJoint3;

pathogen_biomass_rxn = 'pathogenbiomass525';
pathogen_bio_idx = find(strcmp(modelJoint4.rxns, pathogen_biomass_rxn));
temp_model = changeObjective(modelJoint4, pathogen_biomass_rxn);
FBAsol2 = optimizeCbModel(temp_model);
biomass_flux = FBAsol2.x(pathogen_bio_idx);
fprintf('The pathogen biomass flux is: %f\n', biomass_flux);



human_biomass_rxn = 'humanbiomass_reaction';
human_bio_idx = find(strcmp(modelJoint4.rxns, human_biomass_rxn));
temp_model = changeObjective(modelJoint4, human_biomass_rxn);
FBAsol2 = optimizeCbModel(temp_model);
biomass_flux = FBAsol2.x(human_bio_idx);
fprintf('The human biomass flux is: %f\n', biomass_flux);

%%
modelJoint=modelJoint4;
save('C:/Users/USER/Downloads/NEW_model/mkn28_0_hp_integrated_updated.mat','modelJoint');

%%

modelJoint1=load('C:/Users/USER/Downloads/NEW_model/mkn28_0_hp_integrated_updated.mat').modelJoint;
modelJoint2=load('C:/Users/USER/Downloads/NEW_model/mkn28_2_hp_integrated_updated.mat').modelJoint;
modelJoint3=load('C:/Users/USER/Downloads/NEW_model/mkn28_6_hp_integrated_updated.mat').modelJoint;


%%
[minFlux, maxFlux] = fluxVariability(modelJoint3);
fluxTable = table(modelJoint3.rxns, minFlux, maxFlux,'VariableNames', {'ReactionID', 'minFlux', 'maxFlux'});
csvFileName = 'C:/Users/USER/Downloads/NEW_model/mkn28_6_hp_integrated_flux_variability.csv';
writetable(fluxTable, csvFileName);
disp(['Flux variability results saved to: ', csvFileName]);

%%
[minFlux, maxFlux] = fluxVariability(modelJoint1);
fluxTable = table(modelJoint1.rxns, minFlux, maxFlux,'VariableNames', {'ReactionID', 'minFlux', 'maxFlux'});
csvFileName = 'C:/Users/USER/Downloads/NEW_model/mkn28_0_hp_integrated_flux_variability.csv';
writetable(fluxTable, csvFileName);
disp(['Flux variability results saved to: ', csvFileName]);

%%

[minFlux, maxFlux] = fluxVariability(modelJoint2);
fluxTable = table(modelJoint2.rxns, minFlux, maxFlux,'VariableNames', {'ReactionID', 'minFlux', 'maxFlux'});
csvFileName = 'C:/Users/USER/Downloads/NEW_model/mkn28_2_hp_integrated_flux_variability.csv';
writetable(fluxTable, csvFileName);
disp(['Flux variability results saved to: ', csvFileName]);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exchanges = {
    'EX_glc_D[e]'
    'EX_cl[e]'
    'EX_na1[e]'
    'EX_ca24g[e]'
    'EX_k[e]'
    'EX_so4[e]'
    'EX_ncam[e]'
    'EX_pnto_R[e]'
    'EX_chol[e]'
    'EX_fol[e]'
    'EX_inost[e]'
    'EX_pydxn[e]'
    'EX_ribflv[e]'
    'EX_thm[e]'
    'EX_pi[e]'
    'EX_arg_L[e]'
    'EX_cyst_L[e]'
    'EX_his_L[e]'
    'EX_ile_L[e]'
    'EX_leu_L[e]'
    'EX_lys_L[e]'
    'EX_met_L[e]'
    'EX_phe_L[e]'
    'EX_thr_L[e]'
    'EX_trp_L[e]'
    'EX_tyr_L[e]'
    'EX_val_L[e]'
    'EX_pchol_hs[e]' 
    'EX_pe_hs[e]'    
    'EX_ps_hs[e]'    
    'EX_sphmyln_hs[e]' 
    'EX_chsterol[e]'
    'EX_o2[e]'   
    'EX_hco3[e]'
};


normallimits_0h = [
   -8.793506426214   % EX_glc__D_e
 -199.718966257211   % EX_cl_e
 -187.179075066124   % EX_na1_e
   -2.851948111404   % EX_ca2_e
   -8.441765738635   % EX_k_e
   -1.288292678238   % EX_so4_e
   -0.012974025044   % EX_ncam_e
   -0.003318304133   % EX_pnto__R_e
   -0.011305937105   % EX_chol_e
   -0.003589186368   % EX_fol_e
   -0.017587012220   % EX_inost_e
   -0.007758976398   % EX_pydxn_e
   -0.000420965712   % EX_ribflv_e
   -0.004696828314   % EX_thm_e
   -1.605770770356   % EX_pi_e
   -0.945197730779   % EX_arg__L_e
   -0.158283111556   % EX_cyst__L_e
   -0.316566223113   % EX_his__L_e
   -0.628299382213   % EX_ile__L_e
   -0.628299382213   % EX_leu__L_e
   -0.623185713900   % EX_lys__L_e
   -0.159345412831   % EX_met__L_e
   -0.306973316854   % EX_phe__L_e
   -0.638452927253   % EX_thr__L_e
   -0.077589760815   % EX_trp__L_e
   -0.312955199691   % EX_tyr__L_e
   -0.622309680190   % EX_val__L_e
   -10.0/10 % EX_pchol_hs_e uptake rate
   -10.0/10 % EX_pe_hs_e uptake rate
   -10.0/10 % EX_ps_hs_e uptake rate
   -10.0/10 % EX_sphmyln_hs_e uptake rate
   -10.0/10 % EX_chsterol_e uptake rate
   -10.0/10 % EX_o2_e (typically a high negative number for unlimited air)
   -10.0/10 % EX_hco3_e
];

normallimits_2h = [
   -8.396342603586   % EX_glc__D_e
 -190.698543203487   % EX_cl_e
 -178.725023477820   % EX_na1_e
   -2.723138219313   % EX_ca2_e
   -8.060488488357   % EX_k_e
   -1.230106191534   % EX_so4_e
   -0.012388045671   % EX_ncam_e
   -0.003168431000   % EX_pnto__R_e
   -0.010795297892   % EX_chol_e
   -0.003427078682   % EX_fol_e
   -0.016792684603   % EX_inost_e
   -0.007408537725   % EX_pydxn_e
   -0.000401952551   % EX_ribflv_e
   -0.004484693337   % EX_thm_e
   -1.533245201316   % EX_pi_e
   -0.902507326553   % EX_arg__L_e
   -0.151134162935   % EX_cyst__L_e
   -0.302268325870   % EX_his__L_e
   -0.599921875869   % EX_ile__L_e
   -0.599921875869   % EX_leu__L_e
   -0.595039169354   % EX_lys__L_e
   -0.152148484756   % EX_met__L_e
   -0.293108688791   % EX_phe__L_e
   -0.609616830153   % EX_thr__L_e
   -0.074085374225   % EX_trp__L_e
   -0.298820396417   % EX_tyr__L_e
   -0.594202702215   % EX_val__L_e
   -10.0/10 % EX_pchol_hs_e uptake rate
   -10.0/10 % EX_pe_hs_e uptake rate
   -10.0/10 % EX_ps_hs_e uptake rate
   -10.0/10 % EX_sphmyln_hs_e uptake rate
   -10.0/10 % EX_chsterol_e uptake rate
   -10.0/10 % EX_o2_e (typically a high negative number for unlimited air)
   -10.0/10 % EX_hco3_e
];

normallimits_6h = [
   -2.551711457776   % EX_glc__D_e
  -57.954716791303   % EX_cl_e
  -54.315874390941   % EX_na1_e
   -0.827582117999   % EX_ca2_e
   -2.449642874533   % EX_k_e
   -0.373838492712   % EX_so4_e
   -0.003764819943   % EX_ncam_e
   -0.000962909932   % EX_pnto__R_e
   -0.003280771953   % EX_chol_e
   -0.001041514901   % EX_fol_e
   -0.005103422732   % EX_inost_e
   -0.002251510150   % EX_pydxn_e
   -0.000122156393   % EX_ribflv_e
   -0.001362931923   % EX_thm_e
   -0.465964710171   % EX_pi_e
   -0.274278741902   % EX_arg__L_e
   -0.045930805046   % EX_cyst__L_e
   -0.091861610092   % EX_his__L_e
   -0.182320755203   % EX_ile__L_e
   -0.182320755203   % EX_leu__L_e
   -0.180836864091   % EX_lys__L_e
   -0.046239065051   % EX_met__L_e
   -0.089077927721   % EX_phe__L_e
   -0.185267124485   % EX_thr__L_e
   -0.022515100585   % EX_trp__L_e
   -0.090813758484   % EX_tyr__L_e
   -0.180582655457   % EX_val__L_e
   -10.0/10 % EX_pchol_hs_e uptake rate
   -10.0/10 % EX_pe_hs_e uptake rate
   -10.0/10 % EX_ps_hs_e uptake rate
   -10.0/10 % EX_sphmyln_hs_e uptake rate
   -10.0/10 % EX_chsterol_e uptake rate
   -10.0/10 % EX_o2_e (typically a high negative number for unlimited air)
   -10.0/10 % EX_hco3_e
];


clc;
exchanges = {
    'EX_glc_D[e]'
    'EX_cl[e]'
    'EX_na1[e]'
    'EX_ca2[e]'
    'EX_k[e]'
    'EX_so4[e]'
    'EX_ncam[e]'
    'EX_pnto_R[e]'
    'EX_chol[e]'
    'EX_fol[e]'
    'EX_inost[e]'
    'EX_pydxn[e]'
    'EX_ribflv[e]'
    'EX_thm[e]'
    'EX_pi[e]'
    'EX_arg_L[e]'
    'EX_cyst_L[e]'
    'EX_his_L[e]'
    'EX_ile_L[e]'
    'EX_leu_L[e]'
    'EX_lys_L[e]'
    'EX_met_L[e]'
    'EX_phe_L[e]'
    'EX_thr_L[e]'
    'EX_trp_L[e]'
    'EX_tyr_L[e]'
    'EX_val_L[e]'
    'EX_pchol_hs[e]' 
    'EX_pe_hs[e]'    
    'EX_ps_hs[e]'    
    'EX_sphmyln_hs[e]' 
    'EX_chsterol[e]'
    'EX_o2[e]'   
    'EX_hco3[e]'
};

normallimits_0h = [
    -3.517064763
-32.7450842
-32.7450842
-0.134125355
-1.688191156
-0.128830084
-0.002594556
-0.000165899
-0.01
-0.000717768
-0.061548636
-0.001536582
-0.00016837
-0.000939275
-1.783300519
-0.363834293
-0.065734281
-0.030632499
-0.120815206
-0.120815206
-0.069188158
-0.031866023
-0.028775984
-0.0531993
-0.007758231
-0.035170648
-0.054108692
-0.01
-0.01
-0.01
-0.01
-0.01
-0.01
-7.536567802
    
];

normallimits_2h = [
    -3.359183155
-31.27515205
-31.27515205
-0.128104446
-1.612407981
-0.12304688
-0.002478086
-0.000158452
-0.009551098
-0.000685548
-0.058785707
-0.001467604
-0.000160812
-0.000897111
-1.703247869
-0.347501713
-0.062783458
-0.029257401
-0.115391792
-0.115391792
-0.066082291
-0.030435552
-0.027484226
-0.050811175
-0.007409963
-0.033591832
-0.051679744
-0.009551098
-0.009551098
-0.009551098
-0.009551098
-0.009551098
-0.009551098
-7.19825005

];

normallimits_6h = [
    -3.063645264
-28.52359251
-28.52359251
-0.116833933
-1.470549787
-0.112221327
-0.002260066
-0.000144512
-0.008710801
-0.000625234
-0.053613794
-0.001338486
-0.000146664
-0.000818184
-1.553397665
-0.316928827
-0.057259827
-0.026683361
-0.105239727
-0.105239727
-0.060268431
-0.02775786
-0.025066189
-0.046340854
-0.006758041
-0.030636453
-0.047133007
-0.008710801
-0.008710801
-0.008710801
-0.008710801
-0.008710801
-0.008710801
-6.564954531

];

normallimits_2h = [
-0.020619989
-0.01546499
-0.177758521
-0.058579512
-0.023255624
-0.032115798
-0.021040805
-0.936373433
-0.01496612
-0.023610672
-0.059026685
-0.059026685
-0.033803258
-0.015568783
-0.014059083
-0.026895638
-0.044185691
-0.025991582
-0.003790439
-0.017183323
-0.026435884
-0.001335307
-0.000126762
-0.003313927
-8.10533E-05
-0.00035068
-0.000750728
-8.22605E-05
-0.000458902
-0.030070817
-1.718332331
-0.003012455
-0.000503745
-20.60808047
-20.89534204
-0.065529451
-0.825150759
-0.062966902
-3.682146891
-0.48857
    
];

normallimits_6h = [
-0.006266566
-0.004699924
-0.054022116
-0.017802743
-0.007067554
-0.009760226
-0.006394455
-0.284570742
-0.004548313
-0.007175456
-0.017938642
-0.017938642
-0.010273058
-0.004731467
-0.004272658
-0.008173781
-0.013428355
-0.007899032
-0.001151942
-0.005222138
-0.008034059
-0.00040581
-3.8524E-05
-0.001007127
-2.46327E-05
-0.000106574
-0.000228152
-2.49996E-05
-0.000139464
-0.009138741
-0.522213776
-0.000915507
-0.000153092
-6.262946534
-6.350247429
-0.01991488
-0.250769357
-0.019136102
-1.119031398
-0.14848
    
];

%%
bioIdx = find(strcmp(model.rxns, 'biomass_reaction'));

% Get all metabolites participating in biomass reaction
bio_met_indices = find(model.S(:, bioIdx) ~= 0);
bio_mets = model.mets(bio_met_indices);
bio_stoich = full(model.S(bio_met_indices, bioIdx));

% Display them with stoichiometry
disp('Metabolites in biomass reaction:');
for i = 1:numel(bio_mets)
    fprintf('  %+.4f  %s\n', bio_stoich(i), bio_mets{i});
end

% Specifically check if ca2 or cl are among them
has_ca2 = any(contains(bio_mets, 'ca2'));
has_cl  = any(contains(bio_mets, 'cl'));
fprintf('\nCalcium (ca2) in biomass: %d\n', has_ca2);
fprintf('Chloride (cl) in biomass: %d\n', has_cl);

%%
% Extract all suffixes from model.genes
gene_strings = string(model.genes);
suffixes = extractAfter(gene_strings, '.');
suffix_nums = str2double(suffixes);

% Summary
max_isoform = max(suffix_nums);
fprintf('Maximum isoform number: %d\n', max_isoform);

% Distribution of isoform numbers
fprintf('\nIsoform distribution:\n');
for i = 1:max_isoform
    count = sum(suffix_nums == i);
    if count > 0
        fprintf('  .%d  →  %d genes\n', i, count);
    end
end

% Total genes per category
fprintf('\nTotal model genes: %d\n', numel(model.genes));
fprintf('Unique base Entrez IDs: %d\n', numel(unique(extractBefore(gene_strings, '.'))));


%%
clc;
model=model_new1;
metID = 'xol7ah[c]';
metIdx = find(strcmp(model.mets, metID));
rxnIdx = find(model.S(metIdx,:) ~= 0);
for i = 1:length(rxnIdx)

    fprintf('Reaction: %s\n', model.rxns{rxnIdx(i)});

    sub = model.subSystems{rxnIdx(i)};

    if iscell(sub)
        sub = strjoin(sub, ', ');
    end

    fprintf('Subsystem: %s\n\n', sub);

end
%%