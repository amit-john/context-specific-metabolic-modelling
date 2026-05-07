%%  Loading Cobratool box
clc;
addpath(genpath('C:/Users/USER/Documents/MATLAB/cobratoolbox'));
initCobraToolbox(0)
changeCobraSolver('gurobi','all')
changeCobraSolverParams('LP', 'feasTol', 1e-8);
addpath(genpath('C:/Users/USER/Desktop/Desktop/personal_docs/IITM_academics/5th_year/Final year Project/Making the MKN GSM model/SprintGapFiller-main/SprintGapFiller-main'));

%% Loading the Human model
model=load('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Datasets/Recon3D_301/Recon3D_301/Consistant_Recon3D.mat').model;

%% Loading the transcriptomic data 

t_data=readtable("C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/MKN28/rnaseq_tpm_MKN28.csv","VariableNamingRule","preserve");
gene_names=readtable("C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Datasets/Expression data/symbol2entrez.csv");

%%
td = t_data;
td.gene_symbol = upper(string(td.gene_symbol));
td.expression  = double(td.expression);

gn = gene_names;
gn.SYMBOL   = upper(string(gn.SYMBOL));
% make ENTREZID text
gn.ENTREZID = string(gn.ENTREZID);

%% 1) Collapse duplicates in t_data at SYMBOL level (choose @mean or @max)
[Gsym, sym_ids] = findgroups(td.gene_symbol);
sym_expr = splitapply(@mean, td.expression, Gsym);  % or @max
expr_by_symbol = table(sym_ids, sym_expr, ...
    'VariableNames', {'SYMBOL','Expr'});

%% 2) Convert model.genes (e.g., '55293.1') -> ENTREZ (strip version after '.')
mg_raw = string(model.genes(:));
mg_entrez = extractBefore(mg_raw, '.');          % '55293'
% If some genes have no '.', extractBefore returns <missing>, fix them:
noDot = ismissing(mg_entrez);
mg_entrez(noDot) = mg_raw(noDot);               % keep as-is if no dot

%% 3) Map ENTREZ -> SYMBOL using gene_names
[tf_map, idx_map] = ismember(mg_entrez, gn.ENTREZID);
model_symbol = strings(numel(mg_entrez),1);
model_symbol(:) = missing;
model_symbol(tf_map) = gn.SYMBOL(idx_map(tf_map));

%% 4) Join SYMBOL -> expression (from MKN28)
[tf_expr, idx_expr] = ismember(model_symbol, expr_by_symbol.SYMBOL);

% 5) Build gene_expression struct aligned to model order, drop NaNs
vals = nan(numel(mg_raw),1);
vals(tf_expr) = expr_by_symbol.Expr(idx_expr(tf_expr));

keep = ~isnan(vals);                    % only genes we could map & found expr
gene_expression = struct();
gene_expression.genes   = model.genes(keep);      % keep original model IDs
gene_expression.value   = vals(keep);
gene_expression.context = {'MKN28'};
%%

bio_id = find(ismember(model.rxns,'biomass_reaction'));
ATP_id = find(ismember(model.rxns,'DM_atp_c_'));
% If you specifically want extracellular [e] exchanges only:
%isExtracellular = contains(model.rxns,'[e]') & contains(model.rxns, 'EX_');
exchanges={'EX_glc_D[e]','EX_hco3[e]','EX_k[e]','EX_so4[e]','EX_pnto_R[e]','EX_fol[e]','EX_inost[e]','EX_pydxn[e]','EX_ribflv[e]','EX_thm[e]','EX_pi[e]','EX_chol[e]','EX_arg_L[e]','EX_Lcystin[e]','EX_gln_L[e]','EX_gly[e]','EX_his_L[e]','EX_ile_L[e]','EX_leu_L[e]','EX_lys_L[e]','EX_met_L[e]','EX_phe_L[e]','EX_ser_L[e]','EX_thr_L[e]','EX_trp_L[e]','EX_tyr_L[e]','EX_val_L[e]','EX_asn_L[e]','EX_asp_L[e]','EX_glu_L[e]','EX_pro_L[e]','EX_btn[e]','EX_aqcobal[e]'};

indi=[];
for i=1:length(exchanges)
    idx_logical = ismember(model.rxns, exchanges{i});
    idx = find(idx_logical);
    
    % Append the single index to the array
    if ~isempty(idx)
        indi(end+1) = idx;
    else
        fprintf('Warning: Reaction %s not found in the model.\n', exchanges{i});
    end
end

%%

coreRxn = [bio_id;ATP_id];
for i=1:length(indi)
    coreRxn(end+1)=indi(i);
end
[RxnImp,Contexts] = GiniReactionImportance(gene_expression,model,90,10,coreRxn);
RxnImp(isnan(RxnImp))=0;
con_id=find(strcmp(Contexts, 'MKN28'));
curr_core = RxnImp(:,con_id);
% weights are for those non-core reactions more the weights, less the
% probabilty of getting added
weights = zeros(numel(curr_core),1);
weights(curr_core<1)=1-curr_core(curr_core<1);
% core reactions are those with RxnImp more than 1
core = find(curr_core>=1);

%%

clc;
model_new = sprintcore(model_new,core,1e-4,[],weights);


%% Constrainst
%id_id=find(strcmp(model_new.rxns,'EX_ncam[e]'));

limits=[420.4458,901.4964,203.215,31.4728,0.153,0.08583,7.348332,0.2227,0.002,0.112497,403.4007,0.814377,43.5597,7.878624,151.512,5.0377,3.7006,14.4693,14.4693,7.38621,3.8067,3.4393,10.79523,6.3635,0.928,4.0908,6.465,14.317884,5.7195,5.1574,6.591,31.059,0.0001397];


exIdx = find(findExcRxns(model_new));                 % indices into mp.rxns
[~, posInExIdx, posInDMEM] = intersect(model_new.rxns(exIdx), exchanges, 'stable');

model_new.lb(exIdx) = 0; 
model_new.ub(exIdx) = 1000; 
model_new.lb(exIdx(posInExIdx)) = -limits(posInDMEM);  % overwrite with DMEM caps


%%
exchangerxns={'EX_glc_D[e]'
'EX_cl[e]'
'EX_h2o[e]'
'EX_hco3[e]'
'EX_h[e]'
'EX_na1[e]'
'EX_ca2[e]'
'EX_fe3[e]'
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
'EX_pyr[e]'
'EX_creat[e]'
'EX_bilirub[e]'
'EX_pi[e]'
'EX_urea[e]'
'EX_chsterol[e]'
'EX_cortsn[e]'
'EX_tststerone[e]'
'EX_prgstrn[e]'
'EX_prostge1[e]'
'EX_prostge2[e]'
'EX_prostgf2[e]'
'EX_retinol[e]'
'EX_avite1[e]'
'EX_M02050[e]'
'EX_arg_L[e]'
'EX_Lcystin[e]'
'EX_gln_L[e]'
'EX_gly[e]'
'EX_his_L[e]'
'EX_ile_L[e]'
'EX_leu_L[e]'
'EX_lys_L[e]'
'EX_met_L[e]'
'EX_phe_L[e]'
'EX_ser_L[e]'
'EX_thr_L[e]'
'EX_trp_L[e]'
'EX_tyr_L[e]'
'EX_val_L[e]'
'EX_o2[e]'};

cancerlimits=[-16.04515122
-76.42102832
-1.207674955
-11.5168974
-1.143436078
-1.121735968
-1.162371598
-0.000159655
-3.460764588
-0.523523303
-0.020471667
-0.00524659
-0.017905744
-0.005663797
-0.024978352
-0.012157168
-0.000664187
-0.0074118
-0.625
-0.014775414
-0.000427597
-0.065809563
-0.166500167
-5.00926E-05
-8.66972E-08
-8.66791E-08
-1.59003E-05
-1.06394E-06
-1.06394E-06
-2.11573E-06
-1.96368E-05
-0.000145109
-0.004868339
-0.249216747
-0.124920179
-2.497605036
-0.249766884
-0.125220627
-0.500304948
-0.500304948
-0.499589379
-0.125660977
-0.249712452
-0.249785898
-0.498446944
-0.048966801
-0.248367984
-0.501493811
-1000];
%%
clc;
exIdx  = find(findExcRxns(model_new));
exRxns = model_new.rxns(exIdx);
commonRxns = intersect(exRxns, exchangerxns);

dmem_expression_indices=[];
for i=1:length(commonRxns)
    j = find(strcmp(commonRxns{i}, exchangerxns));   % or strcmpi for case-insensitive
    dmem_expression_indices = [dmem_expression_indices; j(:)];
end

dmem_normallb_constraints=[];
for i=1:length(dmem_expression_indices)
    dmem_normallb_constraints(end+1)=cancerlimits(dmem_expression_indices(i));
end
dmem_normallb_constraints=dmem_normallb_constraints';

model_new.lb(dmem_expression_indices) = dmem_normallb_constraints;
%%

model_new=load('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/MKN28/MKN28_model.mat').model_new;
%%
biomass_idx = find(strcmp(model_new.rxns,'biomass_reaction'));
% Get the flux value from the solution vector
temp_model = changeObjective(model_new, 'biomass_reaction');
FBAsol = optimizeCbModel(temp_model);
biomass_flux = FBAsol.x(biomass_idx);
fprintf('The biomass flux of context specific model is: %f\n', biomass_flux);


%%
save("C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/MKN28/MKN28_model.mat",'model_new')

%%
[minFlux, maxFlux] = fluxVariability(model_new);

%%
fluxTable = table(model_new.rxns, minFlux, maxFlux,'VariableNames', {'ReactionID', 'minFlux', 'maxFlux'});

%% Save the Table to a CSV File
% Specify the desired file name. The file will be saved in your current MATLAB folder.
csvFileName = 'C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/MKN28/MKN28_flux_variability.csv';
writetable(fluxTable, csvFileName);
disp(['Flux variability results saved to: ', csvFileName]);
%%
