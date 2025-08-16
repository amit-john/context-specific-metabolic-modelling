%%
clc;
addpath('C:\Users\USER\cobratoolbox')
initCobraToolbox()

%%
t_data=readtable("C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Datasets/MKN_Control/NT_epithelial_CPM.csv","VariableNamingRule","preserve");

%%
model=load('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Datasets/Recon3D_301/Recon3D_301/Consistant_Recon3D.mat').model;

%%
gene_names=readtable("C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Datasets/Expression data/symbol2entrez.csv");

%%
genes_integer = extractBefore(model.genes, '.');
%%
gene_names.ENTREZID = arrayfun(@num2str, gene_names.ENTREZID, 'UniformOutput', false);

%%
clc;
symbols = {};
clean_ids_num = genes_integer;  % This is a cell array of strings

for i = 1:length(clean_ids_num)
    current_id = clean_ids_num{i};
    found = false;
    disp(i)
    for j = 1:height(gene_names)
        entrez_id = gene_names.ENTREZID{j};  % safely extract from cell
        if isnumeric(entrez_id)
            entrez_id_str = num2str(entrez_id);
        else
            entrez_id_str = string(entrez_id);
        end

        if strcmp(current_id, entrez_id_str)
            val = gene_names.SYMBOL{j};
            found = true;
            break;
        end
    end

    if ~found
        val = 'None';
    end

    symbols{end+1} = val;
end

symbols = symbols';

%%

external_names = t_data.Gene; 
symbols = symbols(:); 
[isMatch, matchIdx] = ismember(symbols, external_names);
matched_row_indices = matchIdx(isMatch);


row_indices = cell(size(symbols));

for i = 1:length(symbols)
    disp(i)
    if matchIdx(i) ~= 0
        row_indices{i} = matchIdx(i);  
    else
        row_indices{i} = 'None';  
    end
end

%%
matched_data = t_data(matchIdx(isMatch), :);




%%
idx_model_matched = find(isMatch);                  % indices in model.genes that matched
genes_t_data = model.genes(idx_model_matched);      % exactly same count as matched_data
matched_data.genes = genes_t_data; 

%%

clc;
control_data = matched_data(:, {'genes'});
control_data.expression=log2(matched_data.NormalizedExpression + 1);


%%  Trying out the Sprint core algorihtm
initCobraToolbox(0)
changeCobraSolver('gurobi','all')
changeCobraSolverParams('LP', 'feasTol', 1e-8);
%%
clc;
addpath(genpath('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Making the MKN GSM model/SprintGapFiller-main/SprintGapFiller-main'));
%%

%Need to convert the omics data into a struct
clc;
gene_expression.value=control_data.expression;
gene_expression.genes=control_data.genes;
gene_expression.context={'ACH-000758'};

%%
clc;
context = gene_expression.context;
bio_id = find(ismember(model.rxns,'biomass_reaction'));
ATP_id = find(ismember(model.rxns,'DM_atp_c_'));
coreRxn = [bio_id;ATP_id];
[RxnImp,Contexts] = GiniReactionImportance(gene_expression,model,90,10,coreRxn);
RxnImp(isnan(RxnImp))=0;
model_cleansed = rmfield(model,{'metCharges','metFormulas','metSmiles','metInChIString',...
    'metKEGGID','metPubChemID','rxnNotes','rxnECNumbers','rxnReferences',...
    'rxnKEGGID','metCHEBIID','metPdMap','rxnCOG','rxnKeggOrthology'});

%%
con_id=find(strcmp(Contexts, 'ACH-000758'));

curr_core = RxnImp(:,con_id);
% weights are for those non-core reactions more the weights, less the
% probabilty of getting added
weights = zeros(numel(curr_core),1);
weights(curr_core<1)=1-curr_core(curr_core<1);
% core reactions are those with RxnImp more than 1
core = find(curr_core>=1);


%%
clc;
mf = fastcore(model_cleansed,core,1e-4);
save('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/Gastric_epithelial_control/recon_control_f.mat',"mf")

clc;
mp = sprintcore(model_cleansed,core,1e-4,[],weights);
save('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/Gastric_epithelial_control/recon_control_sp.mat',"mp")

clc;
msw = swiftcore(model_cleansed,core,weights,1e-10,1);
save('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/Gastric_epithelial_control/recon_control_sw.mat',"msw")

%%
taskfile='C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Datasets/metabolic_tasks_richelle_2019.xlsx';
%%
clc;
[taskReport_sprint,~,~]=checkMetabolicTasks(mp, taskfile);
%%
[taskReport_swift,~,~]=checkMetabolicTasks(msw, taskfile);
%%
[taskReport_fast,~,~]=checkMetabolicTasks(mf, taskfile);
%%
clc;
T=table();
METABOLISM_TYPE=taskReport_sprint(:,2);
REACTION_NAME=taskReport_sprint(:,3);
REACTION_DESCRIPTION=taskReport_sprint(:,4);
T=table(METABOLISM_TYPE,REACTION_NAME,REACTION_DESCRIPTION);
model_suffix={'sprintcore','swiftcore','fastcore'};
task={taskReport_sprint,taskReport_swift,taskReport_fast};
for i=1:length(model_suffix)
    t=task{i};
    T.(model_suffix{i}) = t(:,5);
end

writetable(T, 'C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/Metabolic_Control_epicells_Task_Scores.csv');

%% Biomass

models={mp,msw,mf};
clc;
for i=1:length(models)
    old_mod=models{i};
    biomass_ind=findRxnIDs(old_mod, 'biomass_reaction');
    new_mod=changeObjective(old_mod,old_mod.rxns{biomass_ind});
    FBA=optimizeCbModel(new_mod);
    disp(['Biomass flux of model',model_suffix{i},': ',num2str(FBA.f)]);
end

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
normallimits=[-11.66906086
-53.88109497
-0.865090088
-8.116109372
-0.831579896
-0.815798187
-0.819138584
-0.000112511
-2.438846412
-0.368933771
-0.014888307
-0.003815656
-0.013022203
-0.004119076
-0.018165856
-0.008841471
-0.000483039
-0.005390335
-0.45454
-0.010745626
-0.000310976
-0.047860926
-0.121089577
-3.64305E-05
-6.30517E-08
-6.30386E-08
-1.15637E-05
-7.73762E-07
-7.73762E-07
-1.5387E-06
-1.42811E-05
-0.000105533
-0.003540568
-0.181246369
-0.090849949
-1.816418229
-0.181646463
-0.091068454
-0.363853778
-0.363853778
-0.36333337
-0.091388705
-0.181606877
-0.181660291
-0.362502518
-0.035611791
-0.180629093
-0.364718395
-1000];

%% Trying to constrain the model


clc;
exIdx  = find(findExcRxns(mp));
exRxns = mp.rxns(exIdx);
clc;
commonRxns = intersect(exRxns, exchangerxns);

dmem_expression_indices=[];
for i=1:length(commonRxns)
    j = find(strcmp(commonRxns{i}, exchangerxns));   % or strcmpi for case-insensitive
    dmem_expression_indices = [dmem_expression_indices; j(:)];
end

dmem_normallb_constraints=[];
for i=1:length(dmem_expression_indices)
    dmem_normallb_constraints(end+1)=normallimits(dmem_expression_indices(i));
end
dmem_normallb_constraints=dmem_normallb_constraints';

mp.lb(exIdx) = -10;
mp.lb(dmem_expression_indices) = dmem_normallb_constraints;



%%


% 1) Get exchange indices for this model
exIdx = find(findExcRxns(mp));                 % indices into mp.rxns
[~, posInExIdx, posInDMEM] = intersect(mp.rxns(exIdx), exchangerxns, 'stable');
mp.lb(exIdx) = -10;                            % generic rich-media default
mp.lb(exIdx(posInExIdx)) = normallimits(posInDMEM);  % overwrite with DMEM caps

%%
clc;
models={mp,msw,mf};
for num=1:length(models)
    exIdx = find(findExcRxns(models{num}));                 % indices into mp.rxns
    [~, posInExIdx, posInDMEM] = intersect(models{num}.rxns(exIdx), exchangerxns, 'stable');
    models{num}.lb(exIdx) = -10;                            % generic rich-media default
    models{num}.lb(exIdx(posInExIdx)) = normallimits(posInDMEM);  % overwrite with DMEM caps
end

%%
clc;
for i=1:length(models)
    old_mod=models{i};
    biomass_ind=findRxnIDs(old_mod, 'biomass_reaction');
    new_mod=changeObjective(old_mod,old_mod.rxns{biomass_ind});
    FBA=optimizeCbModel(new_mod);
    disp(['Biomass flux of model (Constrained condition)',model_suffix{i},': ',num2str(FBA.f)]);
end

%%
nsp=models{1};
nsw=models{2};
nf=models{3};

%%
clc;
save('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/Gastric_epithelial_control/recon_control_constrained_f.mat',"nf")
save('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/Gastric_epithelial_control/recon_control_constrained_sw.mat',"nsw")
save('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/Gastric_epithelial_control/recon_control_constrained_sp.mat',"nsp")

%% Constraining the cancer models

sp0=load("C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/hp_0hour/recon_uninfected_sp.mat").mp;
sp2=load("C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/hp_2hour/recon_wt2h_sp.mat").mp;
sp6=load("C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/hp_6hour/recon_wt6h_sp.mat").mp;

sw0=load("C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/hp_0hour/recon_uninfected_sw.mat").msw;
sw2=load("C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/hp_2hour/recon_wt2h_sw.mat").msw;
sw6=load("C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/hp_6hour/recon_wt6h_sw.mat").msw;

f0=load("C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/hp_0hour/recon_uninfected_f.mat").mf;
f2=load("C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/hp_2hour/recon_wt2h_f.mat").mf;
f6=load("C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/hp_6hour/recon_wt6h_f.mat").mf;

%%
clc;
inf_models={sp0,sp2,sp6,sw0,sw2,sw6,f0,f2,f6};
inf_suf={'sp0','sp2','sp6','sw0','sw2','sw6','f0','f2','f6'};
clc;
for i=1:length(inf_models)
    old_mod=inf_models{i};
    biomass_ind=findRxnIDs(old_mod, 'biomass_reaction');
    new_mod=changeObjective(old_mod,old_mod.rxns{biomass_ind});
    FBA=optimizeCbModel(new_mod);
    disp(['Biomass flux of model (unconstrained condition)',inf_suf{i},': ',num2str(FBA.f)]);
end


%% Constraining cancer models

for num=1:length(inf_models)
    exIdx = find(findExcRxns(inf_models{num}));                 % indices into mp.rxns
    [~, posInExIdx, posInDMEM] = intersect(inf_models{num}.rxns(exIdx), exchangerxns, 'stable');
    inf_models{num}.lb(exIdx) = -5;                            % generic rich-media default
    inf_models{num}.lb(exIdx(posInExIdx)) = cancerlimits(posInDMEM);  % overwrite with DMEM caps
end
%%
uncommon=setdiff(exchangerxns,inf_models{1}.rxns(exIdx));
idx = find(contains(inf_models{1}.rxns, 'glc'));
%%

plot(,)
%%  FBA after constaining
disp("after")

for i=1:length(inf_models)
    old_mod=inf_models{i};
    biomass_ind=findRxnIDs(old_mod, 'biomass_reaction');
    new_mod=changeObjective(old_mod,old_mod.rxns{biomass_ind});
    FBA=optimizeCbModel(new_mod);
    disp(['Biomass flux of model (Constrained condition)',inf_suf{i},': ',num2str(FBA.f)]);
end

%%
m1=inf_models{1};
m2=inf_models{2};
m3=inf_models{3};
m4=inf_models{4};
m5=inf_models{5};
m6=inf_models{6};
m7=inf_models{7};
m8=inf_models{8};
m9=inf_models{9};

%%

save('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/hp_0hour/recon_uninfected_constrained_sp.mat',"m1")
save('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/hp_2hour/recon_wt2_constrained_sp.mat',"m2")
save('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/hp_6hour/recon_wt6_constrained_sp.mat',"m3")

save('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/hp_0hour/recon_uninfected_constrained_sw.mat',"m4")
save('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/hp_2hour/recon_wt2_constrained_sw.mat',"m5")
save('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/hp_6hour/recon_wt6_constrained_sw.mat',"m6")

save('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/hp_0hour/recon_uninfected_constrained_f.mat',"m7")
save('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/hp_2hour/recon_wt2_constrained_f.mat',"m8")
save('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/hp_6hour/recon_wt6_constrained_f.mat',"m9")

%%



%%
normal_f=load('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/Gastric_epithelial_control/recon_control_f.mat').mf;
normal_sp=load('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/Gastric_epithelial_control/recon_control_sp.mat').mp;
normal_sw=load('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/Gastric_epithelial_control/recon_control_sw.mat').msw;

n_models={normal_f,normal_sw,normal_sp};
n_suf={'normal_f','normal_sw','normal_sp'};
%%
for i=1:length(n_models)
    old_mod=n_models{i};
    biomass_ind=findRxnIDs(old_mod, 'biomass_reaction');
    new_mod=changeObjective(old_mod,old_mod.rxns{biomass_ind});
    FBA=optimizeCbModel(new_mod);
    disp(['Biomass flux of model (Constrained condition)',n_suf{i},': ',num2str(FBA.f)]);
end

%%
for num=1:length(n_models)
    exIdx = find(findExcRxns(n_models{num}));                 % indices into mp.rxns
    [~, posInExIdx, posInDMEM] = intersect(n_models{num}.rxns(exIdx), exchangerxns, 'stable');
    n_models{num}.lb(exIdx) = -5;                            % generic rich-media default
    n_models{num}.lb(exIdx(posInExIdx)) = cancerlimits(posInDMEM);  % overwrite with DMEM caps
end