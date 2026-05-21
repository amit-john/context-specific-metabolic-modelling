%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clean code for MKN28
% Load cobratoolbox
clc;
addpath(genpath('C:/Users/USER/Documents/MATLAB/cobratoolbox'));
initCobraToolbox(0)
changeCobraSolver('gurobi','all')
changeCobraSolverParams('LP', 'feasTol', 1e-8);
addpath(genpath('C:/Users/USER/Desktop/Desktop/personal_docs/IITM_academics/5th_year/Final year Project/Making the MKN GSM model/SprintGapFiller-main/SprintGapFiller-main'));
%%
model=load('C:/Users/USER/Downloads/Recon3D (1).mat').Recon3D;

%Changing the objective function
biomassTargetID = 'BIOMASS_reaction'; 
bioIdx = find(strcmp(model.rxns, biomassTargetID)); % Use strcmp for exact match
initial_obj_index=find(model.c==1);
model.c(initial_obj_index)=0;
model.c(bioIdx)=1;


exchanges = {
    'EX_glc__D_e'
    'EX_cl_e'
    'EX_na1_e'
    'EX_ca2_e'
    'EX_k_e'
    'EX_so4_e'
    'EX_ncam_e'
    'EX_pnto__R_e'
    'EX_chol_e'
    'EX_fol_e'
    'EX_inost_e'
    'EX_pydxn_e'
    'EX_ribflv_e'
    'EX_thm_e'
    'EX_pi_e'
    'EX_arg__L_e'
    'EX_cyst__L_e'
    'EX_his__L_e'
    'EX_ile__L_e'
    'EX_leu__L_e'
    'EX_lys__L_e'
    'EX_met__L_e'
    'EX_phe__L_e'
    'EX_thr__L_e'
    'EX_trp__L_e'
    'EX_tyr__L_e'
    'EX_val__L_e'
    'EX_pchol_hs_e' 
    'EX_pe_hs_e'    
    'EX_ps_hs_e'    
    'EX_sphmyln_hs_e' 
    'EX_chsterol_e'
    'EX_o2_e'   
    'EX_hco3_e'
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
   -10.0 % EX_pchol_hs_e uptake rate
   -10.0 % EX_pe_hs_e uptake rate
   -10.0 % EX_ps_hs_e uptake rate
   -10.0 % EX_sphmyln_hs_e uptake rate
   -10.0 % EX_chsterol_e uptake rate
   -10.0 % EX_o2_e (typically a high negative number for unlimited air)
   -100.0 % EX_hco3_e
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
   -10.0 % EX_pchol_hs_e uptake rate
   -10.0 % EX_pe_hs_e uptake rate
   -10.0 % EX_ps_hs_e uptake rate
   -10.0 % EX_sphmyln_hs_e uptake rate
   -10.0 % EX_chsterol_e uptake rate
   -10.0 % EX_o2_e (typically a high negative number for unlimited air)
   -100.0 % EX_hco3_e
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
   -10.0 % EX_pchol_hs_e uptake rate
   -10.0 % EX_pe_hs_e uptake rate
   -10.0 % EX_ps_hs_e uptake rate
   -10.0 % EX_sphmyln_hs_e uptake rate
   -10.0 % EX_chsterol_e uptake rate
   -10.0 % EX_o2_e (typically a high negative number for unlimited air)
   -100.0 % EX_hco3_e
];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transcriptomic Data
t_data=readtable("C:/Users/USER/Downloads/GES1_control_ENTREZ_log2.csv","VariableNamingRule","preserve");  
%gene_identifer



% Shortening the t_dataset
mg_entrez = str2double(extractBefore(string(model.genes(:)), '_AT')); 
%Finding the matching indices of t_data with the model genes

[tf1, loc] = ismember(mg_entrez, t_data.GENE);  
% Shortening the t_dataset
idx_model = find(tf1);    % Just for checking      
idx_tdata = loc(tf1);

% Generating the expression structure
gene_expression = struct();
gene_expression.genes   = model.genes(idx_model);     
gene_expression.value   = t_data.GSM1923245(idx_tdata);
gene_expression.context = {'GES-1'};


biomassIdx = find(contains(model.rxns,'biomass','IgnoreCase',true));
%model.rxnNames(biomassIdx)
atpIdx=find(contains(model.rxns,'ATP','IgnoreCase',true));
%model.rxnNames(atpIdx)
context = gene_expression.context;


%
coreRxn = [biomassIdx;atpIdx];
[tf1, loc1] = ismember(model.rxns,exchanges);
idx_model_ids = find(tf1);    % Just for checking      
idx_exchanges = loc1(tf1);
for i=1:length(idx_model_ids)
    coreRxn(end+1)=idx_model_ids(i);
end

%

if ~isfield(model,'rules')
    model.rules = cell(size(model.grRules));
    for i = 1:numel(model.grRules)
        if isempty(model.grRules{i})
            model.rules{i} = '';
        else
            % Convert gene rules from grRules syntax to rules syntax
            tmp = model.grRules{i};
            tmp = strrep(tmp,' and ',' & ');
            tmp = strrep(tmp,' or ',' | ');
            model.rules{i} = tmp;
        end
    end
end


[RxnImp, Contexts] = GiniReactionImportance_rules(gene_expression, model, 90, 10, coreRxn);
RxnImp(isnan(RxnImp)) = 0;
con_id    = find(strcmp(Contexts, 'GES-1'));
curr_core = RxnImp(:, con_id);
weights   = zeros(numel(curr_core),1);
weights(curr_core < 1) = 1 - curr_core(curr_core < 1);
core      = find(curr_core >= 1);

% Making the model using sprintcore
model_new = sprintcore(model,core,1e-4,[],weights);

% Sanity check

% Addition
[model_new, exClID] = addExchangeRxn(model_new, {'cl_e'});
% Change bounds after creation using changeRxnBounds if default is not -1000 to 1000
model_new = changeRxnBounds(model_new, exClID, -1000, 'l'); % Set lower bound
model_new = changeRxnBounds(model_new, exClID, 1000, 'u'); % Set upper bound

% exClID will be a cell array, use char() to print nicely in fprintf
fprintf('Added chloride exchange: %s\n', char(exClID{1})); 

% Find index of new reaction (Optional)
Idx = find(contains(model_new.rxns, exClID{1}, 'IgnoreCase', true));

% Add calcium exchange reaction
[model_new, exCaID] = addExchangeRxn(model_new, {'ca2_e'});
model_new = changeRxnBounds(model_new, exCaID, -1000, 'l'); 
model_new = changeRxnBounds(model_new, exCaID, 1000, 'u'); 

fprintf('Added calcium exchange: %s\n', char(exCaID{1}));

% Verification assertions (These are correct)
assert(any(strcmp(model_new.mets,'cl_e')), 'cl_e metabolite not found');
assert(any(strcmp(model_new.mets,'ca2_e')), 'ca2_e metabolite not found');

%
FBAsol_test = optimizeCbModel(model_new);% biomass value
FBAsol_test.stat    
fprintf('biomass flux %f',FBAsol_test.x(bioIdx))



[model_new, addedRxns] = addTransportsForExchanges(model_new, exchanges, '_e', '_c', 'MEM_transports');

disp('Newly added transport reactions:');
disp(addedRxns);
%
bioIdx = find(strcmp(model_new.rxns, biomassTargetID));
model_new.rxns(find(model_new.c==1))
initial_obj_index=find(model_new.c==1);
model_new.c(initial_obj_index)=0;
model_new.c(bioIdx)=1;
%
FBAsol_test = optimizeCbModel(model_new);
FBAsol_test.f
FBAsol_test.stat

%
clc;
exIdx = find(findExcRxns(model_new));
model_new.lb(exIdx) = -1;        % no uptake by default
model_new.ub(exIdx) = 1000;

[tf2, loc2] = ismember(model_new.rxns,exchanges);
idx_model_ids_new = find(tf2);    % Just for checking      
idx_exchanges_new = loc2(tf2); 

model_new.lb(idx_model_ids_new) = normallimits_0h(idx_exchanges_new); 

%

biomassTargetID = 'BIOMASS_reaction';
bioIdx = find(strcmp(model_new.rxns, biomassTargetID));
FBAsol = optimizeCbModel(model_new);% biomass value
FBAsol.stat 
FBAsol.f
fprintf('biomass flux after constraining %f',FBAsol.x(bioIdx))
%
save('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/DDP_models/GES_1.mat','model_new')
%%
%
pathogen=load("C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Datasets/Helicobacter Pylori/iIT341.mat").iIT341;
pathogen_biomass_ind=findRxnIDs(pathogen, 'BIOMASS_HP_published');
pathogen_biomass_name = pathogen.rxns{pathogen_biomass_ind};
%


temp_model = changeObjective(pathogen, 'BIOMASS_HP_published');
FBAsol = optimizeCbModel(temp_model);
FBAsol.stat 
biomass_flux = FBAsol.x(pathogen_biomass_ind);
fprintf('The pathogen biomass flux in joint model is: %f\n', biomass_flux);
%%


[modelJoint] = createMultipleSpeciesModel( ...
    {pathogen}, ...                      % The pathogen model in a cell array
    {pathogen_biomass_name}, ...         % Its biomass reaction name
    'modelHost', model_new, ...          % Your human model as the host
    'nameTagsModels', {'pathogen'}, ...  % A tag for the pathogen's reactions
    'nameTagHost', 'human', ...          % A tag for the human's reactions
    'mergeGenesFlag', true ...           % Explicitly set to true to merge genes and GPRs
);

MOI=100;
pathogen_uptake_rxns = modelJoint.rxns(startsWith(modelJoint.rxns,'pathogenIEX_'));
upperbounds = modelJoint.ub(ismember(modelJoint.rxns, pathogen_uptake_rxns));
moi_100_ub = upperbounds * MOI;
% Now, apply these new bounds to your joint model
modelJoint = changeRxnBounds(modelJoint, pathogen_uptake_rxns, moi_100_ub, 'u');

%

clc;
exIdx = find(startsWith(modelJoint.rxns, 'EX_'));
modelJoint.lb(exIdx) = -1;        % no uptake by default
modelJoint.ub(exIdx) = 1000;

[tf2, loc2] = ismember(modelJoint.rxns,exchanges);
idx_model_ids_new = find(tf2);    % Just for checking      
idx_exchanges_new = loc2(tf2); 

modelJoint.lb(idx_model_ids_new) = normallimits_0h(idx_exchanges_new); 

pathogen_biomass_idx = find(strcmp(modelJoint.rxns, 'pathogenBIOMASS_HP_published'));
temp_model = changeObjective(modelJoint, 'pathogenBIOMASS_HP_published');
FBAsol = optimizeCbModel(temp_model);
FBAsol.stat 
biomass_flux = FBAsol.x(pathogen_biomass_idx);
fprintf('The pathogen biomass flux in joint model is: %f\n', biomass_flux);

biomass_idx = find(strcmp(modelJoint.rxns,'humanBIOMASS_reaction'));
temp_model = changeObjective(modelJoint, 'humanBIOMASS_reaction');
FBAsol = optimizeCbModel(temp_model);
FBAsol.stat 
biomass_flux = FBAsol.x(biomass_idx);
fprintf('The human biomass flux in joint model: %f\n', biomass_flux);

%
save('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/DDP_models/GES_1_0_intergrated.mat','modelJoint')



