%%
clc;
addpath('C:\Users\USER\cobratoolbox')
initCobraToolbox()

%%
t_data=readtable("C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Datasets/HP_time_data/GSE202165_countdata_norm.xls","VariableNamingRule","preserve");
%%
model=load('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Datasets/Recon3D_301/Recon3D_301/Consistant_Recon3D.mat').model;
%model=load('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Making the MKN GSM model/model_final.mat').model;

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

external_names = t_data.external_gene_name; 
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
clc;
genes_t_data={};
for i=1:length(row_indices)
    if row_indices{i}~='None'
        genes_t_data{end+1}=model.genes{i};
    end
end
genes_t_data = genes_t_data';
%%
matched_data.genes = genes_t_data;
%%
clc;
hpw1_data = matched_data(:, {'genes'});
%%
hpw1_data.expression=log2(matched_data.Mock1 + 1);

%%

transcriptome=matched_data(:, {'genes'});
transcriptome.zero=log2(matched_data.Mock1 + 1);
transcriptome.two=log2(matched_data.Hpwt_2h_1 + 1);
transcriptome.six=log2(matched_data.Hpwt_6h_1 + 1);





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
gene_expression.value=hpw1_data.expression;
gene_expression.genes=hpw1_data.genes;
gene_expression.context={'ACH-000758'};
%%

clc;
context = gene_expression.context;
bio_id = find(ismember(model.rxns,'biomass_reaction'));
ATP_id = find(ismember(model.rxns,'DM_atp_c_'));
coreRxn = [bio_id;ATP_id];
[RxnImp,Contexts] = GiniReactionImportance(gene_expression,model,90,10,coreRxn);
RxnImp(isnan(RxnImp))=0;
%%
model_cleansed = rmfield(model,{'metCharges','metFormulas','metSmiles','metInChIString',...
    'metKEGGID','metPubChemID','rxnNotes','rxnECNumbers','rxnReferences',...
    'rxnKEGGID','metCHEBIID','metPdMap','rxnCOG','rxnKeggOrthology'});

%%
con_id=find(strcmp(Contexts, 'ACH-000758'));
%%
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
save('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/hp_2hour/recon_mut2h_f.mat',"mf")
%%
clc;
mp = sprintcore(model_cleansed,core,1e-4,[],weights);
save('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/hp_2hour/recon_mut2h_sp.mat',"mp")
%%
clc;
msw = swiftcore(model_cleansed,core,weights,1e-10,1);
save('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/hp_2hour/recon_mut2h_sw.mat',"msw")

%%  %%%%%%%%%%%%%%%

%model_rec=load('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Datasets/Recon3D_301/Recon3D_301/Recon3D_301.mat').Recon3D;

%%
%biomass_rec_idx = find(contains(lower(model_rec.rxns), 'biomass'));
%ATP_rec_idx = find(contains(lower(model_rec.rxns), 'atp'));

% Display matched reaction IDs and names
%model_rec.rxns(biomass_idx)

%bio_rec_id = find(ismember(model_rec.rxns,'biomass_reaction'));
%ATP_rec_id = find(ismember(model_rec.rxns,'DM_atp_c_'));
%%
%clc;
%b_rxn = findBlockedReaction(model_rec,"FVA");
%%
%important_rxn={};
%for i=1:length(biomass_rec_idx)
%    important_rxn{end+1}=model_rec.rxns(biomass_rec_idx(i));
%end
%for i=1:length(ATP_rec_idx)
%    important_rxn{end+1}=model_rec.rxns(ATP_rec_idx(i));
%end 
%%
%important_rxn = cellfun(@char, important_rxn, 'UniformOutput', false);
%%
%clc;
%reactionsToRemove = setdiff(b_rxn,important_rxn);
%%
% Optionally remove them
%model_removed = removeRxns(model_rec, reactionsToRemove);
%%
%model=model_removed;
%% Generated consistant Recon3D
%save('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Datasets/Recon3D_301/Recon3D_301/Consistant_Recon3D.mat',"model")


%%
clc;
[contextModel, rxns_kept] = MEMETO(model,gene_expression,'ThH', 90,'ThL', 10,'taskFile', 'C:/Users/USER/Downloads/metabolic_tasks_richelle_2019.xlsx','tol', 1e-6);