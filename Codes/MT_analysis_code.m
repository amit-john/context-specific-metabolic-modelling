%%
clc;
addpath('C:\Users\USER\cobratoolbox')
initCobraToolbox()
%%
save('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/transcriptome.mat', 'transcriptome');
%%

clc;
model_sp0=load("C:\Users\USER\Desktop\Desktop\personal docs\IITM academics\5th year\Final year Project\Generated_models\hp_0hour\recon_uninfected_sp.mat").mp;
model_sp2=load("C:\Users\USER\Desktop\Desktop\personal docs\IITM academics\5th year\Final year Project\Generated_models\hp_2hour\recon_wt2h_sp.mat").mp;
model_sp6=load("C:\Users\USER\Desktop\Desktop\personal docs\IITM academics\5th year\Final year Project\Generated_models\hp_6hour\recon_wt6h_sp.mat").mp;

model_sw0=load("C:\Users\USER\Desktop\Desktop\personal docs\IITM academics\5th year\Final year Project\Generated_models\hp_0hour\recon_uninfected_sw.mat").msw;
model_sw2=load("C:\Users\USER\Desktop\Desktop\personal docs\IITM academics\5th year\Final year Project\Generated_models\hp_2hour\recon_wt2h_sw.mat").msw;
model_sw6=load("C:\Users\USER\Desktop\Desktop\personal docs\IITM academics\5th year\Final year Project\Generated_models\hp_6hour\recon_wt6h_sw.mat").msw;

model_f0=load("C:\Users\USER\Desktop\Desktop\personal docs\IITM academics\5th year\Final year Project\Generated_models\hp_0hour\recon_uninfected_f.mat").mf;
model_f2=load("C:\Users\USER\Desktop\Desktop\personal docs\IITM academics\5th year\Final year Project\Generated_models\hp_2hour\recon_wt2h_f.mat").mf;
model_f6=load("C:\Users\USER\Desktop\Desktop\personal docs\IITM academics\5th year\Final year Project\Generated_models\hp_6hour\recon_wt6h_f.mat").mf;
%%
clc;
genes=unique([model_sp0.genes;model_sp2.genes;model_sp6.genes;model_sw0.genes;model_sw2.genes;model_sw6.genes;model_f0.genes;model_f2.genes;model_f6.genes]);

%%
clc;
model_names = {model_sp0, model_sp2, model_sp6,model_sw0, model_sw2, model_sw6,model_f0,  model_f2,  model_f6};
genePresence = zeros(length(genes), length(model_names));

% Fill matrix: 1 if gene present in model, else 0
for i = 1:length(model_names)
    geneList = model_names{i}.genes;
    genePresence(:, i) = ismember(genes, geneList);
end
%%
clc;
modelNames = {'sp0','sp2','sp6','sw0','sw2','sw6','f0','f2','f6'};
geneTable = array2table(genePresence, 'VariableNames', modelNames, 'RowNames', genes);

%%
addpath(genpath("C:/Users/USER/AppData/Roaming/MathWorks/MATLAB Add-Ons/Collections/Uniform Manifold Approximation and Projection (UMAP)"))
%%
clc;
dataMatrix = table2array(geneTable(:,1:end));
[reduction, umapStruct] = run_umap(dataMatrix, 'n_neighbors', 15, 'min_dist', 0.1, 'metric', 'hamming');
scatter(reduction(:,1), reduction(:,2), 20, 'filled');
title('UMAP projection of gene presence matrix');
xlabel('UMAP 1');
ylabel('UMAP 2');
%%

% Now use DBSCAN on reduced data
labels = dbscan(reduction, 0.5, 5);  % tune eps and minPts


% Visualize with colors
gscatter(reduction(:,1), reduction(:,2), labels);
title('UMAP + DBSCAN Clustering of Genes');
%%
transcriptome=load("C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/transcriptome.mat").transcriptome;
Names_gene = transcriptome.genes;
exprMatrix = table2array(transcriptome(:, 2:end));
%%
exprMatrixNorm = (exprMatrix - mean(exprMatrix, 2)) ./ std(exprMatrix, 0, 2);
%%
[umapEmbedding, umapStruct] = run_umap(exprMatrixNorm, 'n_neighbors', 15, 'min_dist', 0.1, 'metric', 'euclidean');
%%
clusterLabels = dbscan(umapEmbedding, 0.5, 5);
%%
gscatter(umapEmbedding(:,1), umapEmbedding(:,2), clusterLabels);
title('Gene Clusters from Transcriptomic Data (UMAP + DBSCAN)');
xlabel('UMAP 1');
ylabel('UMAP 2');
%%
clc;
genes_model = geneTable.Properties.RowNames;
genes_expr = transcriptome.genes;

% Find common genes
[commonGenes, idxModel, idxExpr] = intersect(genes_model, genes_expr);

% Extract aligned cluster labels
labels_model = labels(idxModel);
labels_expr = clusterLabels(idxExpr);

n = length(labels_model);
jaccardVals = zeros(n,1);

for i = 1:n
    a = labels_model == labels_model(i);
    b = labels_expr == labels_expr(i);
    jaccardVals(i) = sum(a & b) / sum(a | b);
end

meanJaccard = mean(jaccardVals);
disp(['Mean Jaccard Similarity (on aligned genes): ', num2str(meanJaccard)]);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









%% Metabolic tasks
[taskReport, essentialReactions, taskStructure]=checkMetabolicTasks(model_sp0, 'C:/Users/USER/Downloads/metabolic_tasks_richelle_2019.xlsx');
%%
clc;
models={model_sp0,model_sp2,model_sp6,model_sw0,model_sw2,model_sw6,model_f0,model_f2,model_f6};
model_suffix={'sp0','sp2','sp6','sw0','sw2','sw6','f0','f2','f6'};
taskfile='C:/Users/USER/Downloads/metabolic_tasks_richelle_2019.xlsx';
T = table();
for i=1:length(models)
    disp("iteration")
    disp(i)
    [taskReport,~,~]=checkMetabolicTasks(models{i}, taskfile);
    if i==1
        METABOLISM_TYPE=taskReport(:,2);
        REACTION_NAME=taskReport(:,3);
        REACTION_DESCRIPTION=taskReport(:,4);
        T=table(METABOLISM_TYPE,REACTION_NAME,REACTION_DESCRIPTION);
    end
    T.(model_suffix{i}) = taskReport(:,5);
end
writetable(T, 'Metabolic_Task_Scores.csv');

%%
clc;
[taskReport1,~,~]=checkMetabolicTasks(models{1}, taskfile);
%%
[taskReport2,~,~]=checkMetabolicTasks(models{2}, taskfile);
%%
[taskReport3,~,~]=checkMetabolicTasks(models{3}, taskfile);
%%
[taskReport4,~,~]=checkMetabolicTasks(models{4}, taskfile);
%%
[taskReport5,~,~]=checkMetabolicTasks(models{5}, taskfile);
%%
[taskReport6,~,~]=checkMetabolicTasks(models{6}, taskfile);
%%
[taskReport7,~,~]=checkMetabolicTasks(models{7}, taskfile);
%%
[taskReport8,~,~]=checkMetabolicTasks(models{8}, taskfile);
%%
[taskReport9,~,~]=checkMetabolicTasks(models{9}, taskfile);
%%
clc;
T=table();
METABOLISM_TYPE=taskReport1(:,2);
REACTION_NAME=taskReport1(:,3);
REACTION_DESCRIPTION=taskReport1(:,4);
T=table(METABOLISM_TYPE,REACTION_NAME,REACTION_DESCRIPTION);
model_suffix={'sp0','sp2','sp6','sw0','sw2','sw6','f0','f2','f6'};
task={taskReport1,taskReport2,taskReport3,taskReport4,taskReport5,taskReport6,taskReport7,taskReport8,taskReport9};
for i=1:length(model_suffix)
    t=task{i};
    T.(model_suffix{i}) = t(:,5);
end

writetable(T, 'C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/Metabolic_Task_Scores.csv');

%% Obtaining the Biomass Objective Function after performing FBA

clc;
for i=1:length(models)
    old_mod=models{i};
    biomass_ind=findRxnIDs(old_mod, 'biomass_reaction');
    new_mod=changeObjective(old_mod,old_mod.rxns{biomass_ind});
    FBA=optimizeCbModel(new_mod);
    disp(['Biomass flux of model',model_suffix{i},': ',num2str(FBA.f)]);
end


%%  FLUX VARIABILITY Analysis

models = {model_sp0, model_sp2, model_sp6};
names = {'sp0', 'sp2', 'sp6'};

for i = 1:length(models)
    model = models{i};
    [minFlux, maxFlux] = fluxVariability(model, 100);
    T = table(model.rxns, minFlux, maxFlux);
    writetable(T, [names{i} '_FVA.csv']);
end

