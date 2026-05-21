%%  Loading Cobratool box
clc;
addpath(genpath('C:/Users/USER/Documents/MATLAB/cobratoolbox'));
initCobraToolbox(0)
changeCobraSolver('gurobi','all')
changeCobraSolverParams('LP', 'feasTol', 1e-8);
addpath(genpath('C:/Users/USER/Desktop/Desktop/personal_docs/IITM_academics/5th_year/Final year Project/Making the MKN GSM model/SprintGapFiller-main/SprintGapFiller-main'));

%%
%model1=load('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Models/Host-pathogen-interaction-models-gastric/Host-pathogen-interaction-models-gastric/Infected host model/Gastric/tissueModel_hp_0.mat').tissueModel_hp_0;

model=load('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/Normal_GES_pylori/ges1_control.mat').model_new;
model1=load('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/Normal_GES_pylori/ges1_hp.mat').model_new;
%%
filepath = fullfile('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/Normal_GES_pylori/', 'ges1_control_sbml.xml');
writeCbModel(model, 'format', 'sbml', 'fileName', filepath);
%%
filepath = fullfile('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/Final year Project/Generated_models/Normal_GES_pylori/', 'ges1_hp_sbml.xml');
writeCbModel(model1, 'format', 'sbml', 'fileName',filepath);
%%
model3=load('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/DDP_models/MKN28_hp_6_non_intergrated.mat').model_new;
filepath = fullfile('C:/Users/USER/Desktop/Desktop/personal docs/IITM academics/5th year/DDP_models/', 'MKN28_hp_6_non_intergrated.xml');
writeCbModel(model3, 'format', 'sbml', 'fileName',filepath);