clc
close
clear
%% Preparation of the parallel computations. Comment out this section if the workers are ready.
% p = gcp();
% parfor i = 1 : p.NumWorkers
%     changeCobraSolver('ibm_cplex');
% end

%% Loading models prepared for the minimal medium / low oxygen condition
load iYL1228_ML.mat     %#1
model1_ML = model;
load iJN1463_ML.mat     %#2
model2_ML = model;
load iML1515_ML.mat     %#3
model3_ML = model;
load STM_v1_0_ML.mat    %#4
model4_ML = model;
load iSDY_1059_ML.mat   %#5
model5_ML = model;
load iPC815_ML.mat      %#6
model6_ML = model;
% ML_models contains all loaded ML models as a single variable.
ML_models = {model1_ML, model2_ML, model3_ML, model4_ML, model5_ML, model6_ML};

%% Loading models prepared for the rich medium / high oxygen condition
load iYL1228_RH.mat     %#1
model1_RH = model;
load iJN1463_RH.mat     %#2
model2_RH = model;
load iML1515_RH.mat     %#3
model3_RH = model;
load STM_v1_0_RH.mat    %#4
model4_RH = model;
load iSDY_1059_RH.mat   %#5
model5_RH = model;
load iPC815_RH.mat      %#6
model6_RH = model;
% RH_models contains all loaded RH models as a single variable.
RH_models = {model1_RH, model2_RH, model3_RH, model4_RH, model5_RH, model6_RH};
% modelNames contains the names of the models
modelNames = {'iYL1228', 'iJN1463', 'iML1515', 'STM_v1_0', 'iSDY_1059', 'iPC815'};

Counter = 1;
for i = 1 : length(modelNames)
    % LetalSetsName is a character variable representing the name of the
    % file that contains the lethal reaction sets of the corresponding
    % model.
    LetalSetsName = [modelNames{i}, '_RH_Results_4_curated_geneAssoc'];
   
    MainTarget = load([modelNames{i}, '_RH']);
    MainTarget = MainTarget.model;
    load(LetalSetsName) % Loading the lethal sets corresponding to the MainTarget
   
    for j = 1 : length(setdiff(1 : length(modelNames), i))
        C = combntns( setdiff(1 : length(modelNames), i) , j ); % Computes all possible combinations of a set of values
        for k = 1 : length(C(:, 1))
            PathogenesID = logical(dec2bin(0 : 2^(length(C(k, :)))-1)-'0');
            for m = 1 : length(PathogenesID(:, 1))
                for n = 1 : length(PathogenesID(m, :))
                    if PathogenesID(m, n)
                        Models{n} = RH_models{C(k,n)};
                    elseif ~PathogenesID(m, n)
                        Models{n} = ML_models{C(k,n)};
                    end
                end
                ModelsToPrsv = Models(~PathogenesID(m, :));
                rxns = MainTarget.rxns;
                for o = 1 : length(ModelsToPrsv)
                    mod = ModelsToPrsv{o};
                    rxns = intersect(rxns, mod.rxns);
                end
                ISave = [];
                for p = 1 : length(LethalSets)
                    I = FindSelectiveCases(Models, LethalSets{p}, rxns, PathogenesID(m, :));
                    nSelective = sum(I == 1);
                    nNonApplicable = sum(I == -1);
                    nCommon = sum(I == 0);
                    Results1{Counter, p} = [nSelective, nNonApplicable, nCommon];
                    ISave = {ISave, I};
                end
                Results2{Counter, 1} = C(k, :);
                Results3{Counter, 1} = PathogenesID(m, :);
                Results4{Counter, 1} = ISave;
                Counter = Counter + 1
            end
        end
    end
end

