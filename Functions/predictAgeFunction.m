function [age,pred_weeks,pI,vars,cl_train] = predictAgeFunction(nw_array,vars2use,age_range,inclusion_array,exclusion_array)
%% Training selection
MIN_AGE = age_range(1);
MAX_AGE = age_range(2);
N_groups = max([length(inclusion_array),length(exclusion_array)]);
MET = 'median';
TH = 3;
% vars2use(vars2use==2) = [2,4]; %1:nw,2:wf,3:mc,4:act

batches = unique([nw_array.PlatingDate]);

inclusion = inclusion_array{1};%{{'PlatingDate',batches(batchIDs)},{'Mutation','a53t'}}; % Cell array of cells; first value property, rest values
exclusion = exclusion_array{1};%{{'DIV',28,9}};%,{'Mutation','a53t'}};,{'Treatment','LNA7'}

norm_inclusion = {{'PlatingDate',batches}}; norm_exclusion = {};
norm_array = filterObjectArray(nw_array,norm_inclusion,norm_exclusion);
[norm_rec,tp] = groupRecordingsByChip(norm_array,MIN_AGE,MAX_AGE);
[norm_table,~,vars_idx,~,~,batch_norm] = generateInputMatrix(norm_rec,MET,TH);
norm_table = norm_table(:,ismember(vars_idx,vars2use));
batch_mu = nan(max(batch_norm),size(norm_table,2));
batch_sd = nan(max(batch_norm),size(norm_table,2));

nw_train = filterObjectArray(nw_array,inclusion,exclusion);

[rec,tp] = groupRecordingsByChip(nw_train,MIN_AGE,MAX_AGE);
[train_table,~,vars_idx,age,cl_train,batch_train] = generateInputMatrix(rec,MET,TH);
train_table = train_table(:,ismember(vars_idx,vars2use));
vars = train_table.Properties.VariableNames;
for i = unique(batch_train)
    batch_mu(i,:) = mean(norm_table{batch_train==i,:});
    batch_sd(i,:) = std(norm_table{batch_train==i,:});
    train_table{batch_train==i,:} = (train_table{batch_train==i,:}-batch_mu(i,:))./batch_sd(i,:);
end

%% Regression training
t = templateTree('Surrogate','on','MinLeafSize',1,'NumVariablesToSample','all');
clf = fitrensemble(train_table,age,'Method','Bag','NumLearningCycles',250,'Learners',t);
pI = clf.oobPermutedPredictorImportance;
cvmodel = crossval(clf,'leaveout','on');
wt_pred = cvmodel.kfoldPredict;
wt_y = cvmodel.Y;

if N_groups > 1
    %% Test prediction
    % Normal comparison
    inclusion = inclusion_array{2};%{{'PlatingDate',batches(batchIDs)},{'Mutation','a53t'}};
    exclusion = exclusion_array{2};%{{'DIV',28,9},{'Treatment','LNA7'},{'Mutation','wt'}};
    nw_test = filterObjectArray(nw_array,inclusion,exclusion);
    [rec,~] = groupRecordingsByChip(nw_test,MIN_AGE,MAX_AGE);
    [test_table,~,vars_idx,mut_y,cl,batch_test] = generateInputMatrix(rec,MET,TH);
    test_table = test_table(:,ismember(vars_idx,vars2use));
    for i = unique(batch_test)
        test_table{batch_test==i,:} = (test_table{batch_test==i,:}-batch_mu(i,:))./batch_sd(i,:);
    end
    mut_pred = predict(clf,test_table);
else
    mut_pred = [];
end
%% Treatment prediction
if N_groups == 3
    % Normal prediction
    inclusion = inclusion_array{3};%{{'PlatingDate',batches(3)}};
    exclusion = exclusion_array{3};%{{'DIV',28,9},{'Treatment','-','ntLNA'}};
    nw_treat = filterObjectArray(nw_array,inclusion,exclusion);
    [rec,~] = groupRecordingsByChip(nw_treat,MIN_AGE,MAX_AGE);
    [treatment_table,vars,vars_idx,treatment_y,cl,batch_treat] = generateInputMatrix(rec,MET,TH);
    treatment_table = treatment_table(:,ismember(vars_idx,vars2use));
    for i = unique(batch_treat)
        treatment_table{batch_treat==i,:} = (treatment_table{batch_treat==i,:}-batch_mu(i,:))./batch_sd(i,:);
    end
    
    treatment_pred = predict(clf,treatment_table);
else
    treatment_pred = [];
    
end
%% Calculate Deviations/Distributions
pred_weeks = cell(1,tp);
for i = 1:tp
    wt_pred_week = wt_pred(wt_y+1==i*7);
    if N_groups>1
        mut_pred_week = mut_pred(mut_y+1==i*7)';
    else
        mut_pred_week = [];
    end
    if N_groups == 3
        treat_pred_week = treatment_pred(treatment_y+1==i*7)';
    else
        treat_pred_week = [];
    end
    pred_weeks{i} = {wt_pred_week mut_pred_week treat_pred_week}; pred_weeks{i} = pred_weeks{i}(1:N_groups);
end

