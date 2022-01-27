function [train_err,sf_acc,feat_imp,input_features] = single_feature_classification(nw_array,age_range,mat_idx,kf,TH,N_hyper)
arguments
    nw_array
    age_range (1,2) {isnumeric}
    mat_idx {isnumeric} = 1:2
    kf (1,1) {isnumeric} = 5
    TH (1,1) {isnumeric} = 3
    N_hyper (1,1) {isnumeric} = 100
end
nw_sel = []; sc_sel = [];
[input_mat,feature_names,mutation,batch,~] = generate_classifier_input(nw_array,age_range,mat_idx,nw_sel,sc_sel,TH);
input_features = vertcat(feature_names{mat_idx});

%% Single-feature classification
cv_split = cvpartition(mutation,'KFold',kf);
sf_acc = nan(kf,length(input_features));
train_err = nan(kf,length(input_features));
feat_imp = nan(length(input_features),5,kf);
tp = size(input_mat,2)/length(input_features);

%%
parfor f = 1:length(input_features)
    feature_input = input_mat(:,(f-1)*tp+1:f*tp);
    for k = 1:kf
        X_train = nan(sum(cv_split.training(k)),size(feature_input,2));%input(c.training(k),:);
        Y_train = mutation(cv_split.training(k));
        X_test = nan(sum(cv_split.test(k)),size(feature_input,2));%input(c.test(k),:);
        Y_test = mutation(cv_split.test(k));
        
        %Batch normalization
        b_train = 0; b_test = 0;
        for b = unique(batch)
            train_idx = batch==b & cv_split.training(k)';
            test_idx = batch==b & cv_split.test(k)';
            batch_mat = feature_input(train_idx,:);
            mu = mean(batch_mat);
            sd = std(batch_mat);
            X_train(b_train+1:sum(train_idx)+b_train,:) = (batch_mat-mu)./sd;
            X_test(b_test+1:sum(test_idx)+b_test,:) = (feature_input(test_idx,:)-mu)./sd;
            b_train = b_train+sum(train_idx);
            b_test = b_test+sum(test_idx);
        end
        if N_hyper > 0
            hyperparams = {'NumLearningCycles','MinLeafSize','MaxNumSplits','SplitCriterion','NumVariablesToSample'};
            t = templateTree('Reproducible',true);
            clf = fitcensemble(X_train,Y_train,'Method','Bag','OptimizeHyperparameters',hyperparams,'Learners',t, ...
                'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',N_hyper,'ShowPlots',false));
        else
            t = templateTree('Surrogate','on','MinLeafSize',1,'NumVariablesToSample','all','Reproducible',true);
            clf = fitcensemble(X_train,Y_train,'Method','Bag','NumLearningCycles',500,'Learners',t);
        end
        
        try
            feat_imp(f,:,k) = clf.oobPermutedPredictorImportance;
        catch
            feat_imp(f,:,k) = clf.predictorImportance;
        end
        train_err(k,f) = resubLoss(clf,'LossFun','classiferror');
        pred = predict(clf,X_test);
        sf_acc(k,f) = sum(pred'==Y_test)/length(Y_test);
    end
%     fprintf('Finished feature %i/%i',f,length(input_features))
end
