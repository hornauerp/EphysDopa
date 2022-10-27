function [clf,train_acc] = create_classifier(X_train,Y_train,alg,N_hyper)
if N_hyper > 0
    switch alg
        case 'svm'
            clf = fitcsvm(X_train,Y_train,'OptimizeHyperparameters','all','HyperparameterOptimizationOptions',...
                struct('AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',N_hyper,'ShowPlots',false,...
                'Verbose',0));
        case 'cnb'
            clf = fitcnb(X_train,Y_train,'OptimizeHyperparameters','all','HyperparameterOptimizationOptions',...
                struct('AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',N_hyper,'ShowPlots',false,...
                'Verbose',0));
        case 'knn'
            clf = fitcknn(X_train,Y_train,'OptimizeHyperparameters','all','HyperparameterOptimizationOptions',...
                struct('AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',N_hyper,'ShowPlots',false,...
                'Verbose',0));
        case 'rf'
            hyperparams = {'NumLearningCycles','MinLeafSize','MaxNumSplits','SplitCriterion','NumVariablesToSample'};
            t = templateTree('Reproducible',true);
            clf = fitcensemble(X_train,Y_train,'Method','Bag','OptimizeHyperparameters',hyperparams,'Learners',t, ...
                'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',N_hyper,'ShowPlots',false,...
                'Verbose',0));
    end
    train_acc = clf.HyperparameterOptimizationResults.MinObjective;
else
    switch alg
        case 'svm'
            clf = fitcsvm(X_train,Y_train);
        case 'cnb'
            clf = fitcnb(X_train,Y_train);
        case 'knn'
            clf = fitcknn(X_train,Y_train);
        case 'rf'
            t = templateTree('Surrogate','on','MinLeafSize',1,'NumVariablesToSample','all','Reproducible',true);
            clf = fitcensemble(X_train,Y_train,'Method','Bag','NumLearningCycles',500,'Learners',t);
    end
    train_acc = 1-resubLoss(clf,'LossFun','classiferror');
end