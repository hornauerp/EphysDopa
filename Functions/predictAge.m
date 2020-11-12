function [mdl, pred_age, Y_test, score_matrix, weeks, batches, pheno, ids,fi] =  predictAge(rec,mat,batch,cellline,N_ITER,AGE_PRED_MET,...
    AGE_PCA,AGE_PCA_IND,EXC,INC,FEAT_IMP)
assert(isempty(EXC)|isempty(INC),'Either specify features to be included or excluded, not both!')
fi = [];
tp = length(rec{1}); %Recorded timepoints

if isempty(EXC) && isempty(INC)
    XI = 0;
    II = 0;
elseif isempty(EXC)
    if length(AGE_PCA_IND) > length(INC)
        AGE_PCA_IND = 1:length(INC);
    end
    II = 1;
    XI = 0;
elseif isempty(INC)
    if length(AGE_PCA_IND) > size(mat,2)/tp-length(EXC)
        AGE_PCA_IND = size(mat,2)/tp-length(EXC);
    end
    XI = 1;
    II = 0;
end

if FEAT_IMP
    t = templateTree('Surrogate','on','MinLeafSize',1,'NumVariablesToSample','all');
    fi = [];
    AGE_PCA = 0;
    N_CYCLES = 100;
end
num_rec = length(rec);
Y_test= zeros(N_ITER,num_rec);
pred_age = zeros(N_ITER,num_rec);
score_matrix = []; 
weeks = [];
batches = [];
pheno = [];
ids = [];

    
for i = 1:N_ITER
    rmd = rem(length(rec),tp);
    sampling = [ceil(randperm(num_rec-rmd)/floor(num_rec/tp)),randi(tp,[1 rmd])];
%     ind = (sampling-1).*(size(mat,2)/tp)+1;
    mat_in = zeros(size(mat,1),size(mat,2)/tp);
    for n = 1:size(mat,1)
%         mat_in(n,:) = mat(n,ind(n):ind(n)+(size(mat,2)/tp)-1);
          mat_in(n,:) = mat(n,linspace(sampling(n),size(mat,2)+sampling(n)-tp,(size(mat,2)/tp)));
    end
    if XI == 1
        mat_in(:,EXC) = []; 
    elseif II == 1
        mat_in = mat_in(:,INC);
    end
    for ii = 1:size(mat,1)
        X_train = mat_in; 
        X_test = X_train(ii,:);
        X_train(ii,:) = [];
        if AGE_PCA == 1
            mu = mean(X_train);
            sd = std(X_train);
            X_train = (X_train-mu)./sd;
            [age_coeff,age_score,~,~,~] = pca(X_train,...
            'VariableWeights','variance');
            X_train = age_score(:,AGE_PCA_IND);
            X_test = (X_test-mu)./sd;
            X_test = X_test*age_coeff(:,AGE_PCA_IND);
        else
            mu = mean(X_train);
            sd = std(X_train);
            X_train = (X_train-mu)./sd;
            X_test = (X_test-mu)./sd;
        end
        Y_train = sampling*7;
        Y_test(i,ii) = sampling(ii)*7;
        Y_train(ii) = [];
        if FEAT_IMP
            if AGE_PRED_MET == "REG"
                mdl = fitrensemble(X_train,Y_train,'Method','Bag','NumLearningCycles',N_CYCLES,'Learners',t);
            elseif AGE_PRED_MET == "CLF"
                mdl = fitcensemble(X_train,Y_train,'Method','Bag','NumLearningCycles',N_CYCLES,'Learners',t);
            end
            fi = [fi; mdl.oobPermutedPredictorImportance];
        else
            if AGE_PRED_MET == "REG"
                mdl = fitrgp(X_train,Y_train,'BasisFunction','constant','FitMethod','exact','PredictMethod','exact');
            elseif AGE_PRED_MET == "CLF"
                mdl = fitcknn(X_train,Y_train,'Distance','cityblock','IncludeTies',true,'NumNeighbors',1);
            end
        end
%         mdl = fitcnb(X_train,Y_train,'DistributionNames','kernel','Width',0.31);
%         mdl = fitcdiscr(X_train,Y_train,'Delta',0.009393,'Gamma',0.61216);
%         mdl = stepwiseglm(X_train,Y_train);
%                 mdl = fitrsvm(X_train,Y_train,'KernelFunction','linear');
        pred_age(i,ii) = predict(mdl,X_test);
        score_matrix = [score_matrix; X_test];%X_train; 
        weeks = [weeks; sampling(ii)*7];%Y_train'; 
        b_train = batch;
        b_test = batch(ii);
        b_train(ii) = [];
        batches = [batches; b_test];
        pheno = [pheno cellline(ii)];
        ids = [ids string([rec{ii}(1).ChipID])];
    end
end
end