function [train_acc,feature_names,test_acc,mutation] = classify_genotypes(nw_array,age_range,mat_idx,nw_sel,sc_sel,alg,kf,TH,N_hyper)
arguments
   nw_array
   age_range (1,2) {isnumeric}
   mat_idx {isnumeric} = 1:2
   nw_sel = []
   sc_sel = []
   alg = 'rf'
   kf (1,1) {isnumeric} = 5
   TH (1,1) {isnumeric} = 3
   N_hyper (1,1) {isnumeric} = 100
end

[input_mat,feature_names,mutation,batch] = generate_classifier_input(nw_array,age_range,mat_idx,nw_sel,sc_sel,TH);

%% Main loop
test_acc = nan(1,kf);
train_acc = nan(1,kf);
cv_split = cvpartition(mutation,'KFold',kf);
for k = 1:kf
    X_train = nan(sum(cv_split.training(k)),size(input_mat,2));%input(c.training(k),:);
    Y_train = mutation(cv_split.training(k));
    X_test = nan(sum(cv_split.test(k)),size(input_mat,2));%input(c.test(k),:);
    Y_test = mutation(cv_split.test(k));
    
    %Batch normalization
    b_train = 0; b_test = 0;
    for b = unique(batch)
        train_idx = batch==b & cv_split.training(k)';
        test_idx = batch==b & cv_split.test(k)';
        batch_mat = input_mat(train_idx,:);
        mu = mean(batch_mat);
        sd = std(batch_mat);
        X_train(b_train+1:sum(train_idx)+b_train,:) = (batch_mat-mu)./sd;
        X_test(b_test+1:sum(test_idx)+b_test,:) = (input_mat(test_idx,:)-mu)./sd;
        b_train = b_train+sum(train_idx);
        b_test = b_test+sum(test_idx);
    end
    [clf,train_acc_k] = create_classifier(X_train,Y_train,alg,N_hyper);
    train_acc(k) = train_acc_k;
    pred = predict(clf,X_test);
    test_acc(k) = sum(pred'==Y_test)/length(Y_test);
end