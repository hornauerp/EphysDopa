function [clf,train_acc,pred,test_mutation,test_treatment] = evaluate_intervention(train_array,test_array,age_range,mat_idx,nw_sel,sc_sel,alg,TH,N_hyper)
arguments
   train_array
   test_array
   age_range (1,2) {isnumeric}
   mat_idx {isnumeric} = 1:2
   nw_sel = []
   sc_sel = []
   alg = 'rf'
   TH (1,1) {isnumeric} = 3
   N_hyper (1,1) {isnumeric} = 100
end

[train_mat,feature_names,train_mutation,train_batch,~] = generate_classifier_input(train_array,age_range,mat_idx,nw_sel,sc_sel,TH);
norm_train = normalizeByBatch(train_mat,train_batch);

[test_mat,feature_names,test_mutation,test_batch,test_treatment] = generate_classifier_input(test_array,age_range,mat_idx,nw_sel,sc_sel,TH);
norm_test = normalizeByBatch(test_mat,test_batch);

[clf,train_acc] = create_classifier(norm_train,train_mutation,alg,N_hyper);

pred = predict(clf,norm_test);