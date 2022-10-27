function [input_mat,feature_names,mutation,batch,treatment] = generate_classifier_input(nw_array,age_range,mat_idx,nw_sel,sc_sel,TH)

MIN_AGE = age_range(1);
MAX_AGE = age_range(2);
MET = 'median';
PCA_IND = 0; %If raw feature values (0) or PCA values (1) should be used
[rec,~] = groupRecordingsByChip(nw_array,MIN_AGE,MAX_AGE);

[nw_mat,nw_prop,mutation,batch,treatment] = generateNWMat(rec,PCA_IND,nw_sel);

[sc_mat,wf_mat,~,act_mat,sc_prop,~,~,~,~,~] = generateSCMat(rec,MET,TH,PCA_IND,sc_sel);

% sc_mat = [wf_mat act_mat];
mats = {nw_mat sc_mat};
input_mat = [mats{mat_idx}]; input_mat(isnan(input_mat)) = 0;
feature_names = {nw_prop; sc_prop}; 