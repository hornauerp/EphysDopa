function [color,mat,batch] = preparePCAinput(nw_array,MIN_AGE,MAX_AGE,featureGroup,TH)
[rec,~] = groupRecordingsByChip(nw_array,MIN_AGE,MAX_AGE);
[nw_mat,nw_prop,mutation,batch,treatment] = generateNWMat(rec);
[sc_mat,wf_mat,mc_mat,act_mat,sc_prop,wf_prop,mc_prop,act_prop,~,~] = generateSCMat(rec,'median',TH);
mat = eval(featureGroup);
mat = normalizeByBatch(mat,batch);
color = mutation*2+treatment+1;