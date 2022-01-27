function [mat,vars,vars_idx,age,cl,batch] = generateInputMatrix(rec,MET,TH)
nw_array = reshape([rec{:}],1,[]);
nw_idx = [17:18,20:31];%11,12,
nw_vars = properties(WholeNetwork);
nw_prop = nw_vars(nw_idx);
tmp_vars = properties(Template);
tmp_idx = [4,11:15,17,19,21:27,41:48];
wf_prop = [1:8,10:11];
mc_prop = [9,16:20];
act_prop = [12:15,21:23];
sc_prop = tmp_vars(tmp_idx);
vars = [nw_prop' sc_prop'];
nw_mat = zeros(numel(nw_array),length(nw_prop));
sc_mat = zeros(numel(nw_array),length(sc_prop));
for p = 1:length(nw_prop)
    nw_mat(:,p) = [nw_array.(nw_prop{p})];
end
for j = 1:numel(sc_prop)
    sc_mat(:,j) = arrayfun(@(x) mean(rmoutliers([x.Templates.(sc_prop{j})],MET,'ThresholdFactor',TH),'omitnan'),nw_array)';
end
batches = unique([nw_array.PlatingDate]);
batch = arrayfun(@(x) find(x.PlatingDate==batches),nw_array);
% sc_mat = normalizeByBatch(sc_mat,batch);
mat = [nw_mat sc_mat];
% mat = normalizeByBatch(mat,batch);
age = [nw_array.DIV];
cl = [nw_array.CellLine];
vars_idx = nan(1,size(mat,2));
vars_idx(1:length(nw_prop)) = 1;
vars_idx(length(nw_prop)+wf_prop) = 2;
vars_idx(length(nw_prop)+mc_prop) = 3;
vars_idx(length(nw_prop)+act_prop) = 4;
mat = array2table(mat,'VariableNames',vars);
end