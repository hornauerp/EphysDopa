function [mat,vars,vars_idx,age,cl,batch] = generateInputMatrix(rec,MET,TH)
nw_array = reshape([rec{:}],1,[]);
nw_idx = [16:17,19:30];%11,12,
nw_vars = properties(WholeNetwork);
nw_prop = nw_vars(nw_idx);
tmp_vars = properties(Unit);
tmp_idx = [4,11:26,45:47]; tmp_idx([7,9,11]) = [];
wf_prop = [1:10];
mc_prop = [];
act_prop = [11:17];
sc_prop = tmp_vars(tmp_idx);
vars = [nw_prop' sc_prop'];
nw_mat = zeros(numel(nw_array),length(nw_prop));
sc_mat = zeros(numel(nw_array),length(sc_prop));
for p = 1:length(nw_prop)
    nw_mat(:,p) = [nw_array.(nw_prop{p})];
end
for j = 1:numel(sc_prop)
    sc_mat(:,j) = arrayfun(@(x) mean(rmoutliers([x.Units.(sc_prop{j})],MET,'ThresholdFactor',TH),'omitnan'),nw_array)';
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