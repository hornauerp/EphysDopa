function norm_mat = normalizeByBatch(mat,batch)
N_BATCH = unique(batch);
norm_mat = zeros(size(mat));
for b = 1:length(N_BATCH)
    batch_mat = mat(batch==b,:);
    norm_mat(batch==b,:) = normalize(batch_mat);%,'norm'
end