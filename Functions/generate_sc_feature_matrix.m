function feat_mat = generate_sc_feature_matrix(Unit_array,sc_prop,norm_idx)
feat_mat = zeros(length(sc_prop),length([Unit_array.ID]));
for f = 1:length(sc_prop)
    feat_vals = [Unit_array.(sc_prop(f))];
    feat_vals(isnan(feat_vals)) = 0;
    feat_mat(f,:) = feat_vals;
end
if norm_idx
    feat_mat = normalize(feat_mat');
end