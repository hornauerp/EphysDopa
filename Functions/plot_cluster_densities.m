function value_array = plot_cluster_densities(reduction,cluster_idx_array,cluster_name_array,n_bins)
% Plot

id_arrays = cellfun(@(x) unique(x),cluster_idx_array,'un',0);
group_comb = combvec(id_arrays{:});
id_matrix = vertcat(cluster_idx_array{:});
value_array = zeros(n_bins,n_bins,size(group_comb,2));

figure('Color','w');
tiledlayout('flow','TileSpacing','compact')
for i = 1:size(group_comb,2)
    plot_idx = all(id_matrix == group_comb(:,i));
    nexttile
    data = [reduction(plot_idx,1),reduction(plot_idx,2)];
    values = hist3(data,[n_bins n_bins]);
    value_array(:,:,i) = values;
%     values = values.'./max(max(values));
    imagesc(values)
    title(cluster_name_array + " " + group_comb(:,i)')
    xticks([])
    yticks([])
end