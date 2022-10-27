function [group_idx, group_count, group_comb] = groupby_multiple_idx(cluster_idx_array)
% Group by several idx values that are concatenated in a cell array
id_arrays = cellfun(@(x) unique(x),cluster_idx_array,'un',0);
group_comb = combvec(id_arrays{:});
id_matrix = vertcat(cluster_idx_array{:});
group_count = arrayfun(@(x) sum(all(id_matrix == group_comb(:,x))),1:size(group_comb,2));
group_idx = arrayfun(@(x) all(id_matrix == group_comb(:,x))*x,1:size(group_comb,2),'un',0);
group_idx = sum(vertcat(group_idx{:}));