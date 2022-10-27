function new_cluster_assignment = reassign_cluster_ids(cluster_assignment,new_cluster_groups)
% Provide old cluster assignment and reassign them based on cluter ID
% groups
arguments
   cluster_assignment = isnumeric %Array of cluster assignments
   new_cluster_groups = iscell %Cell array of cluster id lists
end
bin_clust = arrayfun(@(x) ismember(cluster_assignment,new_cluster_groups{x})*x,1:length(new_cluster_groups),'un',0);
new_cluster_assignment = sum(vertcat(bin_clust{:}));