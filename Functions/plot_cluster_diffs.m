function plot_cluster_diffs(value_array,ids)
figure('Color','w');
diff_mat = value_array(:,:,ids(1))./sum(sum(value_array(:,:,ids(1)))) - value_array(:,:,ids(2))./sum(sum(value_array(:,:,ids(2))));
imagesc(diff_mat)
max_val = max(max(abs(diff_mat)));
caxis([-max_val max_val])
colormap(othercolor('RdBu9',100))
colorbar