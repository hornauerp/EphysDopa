function plot_cluster_scatter(mapped_data,idx)
figure('Color','w');
if size(mapped_data,2) == 2
    scatter(mapped_data(:,1),mapped_data(:,2),5,idx,'filled')
else
    scatter3(mapped_data(:,1),mapped_data(:,2),mapped_data(:,3),5,idx,'filled')
end
N_clust = double(max(idx));
c_map = othercolor('Set19',N_clust);
colormap(c_map)
cb = colorbar; cb.Title.String = "Cluster Index";
cb.Ticks = linspace(cb.Limits(1) + (N_clust-1)/(2*N_clust), cb.Limits(2) - (N_clust-1)/(2*N_clust), N_clust); cb.TickLabels = 1:N_clust;
end