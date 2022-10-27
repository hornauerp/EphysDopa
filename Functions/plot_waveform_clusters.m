function plot_waveform_clusters(norm_mat,idx)
N_clust = max(idx);
c_map = othercolor('Set19',N_clust);
figure;
tiledlayout('flow')
for i = 1:max(idx)
    nexttile
    plot(norm_mat(idx==i,:)','Color',[c_map(i,:) 0.1],'LineWidth',0.01)
    hold on
    plot(mean(norm_mat(idx==i,:)),'k')
    ylim([-1 1])
    title("Cluster " + i)
end
nexttile
for i = 1:max(idx)
   plot(mean(norm_mat(idx==i,:)),'Color',c_map(i,:))
   hold on
   ylim([-1 1])
end
end