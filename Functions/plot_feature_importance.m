function plot_feature_importance(sf_acc,feat_imp,input_features)
fi = mean(feat_imp,3); fi(fi<=0)=nan; %Remove negative feature importances to make them scalable for the scatter size
mean_sf_acc = mean(sf_acc);
sd_sf_acc = std(sf_acc);
[sorted_sf_acc,sorted_idx] = sort(mean_sf_acc,'descend');
sorted_sf_sd = sd_sf_acc(sorted_idx);
% [Y,X] = ndgrid(1:size(fi,1),1:size(fi,2));
fi = fi(sorted_idx,:);

figure('Color','w');
s1 = subplot('Position',[0.3 0.1 0.4 0.8]);
imagesc(fi)
xlabel('Week')
xticks(1:5)
yticks(1:length(fi))
yticklabels(input_features(sorted_idx))
ylim([0.5 length(fi)+0.5])
ax1 = get(gca);
% ax1.YAxis.Visible = 'off';
% ax1.XRuler.TickLength = 0;
cb = colorbar('Location','northoutside'); cb.Title.String = "Feature importance"; cb.Position = [0.3 0.91 0.4 0.02];
set(s1,'FontSize',7)

s2 = subplot('Position',[0.71 0.1 0.28 0.8]);
b = barh(1:length(sorted_sf_acc),sorted_sf_acc(end:-1:1),1,'k','EdgeColor','w');
hold on
er = errorbar(b.YData,b.XData,zeros(1,length(sorted_sf_sd)),sorted_sf_sd(end:-1:1),'horizontal','Color',[0.5 0.5 0.5],...
    'linestyle','none','CapSize',1);
axis tight
yticks([])
xlabel("Accuracy")
set(gca,'FontSize',7)
box off
axis tight
set(gca,'XLim',[0.5 1.1])
set(gca,'XTick',[0.6 0.8 1.0])
