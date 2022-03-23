function ax = predictAgeBoxplot(age,pred_weeks)
N_groups = max(cellfun(@length, pred_weeks));
X = unique(age);
marker_size = 1;
preds = [pred_weeks{:}];
data = arrayfun(@(x) reshape(preds{x},1,[]),1:length(preds),'un',0); data = [data{:}];
groups = arrayfun(@(x) ones(1,length(preds{x}))*x,1:length(preds),'un',0);
groups = [groups{:}];

% labels = arrayfun(@(x) ['Week ' num2str(x)],1:length(X),'un',0);
labels = 1:length(X);

positions = [];
x_wt = [];
offset = 0.75;
for i = 1:length(preds)
    group_idx = rem(i-1,N_groups);
   if group_idx == 0
        pos = ceil(i/N_groups)*N_groups;
        x_wt = [x_wt pos];
   else
      pos = offset*group_idx+ceil(i/N_groups)*N_groups;
   end
    positions = [positions repelem(pos,length(preds{i}))];
end
xdata = unique(positions);
xticks = arrayfun(@(x) mean(xdata(find(xdata==x):find(xdata==x)+N_groups-1)),x_wt);
%%
leg = ["WT", "A53T", "LNA"];
if N_groups ==4
    leg = ["WT", "A53T", "LNA WT", "LNA A53T"];
end
% figure('color','w');
ax = axes('NextPlot','add','FontSize',6.5);

arrayfun(@(i) boxchart(ax,groups(ismember(groups,i:N_groups:max(groups))),data(ismember(groups,i:N_groups:max(groups))),...
    'LineWidth',0.5,'MarkerSize',marker_size,'XData',positions(ismember(groups,i:N_groups:max(groups)))),1:N_groups)
plot(ax,[0 max(xticks)*1.1],[0 length(labels)*7*1.2],'k--'); p.Color(4) = 0.3;
set(gca,'XTick',xticks)
set(gca,'XTickLabels',labels)
set(gca,'YTick',X+1)
set(gca,'YTickLabels',labels)
% set(ax,'FontSize',7)
xlabel('Culture age [weeks]')
ylabel('Predicted age [weeks]')
xlim([0 max(xticks)*1.1])
ylim([0 length(labels)*7*1.2])
% legend(leg(1:N_groups),'Location','northwest','Box','off')
% p = plot(get(gca,'XLim'),get(gca,'YLim'),'k--'); p.Color(4) = 0.3;
% set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');