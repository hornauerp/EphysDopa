function [ax,leg] = subsamplingPlot(ax,mat_1,mat_2,y_lim,y_label)

ss = [16,32,64,128,256,512,1];
x_cell = arrayfun(@(x) x:1:x+length(ss)-2,linspace(0.55,1.45,4),'un',0);
markerstyles = {'o', '^'};
marker_sz = 3;
colors = othercolor('Set13',10);

if size(mat_1,3)>1
    data_mean = cellfun(@(x) mean(x,[1 3]),{mat_1 mat_2},'un',0);
    data_sd = cellfun(@(x) std(mean(x,3)),{mat_1 mat_2},'un',0);
    
else
    data_mean = cellfun(@mean,{mat_1 mat_2},'un',0);
    data_sd = cellfun(@std,{mat_1 mat_2},'un',0);
end

% ax = axes('NextPlot','add');
hold on
cellfun(@(x,m,s,mkr,c) errorbar(ax, x, m, s,'Marker',mkr,'MarkerFaceColor',c,'MarkerSize',marker_sz,'Color',c,'CapSize',0,'LineWidth',1),...
    x_cell([2,3]),data_mean,data_sd,markerstyles,{colors(3,:) colors(9,:)})
ax.XTick = 1:6;
ax.XLim = [0.5 length(ss)-0.5];
ax.XTickLabel = ss;
ax.XLabel.String = '# neurons';
ax.YLim = y_lim;
ax.YLabel.String = y_label;
ax.XTickLabelRotation = 0;
leg = legend({'SC','NW'},'Location','northeast','box','off');
leg.ItemTokenSize = [10 10];
set(gca,'FontSize',7)