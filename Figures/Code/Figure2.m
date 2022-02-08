addpath(genpath('H:\Data'))
%% Percentage change
mult_fact = 12357834; n_groups = 2;
cm_gray = flipud(gray(6));
cm_gray = cm_gray(1:4,:);
% cm_gray = cm_gray(1:end-70,:);
pVal1 = 0.0001; pVal2 = 0.001; pVal3 = 0.01; pVal4 = 0.05; tps = 6:7:40;
% load('H:\Data\sf_mat_500.mat')
load('H:\Data\AnalysisResults\Batch_12\Hyperparameter\sf_classification.mat')
sf_acc_mat = sf_acc;
sf_sd = std(sf_acc);
sf_acc = mean(sf_acc);

[sorted,idx] = sort(sf_acc,'descend');
[nw_sorted,nw_idx] = sort(sf_acc(1:14),'descend'); nw_sorted_sd = sf_sd(nw_idx); nw_idx = nw_idx+2;
[tmp_sorted,tmp_idx] = sort(sf_acc(15:end),'descend'); tmp_sorted_sd = sf_sd(:,15:end); tmp_sorted_sd = tmp_sorted_sd(tmp_idx);
[wf_sorted,wf_idx] = sort(sf_acc(15:24),'descend');
[act_sorted,act_idx] = sort(sf_acc(25:end),'descend');

feature_group = "nw";
% load('feat_imp_over_time.mat')
if feature_group =="tmp"
%     load('H:\Data\Prism\CDI\TemplateFeatures\tmp_data.mat')
    load('H:\Data\AnalysisResults\Batch_12\tmp_matrix.mat');tmp_matrix(:,:,2) = []; tmp_sd(:,:,2) = [];
    color_mat = load('H:\Data\AnalysisResults\Batch_12\sc_color_mat.mat','sc_color_mat');color_mat = color_mat.sc_color_mat;
    color_mat(2,:) = [];
    color_mat = color_mat(tmp_idx,:);
    sorted = tmp_sorted;
    err = tmp_sorted_sd;
    sf_acc_mat = sf_acc_mat(:,15:end);
    sel_feat_imp = mean(feat_imp(15:end,:,:),3);
    idx = tmp_idx;
    sel_feat_imp = sel_feat_imp(idx,:);
    whole_mat = tmp_matrix(:,:,tmp_idx);%cat(3,nw_matrix,tmp_matrix);
    sd = tmp_sd(:,:,tmp_idx);
    tmp_vars = {'\it AMPL','\it HLFW','\it ASYM','\it T2PR','\it T2PD','\it AUCP1','\it AUCP2','\it AUCT','\it RISE','\it DECAY',...
        'ISIM','ISIV','ISICV','PACF','SCRF','SCRM','RFIT'};
    labels = tmp_vars(tmp_idx)'; labels(ismember(labels,'Fano')) = {'ISICV'};
    scaling_factor = 1;
    feature_idx = [7 3 2]; %For templates
    ylabels = ["RFIT","ASYM","SCRF [Hz]"];
% ylabels = ["AUCpeak2", "Reg" + newline + "Freq [Hz]","PACF"];
    fig_title = 'Single-cell features';
    
elseif feature_group=="nw"
%     load('H:\Data\Prism\CDI\NetworkFeatures\nw_data.mat')
load('H:\Data\AnalysisResults\Batch_12\nw_matrix.mat')
    color_mat = load('H:\Data\AnalysisResults\Batch_12\nw_color_mat.mat','nw_color_mat');color_mat = color_mat.nw_color_mat;
    color_mat = color_mat(nw_idx-2,:);
    sf_acc_mat = sf_acc_mat(:,1:14);
    sel_feat_imp = mean(feat_imp(1:14,:,:),3);
    idx = nw_idx;
    sel_feat_imp = sel_feat_imp(idx-2,:);
    sorted = nw_sorted;
    err = nw_sorted_sd;
    whole_mat = nw_matrix(:,:,nw_idx); 
    sd = nw_sd(:,:,nw_idx); 
    nw_vars(3:end) = {'IBIV','IBIM','MBD','VBD','INTRABF','INTERBF','BRT','BRV','BDT','BDV','SYNC','NRF','NRM','NRFIT'};
    labels = nw_vars(nw_idx)';%labels{1} = 'BaselineFreq';
    scaling_factor = 1;
    feature_idx = [7 2 1];
%     ylabels = ["Burst" + newline + "duration [s]","Reg" + newline + "Freq [Hz]","Synchronicity"];
    ylabels = [ "BDT [s]","IBIM [s]","NRF [Hz]"];
    fig_title = 'Network features';
end


boxsize = 0.02; x_loc = 0.15; h_height = 0.3; diff_heatmaps = 0.01; diff_colorbar = 0.01;cb_height = 0.015;
fontsz = 7;
ref_idx = 2;
comp_idx = 1;
comp_mat = squeeze(whole_mat(ref_idx,:,:)./whole_mat(comp_idx,:,:));
comp_mat(comp_mat<1) = -1./comp_mat(comp_mat<1)+1;
comp_mat(comp_mat>0) = comp_mat(comp_mat>0)-1;
c = colormap(cm_gray);
figure('Color','w');
subplot('Position',[x_loc,...
                    h_height,...
                    size(comp_mat,1)*boxsize,...
                    size(comp_mat,2)*boxsize])
h = heatmap(comp_mat','ColorLimits',[-1 1],'FontSize',fontsz);
cm = othercolor('RdBu9',n_groups*20); cm((size(cm,1)/2)-2:size(cm,1)/2+2,:) = [];
h.Colormap = cm;
colorbar off
h.XLabel = 'Week';
for l = 1:length(labels)
   if ismember(l,feature_idx)
       labels{l} = ['\bf ' labels{l}];
   end
end
h.YData = labels;
h_struct = struct(h);
h_struct.YAxis.FontSize = 6;
ax = axes;
% [imp_y, imp_x] = find(sel_feat_imp(end:-1:1,:)==max(sel_feat_imp(end:-1:1,:),[],2));
% scatter(imp_x,imp_y,'filled','MarkerFaceColor','k','SizeData',5,'MarkerFaceAlpha',0.8)
[Y,X] = ndgrid(1:size(sel_feat_imp,1),1:size(sel_feat_imp,2));
c_scatter = sel_feat_imp(end:-1:1,:); c_scatter(c_scatter<=0) = nan;
scatter(X(:),Y(:),c_scatter(:)*10,'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.8)

ax.Color = 'none';
ax.Position = h.Position;
ax.XTick = [];
ax.YTick = [];
xlim([0.5 5.5])
ylim([0.5 length(labels)+0.5])



% Feature importances
subplot('Position',[x_loc+size(comp_mat,1)*boxsize+diff_heatmaps,...
    h_height,...
    size(color_mat,2)*boxsize,...
    size(color_mat,1)*boxsize])
% b = barh(1:size(imp_mat,2),mean(imp_mat),1,'k','EdgeColor','w');
b = barh(1:length(sorted),sorted(end:-1:1),1,'k','EdgeColor','w');
hold on
capsz = 1;
alpha = 0.5;
er = errorbar(b.YData,b.XData,zeros(1,length(err)),err(end:-1:1),'horizontal','Color',[alpha alpha alpha],...
    'linestyle','none','CapSize',capsz);
% b = boxchart(sf_acc_mat(:,idx(end:-1:1)),'Orientation','horizontal','MarkerSize',1);
% for i = 1:length(idx)
%    boxchart(repelem(i,size(sf_acc_mat,1)),sf_acc_mat(:,idx(length(idx)-i+1)),'Orientation','horizontal',...
%        'MarkerSize',1,'BoxFaceColor','k','MarkerColor','k')
%    hold on
% end
axis tight
yticks([])
% xlabel("Feature" + newline + "Importance")
xlabel("Accuracy")
set(gca,'FontSize',7)
box off
axis tight
set(gca,'XLim',[0.5 1.1])
set(gca,'XTick',[0.6 0.8 1.0])

subplot('Position',[x_loc+size(comp_mat,1)*boxsize+diff_heatmaps+0.1,...
                    h_height,...
                    size(color_mat,2)*boxsize,...
                    size(color_mat,1)*boxsize])
h2 = heatmap(color_mat*mult_fact,'FontSize',fontsz);

h2.Colormap = c;
h2.XDisplayLabels = ["\bf G", "\bf G x T", "\bf T-WT" "\bf T-A53T"];
h2.YDisplayLabels = nan(1,size(h2.ColorData,1));
axs = struct(gca);
cb = axs.Colorbar;
cb_ticks = cb.Limits(2)/8;
cb.Ticks = cb_ticks*[1 3 5 7];
cb.TickLabels = ["\it ns" pVal4 pVal3 "< 0.001"];
cb.Title.String = ["\it p-\itvalues" + newline + "(adj.)"];
cb.Title.FontWeight = 'bold';

% c = othercolor('RdBu4',size(whole_mat,1));
c = h.Colormap(linspace(1,size(h.Colormap,1),size(whole_mat,1)),:);
r = [];plots = [];
for p = 1:length(feature_idx)
    x_plot = 0.15; y_plot = 0.1;
    plots{p} = subplot('Position',[x_loc+0.45 h.Position(2)+0.5*h.Position(4)-0.5*y_plot+(p-2)*0.15 x_plot y_plot]); %0.1+0.15*p
    jitter = linspace(-1,1,size(whole_mat,1));
    for i = 1:size(whole_mat,1)
        x = tps+1+jitter(i);
        y = whole_mat(i,:,feature_idx(p))/scaling_factor;
        xx = linspace(min(x),max(x),100); yy = spline(x,y,xx);
        plot(xx,yy,'Color',c(i,:))
        hold on
        errorbar(tps+1+jitter(i),whole_mat(i,:,feature_idx(p))/scaling_factor,sd(i,:,feature_idx(p))/scaling_factor,...
            'LineWidth',1,'Color',c(i,:),'CapSize',0,'LineStyle','none','Marker','o','MarkerSize',2,...
            'MarkerFaceColor',c(i,:));
        set(gca,'FontSize',fontsz)
        marg = get(gca,'ylabel');
        set(marg,'Margin',3)
    end
    box off
    if p ==1
    xlabel('Week')
    end
%     ylabel(labels(feature_idx(p)))
% ylabel([labels{feature_idx(p)} ' ' char(ylabels(p))])
ylabel(ylabels(p))
xlim([3 max(xlim)])
    xticks(7:7:35)
    xticklabels(1:5)
    yt = get(gca,'YTick');
    yticks(linspace(min(yt),max(yt),2)) %Change here to insert more ticks
    yl = get(gca,'YLim');
    ylim([yl(1)-max(abs(yl))*0.1 yl(2)])
    set(gca,'Fontsize',fontsz-1)
%     set(plots{p}.YLabel,'Rotation',0)
%     set(plots{p}.YLabel,'VerticalAlignment','middle')
%     set(plots{p}.YLabel,'HorizontalAlignment','right')
%     subplot('Position',[x_loc+boxsize,h_height+boxsize*feature_idx(p),size(color_mat,2)*boxsize,boxsize])
r{p} = subplot('Position',[0.8 0.1 0.1 0.1]);
axis off
     rectangle('Position',[1 1 1 1],'FaceColor','none','LineWidth',1.5);
    set(gca,'Position',[x_loc,h_height+boxsize*(size(comp_mat,2)-feature_idx(p)),size(comp_mat,1)*boxsize,boxsize])
    
%     a1 = annotation('line',[r{p}.Position(1)+r{p}.Position(3) plots{p}.Position(1)],...
%         [r{p}.Position(2)+r{p}.Position(4) plots{p}.Position(2)+plots{p}.Position(4)],'LineStyle','--');
%     a1.Color(4) = 0;
%     uistack(a1,'bottom')
%     a2 = annotation('line',[r{p}.Position(1)+r{p}.Position(3) plots{p}.Position(1)],...
%         [r{p}.Position(2) plots{p}.Position(2)],'LineStyle','--');
%     a2.Color(4) = 0.9;
%     uistack(a2,'bottom')

% subplot('Position',[0.1 0.1, 0.1,0.1])
% plot(1,1,'ok','MarkerSize',8)
% text(1.1,1.3,['\bf' char(string(abs(p-length(feature_idx))+1))], 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', fontsz );
% axis off
% set(gca,'Position',[h.OuterPosition(1)-0.02,r{p}.Position(2),0.01,0.01])
% subplot('Position',[0.1 0.1, 0.1,0.1])
% plot(1,1,'ok','MarkerSize',8)
% text(1,1.3,['\bf' char(string(abs(p-length(feature_idx))+1))], 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', fontsz );
% axis off
% set(gca,'Position',[plots{p}.OuterPosition(1)-0.02,plots{p}.Position(2)+(plots{p}.Position(4)/2)-0.006,0.01,0.01])
end
subplot('Position',[0.9 0.1 0.02 0.02])
axis off
ax = gca;
ax.Colormap = cm;%othercolor('RdBu9',n_groups*20);
cb2 = colorbar('Location','north');
cb2.Ticks = [0 0.5 1];
cb2.TickLength = 0.02;
cb2.TickDirection = 'out';
cb2.TickLabels = {'<0.3','1','>3.0'};%{'A53T>WT','WT>A53T'};
cb2.Title.String = ["\bf WT/A53T"];%"\it Rel expression";
cb2.Position = [x_loc,...
                size(comp_mat,2)*boxsize+h_height+diff_colorbar,...
                size(comp_mat,1)*boxsize,...
                cb_height];
cb2.FontSize = fontsz-1;
% cb2.FontWeight = 'bold';
cb2.Title.Position = [21 14 0];
cb2.Ruler.TickLabelGapOffset = -1;
ttl = subplot('Position',[0.1 0.1 0.1 0.1]);
ttl_txt = text(0.2,1,fig_title,'HorizontalAlignment','center','FontSize',8);
axis off
set(gca,'Position',[(h.OuterPosition(1)+plots{1}.OuterPosition(1)+plots{1}.Position(3))/2 cb2.Position(2) 0.1 0.1])
%%
f = gcf;
exportgraphics(f,fullfile('C:\Users\Philipp\Documents\Promotion\Drafts\2020\Paper\Figures','SCFeatures.tif'),'Resolution',200)
% subplot('Position',[0.1 0.1 0.1 0.1])
% vars = [nw_prop'];%nw_prop' act_prop' 
% b = bar(max_imp,'FaceColor','k','FaceAlpha',0.8);
% box off
% set(gca,'Fontsize',fontsz-1)
% vars{imp_idx(2)} = 'BaselineFreq';
% set(gca,'XTickLabel',vars(imp_idx))
% set(gca,'XTickLabelRotation',45)
% set(gca,'Position',[plots{1}.OuterPosition(1)+plots{1}.OuterPosition(3)+0.05...
%     plots{3}.Position(2) x_plot y_plot])
% title('Top predictors')
% ylabel('Importance')