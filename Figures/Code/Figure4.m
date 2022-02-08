addpath(genpath('H:\Data'))
%% Load data
% Load asyn data
path = 'C:\Users\Philipp\Documents\Promotion\Images\DA_LNA\20200811\';
paths = find_subpaths(path,[],{'Optimized'});
paths = paths(2:end);
% paths = {'C:\Users\Philipp\Documents\Promotion\Images\DA_LNA\20200811\5206_NT_WT',...
%     'C:\Users\Philipp\Documents\Promotion\Images\DA_LNA\20200811\5215_NT_A53T'};

means = cell(1,length(paths));
conditions = cell(1,length(paths));
for i = 1:length(paths)
    condition = strsplit(paths{i},filesep);
    conditions{i} = condition{end};
   tbl =  readtable([paths{i} '\Results2.csv']);
   measures = string(tbl.Label);
   rows = startsWith(measures,"conf");
   means{i} = tbl.Mean(rows);
end

%%
wrapper = @(x) rmoutliers(x,'ThresholdFactor',5);
means = cellfun(wrapper,means,'UniformOutput',false);
results_table = nan(max(cellfun(@length, means)),length(means));
for i = 1:length(means)
   results_table(1:length(means{i}),i) = means{i}; 
end
%% Statistics
wt_idx = find(endsWith(string(conditions),"WT"));
WT_comparison = ranksum(means{wt_idx(1)},means{wt_idx(2)});
a53t_idx = find(endsWith(string(conditions),"A53T"));
A53T_comparison = ranksum(means{a53t_idx(1)},means{a53t_idx(2)});
[p,tbl,stats] = anova1(results_table);

%% PCA
% load('H:\Data\sf_acc.mat')
load('H:\Data\AnalysisResults\Batch_12\single_feature_results.mat')
% sf_acc = mean(sf_acc);
[sorted,idx] = sort(sf_acc,'descend'); idx = idx+2;
[nw_sorted,nw_idx] = sort(sf_acc(1:14),'descend');nw_idx = nw_idx+2;
[tmp_sorted,tmp_idx] = sort(sf_acc(15:end),'descend');
MIN_AGE = 0;
MAX_AGE = 33;
batchIDs = 3;
feat_th = 0.8;
pca_ind = ""; 
TrainFilter = []; 
if pca_ind == "nw"
    FeatureGroup = "[nw_mat]"; feat_idx = nw_idx(1:sum(nw_sorted>feat_th))-2;
elseif pca_ind == "tmp"
    FeatureGroup = "[wf_mat ]"; feat_idx = [];tmp_idx(1:sum(tmp_sorted>feat_th));
else
    FeatureGroup = "[nw_mat wf_mat act_mat]"; feat_idx = idx(1:sum(sorted>feat_th))-2;
end
[nw_train_color,nw_train_mat,~] = preparePCAinput(MIN_AGE,MAX_AGE,nw_raw,batchIDs,TrainFilter,FeatureGroup,feat_idx);
% figure;plotPCA(nw_train_mat,nw_train_color,[],[],[],15);
%%
nw_th = 0.799; tmp_th = 0.8; mult_fact = 12357834;
nw_sel = nw_idx(1:sum(nw_sorted>=nw_th));
tmp_sel = tmp_idx(1:sum(tmp_sorted>=tmp_th));
t_max = 4;
%% Plots
p_vals = nan(4,4);
p_vals(1,2) = 0.0027; p_vals(1,3) = 0.9994; p_vals(1,4) = 0.0001; p_vals(2,3) = 0.0029; p_vals(2,4) = 0.0002; p_vals(3,4) = 0.0001; 
labels = ["A53T","WT"];
data = cellfun(@mean,means([4 1; 2 3]));%[70.9 30.6; 77.1 21.4];
err = cellfun(@std,means([4 1; 2 3]));%[7.5 13.7; 22.1 6.1];
nbars = size(data, 2);
fontsz = 7;
%Parameters
star_sz = 16;
star_lateral = 0.1;
star_vertical = 1.07;
sigbar1_vertical = 1;
sigbar2_vertical = 1.2;
capsz = 3;
bar_position = [0.5 0.7 0.2 0.2];
figure('Color','w');
% subplot('Position',bar_position)
% b = bar(data,'grouped');
% b(2).FaceColor = [1 1 1];
% b(1).FaceColor = 0.2*[1 1 1];
% set(gca,'XTickLabel',labels)
% 
% hold on
% x = [];
% for i = 1:nbars
%     x = [x ; b(i).XEndPoints];
% end
% er = errorbar(x',data,err,err,'k','linestyle','none','CapSize',capsz);
% ylabel("Fluorescence")%" + newline + "intensity"
% xlabel('Cell line')
% y_lim = get(gca,'YLim');
% set(gca,'YLim',y_lim*1.3)
% y_sigbar1 = max(y_lim)*sigbar1_vertical;
% y_sigbar2 = max(y_lim)*sigbar2_vertical;
% sigbar_1 = er(1).XData;
% sigbar_2 = er(2).XData;
% 
% plot([sigbar_1(1),sigbar_2(1)],[y_sigbar1,y_sigbar1],'k')
% plot([sigbar_1(2),sigbar_2(2)],[y_sigbar1,y_sigbar1],'k')
% plot([sigbar_1(1),sigbar_1(2)],[y_sigbar2,y_sigbar2],'k')
% 
% x_scatter1 = [((sigbar_1(1)+sigbar_2(1))/2 -star_lateral), (sigbar_1(1)+sigbar_2(1))/2,((sigbar_1(1)+sigbar_2(1))/2 +star_lateral)];
% % x_scatter2 = [((sigbar_1(2)+sigbar_2(2))/2 -star_lateral), (sigbar_1(2)+sigbar_2(2))/2,((sigbar_1(2)+sigbar_2(2))/2 +star_lateral)];
% x_scatter2 = [0.5*(sigbar_1(2)+sigbar_2(2))-0.5*star_lateral, 0.5*(sigbar_1(2)+sigbar_2(2))+0.5*star_lateral];
% % x_scatter3 = [0.5*(sigbar_1(1)+sigbar_1(2))-0.5*star_lateral, 0.5*(sigbar_1(1)+sigbar_1(2))+0.5*star_lateral];
% x_scatter3 = [((sigbar_1(1)+sigbar_1(2))/2 -star_lateral), (sigbar_1(1)+sigbar_1(2))/2,((sigbar_1(1)+sigbar_1(2))/2 +star_lateral)];
% 
% s1 = scatter(x_scatter1,repelem(y_sigbar1*star_vertical,3),star_sz,'k*');
% s2 = scatter(x_scatter2,repelem(y_sigbar1*star_vertical,2),star_sz,'k*');
% s3 = scatter(x_scatter3,repelem(y_sigbar2*star_vertical,3),star_sz,'k*');
% box off
% set(gca,'FontSize',fontsz)

% PCA plot
subplot('Position',[0.1 0.65 0.15 0.15])
p = plot(nan,nan,'s','Color','k','MarkerFaceColor',0.2*[1 1 1]);
hold on
plot(nan,nan,'s','Color','k','MarkerFaceColor','w')
axis off
lb = legend(["Control","LNA"],'Box','off','NumColumns',2);lb.ItemTokenSize(1) = 5;
set(gca,'TickLength',[0 0])
set(gca,'FontSize',fontsz)
set(lb,'Position',[bar_position(1)+0.08 bar_position(2)+0.23 0.01 0.01])
box off
hold off

subplot('Position',[0.65 0.09 0.2 0.2])
viewAngle = [-68 12]; markerSize = 15;
nw_silhouette = plotPCA(nw_train_mat,nw_train_color,[],[],viewAngle,markerSize,0);
set(gca,'Position',[0.6 0.15 0.3 0.3])
% title('Network features')
set(gca,'FontSize',7)
x = get(gca,'Xlabel');y= get(gca,'Ylabel');z= get(gca,'Zlabel');
% x.Rotation = -7; x.Position = [2 -7.5 -8];y.Rotation = 50; y.Position = [7.5 2 -9];
% x.Rotation = -7; x.Position = [2 -6.5 -6];y.Rotation = 50; y.Position = [6.5 5 -9];
l = legend({'WT','WT+LNA','A53T','A53T+LNA'},'Box','off','NumColumns',2);l.ItemTokenSize(1) = 10;
% set(l,'Position',[0.78 0.31 0.01 0.01])
set(l,'Position',[0.75 0.5 0.01 0.01])
set(gca,'FontSize',fontsz)
% xlim([-6 6]);ylim([-6 6]); zlim([-6 6])

% Timeline plot
load('H:\Data\Prism\CDI\LNATreatment\WT\nw_matrix_3.mat')
wt_lna = load('H:\Data\AnalysisResults\Batch_3\WT_LNA\nw_matrix.mat');
a53t_lna = load('H:\Data\AnalysisResults\Batch_3\A53T_LNA\nw_matrix.mat');
nw_matrix = [a53t_lna.nw_matrix;wt_lna.nw_matrix];
nw_sd = [a53t_lna.nw_sd; wt_lna.nw_sd];
nw_vars = wt_lna.nw_vars;
n_groups = 4;
tl_plot = subplot('Position',[0.6 0.65 0.25 0.25]);
c = othercolor('RdBu4',n_groups); c([1 2 4 3],:) = c;
tps = (1:4)*7-1; scaling_factor = 1;
feature_idx = 14; whole_mat = nw_matrix(:,1:t_max,:); sd = nw_sd(:,1:t_max,:);
jitter = linspace(-1,1,size(whole_mat,1));
for i = 1:size(whole_mat,1)
    x = tps+1+jitter(i);
    y = whole_mat(i,:,feature_idx)/scaling_factor;
%     xx = linspace(min(x),max(x),100); yy = spline(x,y,xx);

    plot(x,y,'Color',c(i,:))
    hold on
        er = errorbar(tps+1+jitter(i),whole_mat(i,:,feature_idx)/scaling_factor,sd(i,:,feature_idx)/scaling_factor,...
        'LineWidth',1,'Color',c(i,:),'CapSize',0,'LineStyle','none','Marker','o','MarkerSize',2,...
        'MarkerFaceColor',c(i,:));
     set(get(get(er,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(gca,'FontSize',fontsz)
end
box off
xlabel('Week')
ylabel("NRF [Hz]")
xlim([3 max(xlim)])
xticks(7:7:35)
xticklabels(1:5)
yt = get(gca,'YTick');
yticks(linspace(min(yt),max(yt),2)) %Change here to insert more ticks
yl = get(gca,'YLim');
set(gca,'YLim',[yl(1)-max(abs(yl))*0.1 yl(2)])
% ylim([yl(1)-max(abs(yl))*0.1 yl(2)])
set(gca,'Fontsize',fontsz)
ax = gca;
lt = legend(ax.Children([4,2,8,6]),{'WT','WT+LNA','A53T','A53T+LNA'},'Box','off','NumColumns',1);lt.ItemTokenSize(1) = 10;
set(lt,'Position',[tl_plot.Position(1)+tl_plot.Position(3)+0.06 tl_plot.Position(2)+0.2 0.01 0.01])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('H:\Data\Prism\CDI\LNATreatment\WT\nw_data.mat')
whole_mat = nw_matrix;
nw_vars(3:end) = {'IBIV','IBIM','MBD','VBD','INTRABF','INTERBF','BRT','BRV','BDT','BDV','SYNC','NRF','NRM','NRFIT'};
labels = nw_vars;
boxsize = 0.025; x_loc = 0.1; h_height = 0.20; diff_heatmaps = 0.01; diff_colorbar = 0.01;cb_height = 0.015;
fontsz = 7;
ref_idx = 3;
comp_idx = 4;
comp_mat = squeeze(whole_mat(ref_idx,1:t_max,:)./whole_mat(comp_idx,1:t_max,:));
comp_mat(comp_mat<1) = -1./comp_mat(comp_mat<1)+1;
comp_mat(comp_mat>0) = comp_mat(comp_mat>0)-1;
comp_mat = comp_mat(:,nw_sel);
subplot('Position',[x_loc,...
                    h_height,...
                    size(comp_mat,1)*boxsize,...
                    size(comp_mat,2)*boxsize])
h1 = heatmap(comp_mat','ColorLimits',[-1 1],'FontSize',fontsz);
cm = othercolor('RdBu9',2*20); cm((size(cm,1)/2)-2:size(cm,1)/2+2,:) = [];
h1.Colormap = cm;
colorbar off
% h1.Title = 'WT';
h1.XLabel = 'Week';
labels = labels(nw_sel);
% labels{end} = 'ResonanceFit';
% labels{2} = 'RegMagnitude';
% labels{3} = 'RegFrequency';
for l = 1:length(labels)
   if ismember(l,feature_idx)
       labels{l} = ['\bf ' labels{l}];
   end
end
h1.YData = labels;
% h_struct = struct(h);
% h_struct.YAxis.FontSize = fontsz;

% NW WT p-values
load('H:\Data\AnalysisResults\Batch_3\WT_LNA\nw_color_mat.mat')
color_mat = nw_color_mat;
subplot('Position',[0.1 0.9 0.01 0.01])
color_mat = color_mat(nw_sel-2,3:4);
h1p = heatmap(color_mat.*mult_fact);
h1p.XDisplayLabels = ["\bf L", "\bf L x T"];
set(struct(h1p).NodeChildren(3), 'XTickLabelRotation', 90);
h1p.YDisplayLabels = nan(1,size(h1p.ColorData,1));
cm_gray = flipud(gray);
cm_gray = cm_gray(1:end-70,:);
h1p.Colormap = cm_gray;
h1p.FontSize = fontsz;
colorbar off
set(h1p,'Position',[h1.Position(1)+h1.Position(3)+0.01,h1.Position(2),2*boxsize,h1.Position(4)])

% Heatmap SC WT
load('H:\Data\Prism\CDI\LNATreatment\WT\tmp_matrix_3.mat')
load('H:\Data\Prism\CDI\LNATreatment\WT\tmp_data.mat')
wt_lna = load('H:\Data\AnalysisResults\Batch_3\WT_LNA\tmp_matrix.mat');
a53t_lna = load('H:\Data\AnalysisResults\Batch_3\A53T_LNA\tmp_matrix.mat');
tmp_matrix = [a53t_lna.tmp_matrix;wt_lna.tmp_matrix];
tmp_sd = [a53t_lna.tmp_sd; wt_lna.tmp_sd];
tmp_matrix(:,:,2) = []; tmp_sd(:,:,2) = []; 
no_wf_idx = [1,3];
tmp_vars = {'AMPL','HLFW','ASYM','T2PR','T2PD','AUCP1','AUCP2','AUCT','RISE','DECAY','ISIM','ISIV','ISICV','PACF','SCRF','SCRM','RFIT'};
labels = tmp_vars;
whole_mat = tmp_matrix;
fontsz = 7;
ref_idx = 3;
comp_idx = 4;
comp_mat = squeeze(whole_mat(ref_idx,1:t_max,:)./whole_mat(comp_idx,1:t_max,:));
comp_mat(comp_mat<1) = -1./comp_mat(comp_mat<1)+1;
comp_mat(comp_mat>0) = comp_mat(comp_mat>0)-1;
comp_mat = comp_mat(:,tmp_sel);
comp_mat = comp_mat(:,no_wf_idx); %Remove wf features
subplot('Position',[x_loc,...
                    h1.Position(2)+h1.Position(4)+diff_heatmaps,...
                    size(comp_mat,1)*boxsize,...
                    size(comp_mat,2)*boxsize])
h1b = heatmap(comp_mat','ColorLimits',[-1 1],'FontSize',fontsz);
h1b.Colormap = cm;
colorbar off
h1b.Title = 'WT';
% h1b.XLabel = 'Week';
labels = labels(tmp_sel); 
labels = labels(no_wf_idx);
Ax = gca; Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% labels{end} = 'ResonanceFit';
% labels{end-1} = 'RegMag';
% labels{end-2} = 'RegFreq';
% for l = 1:length(labels)
%    if ismember(l,feature_idx)
%        labels{l} = ['\bf ' labels{l}];
%    end
% end
h1b.YData = labels;

% SC WT p-values
color_mat = load('H:\Data\AnalysisResults\Batch_3\WT_LNA\sc_color_mat.mat');
color_mat = color_mat.sc_color_mat; color_mat(2,:) = [];
subplot('Position',[0.1 0.9 0.01 0.01])
color_mat = color_mat(tmp_sel,3:4);
color_mat = color_mat(no_wf_idx,:); %Remove waveform features
h1bp = heatmap(color_mat.*mult_fact);
h1bp.XDisplayLabels = nan(1,size(h1bp.ColorData,2));
h1bp.YDisplayLabels = nan(1,size(h1bp.ColorData,1));
if sum(color_mat)==0
    h1bp.Colormap = [1 1 1];
else
    h1bp.Colormap = cm_gray;
end
h1bp.FontSize = fontsz;
set(h1bp,'Position',[h1b.Position(1)+h1b.Position(3)+0.01,h1b.Position(2),2*boxsize,h1b.Position(4)])
colorbar off

% Heatmap NW A53T
% load('H:\Data\Prism\CDI\LNATreatment\A53T\nw_matrix_3.mat')
% load('H:\Data\Prism\CDI\LNATreatment\A53T\nw_data.mat')
% load('H:\Data\AnalysisResults\Batch_3\A53T_LNA\nw_matrix.mat')
% load('H:\Data\AnalysisResults\Batch_3\A53T_LNA\nw_color_mat.mat')
whole_mat = nw_matrix;
ref_idx = 1;
comp_idx = 2;
comp_mat = squeeze(whole_mat(ref_idx,1:t_max,:)./whole_mat(comp_idx,1:t_max,:));
comp_mat(comp_mat<1) = -1./comp_mat(comp_mat<1)+1;
comp_mat(comp_mat>0) = comp_mat(comp_mat>0)-1;
comp_mat = comp_mat(:,nw_sel);
subplot('Position',[x_loc+h1.Position(3)+h1p.Position(3)+0.03,...
                    h_height,...
                    size(comp_mat,1)*boxsize,...
                    size(comp_mat,2)*boxsize])
h2 = heatmap(comp_mat','ColorLimits',[-1 1],'FontSize',fontsz);

h2.Colormap = cm;
colorbar off
% h2.Title = 'A53T';
Ax = gca; Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
h2.XLabel = 'Week';

% NW A53T p-values
load('H:\Data\AnalysisResults\Batch_3\A53T_LNA\nw_color_mat.mat')
color_mat = nw_color_mat;
subplot('Position',[0.1 0.9 0.01 0.01])
color_mat = color_mat(nw_sel-2,3:4);
h2p = heatmap(color_mat.*mult_fact);
h2p.XDisplayLabels = ["\bf L", "\bf L x T"];
set(struct(h2p).NodeChildren(3), 'XTickLabelRotation', 90);
h2p.YDisplayLabels = nan(1,size(h2p.ColorData,1));
h2p.Colormap = cm_gray;
h2p.FontSize = fontsz;
set(h2p,'Position',[h2.Position(1)+h2.Position(3)+0.01,h2.Position(2),2*boxsize,h2.Position(4)])
colorbar off


% Heatmap SC A53T
% load('H:\Data\Prism\CDI\LNATreatment\A53T\tmp_matrix_3.mat')
% load('H:\Data\Prism\CDI\LNATreatment\A53T\tmp_data.mat')
whole_mat = tmp_matrix;
ref_idx = 1;
comp_idx = 2;
comp_mat = squeeze(whole_mat(ref_idx,1:t_max,:)./whole_mat(comp_idx,1:t_max,:));
comp_mat(comp_mat<1) = -1./comp_mat(comp_mat<1)+1;
comp_mat(comp_mat>0) = comp_mat(comp_mat>0)-1;
comp_mat = comp_mat(:,tmp_sel);
comp_mat = comp_mat(:,no_wf_idx); % Remove wf features
subplot('Position',[x_loc+h1.Position(3)+h1p.Position(3)+0.03,...
                    h2.Position(2)+h2.Position(4)+diff_heatmaps,...
                    size(comp_mat,1)*boxsize,...
                    size(comp_mat,2)*boxsize])
h2b = heatmap(comp_mat','ColorLimits',[-1 1],'FontSize',fontsz);
cm = othercolor('RdBu9',2*20); cm((size(cm,1)/2)-2:size(cm,1)/2+2,:) = [];
h2b.Colormap = cm;
colorbar off
h2b.Title = 'A53T';
Ax = gca; Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% h2b.XLabel = 'Week';

% SC A53T p-values
color_mat = load('H:\Data\AnalysisResults\Batch_3\A53T_LNA\sc_color_mat.mat');
color_mat = color_mat.sc_color_mat; color_mat(2,:) = [];
subplot('Position',[0.1 0.9 0.01 0.01])
color_mat = color_mat(tmp_sel,3:4);
color_mat = color_mat(no_wf_idx,:); %Remove waveform features
h2bp = heatmap(color_mat.*mult_fact);
h2bp.XDisplayLabels = nan(1,size(h2bp.ColorData,2));
h2bp.YDisplayLabels = nan(1,size(h2bp.ColorData,1));
if sum(color_mat)==0
    h2bp.Colormap = [1 1 1];
else
    h2bp.Colormap = cm_gray;
end
h2bp.FontSize = fontsz;
set(h2bp,'Position',[h2b.Position(1)+h2b.Position(3)+0.01,h2b.Position(2),2*boxsize,h2b.Position(4)])
colorbar off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Colorbar color
cbc = subplot('Position',[0.9 0.1 0.02 0.02]);
axis off
ax = gca;
cbc.Colormap = cm;%othercolor('RdBu9',n_groups*20);
cb2 = colorbar('Location','southoutside','FontSize',fontsz);
cb2.Ticks = [0 0.5 1];
cb2.TickLength = 0.01;
cb2.TickDirection = 'out';
% cb2.Direction = 'reverse';
% cb2.TickLabels = {['LNA > Control' repelem(' ',22)],[repelem(' ',22) 'Control > LNA']};
cb2.TickLabels = {'<0.3','1.0','>3.0'};
cb2.Title.String = "\bf Control/LNA";%"\it Rel expression";
cb2.Title.Position = [37.8 -20 0];
cb2.Ruler.TickLabelGapOffset = -1;
% cb2.Position = [x_loc,...
%                 size(comp_mat,2)*boxsize+h_height+diff_colorbar,...
%                 h2.Position(1)+h2.Position(3)-x_loc,...
%                 cb_height];
cb2.Position = [x_loc+0.5*h1.Position(3),...
                h1.OuterPosition(2)-0.03,...
                h2.Position(1)+h2.Position(3)-(x_loc+0.5*h1.Position(3)+0.5*h2.Position(3)),...
                cb_height];
            
% Colorbar p-vals
pc = subplot('Position',[0.2 0.7 0.2 0.2]);

axis off

cb = colorbar(pc,'Location','eastoutside','FontSize',fontsz);
cm_gray = flipud(gray(6));
cm_gray = cm_gray(1:4,:);
set(pc,'Colormap',cm_gray)
cb_ticks = cb.Limits(2)/8;
cb.Ticks = cb_ticks*[1 3 5 7];
cb.TickDirection = 'out';
cb.TickLabels = ["\it ns" 0.05 0.01 "< 0.001"];
cb.Title.String = "\it p-\itvalues" + newline + "(adj.)";
cb.Title.FontWeight = 'bold';
% cb.Limits = [min(cb.Ticks) max(cb.Ticks)];
cb.Position = [h2p.Position(1)+h2p.Position(3)+diff_heatmaps,...
               h2p.Position(2),...
               0.01,...
               h2bp.Position(2)+h2bp.Position(4)-h2p.Position(2)];
            
% subplot('Position',[0.1 0.8 0.1 0.1])
% scatter(nan,nan,'ks','filled');hold on;scatter(nan,nan,'ws','filled','MarkerEdgeColor','k')
% axis off
% l = legend({'Untreated','LNA'},'NumColumns',2,'Box','off','Position',[bar_position(1)+0.04,bar_position(2)+bar_position(4),0.1,0.1]);
% l.ItemTokenSize(1) = 10;

%%
f = gcf;
exportgraphics(f,fullfile('C:\Users\Philipp\Documents\Promotion\Drafts\2020\Paper\Figures','Figure_4_new.tif'),'Resolution',400)