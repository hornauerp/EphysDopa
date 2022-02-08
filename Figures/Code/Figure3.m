addpath(genpath('/cluster/project/bsse_sdsc/BELUB/Scripts/'))
addpath(genpath('H:\Data'))
%% Load Data
cls = [1:4,9:14];
path = 'H:\Data\NewAnalysis';
array_name = 'nw_array.mat';
nw_raw = LoadNW(path,array_name,cls);

%% PCA data
MIN_AGE = 0;
MAX_AGE = 40;
batchIDs = 1:2;
TrainFilter = [];
TestFilter = {{'Treatment','-','ntLNA'}};
markerSize = 15;
feat_th = 0.8;
FeatureGroup = "nw_mat"; feat_idx = [];%nw_idx(1:sum(nw_sorted>feat_th))-2;
[nw_train_color,nw_train_mat,~] = preparePCAinput(MIN_AGE,MAX_AGE,nw_raw,batchIDs,TrainFilter,FeatureGroup,feat_idx);
FeatureGroup = "[wf_mat act_mat]"; feat_idx = [];%tmp_idx(1:sum(tmp_sorted>feat_th));
[sc_train_color,sc_train_mat,~] = preparePCAinput(MIN_AGE,MAX_AGE,nw_raw,batchIDs,TrainFilter,FeatureGroup,feat_idx);
FeatureGroup = "[nw_mat wf_mat act_mat]"; feat_idx = [];%idx(1:sum(sorted>feat_th))-2;
[comb_train_color,comb_train_mat,~] = preparePCAinput(MIN_AGE,MAX_AGE,nw_raw,batchIDs,TrainFilter,FeatureGroup,feat_idx);
% [test_color,test_mat,batch] = preparePCAinput(MIN_AGE,MAX_AGE,nw_raw,batchIDs,TestFilter,FeatureGroup);
test_color = [];test_mat = [];

%% Prediction accuracies
tps = 6;
mat_idx = 2; sc_acc = nan(5,tps); sc_imp_vec = [];
for week = 1:tps
    [sc_imp,sc_avg_acc,mutation] = iterateClassification(nw_raw,week,mat_idx);%mat_idx=1 for nw, mat_idx=2 for sc
    sc_acc(:,week) = sc_avg_acc;
    sc_imp_vec(week,:) = sc_imp;
end

mat_idx = 1; nw_acc = nan(5,tps); nw_imp_vec = [];
for week = 1:tps
    [nw_imp,nw_avg_acc,~] = iterateClassification(nw_raw,week,mat_idx);%mat_idx=1 for nw, mat_idx=2 for sc
    nw_acc(:,week) = nw_avg_acc;
    nw_imp_vec(week,:) = nw_imp;
end
%%
mat_idx = 1:2; comb_acc = nan(5,tps); comb_imp_vec = [];
for week = 1:tps
    [comb_imp,comb_avg_acc,~] = iterateClassification(nw_raw,week,mat_idx);%mat_idx=1 for nw, mat_idx=2 for sc
    comb_acc(:,week) = comb_avg_acc;
    comb_imp_vec(week,:) = comb_imp;
end
%% Age prediction data
batches = unique([nw_raw.PlatingDate]);
batchIDs = 1:2; age_range = [0 40]; N_groups = 2;%3 includes treatment
inclusion_array = {{{'PlatingDate',batches(batchIDs)}},...%,{'Mutation','wt'}
    {{'PlatingDate',batches(batchIDs)}}};%,{'Treatment','LNA7'}
exclusion_array = {{{'DIV',28,9},{'Treatment','LNA7'}},...
    {{'DIV',28,9}}};
% vars2use = [1]; %1:nw,2:wf,3:mc,4:act
% [age,nw_pred_weeks,pI_nw,vars_nw,cl_nw] = predictAgeFunction(nw_raw,vars2use,N_groups,age_range,inclusion_array,exclusion_array);
% vars2use = [2,4];
% [~,sc_pred_weeks,pI_sc,vars_sc,cl_sc] = predictAgeFunction(nw_raw,vars2use,N_groups,age_range,inclusion_array,exclusion_array);
vars2use = [1,2,4];
[age,comb_pred_weeks,pI_comb,vars_comb,cl_comb] = predictAgeFunction(nw_raw,vars2use,N_groups,age_range,inclusion_array,exclusion_array);
%%
sc_imp = pI_comb(15:end); nw_imp = pI_comb(1:14); sc_vars = vars_comb(15:end); nw_vars = vars_comb(1:14);
[nw_imp,nw_idx] = maxk(nw_imp,3); [sc_imp,sc_idx] = maxk(sc_imp,3);

figure('Color','w')
subplot(1,2,1)
bar(nw_imp,'k','FaceAlpha',0.8)
box off
xticklabels(nw_vars(nw_idx))
xtickangle(45)
title('Network features')
set(gca,'FontSize',7)
subplot(1,2,2)
bar(sc_imp,'k','FaceAlpha',0.8)
box off
xtl = sc_vars(sc_idx); xtl{2} = 'RegMagnitude';
xticklabels(xtl)
xtickangle(45)
title('Single cell features')
set(gca,'FontSize',7)
%% PCA plots
ss = [16,32,64,128,256,512,1];
x_cell = arrayfun(@(x) x:1:x+length(ss)-2,linspace(0.55,1.45,4),'un',0);
load('H:\Data\AnalysisResults\Batch_12\pca_input.mat')
x_sz = 0.2; y_sz = 0.2*3/4;x_offset = 0.1;y_offset = 0.1;x_dist = 0.3; y_dist = 0.32*3/4;
figure('Color','w','Position',[400 400 560 420*4/3]);
markerSize = 15; test_mat = []; test_color = [];
subplot('Position',[x_offset y_offset+3*y_dist x_sz y_sz])
viewAngle = [-43 27];%Single Cell
[sc_silhouette_coef,sc_sum_explained] = plotPCA(sc_train_mat,sc_train_color,test_mat,test_color,viewAngle,markerSize);
title('Single-cell features','FontWeight','normal')
set(gca,'FontSize',7)
x = get(gca,'Xlabel');y= get(gca,'Ylabel');z= get(gca,'Zlabel');
x.Rotation = 20;y.Rotation = -20; x.Position = [-5 -11 -10]; y.Position = [-21 -5 -9];

subplot('Position',[x_offset+x_dist y_offset+3*y_dist x_sz y_sz])
viewAngle = [35 22];%Network
[nw_silhouette_coef,nw_sum_explained] = plotPCA(nw_train_mat,nw_train_color,test_mat,test_color,viewAngle,markerSize);
title('Network features','FontWeight','normal')
set(gca,'FontSize',7)
x = get(gca,'Xlabel');y= get(gca,'Ylabel');z= get(gca,'Zlabel');
x.Rotation = -10; x.Position = [3 -8 -7];y.Rotation = 25; y.Position = [8 6 -9];
zticks([-5 0 5]);zlim([-5 6])

subplot('Position',[x_offset+2*x_dist y_offset+3*y_dist x_sz y_sz])
viewAngle = [16 37];%Combined
[comb_silhouette_coef,comb_sum_explained] = plotPCA(comb_train_mat,comb_train_color,test_mat,test_color,viewAngle,markerSize);
title('Combined features','FontWeight','normal')
set(gca,'FontSize',7)
x = get(gca,'Xlabel');y= get(gca,'Ylabel');z= get(gca,'Zlabel');
x.Rotation = -7; x.Position = [3 -12 -17.5]; y.Rotation = 57; y.Position = [14 2 -18];
yticks([0 10])

l = legend({'WT','A53T'},'Position',[0.05+x_dist y_offset+3*y_dist+y_sz 0.01 0.01],'Box','off');l.ItemTokenSize(1) = 10;


% Prediction accuracies
capsz = 3;
% load('H:\Data\accuracies_week_500.mat')
% load('H:\Data\AnalysisResults\Batch_12\weekly_accuracies.mat')
load('H:\Data\AnalysisResults\Batch_12\Hyperparameter\sc_weekly_classification.mat')
sc_acc = test_acc_mat';
chance_line = 14/26;%max([sum(mutation)/length(mutation) 1-sum(mutation)/length(mutation)]);
subplot('Position',[x_offset y_offset+2*y_dist x_sz y_sz])
bar([mean(sc_acc(:,1:5)) 0 mean(sc_acc(:,6))],'k','FaceAlpha',0.8)
hold on
err = [std(sc_acc(:,1:5)) 0 std(sc_acc(:,6))];
er = errorbar(1:7,[mean(sc_acc(:,1:5)) 0 mean(sc_acc(:,6))],err,err,'k','linestyle','none','CapSize',capsz);
yline(chance_line,'k:','LineWidth',2,'Alpha',0.3);
xlabel('Week')
set(gca,'XTickLabel',{1:5,'All'})
set(gca,'XTick',[1:5,7])
set(gca,'YLim',[0.4 1.05])
set(gca,'YTick',[0.4 0.6 0.8 1])
ylabel('Accuracy')
box off
set(gca,'FontSize',7)

load('H:\Data\AnalysisResults\Batch_12\Hyperparameter\nw_weekly_classification.mat')
nw_acc = test_acc_mat';
subplot('Position',[x_offset+x_dist y_offset+2*y_dist x_sz y_sz])
bar([mean(nw_acc(:,1:5)) 0 mean(nw_acc(:,6))],'k','FaceAlpha',0.8)
hold on
err = [std(nw_acc(:,1:5)) 0 std(nw_acc(:,6))];
er = errorbar(1:7,[mean(nw_acc(:,1:5)) 0 mean(nw_acc(:,6))],err,err,'k','linestyle','none','CapSize',capsz);
yline(chance_line,'k:','LineWidth',2,'Alpha',0.3);
xlabel('Week')
set(gca,'XTickLabel',{1:5,'All'})
set(gca,'XTick',[1:5,7])
set(gca,'YLim',[0.4 1.05])
set(gca,'YTick',[0.4 0.6 0.8 1])
ylabel('Accuracy')
box off
set(gca,'FontSize',7)

load('H:\Data\AnalysisResults\Batch_12\Hyperparameter\comb_weekly_classification.mat')
comb_acc = test_acc_mat';
subplot('Position',[x_offset+2*x_dist y_offset+2*y_dist x_sz y_sz])
bar([mean(comb_acc(:,1:5)) 0 mean(comb_acc(:,6))],'k','FaceAlpha',0.8)
hold on
err = [std(comb_acc(:,1:5)) 0 std(comb_acc(:,6))];
er = errorbar(1:7,[mean(comb_acc(:,1:5)) 0 mean(comb_acc(:,6))],err,err,'k','linestyle','none','CapSize',capsz);
yl = yline(chance_line,'k:','LineWidth',2,'Alpha',0.3);
xlabel('Week')
set(gca,'XTickLabel',{1:5,'All'})
set(gca,'XTick',[1:5,7])
set(gca,'YLim',[0.4 1.05])
set(gca,'YTick',[0.4 0.6 0.8 1])
ylabel('Accuracy')
box off
set(gca,'FontSize',7)
txt = text(8.5,chance_line,'Chance','HorizontalAlignment','center','Rotation',-90,'FontSize',7);

% Axis breaks
annotation('line',[0.845 0.85],[0.57 0.59])
annotation('line',[0.84 0.845],[0.57 0.59])
annotation('line',[0.845-x_dist 0.85-x_dist],[0.57 0.59])
annotation('line',[0.84-x_dist 0.845-x_dist],[0.57 0.59])
annotation('line',[0.845-2*x_dist 0.85-2*x_dist],[0.57 0.59])
annotation('line',[0.84-2*x_dist 0.845-2*x_dist],[0.57 0.59])

% Age prediction
% load('H:\Data\age_pred_12_250.mat')
load('H:\Data\AnalysisResults\Batch_12\age_pred.mat')
apcm = othercolor('RdBu9',2);
N_groups = 2;
% sc_train = arrayfun(@(x) sc_pred_weeks{x}(1),1:5,'un',0);
% sc_train = arrayfun(@(x) {sc_train{x}{:}(ismember(cl_sc(1:5:end),[1 3])),sc_train{x}{:}(ismember(cl_sc(1:5:end),[2 4]))},1:5,'un',0);
% predictAgeBoxplot(age,N_groups,sc_train)
% groups = {'WT','A53T','LNA'};
% % legend(groups(1:N_groups),'Location','southeast','box','off');
% set(gca,'Position',[x_offset y_offset x_sz y_sz])
% set(gca,'FontSize',7)
% set(gca,'ColorOrder',flipud(apcm))
% 
% nw_train = arrayfun(@(x) nw_pred_weeks{x}(1),1:5,'un',0);
% nw_train = arrayfun(@(x) {nw_train{x}{:}(ismember(cl_nw(1:5:end),[1 3])),nw_train{x}{:}(ismember(cl_nw(1:5:end),[2 4]))},1:5,'un',0);
% predictAgeBoxplot(age,N_groups,nw_train)
% set(gca,'Position',[x_offset+x_dist y_offset x_sz y_sz])
% set(gca,'FontSize',7)
% set(gca,'ColorOrder',flipud(apcm))

% subplot('Position',[x_offset+2*x_dist y_offset x_sz y_sz])
% predictAgePlot(age,N_groups,comb_pred_weeks)
comb_train = arrayfun(@(x) comb_pred_weeks{x}(1),1:5,'un',0);
comb_train = arrayfun(@(x) {comb_train{x}{:}(ismember(cl_comb(1:5:end),[1 3])),comb_train{x}{:}(ismember(cl_comb(1:5:end),[2 4]))},1:5,'un',0);
predictAgeBoxplot(age,N_groups,comb_train)
hold on
x_lim = get(gca,'XLim');
y_lim = get(gca,'YLim');
p = plot(x_lim,y_lim,'k--'); p.Color(4) = 0.5;
% hold on
% weekly_rmse = arrayfun(@(x) sqrt(immse(comb_pred_weeks{x}{1},ones(size(comb_pred_weeks{x}{1}))*age(x))),1:5);
% xt = get(gca,'XTick');
% arrayfun(@(x) text(xt(x),40,num2str(weekly_rmse(x),2),'FontSize',fontsz,'HorizontalAlignment','center'),1:5)

% yyaxis right
% wp = plot(get(gca,'Xtick'),weekly_rmse,'k:');
% wp.Color(4) = 0.5;
% ylabel('RMSE')
set(gca,'Position',[x_offset y_offset+y_dist x_sz y_sz])
% set(gca,'Position',[x_offset+2*x_dist y_offset x_sz y_sz])
set(gca,'FontSize',7)
set(gca,'ColorOrder',flipud(apcm))
% set(gca,'ColorOrder',zeros(7,3))

subplot('Position',[0.9 0.1 0.1 0.1])
[imp_vals,imp_idx] = maxk(pI_comb,5);
bar(imp_vals,'FaceColor','k','FaceAlpha',0.8)
ylabel('Feature importance')
% xticklabels(vars_comb(imp_idx))
% xticklabels({'INTRABF','BRT','INTERBF','ISIV','T2PD'})
xticklabels({'AUCP1','INTRABF','T2PD','INTERBF','NRM'})
xtickangle(45)
box off
set(gca,'FontSize',7)
set(gca,'Position',[x_offset+x_dist y_offset+y_dist x_sz y_sz])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
markerstyles = {'o', '^'};
marker_sz = 3;
colors = othercolor('Set13',10);

% subplot('Position',[0.9 0.1 0.1 0.1])
% load('H:\Data\subsampling_age_template.mat')
load('H:\Data\AnalysisResults\Batch_12\subsampling_age_sc.mat')
tmp_rmse = rmse_array;
% load('H:\Data\subsampling_age_network.mat')
load('H:\Data\AnalysisResults\Batch_12\subsampling_age_nw.mat')
rmse_mean = cellfun(@mean,{tmp_rmse rmse_array},'un',0);
rmse_sd = cellfun(@std,{tmp_rmse rmse_array},'un',0);

ax = axes('NextPlot','add');
cellfun(@(x,m,s,mkr,c) errorbar(x, m, s,'Marker',mkr,'MarkerFaceColor',c,'MarkerSize',marker_sz,'Color',c,'CapSize',0,'LineWidth',1),...
    x_cell([2,3]),rmse_mean,rmse_sd,markerstyles,{colors(3,:) colors(9,:)})
ax.XTick = 1:6;
ax.XLim = [0.5 length(ss)-0.5];
ax.XTickLabel = ss;
ax.XLabel.String = '# neurons';
ax.YLim = [5.5 9];
ax.YLabel.String = "Age prediction" + newline + "RMSE [days]";
leg = legend({'SC','NW'},'Location','northeast','box','off');
leg.ItemTokenSize = [10 10];
set(gca,'FontSize',7)
set(gca,'Position',[x_offset+2*x_dist y_offset+y_dist x_sz y_sz])

% Last row
last_row_offset = 0.02;
load('H:\Data\AnalysisResults\Batch_12\purity.mat')
ax = axes('NextPlot','add');
silh_mean = cellfun(@mean,{purity_sc_pca(:,1:6),purity_nw_pca(:,1:6)},'un',0);
silh_sd = cellfun(@std,{purity_sc_pca(:,1:6),purity_nw_pca(:,1:6)},'un',0);
cellfun(@(x,m,s,mkr,c) errorbar(x, m, s,'Marker',mkr,'MarkerFaceColor',c,'MarkerSize',marker_sz,'Color',c,'CapSize',0,'LineWidth',1),...
    x_cell([2,3]),silh_mean,silh_sd,markerstyles,{colors(3,:) colors(9,:)})
ylabel('Cluster purity')
xticks(1:6)
xticklabels(ss)
xlabel('# neurons')
xlim([0.5 6.5])
set(gca,'FontSize',7)
set(gca,'Position',[x_offset y_offset-last_row_offset x_sz y_sz])


load('H:\Data\AnalysisResults\Batch_12\subsampling_significances.mat')
combined_mats{1} = (sensitivity_mats{1}*sum(sig_idx{1})+sensitivity_mats{2}*sum(sig_idx{2}))/(sum(sig_idx{1})+sum(sig_idx{2}));
combined_mats{2} = (sensitivity_mats{3}*sum(sig_idx{3})+sensitivity_mats{4}*sum(sig_idx{4}))/(sum(sig_idx{3})+sum(sig_idx{4}));
comb_mean = cellfun(@mean,combined_mats,'un',0);
comb_sd = cellfun(@std,combined_mats,'un',0);
ax = axes('NextPlot','add');
cellfun(@(x,m,s,mkr,c) errorbar(x, m, s,'Marker',mkr,'MarkerFaceColor',c,'MarkerSize',marker_sz,'Color',c,'CapSize',0,'LineWidth',1),...
    x_cell([2,3]),comb_mean,comb_sd,markerstyles,{colors(3,:) colors(9,:)})
ax.XTick = 1:6;
ax.XLim = [0.5 length(ss)-0.5];
ax.XTickLabel = ss;
ax.XLabel.String = '# neurons';
ax.YLim = [0 1.1];
ax.YTick = [0:0.2:1];
ax.YTickLabel = [0:20:100];
ax.YLabel.String = "Detected significant" + newline + "features [%]";

set(gca,'FontSize',7)
set(gca,'Position',[x_offset+x_dist y_offset-last_row_offset x_sz y_sz])


% load('H:\Data\AnalysisResults\Batch_12\subsampling_rf_sc.mat')
load('H:\Data\AnalysisResults\Batch_12\Hyperparameter\sc_subsampling_classification.mat')
temp_acc = test_acc_mat;
% load('H:\Data\AnalysisResults\Batch_12\subsampling_rf_nw.mat')
load('H:\Data\AnalysisResults\Batch_12\Hyperparameter\nw_subsampling_classification.mat')
mean_acc = cellfun(@(x) mean(x,[1 3]),{temp_acc test_acc_mat},'un',0);
sd_acc = cellfun(@(x) std(mean(x,3)),{temp_acc test_acc_mat},'un',0);
ax = axes('NextPlot','add');
cellfun(@(x,m,s,mkr,c) errorbar(x, m, s,'Marker',mkr,'MarkerFaceColor',c,'MarkerSize',marker_sz,'Color',c,'CapSize',0,'LineWidth',1),...
    x_cell([2,3]),mean_acc,sd_acc,markerstyles,{colors(3,:) colors(9,:)})
box off
ylim([0.5 1])
ylabel('Classification accuracy')
xticks(1:6)
xticklabels(ss)
xlabel('# neurons')
xlim([0.5 6.5])
set(gca,'FontSize',7)
set(gca,'Position',[x_offset+2*x_dist y_offset-last_row_offset x_sz y_sz])


%%
f = gcf;
exportgraphics(f,fullfile('C:\Users\Philipp\Documents\Promotion\Drafts\2020\Paper\Figures','Figure3.tif'),'Resolution',300)