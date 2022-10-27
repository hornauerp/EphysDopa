addpath(genpath('/cluster/home/phornauer/Git/EphysDopa'))
addpath(genpath('/cluster/home/phornauer/Git/drtoolbox'))
addpath(genpath('F:\Promotion\Git\EphysDopa'))
addpath(genpath('C:\Users\Philipp\Downloads\umapAndEppFileExchange (4.1)'))

%% Load full dataset
load('F:\Promotion\Git\EphysDopa\Data\full_dataset_v2.mat')

%% Parameters
% Select features to be included in the single-cell clustering

wf_prop = [
    "Halfwidth",...
    "Asymmetry",...
    "T2Pratio",...
    "T2Pdelay",...
    "AUCpeak1",...
    "AUCpeak2",...
    "AUCtrough",...
    "Rise",...
    "Decay"
    ];

sc_prop = ["ISIMean",...
    "ISIVar",...
    "Fano",...
    "PACF",...
    "ResonanceFrequency",...
    "ResonanceMagnitude",...
    "ResonanceFit"
    ];

norm_idx = 1; %Normalize wf to max amplitude
rm_axon = 1; %Remove axonal signals (positive peak)
noise_cutoff = 14; %Cutoff for determining noise templates
plot_noisy = true; %Plot noisy templates
batches = unique([nw_array.PlatingDate]); %Refer to different batches via their PlatingDate property

%% Generate and clean input data
inclusion = {{'PlatingDate',batches([1])},{'DIV',27,34}}; %e.g. specify the batches to analyze by indexing the batches array // here we use the first two
exclusion = {{'Treatment','LNA7'}};%{'Mutation','wt'},

[feat_mat,wf_full,mut_idx,age_idx,treat_idx,unit_array] = generate_single_cell_clustering_input(nw_array,inclusion,exclusion,sc_prop,norm_idx,rm_axon,noise_cutoff,plot_noisy);

%% UMAP dimensionality reduction
rescaled_feat = normalize(feat_mat,'range',[-1 1]);
full = [wf_full rescaled_feat];
[reduction, umap, clusterIdentifiers, extras]=run_umap(full,'n_components',2,'n_neighbors',100,'min_dist',0.1,'cluster_detail','adaptive','spread',1,'sgd_tasks',20,...
    'verbose','none','color_file','C:\Users\Philipp\Downloads\umapAndEppFileExchange (4.1)\colorsByName.properties');

%% DBSCAN clustering
minpts = 50;
epsilon = 0.4;
db_labels = dbscan(reduction,epsilon,minpts);
fprintf('Detected %i clusters\n', max(db_labels))

%% Cluster plots
cluster_assignment = clusterIdentifiers; %Select array to determine cluster coloring
plot_cluster_scatter(reduction,cluster_assignment)
plot_waveform_clusters(wf_full,cluster_assignment)

%% Reassign clusters
new_cluster_groups = {[1],[2],[3:6]}; %IDs by to group the clusters to larger ones
new_cluster_assignment = reassign_cluster_ids(cluster_assignment,new_cluster_groups);
cluster_idx_array = {age_idx,mut_idx,new_cluster_assignment}; %1.index for grouping, 2.index directly compared, bin_clust must be last
[group_idx, group_count, group_comb] = groupby_multiple_idx(cluster_idx_array);
N_bars = length(group_count)/length(unique(new_cluster_assignment));
total_units = arrayfun(@(x) sum(group_count(x:N_bars:end)),1:N_bars);
ratios = arrayfun(@(x) group_count(x:N_bars:end)./total_units(x),1:N_bars,'un',0);

%%
bar_plot_data = vertcat(ratios{:});
figure('Color','w');
bar(bar_plot_data*100,'stacked')
box off

ylabel('Putative DA neurons [%]')
xlabel('Age')
legend(["WT","A53T"],'Location','northwest','Box','off')
xticklabels(unique(age_idx))

%% Plot single cell densities
n_bins = 50;
cluster_idx_array = {age_idx,mut_idx};
cluster_name_array = ["Age","Mutation"];
cluster_density_array = plot_cluster_densities(reduction,cluster_idx_array,cluster_name_array,n_bins);

%% Plot cluster difference
ids = [1,2];
plot_cluster_diffs(cluster_density_array,ids)

%% Feature comparisons
feat = 'ResonanceFrequency';
cluster_feature = arrayfun(@(x) nanmean(rmoutliers([unit_array(group_idx==x).(feat)])),unique(group_idx));
cluster_std = arrayfun(@(x) nanstd(rmoutliers([unit_array(group_idx==x).(feat)])),unique(group_idx));
plot_data = reshape(cluster_feature,length(unique(age_idx)),[]);
plot_std = reshape(cluster_std,length(unique(age_idx)),[]);
figure;
% plot(plot_data); 
errorbar(plot_data,plot_std)
legend(["WT 1","A53T 1","WT 2","A53T 2"])

%% 
feat = 'Halfwidth';
feature_sc = arrayfun(@(x) rmoutliers([unit_array(cluster_idx==x).(feat)]),unique(cluster_idx),'un',0);

feat_sel = feature_sc(4:N_bars:end);
x = arrayfun(@(x) ones(1,length(feat_sel{x}))*x,1:length(feat_sel),'un',0);
y = [feat_sel{:}];
figure;
swarmchart([x{:}],y,'.')

%% 
y1 = [feat_sel{[1]}];
y2 = [feat_sel{[3]}];
x = [ones(1,length(y1)),ones(1,length(y2))*2];
figure;
boxchart(x,[y1,y2])
% [h,p] = ttest2(y1,y2)
%%
cluster_idx_array = {mut_idx,new_cluster_assignment}; %1.index for grouping, 2.index directly compared, bin_clust must be last
[group_idx, group_count, group_comb] = groupby_multiple_idx(cluster_idx_array);
features = [sc_prop, wf_prop];
feat_array = nan(length(unique(group_idx)),length(features));

for feat = 1:length([sc_prop, wf_prop])
    feature_sc = arrayfun(@(x) rmoutliers([unit_array(group_idx==x).(features(feat))]),unique(group_idx),'un',0);
    clusters = feature_sc;%arrayfun(@(x) [feature_sc{group_comb(end,:)==x}],unique(group_comb(end,:)),'un',0);
%     norm_data = [clusters{:}]/nanmean([clusters{:}]);
    norm_data = normalize([clusters{:}],'range',[0 1]);
    norm_idx = arrayfun(@(x) ones(1,length(clusters{x}))*x,1:length(clusters),'un',0);
    norm_clusters = arrayfun(@(x) norm_data([norm_idx{:}]==x),1:length(clusters),'un',0);
    feat_array(:,feat) = cellfun(@nanmean,norm_clusters');
end
%%
feat_array = feat_array./feat_array(end,:);
% feat_array(feat_array>3) = 3;
figure('Color','w','NextPlot','add');
for k = 1:size(feat_array,1)
    polarplot(linspace(0,2*pi,size(feat_array,2)+1),[feat_array(k,:) feat_array(k,1)],'Marker','.')
    hold on
end

polarhistogram([0.5 0.5 0.5],'BinLimits',[0 pi*3/4],'FaceAlpha',0.2,'EdgeAlpha',0,'FaceColor','k')
polarhistogram([1,2],'BinLimits',[pi*3/4 pi*2],'FaceAlpha',0.2,'EdgeAlpha',0,'FaceColor','k')
set(gca,'ThetaTick',linspace(0,360,size(feat_array,2)+1),'ThetaTickLabel',features,'RLim',[0 3])