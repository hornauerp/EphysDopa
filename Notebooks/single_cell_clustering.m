addpath(genpath('/cluster/home/phornauer/Git/EphysDopa'))
addpath(genpath('/cluster/home/phornauer/Git/drtoolbox'))

%% Load full dataset
load('F:\Promotion\Git\EphysDopa\Data\full_dataset.mat')

%% Select data
batches = unique([nw_array.PlatingDate]); %Refer to different batches via their PlatingDate property
inclusion = {{'PlatingDate',batches([3])},{'DIV',6,13,20,27,34}}; %e.g. specify the batches to analyze by indexing the batches array // here we use the first two
exclusion = {};
sel_nw = filterObjectArray(nw_array,inclusion,exclusion);
Unit_array = arrayfun(@(x) sel_nw(x).Units(sel_nw(x).ActiveChannels),1:length(sel_nw),'UniformOutput',false);
Unit_array = [Unit_array{:}];
mut_idx = double(arrayfun(@(x) contains(Unit_array(x).Network.Mutation,'a53t'),1:length(Unit_array)));
age_idx = arrayfun(@(x) Unit_array(x).Network.DIV,1:length(Unit_array));
treat_idx = double(arrayfun(@(x) contains(Unit_array(x).Network.Treatment,'LNA7'),1:length(Unit_array)));

%% Select features to be included in the single-cell clustering
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

%% Generate matrix of activity features
% feat_wt = generate_sc_feature_matrix(sel_nw,sc_prop,norm_idx);
% feat_a53t = generate_sc_feature_matrix(nw_a53t,sc_prop,norm_idx);

%% Full matrices
[wf_full,ax_full] = generate_wf_matrix(Unit_array,norm_idx,rm_axon);
feat_full = generate_sc_feature_matrix(Unit_array,sc_prop,norm_idx);
feat_full(ax_full,:) = [];
feat_mat = feat_full;
norm_mat = wf_full;
mut_idx(ax_full) = [];
age_idx(ax_full) = [];
treat_idx(ax_full) = [];
% %% Generate waveform matrix
% [wf_wt,ax_wt] = generate_wf_matrix(nw_wt,norm_idx,rm_axon);
% [wf_a53t,ax_a53t] = generate_wf_matrix(nw_a53t,norm_idx,rm_axon);
% feat_wt(ax_wt,:) = [];
% feat_a53t(ax_a53t,:) = [];
% mut_idx = [zeros(1,size(wf_wt,1)), ones(1,size(wf_a53t,1))];
% norm_mat = [wf_wt; wf_a53t];
% feat_mat = [feat_wt; feat_a53t];

%% Estimate intrinsic dimensionality and reduce accordingly
no_dims = intrinsic_dim(feat_mat);
[mapped_feat, ~] = compute_mapping(feat_mat, 'tSNE', floor(no_dims));

%% Cluster activity features
k = 4;
idx = kmeans(mapped_feat,k);
%%
plot_cluster_scatter(mapped_feat,idx)

%% Estimate intrinsic dimensionality and reduce accordingly
no_dims = intrinsic_dim(norm_mat);
[mapped_wf, ~] = compute_mapping(norm_mat, 'tSNE', round(no_dims));

%% Cluster waveform features
k = 4;
idx = kmeans(mapped_wf,k);

%% Plot first 3 dims
plot_cluster_scatter(mapped_wf,idx)

%% Plot waveform clusters
plot_waveform_clusters(norm_mat,idx)

%% Combine activity and waveform
comb_mat = [mapped_wf mapped_feat];
no_dims = intrinsic_dim(comb_mat);
[mapped_comb, ~] = compute_mapping(comb_mat, 'tSNE', floor(no_dims));

%% Find optimal number of clusters
max_clust = 10;
silh_array = nan(1,max_clust-1);
for nc = 2:10
   idx = kmeans(mapped_comb,nc); 
   s = silhouette(mapped_comb,idx);
   silh_array(nc) = mean(s);
end

%% Cluster combined features 
k = 2;
idx = kmeans(mapped_comb,k);
plot_cluster_scatter(mapped_comb,idx)

%% Plot waveform clusters
plot_waveform_clusters(norm_mat,idx)

%% UMAP clustering
rescaled_feat = normalize(feat_mat,'range',[-1 1]);
full = [norm_mat rescaled_feat];
[reduction, umap, clusterIdentifiers, extras]=run_umap(full,'n_components',2,'n_neighbors',100,'min_dist',0.1,'cluster_detail','adaptive','spread',1,'sgd_tasks',20);

%% UMAP plot
% idx = spectralcluster(reduction,2);
plot_cluster_scatter(reduction,clusterIdentifiers)
% colormap(othercolor('RdYlBu9',100))
colorbar
plot_waveform_clusters(norm_mat,clusterIdentifiers)

%% Compare mutations
figure('Color','w');
subplot(1,2,1)
wt_data = [reduction(mut_idx==0,1),reduction(mut_idx==0,2)];
values_wt = hist3(wt_data,[101 101]);
i1 = imagesc(values_wt.'./max(max(values_wt)));
cb = colorbar('LimitsMode','manual');
% scatter(reduction(mut_idx==0,1),reduction(mut_idx==0,2),1,'b','filled')
subplot(1,2,2)
a53t_data = [reduction(mut_idx==1,1),reduction(mut_idx==1,2)];
values = hist3(a53t_data,[101 101]);
i2 = imagesc(values.'./max(max(values_wt)));
colorbar
% scatter(reduction(mut_idx==1,1),reduction(mut_idx==1,2),1,'r','filled')
%% Apply UMAP Unit to new data
[reduction, umap, clusterIdentifiers, extras]=run_umap(wf_a53t,'Unit_file','umap.mat','n_components',2,'n_neighbors',20,'min_dist',0.001,'cluster_detail','adaptive','spread',1,...
    'set_op_mix_ratio',0.5,'target_weight',0);

%% UMAP plot
plot_cluster_scatter(reduction,clusterIdentifiers)
plot_waveform_clusters(norm_mat,clusterIdentifiers)