addpath(genpath('/cluster/home/phornauer/Git/EphysDopa'))
addpath(genpath('/cluster/home/phornauer/Git/drtoolbox'))
addpath(genpath('/cluster/home/phornauer/Git/UMAP/'))
%% 
params.QC.amp = nan; %Amplitude thresholds to exclude units, either [min_amp max_amp] or nan for no amplitude filtering
params.QC.rate = [0.1 10]; %Firing rate thresholds to exclude units, either [min_rate max_rate] or nan for no firing rate filtering
params.QC.rv = 0.02; %Maximum refractory period violations in [%], default is 0.02
params.Bursts.N = 0.0015; %Parameters for burst detection, for details see Bakkum et al., 2014 (https://doi.org/10.3389/fncom.2013.00193)
params.Bursts.ISI_N = 1.5;
params.Bursts.merge_t = 2;
params.Regularity.binning = 0.1; %Binning for regularity calculations
params.Bursts.binning = 0.1; %Binning for burst quantifications

%%
sorting_path = '/cluster/project/bsse_sdsc/BELUB/Measurements/Cell_line_9/200217/5206/Sorting/traces/traces-merged.GUI';
nw_wt = WholeNetwork(sorting_path,params);
sorting_path = '/cluster/project/bsse_sdsc/BELUB/Measurements/Cell_line_10/200217/5215/Sorting/traces/traces-merged.GUI';
nw_a53t = WholeNetwork(sorting_path,params);

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
feat_wt = generate_sc_feature_matrix(nw_wt,sc_prop,norm_idx);
feat_a53t = generate_sc_feature_matrix(nw_a53t,sc_prop,norm_idx);

%% Generate waveform matrix
[wf_wt,ax_wt] = generate_wf_matrix(nw_wt,norm_idx,rm_axon);
[wf_a53t,ax_a53t] = generate_wf_matrix(nw_a53t,norm_idx,rm_axon);
feat_wt(ax_wt,:) = [];
feat_a53t(ax_a53t,:) = [];
mut_idx = [zeros(1,size(wf_wt,1)), ones(1,size(wf_a53t,1))];
norm_mat = [wf_wt; wf_a53t];
feat_mat = [feat_wt; feat_a53t];

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
[reduction, umap, clusterIdentifiers, extras]=run_umap(wf_wt,'n_components',2,'n_neighbors',20,'min_dist',0.001,'cluster_detail','adaptive','spread',1,...
    'set_op_mix_ratio',0.5,'target_weight',0);

%% UMAP plot
plot_cluster_scatter(reduction,clusterIdentifiers)
plot_waveform_clusters(norm_mat,clusterIdentifiers)

%% Apply UMAP template to new data
[reduction, umap, clusterIdentifiers, extras]=run_umap(wf_a53t,'template_file','umap.mat','n_components',2,'n_neighbors',20,'min_dist',0.001,'cluster_detail','adaptive','spread',1,...
    'set_op_mix_ratio',0.5,'target_weight',0);

%% UMAP plot
plot_cluster_scatter(reduction,clusterIdentifiers)
plot_waveform_clusters(norm_mat,clusterIdentifiers)
%%
function plot_waveform_clusters(norm_mat,idx)
figure;
tiledlayout('flow')
for i = 1:max(idx)
    nexttile
    plot(norm_mat(idx==i,:)','Color',[0.8 0.8 0.8 0.5],'LineWidth',0.01)
    hold on
    plot(mean(norm_mat(idx==i,:)),'k')
end
nexttile
for i = 1:max(idx)
   plot(mean(norm_mat(idx==i,:)))
   hold on
end
end
%%
function plot_cluster_scatter(mapped_data,idx)
figure;
if size(mapped_data,2) == 2
    scatter(mapped_data(:,1),mapped_data(:,2),10,idx,'filled')
else
    scatter3(mapped_data(:,1),mapped_data(:,2),mapped_data(:,3),10,idx,'filled')
end
end
