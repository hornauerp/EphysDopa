addpath(genpath('/home/phornauer/Git/EphysDopa'))
addpath(genpath('/home/phornauer/Git/drtoolbox'))
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
sorting_path = '/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/si_test/sorted/';
nw = WholeNetwork(sorting_path,params);

%% Select features to be included in the single-cell clustering
sc_prop = ["ISIMean",...
    "ISIVar",...
    "Fano",...
    "PACF",...
    "ResonanceFrequency",...
    "ResonanceMagnitude",...
    "ResonanceFit"
    ];

%% Generate matrix of activity features
feat_mat = zeros(length(sc_prop),length([nw.Templates.ID]));
for f = 1:length(sc_prop)
    feat_vals = [nw.Templates.(sc_prop(f))];
    feat_vals(isnan(feat_vals)) = 0;
    feat_mat(f,:) = feat_vals;
end
feat_mat = normalize(feat_mat');

%% Generate matrix of peak waveforms
wf_mat = vertcat(nw.Templates.MaxWf);
wf_max = max(abs(wf_mat),[],2);
norm_mat = wf_mat(:,20:end)./wf_max;

%% Remove axonal signal with positive peak
ax_idx = max(norm_mat,[],2)==1;
norm_mat(ax_idx,:) = [];
feat_mat(ax_idx,:) = [];

%% Estimate intrinsic dimensionality and reduce accordingly
no_dims = intrinsic_dim(feat_mat);
[mapped_feat, ~] = compute_mapping(feat_mat, 'tSNE', round(no_dims));

%% Cluster activity features
k = 4;
idx = kmeans(mapped_feat,k);
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
k = 6;
idx = kmeans(mapped_comb,k);
plot_cluster_scatter(mapped_comb,idx)

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