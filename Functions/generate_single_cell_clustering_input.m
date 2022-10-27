function [feat_mat,wf_full,mut_idx,age_idx,treat_idx,Unit_array] = generate_single_cell_clustering_input(nw_array,inclusion,exclusion,sc_prop,norm_idx,rm_axon,noise_cutoff,plot_noisy)
sel_nw = filterObjectArray(nw_array,inclusion,exclusion);
Unit_array = arrayfun(@(x) sel_nw(x).Units(sel_nw(x).ActiveChannels),1:length(sel_nw),'UniformOutput',false);
Unit_array = [Unit_array{:}];

mut_idx = double(arrayfun(@(x) contains(Unit_array(x).Network.Mutation,'a53t'),1:length(Unit_array)));
age_idx = arrayfun(@(x) Unit_array(x).Network.DIV,1:length(Unit_array));
treat_idx = double(arrayfun(@(x) contains(Unit_array(x).Network.Treatment,'LNA7'),1:length(Unit_array)));

% Full matrices
[wf_full,ax_full] = generate_wf_matrix(Unit_array,norm_idx,rm_axon);
feat_full = generate_sc_feature_matrix(Unit_array,sc_prop,norm_idx);
fprintf('Extracted features from %i units\n',length(Unit_array))

% Axon template removal
feat_full(ax_full,:) = [];
feat_mat = feat_full;
mut_idx(ax_full) = [];
age_idx(ax_full) = [];
treat_idx(ax_full) = [];
Unit_array(ax_full) = [];
fprintf('Removed %i axonal templates\n',sum(ax_full))

% Noise removal
noise_indicator = sum(abs(diff(diff(wf_full(:,15:50)')>0))); %Check for non-monotomies in waveform shapes as noise indicators
noise_idx = noise_indicator>noise_cutoff;
noisy_wf = wf_full(noise_idx,:);
feat_mat(noise_idx,:) = [];
wf_full(noise_idx,:) = [];
mut_idx(noise_idx) = [];
age_idx(noise_idx) = [];
treat_idx(noise_idx) = [];
Unit_array(noise_idx) = [];
fprintf('Removed %i noisy templates\n',sum(noise_idx))
if plot_noisy
    figure('Color','w');
    plot(noisy_wf','Color',[0.9 0.9 0.9])
end