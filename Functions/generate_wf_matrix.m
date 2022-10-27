function [wf_mat,ax_idx] = generate_wf_matrix(Unit_array,norm_idx,rm_axon)
wf_mat = vertcat(Unit_array.MaxWf);
wf_max = max(abs(wf_mat),[],2);
if norm_idx
    wf_mat = wf_mat./wf_max;
    if rm_axon
        ax_idx = max(wf_mat,[],2)==1;
        wf_mat(ax_idx,:) = [];
    end
end