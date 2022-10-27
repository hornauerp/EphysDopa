function [tmp_matrix,tmp_sd,nw_matrix,nw_sd,leg,tps] = generate_feature_matrices(eqLNAs,nw_array,tmp_vars,nw_vars,grouping_conditions,th)
nw_array = convert_ntLNA(eqLNAs,nw_array);
[grouping_values,group_combinations,n_groups] = group_recordings(nw_array,grouping_conditions);
tps = unique([nw_array.DIV]);
tp = length(tps);
%% Pre-allocate output
tmp_matrix = nan(n_groups,tp,length(tmp_vars));
tmp_sd = nan(n_groups,tp,length(tmp_vars));
nw_matrix = nan(n_groups,tp,length(nw_vars));
nw_sd = nan(n_groups,tp,length(nw_vars));
leg = cell(size(group_combinations));
%% Main loop
for n = 1:n_groups
    vals = cellfun(@(x,y) x(y),grouping_values,group_combinations(:,n)','un',0); vals = [vals{:}];
    if length(grouping_conditions)==1
        group_cell_array = {[grouping_conditions,vals]};
    else
        grouping_mat = reshape([grouping_conditions,vals],length(grouping_conditions),[]);
        group_cell_array = arrayfun(@(x) {grouping_mat(x,:)},1:size(grouping_mat,1));
    end
    group_array = filterObjectArray(nw_array,group_cell_array,{});
    leg(:,n) = vals;
    for t = 1:tp
       tp_array = findobj(group_array,'DIV',tps(t));
       for v = 1:length(tmp_vars)
           tmp_matrix(n,t,v) = mean(rmoutliers(arrayfun(@(x) mean([x.Units.(tmp_vars{v})],'omitnan'),tp_array),'ThresholdFactor',th),'omitnan');
           tmp_sd(n,t,v) = std(rmoutliers(arrayfun(@(x) mean([x.Units.(tmp_vars{v})],'omitnan'),tp_array),'ThresholdFactor',th),'omitnan');
       end
       for v = 1:length(nw_vars)
           nw_matrix(n,t,v) = mean(rmoutliers([tp_array.(nw_vars{v})],'ThresholdFactor',th),'omitnan');
           nw_sd(n,t,v) = std(rmoutliers([tp_array.(nw_vars{v})],'ThresholdFactor',th),'omitnan');
       end
    end
end