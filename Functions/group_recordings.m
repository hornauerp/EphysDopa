function [grouping_values,group_combinations,n_groups] = group_recordings(nw_array,grouping_conditions)
grouping_values = cell(1,length(grouping_conditions));
for i = 1:length(grouping_conditions)
    value_array = {nw_array.(grouping_conditions{i})};
    if ischar(value_array{1})
        grouping_values{i} = unique(value_array);
    else
        grouping_values{i} = unique([value_array{:}]);
    end
end

group_vectors = cellfun(@(x) 1:length(x),grouping_values,'un',0);
group_combinations = num2cell(combvec(group_vectors{:}));
n_groups = size(group_combinations,2);