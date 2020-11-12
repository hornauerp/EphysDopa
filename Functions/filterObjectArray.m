function array = filterObjectArray(nw_array,inclusion,exclusion)
% Input as cell arrays
for i = 1:length(inclusion)
    if isnumeric(inclusion{i}{2})
        value_vector = [nw_array.(inclusion{i}{1})];
        idx = ismember(value_vector,[inclusion{i}{2:end}]);
    else
        value_vector = {nw_array.(inclusion{i}{1})};
        idx = ismember(value_vector,inclusion{i}(2:end));
    end
    nw_array = nw_array(idx);
end

for i = 1:length(exclusion)
    if isnumeric(exclusion{i}{2})
        value_vector = [nw_array.(exclusion{i}{1})];
        idx = ismember(value_vector,[exclusion{i}{2:end}]);
    else
        value_vector = {nw_array.(exclusion{i}{1})};
        idx = ismember(value_vector,exclusion{i}(2:end));
    end
    nw_array = nw_array(~idx);
end
array = nw_array;