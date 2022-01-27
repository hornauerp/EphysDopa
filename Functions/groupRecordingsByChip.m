function [rec,tp] = groupRecordingsByChip(nw_array,MIN_AGE,MAX_AGE)
nw = nw_array([nw_array.DIV]>=MIN_AGE&[nw_array.DIV]<=MAX_AGE);
ID = unique({nw.ChipID});
chips = cell(1,length(ID));
for i = 1:length(ID)
    chips{i} = findobj(nw,'ChipID',ID{i})';
end
measurements = cellfun(@(x) length(x),chips);
tp = max(measurements);%floor((MAX_AGE-MIN_AGE)/7);
% ids = cellfun(@(x) {x.ChipID},chips);
% DIV = cellfun(@(x) [x.DIV],chips);
rec = chips(measurements==tp);
