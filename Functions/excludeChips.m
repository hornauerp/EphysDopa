function nw_filtered = excludeChips(nw_array,ids)
for i = 1:length(ids)
   nw_array = findobj(nw_array,'-not','ChipID',string(ids(i)));
end
nw_filtered = nw_array;
end