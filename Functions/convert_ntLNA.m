function nw_array = convert_ntLNA(eqLNAs,nw_array)
if eqLNAs
    for i = 1:length(nw_array)
        if nw_array(i).Treatment == "ntLNA"
            nw_array(i).Treatment = '-';
        end
    end
else
    for i = 1:length(nw_array)
        if ismember(nw_array(i).CellLine,[11 12])
            nw_array(i).Treatment = 'ntLNA';
        end
    end
end