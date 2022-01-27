function [rec,cl] = groupRecordingsByDate(nw_array)
    tps = unique([nw_array.DIV]);
    rec = cell(1,length(tps));
    cl = cell(1,length(tps));
    for t = 1:length(tps)
       rec{t} = findobj(nw_array,'DIV',tps(t))';
       cl{t} = [rec{t}.CellLine];
    end
end