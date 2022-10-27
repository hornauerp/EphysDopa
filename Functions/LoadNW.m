function nw = LoadNW(path,array_name,cls,subsampling,iter)

for i = 1:length(cls)
    switch nargin
        case 5
        fname = sprintf('Cell_line_%d/Iteration_%d/%d/%s',cls(i),iter,subsampling,array_name);
        case 4
        fname = sprintf('Cell_line_%d/%d/%s',cls(i),subsampling,array_name);
        case 3
        fname = sprintf('Cell_line_%d/%s',cls(i),array_name);
    end
    C(i)= load(fullfile(path,fname));
end

nw = [C(:).nw_array];
ind = arrayfun(@(x) isempty(x.Units),nw);
nw(ind) = [];
%Transform Fano to CV
for n = 1:length(nw)
   for t = [nw(n).Templates.ID]
       tmp = nw(n).Templates(t);
       tmp.Fano = sqrt(tmp.ISIVar)/tmp.ISIMean;
   end
end
end
