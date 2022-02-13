function [nw_mat,nw_prop,mutation,batch,treatment] = generateNWMat(rec,PCA_IND,var_sel)
arguments
    rec
    PCA_IND (1,1) {isnumeric} = 0
    var_sel = []
end
nw_prop = ["IBIVar",...
    "IBIMean",...
    "BDMean",...
    "BDVar",...
    "IntraBF",...
    "InterBF",...
    "RiseTime",...
    "RiseVelocity",...
    "DecayTime",...
    "DecayVelocity",...
    "Synchronicity",...
    "RegularityFrequency",...
    "RegularityMagnitude",...
    "RegularityFit"];

% nw_idx = [16:17,19:30];%11,12,
% nw_vars = properties(WholeNetwork);
% nw_prop = nw_vars(nw_idx);
if ~isempty(var_sel)
   nw_prop = nw_prop(var_sel); 
end

mutation = zeros(1,length(rec));
batch = zeros(1,length(rec));
treatment = zeros(1,length(rec));

nw = [rec{:}];
plating_dates = unique([nw.PlatingDate]);
if PCA_IND == 0 || length(plating_dates) == 1
    nw_mat = zeros(length(rec),length(rec{1})*length(nw_prop));
for i = 1:length(rec)
    crec = rec{i};
    vec = [];
    for j = 1:numel(nw_prop)
        vec= [vec [crec.(nw_prop{j})]];
    end
    nw_mat(i,:) = vec;

end
else
    nw_mat = zeros(length(rec),length(nw_prop));
    for j = 1:numel(nw_prop)
        vec = zeros(length(rec),length(rec{1}));
        for i = 1:length(rec)
            crec = rec{i};
            vec(i,:)= [crec.(nw_prop{j})];
        end
        [~,score,~,~,explained] = pca(vec);
        nw_mat(:,j) = score(:,1);
    end
end
for i = 1:length(rec)
    crec = rec{i};
    mutation(i) = crec(1).Mutation ~="wt";
    batch(i) = find(crec(1).PlatingDate==plating_dates);
    treatment(i) = isequal(crec(1).Treatment,'LNA7');
end
nw_mat(isnan(nw_mat))=0;
end