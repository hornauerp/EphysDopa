function [sc_mat,wf_mat,mc_mat,act_mat,sc_prop,wf_prop,mc_prop,act_prop,cellline,batch] = generateSCMat(rec,MET,TH,PCA_IND,var_sel)
arguments
   rec
   MET
   TH
   PCA_IND (1,1) {isnumeric} = 0
   var_sel = []
end
sc_prop = ["Amplitude",...
    "Frequency",...
    "Halfwidth",...
    "Asymmetry",...
    "T2Pratio",...
    "T2Pdelay",...
    "AUCpeak1",...
    "AUCpeak2",...
    "AUCtrough",...
    "Rise",...
    "Decay",...
    "ISIMean",...
    "ISIVar",...
    "Fano",...
    "PACF",...
    "ResonanceFrequency",...
    "ResonanceMagnitude",...
    "ResonanceFit"
    ];

% tmp_vars = properties(Template);
% tmp_idx = [4,11:15,17,19,22:27,46:48];
% sc_prop = tmp_vars(tmp_idx);
if ~isempty(var_sel)
wf_idx = find(ismember(1:10,var_sel));
act_idx = find(ismember(11:11+length([11:14,20:22]),var_sel));
else
    wf_idx = 1:10;
    act_idx = 11:17;
end
wf_prop = sc_prop(wf_idx);
% mc_prop = sc_prop([15:19]);%
mc_prop = [];
act_prop = sc_prop(act_idx);
 
wf_mat = []; act_mat = []; mc_mat = [];
cellline = zeros(1,length(rec));
batch = zeros(1,length(rec));
if PCA_IND==0
    sc_mat = zeros(length(rec),length(rec{1})*length(sc_prop));
for i = 1:length(rec)
    
    crec = rec{i};
    vec = [];
    for j = 1:numel(sc_prop)
        v = arrayfun(@(x) mean(rmoutliers([x.Templates.(sc_prop{j})],MET,'ThresholdFactor',TH),'omitnan'),crec);
        vec = [vec v];
    end
    sc_mat(i,:) = vec;
    vec = [];
    if ~isempty(wf_prop)
    for j = 1:numel(wf_prop)
        v = arrayfun(@(x) mean(rmoutliers([x.Templates.(wf_prop{j})],MET,'ThresholdFactor',TH),'omitnan'),crec);
        vec = [vec v];
    end
    wf_mat(i,:) = vec;
    end
    
    vec = [];
    if ~isempty(mc_prop)
    for j = 1:numel(mc_prop)
        v = arrayfun(@(x) mean(rmoutliers([x.Templates.(mc_prop{j})],MET,'ThresholdFactor',TH),'omitnan'),crec);
        vec = [vec v];
    end
    mc_mat(i,:)= vec;
    end
    vec = [];
    if ~isempty(act_prop)
    for j = 1:numel(act_prop)
        v = arrayfun(@(x) mean(rmoutliers([x.Templates.(act_prop{j})],MET,'ThresholdFactor',TH),'omitnan'),crec);
        vec = [vec v];
    end
    act_mat(i,:)= vec;
    end
end
else
    sc_mat = zeros(length(rec),length(sc_prop));
    for j = 1:numel(sc_prop)
        vec = zeros(length(rec),length(rec{1}));
        for i = 1:length(rec)
            crec = rec{i};
            v = arrayfun(@(x) mean(rmoutliers([x.Templates.(sc_prop{j})],MET,'ThresholdFactor',TH),'omitnan'),crec)';
            vec(i,:) = v;
        end
        [~,score,~,~,explained] = pca(vec);
        sc_mat(:,j) = score(:,1);
    end
    wf_mat = zeros(length(rec),length(wf_prop));
    for j = 1:numel(wf_prop)
        vec = zeros(length(rec),length(rec{1}));
        for i = 1:length(rec)
            crec = rec{i};
            v = arrayfun(@(x) mean(rmoutliers([x.Templates.(wf_prop{j})],MET,'ThresholdFactor',TH),'omitnan'),crec)';
            vec(i,:) = v;
        end
        [~,score,~,~,explained] = pca(vec);
        wf_mat(:,j) = score(:,1);
    end
    mc_mat = zeros(length(rec),length(mc_prop));
    for j = 1:numel(mc_prop)
        vec = zeros(length(rec),length(rec{1}));
        for i = 1:length(rec)
            crec = rec{i};
            v = arrayfun(@(x) mean(rmoutliers([x.Templates.(mc_prop{j})],MET,'ThresholdFactor',TH),'omitnan'),crec)';
            vec(i,:) = v;
        end
        [~,score,~,~,explained] = pca(vec);
        mc_mat(:,j) = score(:,1);
    end
    act_mat = zeros(length(rec),length(act_prop));
    for j = 1:numel(act_prop)
        vec = zeros(length(rec),length(rec{1}));
        for i = 1:length(rec)
            crec = rec{i};
            v = arrayfun(@(x) mean(rmoutliers([x.Templates.(act_prop{j})],MET,'ThresholdFactor',TH),'omitnan'),crec)';
            vec(i,:) = v;
        end
        [~,score,~,~,explained] = pca(vec);
        act_mat(:,j) = score(:,1);
    end
end
for i = 1:length(rec)
    crec = rec{i};
    cellline(i) = mod([crec(1).CellLine],2);
    batch(i) = ([crec(1).CellLine]>2)+1;
end
sc_mat(isnan(sc_mat))=0;
wf_mat(isnan(wf_mat))=0;
mc_mat(isnan(mc_mat))=0;
act_mat(isnan(act_mat))=0;
end