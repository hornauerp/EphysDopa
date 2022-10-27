function run_mlm(nw_array,group_names,inclusion_array,exclusion_array,eqLNAs,factor_name,sc_idx,nw_idx,th,F_precision,MSD_precision)
arguments
   nw_array
   group_names (1,2)%Name of the two groups to be compared
   inclusion_array {iscell}
   exclusion_array {iscell}
   eqLNAs = 0
   factor_name = "Separating factor"
   sc_idx = [4,8,11:14,15,17,19,22:27,46:48]
   nw_idx = [17,18,20:31]
   th = 3
   F_precision = 2
   MSD_precision = 3
end
%% Apply inclusion/exclusion criteria
nw_array = convert_ntLNA(eqLNAs,nw_array);
dataset_1 = filterObjectArray(nw_array,inclusion_array{1},exclusion_array{1});
dataset_2 = filterObjectArray(nw_array,inclusion_array{1},exclusion_array{1});
tps = intersect(unique([dataset_1.DIV]),unique([dataset_2.DIV]));
n_tp = length(tps);

%% Select relevant features
sc_vars = properties(Unit);
sc_vars =  sc_vars(sc_idx);
nw_vars = properties(WholeNetwork);
nw_vars = nw_vars(nw_idx);

%% Preallocation
%sep_factor/Treatment
g_mat_sc = nan(1,length(sc_vars));
gt_mat_sc = nan(1,length(sc_vars));
g_mat_nw = nan(1,length(nw_vars));
gt_mat_nw = nan(1,length(nw_vars));
g_anova_sc = cell(1,length(nw_vars));
g_anova_nw = cell(1,length(nw_vars));

%Time
t_wt_nw = nan(1,length(nw_vars));
t_anova_wt_nw = cell(1,length(nw_vars));
t_wt_sc = nan(1,length(nw_vars));
t_anova_wt_sc = cell(1,length(nw_vars));
t_mut_nw = nan(1,length(nw_vars));
t_anova_mut_nw = cell(1,length(nw_vars));
t_mut_sc = nan(1,length(nw_vars));
t_anova_mut_sc = cell(1,length(nw_vars));

%% Single cell variables
for var_idx = 1:length(sc_vars)
    %Dataset 1
    var = sc_vars{var_idx};
    dataset_1_ids = unique({dataset_1.ChipID});
    dataset_1_mat = nan(n_tp,20);
    for t = 1:n_tp
        time_group = findobj(dataset_1,'DIV',tps(t));
        for i = 1:length(dataset_1_ids)
            id_group = findobj(time_group,'ChipID',dataset_1_ids{i});
            if ~isempty(id_group)
                dataset_1_mat(t,i) = rmoutliers(arrayfun(@(x) mean([x.Units.(var)],'omitnan'),id_group),'ThresholdFactor',th);
            end
        end
        dataset_1_mat(t,isoutlier(dataset_1_mat(t,:),'median','ThresholdFactor',3)) = NaN;
    end
    
    %Dataset 1 // temporal model
    week = categorical(repmat([1:n_tp]',[size(dataset_1_mat,2) 1]));
    subject = categorical(repelem(1:size(dataset_1_mat,2),n_tp)');
    sep_factor = categorical([zeros(numel(dataset_1_mat),1)]);
    y = reshape(dataset_1_mat,[],1);
    tbl = table(y,week,subject,sep_factor,'VariableNames',["y","Week","Subject",factor_name]);
    lme = fitlme(tbl,'y~ Week+ (1|Subject)','FitMethod','REML','DummyVarCoding','effects');
    anova_results = anova(lme,'DFMethod','satterthwaite');
    t_wt_sc(var_idx) = anova_results.pValue(2);
    t_anova_wt_sc{var_idx} = anova_results;
    
    %Dataset 2
    dataset_2_ids = unique({dataset_2.ChipID});
    dataset_2_mat = nan(length(unique([dataset_2.DIV])),20);
    for t = 1:n_tp
        time_group = findobj(dataset_2,'DIV',tps(t));
        for i = 1:length(dataset_2_ids)
            id_group = findobj(time_group,'ChipID',dataset_2_ids{i});
            if ~isempty(id_group)
                dataset_2_mat(t,i) = rmoutliers(arrayfun(@(x) mean([x.Units.(var)],'omitnan'),id_group),'ThresholdFactor',th);
            end
        end
        dataset_2_mat(t,isoutlier(dataset_2_mat(t,:),'median','ThresholdFactor',3)) = NaN;
    end
    %Dataset 2 temporal model
    week = categorical(repmat([1:n_tp]',[size(dataset_2_mat,2) 1]));
    subject = categorical(repelem(1:size(dataset_2_mat,2),n_tp)');
    sep_factor = categorical([zeros(numel(dataset_2_mat),1)]);
    y = reshape(dataset_2_mat,[],1);
    tbl = table(y,week,subject,sep_factor,'VariableNames',["y","Week","Subject",factor_name]);
    lme = fitlme(tbl,'y~ Week+ (1|Subject)','FitMethod','REML','DummyVarCoding','effects');
    anova_results = anova(lme,'DFMethod','satterthwaite');
    t_mut_sc(var_idx) = anova_results.pValue(2);
    t_anova_mut_sc{var_idx} = anova_results;
    
    %Mixed model
    mat = [dataset_1_mat dataset_2_mat];
    week = categorical(repmat([1:n_tp]',[size(mat,2) 1]));
    subject = categorical(repelem(1:size(mat,2),n_tp)');
    sep_factor = categorical([zeros(numel(dataset_1_mat),1);ones(numel(dataset_2_mat),1)]);
    y = reshape(mat,[],1);
    tbl = table(y,week,subject,sep_factor,'VariableNames',["y","Week","Subject",factor_name]);
    try
        lme = fitlme(tbl,'y~ Week*sep_factor + (1|Subject)','FitMethod','REML','DummyVarCoding','effects');
        anova_results = anova(lme,'DFMethod','satterthwaite');
        g_mat_sc(var_idx) = anova_results.pValue(3);
        gt_mat_sc(var_idx) = anova_results.pValue(4);
        g_anova_sc{var_idx} = anova_results;
    catch
        sprintf('%s not fittable',var)
    end
end

%% Network features 
for var_idx = 1:length(nw_vars)
    if include_batch == 3
        group = findobj(nw_array,'Treatment','-');
    else
        group = findobj(nw_array,'Mutation','wt');
    end
    var = nw_vars{var_idx};
    dataset_1_ids = unique({group.ChipID});
    dataset_1_mat = nan(length(unique([group.DIV])),20);
    for t = 1:n_tp
        time_group = findobj(group,'DIV',tps(t));
        for i = 1:length(dataset_1_ids)
            id_group = findobj(time_group,'ChipID',dataset_1_ids{i});
            if ~isempty(id_group)
                dataset_1_mat(t,i) = id_group.(var);
            end
        end
        dataset_1_mat(t,isoutlier(dataset_1_mat(t,:),'median','ThresholdFactor',3)) = NaN;
    end
    %WT temporal model
    mat = dataset_1_mat;
    week = categorical(repmat([1:n_tp]',[size(mat,2) 1]));
    subject = categorical(repelem(1:size(mat,2),n_tp)');
    sep_factor = categorical([zeros(numel(mat),1)]);
    y = reshape(mat,[],1);
    tbl = table(y,week,subject,sep_factor,'VariableNames',["y","Week","Subject",factor_name]);
    lme = fitlme(tbl,'y~ Week+ (1|Subject)','FitMethod','REML','DummyVarCoding','effects');
    anova_results = anova(lme,'DFMethod','satterthwaite');
    t_wt_nw(var_idx) = anova_results.pValue(2);
    t_anova_wt_nw{var_idx} = anova_results;
    
    
    dataset_2_ids = unique({dataset_2.ChipID});
    dataset_2_mat = nan(length(unique([group.DIV])),20);
    for t = 1:n_tp
        time_group = findobj(group,'DIV',tps(t));
        for i = 1:length(dataset_2_ids)
            id_group = findobj(time_group,'ChipID',dataset_2_ids{i});
            if ~isempty(id_group)
                dataset_2_mat(t,i) = id_group.(var);
            end
        end
        dataset_2_mat(t,isoutlier(dataset_2_mat(t,:),'median','ThresholdFactor',3)) = NaN;
    end
    %A53T temporal model
    week = categorical(repmat([1:n_tp]',[size(dataset_2_mat,2) 1]));
    subject = categorical(repelem(1:size(dataset_2_mat,2),n_tp)');
    sep_factor = categorical([zeros(numel(dataset_2_mat),1)]);
    y = reshape(dataset_2_mat,[],1);
    tbl = table(y,week,subject,sep_factor,'VariableNames',["y","Week","Subject",factor_name]);
    lme = fitlme(tbl,'y~ Week+ (1|Subject)','FitMethod','REML','DummyVarCoding','effects');
    anova_results = anova(lme,'DFMethod','satterthwaite');
    t_mut_nw(var_idx) = anova_results.pValue(2);
    t_anova_mut_nw{var_idx} = anova_results;
    
    %Mixed model
    mat = [dataset_1_mat dataset_2_mat];
    week = categorical(repmat([1:n_tp]',[size(mat,2) 1]));
    subject = categorical(repelem(1:size(mat,2),n_tp)');
    sep_factor = categorical([zeros(numel(dataset_1_mat),1);ones(numel(dataset_2_mat),1)]);
    y = reshape(mat,[],1);
    tbl = table(y,week,subject,sep_factor,'VariableNames',["y","Week","Subject",factor_name]);
    try
        lme = fitlme(tbl,'y~ Week*sep_factor + (1|Subject)','FitMethod','REML','DummyVarCoding','effects');
        anova_results = anova(lme,'DFMethod','satterthwaite');
        if any(isnan(anova_results.pValue))
            anova_results = anova(lme);
        end
        g_mat_nw(var_idx) = anova_results.pValue(3);
        gt_mat_nw(var_idx) = anova_results.pValue(4);
        g_anova_nw{var_idx} = anova_results;
    catch
        sprintf('%s not fittable',var)
    end
end %Variables
%% Write statistics table
sc_var_names = {'AMPL','HLFW','ASYM','T2PR','T2PD','AUCP1','AUCP2','AUCT','RISE','DECAY','ISIM','ISIV','ISICV','PACF','SCRF','SCRM','RFIT'};
nw_vars = {'IBIV','IBIM','MBD','VBD','INTRABF','INTERBF','BRT','BRV','BDT','BDV','SYNC','NRF','NRM','NRFIT'};
vars = [sc_var_names nw_vars];  labels = vars;
vars_units = [[' [' char(181) 'V]']," [ms]","",""," [ms]","","","",[' [' char(181) 'V/ms]'],...
    [' [' char(181) 'V/ms]']," [s]","","",""," [Hz]","","",""," [s]",...
    " [s]",""," [Hz]"," [Hz]"," [s]"," [coactivity/s]"," [s]"," [coactivity/s]",""," [Hz]",...
    "",""];
scaling = ones(1,length(labels)); scaling([1,7:11]) = 6.2;
stats_table = table('Size',[length(labels)*12,3],'VariableTypes',["string","string","string"]);
var_counter = 1;
p_adj = length(vars);

for n = var_counter:length(labels)
    start_idx = (n-var_counter)*(7+n_tp);
    stats_table(start_idx+1,1).Variables = string(sprintf('%s%s',string(labels{n}),vars_units(n)));
    
    for t = 1:n_tp
        stats_table(start_idx+1+t,1).Variables = string(sprintf('Week %i',t));
    end
    stats_table(start_idx+1,2).Variables = string(group_names{1});
    stats_table(start_idx+1,3).Variables = string(group_names{2});
    stats_table(start_idx+2+n_tp,1).Variables = "Mixed model effects";
    stats_table(start_idx+2+n_tp,2).Variables = "p-Values";
    stats_table(start_idx+2+n_tp,3).Variables = "F statistics";
    stats_table(start_idx+3+n_tp,1).Variables = string(['Time ' group_names{1}]);
    stats_table(start_idx+4+n_tp,1).Variables = string(['Time ' group_names{2}]);
    stats_table(start_idx+5+n_tp,1).Variables = sep_factor;
    stats_table(start_idx+6+n_tp,1).Variables = [sep_factor " x Time"];
    
    if ismember(vars(n),sc_var_names)
        mean_val_wt = arrayfun(@(x) formatMeanSD(tmp_matrix(2,x,n),tmp_sd(2,x,n),MSD_precision,scaling(n)),1:size(tmp_matrix,2),'un',0);
        mean_val_a53t = arrayfun(@(x) formatMeanSD(tmp_matrix(1,x,n),tmp_sd(1,x,n),MSD_precision,scaling(n)),1:size(tmp_matrix,2),'un',0);
        stats_table(start_idx+2:start_idx+1+n_tp,2).Variables = mean_val_wt';
        stats_table(start_idx+2:start_idx+1+n_tp,3).Variables = mean_val_a53t';
        stats_table(start_idx+3+n_tp,2).Variables = convertPval(t_wt_sc(n),p_adj);
        stats_table(start_idx+3+n_tp,3).Variables = formatFstats(t_anova_wt_sc{n},2,F_precision);
        stats_table(start_idx+4+n_tp,2).Variables = convertPval(t_mut_sc(n),p_adj);
        stats_table(start_idx+4+n_tp,3).Variables = formatFstats(t_anova_mut_sc{n},2,F_precision);
        stats_table(start_idx+5+n_tp,2).Variables = convertPval(g_mat_sc(n),p_adj);
        stats_table(start_idx+5+n_tp,3).Variables = formatFstats(g_anova_sc{n},3,F_precision);
        stats_table(start_idx+6+n_tp,2).Variables = convertPval(gt_mat_sc(n),p_adj);
        stats_table(start_idx+6+n_tp,3).Variables = formatFstats(g_anova_sc{n},4,F_precision);
    else
        mapping_idx = n-length(sc_var_names); %Linear model mapping
        val_mapping = mapping_idx+2;
        mean_val_wt = arrayfun(@(x) formatMeanSD(nw_matrix(2,x,val_mapping),nw_sd(2,x,val_mapping),MSD_precision,scaling(n)),1:size(nw_matrix,2),'un',0);
        mean_val_a53t = arrayfun(@(x) formatMeanSD(nw_matrix(1,x,val_mapping),nw_sd(1,x,val_mapping),MSD_precision,scaling(n)),1:size(nw_matrix,2),'un',0);
        stats_table(start_idx+2:start_idx+1+n_tp,2).Variables = mean_val_wt';
        stats_table(start_idx+2:start_idx+1+n_tp,3).Variables = mean_val_a53t';
        stats_table(start_idx+3+n_tp,2).Variables = convertPval(t_wt_nw(mapping_idx),p_adj);
        stats_table(start_idx+3+n_tp,3).Variables = formatFstats(t_anova_wt_nw{mapping_idx},2,F_precision);
        stats_table(start_idx+4+n_tp,2).Variables = convertPval(t_mut_nw(mapping_idx),p_adj);
        stats_table(start_idx+4+n_tp,3).Variables = formatFstats(t_anova_mut_nw{mapping_idx},2,F_precision);
        stats_table(start_idx+5+n_tp,2).Variables = convertPval(g_mat_nw(mapping_idx),p_adj);
        stats_table(start_idx+5+n_tp,3).Variables = formatFstats(g_anova_nw{mapping_idx},3,F_precision);
        stats_table(start_idx+6+n_tp,2).Variables = convertPval(gt_mat_nw(mapping_idx),p_adj);
        stats_table(start_idx+6+n_tp,3).Variables = formatFstats(g_anova_nw{mapping_idx},4,F_precision);
    end
end
if include_batch == 3
    stats_table(1:110,:) = [];
end
% writetable(stats_table,fullfile(save_path,'SupplementaryTable.xlsx'),'WriteVariableNames',0)