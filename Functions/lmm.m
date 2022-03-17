function anova_results = lmm(mat_array)
arguments
    mat_array (1,2) {iscell}
end

mat = [mat_array{:}];
n_tp = size(mat,1);
week = categorical(repmat([1:n_tp]',[size(mat,2) 1]));
subject = categorical(repelem(1:size(mat,2),n_tp)');
sep_factor = categorical([zeros(numel(mat_array{1}),1);ones(numel(mat_array{2}),1)]);
y = reshape(mat,[],1);
tbl = table(y,week,sep_factor,subject,'VariableNames',["y","Week","sep_factor","Subject"]);
model_description = sprintf('%s ~ %s*%s + (1|%s)',string(tbl.Properties.VariableNames));
lme = fitlme(tbl,model_description,'FitMethod','REML','DummyVarCoding','effects');
anova_results = anova(lme,'DFMethod','satterthwaite');