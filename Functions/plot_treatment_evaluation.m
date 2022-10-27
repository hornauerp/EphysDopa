function plot_treatment_evaluation(Y_test,treatment,pred)
mut = {'WT','A53T'};
treat = {'-','LNA'};
mat_out = zeros(length(mut),length(treat));
for m = 0:1
    for tr = 0:1
        cor_pred = sum(Y_test==m&treatment==tr&pred'==m);
        false_pred = sum(Y_test==m&treatment==tr&pred'~=m);
        mat_out(m+1,tr+1) = cor_pred/(cor_pred+false_pred);
    end
end
figure;
h = heatmap(treat,mut,mat_out);
h.XLabel = 'Treatment';
h.YLabel = 'Genotype';
h.Colormap = othercolor('RdYlGn8',100);
h.ColorLimits = [0 1];
title('Classification accuracy')