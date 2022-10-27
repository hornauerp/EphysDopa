function F_stats = formatFstats(anova_model,F_idx,precision)
f = sprintf('0.%df',precision);
F_stats = string(sprintf(['F(%' f ',%' f ') = %' f],anova_model.DF1(F_idx),anova_model.DF2(F_idx),anova_model.FStat(F_idx)));