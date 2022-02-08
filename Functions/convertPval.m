function pval_out = convertPval(pval,p_adj)
pval = pval*p_adj;
if pval < 0.0001
    pval_out = "<0.0001";
elseif pval >= 1
    pval_out = ">0.9999";
else
    pval_out = string(sprintf('%.4f',pval));
end