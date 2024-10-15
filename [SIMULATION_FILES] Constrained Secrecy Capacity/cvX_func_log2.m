function EVAL=cvX_func_log2(Qij,SSPr,PASTQij,PASTSSPr,kappa,kappaPrime)
deltaQij=(Qij-PASTQij+eps)./(PASTQij+eps);
deltaSSPr=(SSPr-PASTSSPr+eps)./(PASTSSPr+eps);
EVAL=kappaPrime*(...
    sum(sum(sum(PASTQij.*(1+kappa.*deltaQij).*log2z(1+kappa.*deltaQij))))...
    +   sum(PASTSSPr.*(1+kappa.*deltaSSPr).*log2z(1+kappa.*deltaSSPr)));
end