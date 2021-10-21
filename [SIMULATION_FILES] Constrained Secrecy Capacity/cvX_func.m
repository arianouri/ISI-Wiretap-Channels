function EVAL=cvX_func(Qij,SSPr,PASTQij,PASTSSPr,kappa,kappaPrime)
deltaQij=(Qij-PASTQij+eps)./(PASTQij+eps);
deltaSSPr=(SSPr-PASTSSPr+eps)./(PASTSSPr+eps);
EVAL=kappaPrime*(...
    sum(sum(sum(PASTQij.*(1+kappa.*deltaQij).*logz(1+kappa.*deltaQij))))...
    +   sum(PASTSSPr.*(1+kappa.*deltaSSPr).*logz(1+kappa.*deltaSSPr)));
end