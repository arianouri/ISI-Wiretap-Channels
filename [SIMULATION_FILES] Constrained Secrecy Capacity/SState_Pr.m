function [SSPr] = SState_Pr(TRANS)
FRACT=[sum(TRANS,3)]';
FRACT=FRACT-diag(ones(size(FRACT,1),1));
FRACT(1,:)=1+FRACT(1,:);
EVAL=zeros(size(FRACT,1),1);
EVAL(1)=1;
SSPr=(FRACT)\(EVAL);
end

