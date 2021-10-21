function [outP] = modif_mod(inP1,inP2)
    outP = mod(inP1,inP2);
    outP(outP==0)=inP2;
end

