function [out]=sim_row(inpMAT,inpVEC)
    log_cmpr=(inpMAT==inpVEC);
    out=zeros(size(inpMAT,1),1);
    for index=1:size(inpMAT,1)
        andCUM=1;
        for col_and=1:size(inpMAT,2)
            andCUM=andCUM&log_cmpr(index,col_and);
        end
        out(index)=andCUM;
    end
end