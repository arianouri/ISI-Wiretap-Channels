function [ TruthTable ] = Truth_Table( Nbit )
Number=0:(2^(Nbit)-1);
Number=dec2bin(Number);
TruthTable=zeros(size(Number));
for rIndex=1:size(Number,1)
    for cIndex=1:size(Number,2)
        TruthTable(rIndex,cIndex)=round(str2double(Number(rIndex,cIndex)));
    end
end
end

