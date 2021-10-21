function [ Capacity ] = bit_rate_R( o_PDF , VarEyAns , n)
Element_Entropy_Y = log2(o_PDF);
Element_Entropy_Y(1) = 0 ;
Entropy_Y = -mean(Element_Entropy_Y);
Entropy_Z = (0.5)*log2(2*pi*exp(1)*VarEyAns);
Capacity = (Entropy_Y/n - Entropy_Z);
end

