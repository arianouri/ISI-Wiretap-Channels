function [ Capacity ] = C_iid_R( o_PDF , VarEyAns )
% This Function Compute The Capacity Of Dicode Partial-Response Channel
% With Identical-Uniform-Distribution Inputs
Element_Entropy_Y = log(o_PDF);
Element_Entropy_Y(1) = 0 ;
%Element_Entropy_Y(size(o_PDF,1)) = -52 ;
%Element_Entropy_Y(size(o_PDF ,1) - 1) = -52 ;
Entropy_Y = -mean(Element_Entropy_Y);
%Element_Entropy_Y = -sum(Element_Entropy_Y);
%Entropy_Y = (1/size(o_PDF,1)+1).*(Element_Entropy_Y) ;
Entropy_Z = (0.5).*log(2.*pi.*exp(1).*VarEyAns);
Capacity = Entropy_Y - Entropy_Z ;
end

