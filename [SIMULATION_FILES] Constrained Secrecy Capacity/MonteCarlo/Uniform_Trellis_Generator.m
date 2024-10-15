function [ jPR_Trellis ] = Uniform_Trellis_Generator( Cardin )
jPR_Trellis=zeros(Cardin);
 for iIndex=1:Cardin
  for jIndex=2*iIndex-1:2*iIndex
   ccIndex=mod(jIndex,Cardin);
    ccIndex(ccIndex==0)=Cardin;
   jPR_Trellis(iIndex,ccIndex)=0.5;
  end
 end
end

