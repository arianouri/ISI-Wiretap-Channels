function [ Trans_Minimal , Trans_Extended , Trellis_Minimal , Trellis_Extended , Trellis_Extended_io , branchIndex_Vec ] = PR_Trellis_EX( FrequencyResponse , Number_of_inputs , State_Cardinality )
isiLength=size(FrequencyResponse,2)-1;
State_Cardinality_MT=2^isiLength;
if State_Cardinality<State_Cardinality_MT
   error("minimum number of states is 2^(isiLength)");
end
Nbit=size(FrequencyResponse,2);
tblInput=2*Truth_Table(Nbit)-1;
Trans_Minimal=State_Cardinality;

sIndex=1;
coIndex=1;
Trellis_Minimal=zeros(2*State_Cardinality,2);
for roIndex=1:State_Cardinality
    if coIndex>State_Cardinality
        coIndex=1;
    end
    for x=1:2
%         [roIndex,coIndex]
        Trans_Minimal(roIndex,coIndex)=sum(inVec(FrequencyResponse).*tblInput(modif_mod(sIndex,2*State_Cardinality_MT),:));
        Trellis_Minimal(sIndex,:)=[roIndex,coIndex];
        sIndex=sIndex+1;
        coIndex=coIndex+1;
    end
end
clear coIndex;
clear roIndex;

Trellis_Extended=zeros(State_Cardinality*2^(Number_of_inputs),Number_of_inputs+1);
for Index=1:size(Trellis_Extended,2)
    odevIndex=1;
 for rIndex=1:size(Trellis_Extended,1)/(State_Cardinality*2^(Index-1)):size(Trellis_Extended,1)
% [Index rIndex odevIndex]
    Trellis_Extended(rIndex:min(rIndex+2^(Number_of_inputs+1-Index),size(Trellis_Extended,1)),Index)...
                                            =Trellis_Minimal(modif_mod(odevIndex,size(Trellis_Minimal,1)),end);
odevIndex=odevIndex+1;
 end
end

aux_io_index=1;
aux_Trellis_Extended_io=Trellis_Extended(:,[1 size(Trellis_Extended,2)]);
numEmBranch = length(SpecialiZe(aux_Trellis_Extended_io(aux_Trellis_Extended_io(:,1)==1,2)));
for ioState=1:State_Cardinality
    Trellis_Extended_io(aux_io_index:aux_io_index+numEmBranch-1,:)...
        =[repmat(ioState,numEmBranch,1),(SpecialiZe(aux_Trellis_Extended_io(aux_Trellis_Extended_io(:,1)==ioState,2)))'];
    aux_io_index = aux_io_index + numEmBranch;
end

TrellisIndex_FULL=zeros(State_Cardinality^2,2);
Frindep=1;
for FrIndex=1:State_Cardinality
    for FcIndex=1:State_Cardinality
        TrellisIndex_FULL(Frindep,:)=[FrIndex,FcIndex];
        Frindep=Frindep+1;
    end
end

branchIndex_Vec=zeros(size(Trellis_Extended,1),1);
for minIndex=1:size(TrellisIndex_FULL,1)
    branchIndex=1;
    for extIndex=1:size(Trellis_Extended,1)
        if [Trellis_Extended(extIndex,1) Trellis_Extended(extIndex,end)]==TrellisIndex_FULL(minIndex,:)
          branchIndex_Vec(extIndex)=branchIndex;
          branchIndex=branchIndex+1;
        end
    end
end



Trans_Extended=zeros(State_Cardinality,State_Cardinality,max(branchIndex_Vec),Number_of_inputs);

for mnrIndex=1:size(Trellis_Extended,1)
    for mncIndex=1:Number_of_inputs
       Trans_Extended(Trellis_Extended(mnrIndex,1),...
                      Trellis_Extended(mnrIndex,end),...
                      branchIndex_Vec(mnrIndex),...
                      mncIndex) = Trans_Minimal(Trellis_Extended(mnrIndex,mncIndex),Trellis_Extended(mnrIndex,mncIndex+1));
    end
end