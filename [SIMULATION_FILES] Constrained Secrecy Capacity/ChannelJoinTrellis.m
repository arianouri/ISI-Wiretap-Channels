function [ Trans_Extended_1   , Trans_Extended_2,                          ...
           Trellis_Extended_1 , Trellis_Extended_2, TRELLIS,               ...
           Trellis_Extended_io_1 , Trellis_Extended_io_2, TRELLIS_io,      ...
           branchIndex_Vec_1  , branchIndex_Vec_2 ]                        ...
                                                                           ...
= ChannelJoinTrellis( Trans_Minimal_1    , Trans_Minimal_2,                ...
                      Trellis_Extended_1 , Trellis_Extended_2,             ...
                      branchIndex_Vec_1  , branchIndex_Vec_2 )
%%
if(size(Trellis_Extended_1,2)~=size(Trellis_Extended_2,2)),error('Er n_INP');end

%%
n_INP=size(Trellis_Extended_1,2)-1;
S_cardinal_1=size(Trellis_Extended_1,1)/2^n_INP;
S_cardinal_2=size(Trellis_Extended_2,1)/2^n_INP;

%% Trellis_Extended

if size(Trellis_Extended_1,1)==max(size(Trellis_Extended_1,1),size(Trellis_Extended_2,1))
    aux_Trellis_Extended_1=Trellis_Extended_1;
    aux_Trellis_Extended_2=Trellis_Extended_1;
        TRELLIS=Trellis_Extended_1;
else
    aux_Trellis_Extended_1=Trellis_Extended_2;
    aux_Trellis_Extended_2=Trellis_Extended_2;
        TRELLIS=Trellis_Extended_2;
end
Trellis_Extended_1=modif_mod(aux_Trellis_Extended_1,S_cardinal_1);
Trellis_Extended_2=modif_mod(aux_Trellis_Extended_2,S_cardinal_2);

%% Trellis_Extended_io

aux_io_index=1;
aux_Trellis_Extended_io_1=Trellis_Extended_1(:,[1 size(Trellis_Extended_1,2)]);
numEmBranch_1 = length(SpecialiZe(aux_Trellis_Extended_io_1(aux_Trellis_Extended_io_1(:,1)==1,2)));
for ioState=1:S_cardinal_1
    Trellis_Extended_io_1(aux_io_index:aux_io_index+numEmBranch_1-1,:)...
        =[repmat(ioState,numEmBranch_1,1),...
         (SpecialiZe(aux_Trellis_Extended_io_1(aux_Trellis_Extended_io_1(:,1)==ioState,2))).'];
    aux_io_index = aux_io_index + numEmBranch_1;
end

aux_io_index=1;
aux_Trellis_Extended_io_2=Trellis_Extended_2(:,[1 size(Trellis_Extended_2,2)]);
numEmBranch_2 = length(SpecialiZe(aux_Trellis_Extended_io_2(aux_Trellis_Extended_io_2(:,1)==1,2)));
for ioState=1:S_cardinal_2
    Trellis_Extended_io_2(aux_io_index:aux_io_index+numEmBranch_2-1,:)...
        =[repmat(ioState,numEmBranch_2,1),...
         (SpecialiZe(aux_Trellis_Extended_io_2(aux_Trellis_Extended_io_2(:,1)==ioState,2))).'];
    aux_io_index = aux_io_index + numEmBranch_2;
end

aux_io_index=1;
aux_TRELLIS_io=TRELLIS(:,[1 size(TRELLIS,2)]);
numEmBranch = length(SpecialiZe(aux_TRELLIS_io(aux_TRELLIS_io(:,1)==1,2)));
for ioState=1:max(S_cardinal_1,S_cardinal_2)
    TRELLIS_io(aux_io_index:aux_io_index+numEmBranch-1,:)...
        =[repmat(ioState,numEmBranch,1),...
         (SpecialiZe(aux_TRELLIS_io(aux_TRELLIS_io(:,1)==ioState,2))).'];
    aux_io_index = aux_io_index + numEmBranch;
end


%% Branch Index

TrellisIndex_FULL=zeros((max(S_cardinal_1,S_cardinal_2))^2,2);
Frindep=1;
for FrIndex=1:max(S_cardinal_1,S_cardinal_2)
    for FcIndex=1:max(S_cardinal_1,S_cardinal_2)
        TrellisIndex_FULL(Frindep,:)=[FrIndex,FcIndex];
        Frindep=Frindep+1;
    end
end

max_branchIndex_1=max(branchIndex_Vec_1);
max_branchIndex_2=max(branchIndex_Vec_2);
branchIndex_Vec_1=zeros(size(Trellis_Extended_1,1),1);
branchIndex_Vec_2=zeros(size(Trellis_Extended_1,1),1);
for minIndex=1:size(TrellisIndex_FULL,1)
    branchIndex_1=1;
    branchIndex_2=1;
    for extIndex=1:size(Trellis_Extended_1,1)
        if [Trellis_Extended_1(extIndex,1) Trellis_Extended_1(extIndex,end)]==TrellisIndex_FULL(minIndex,:)
          branchIndex_Vec_1(extIndex)=branchIndex_1;
          branchIndex_1=branchIndex_1+1;
        end
        if [Trellis_Extended_2(extIndex,1) Trellis_Extended_2(extIndex,end)]==TrellisIndex_FULL(minIndex,:)
          branchIndex_Vec_2(extIndex)=branchIndex_2;
          branchIndex_2=branchIndex_2+1;
        end
    end
end
branchIndex_Vec_1=modif_mod(branchIndex_Vec_1,max_branchIndex_1);
branchIndex_Vec_2=modif_mod(branchIndex_Vec_2,max_branchIndex_2);

%% Trans Extended

if S_cardinal_1>S_cardinal_2
    branch_DIM=max(branchIndex_Vec_1);
    BRANCH_VEC=branchIndex_Vec_1;
else
    branch_DIM=max(branchIndex_Vec_2);
    BRANCH_VEC=branchIndex_Vec_2;
end

Trans_Extended_1=zeros(max(S_cardinal_1,S_cardinal_2),                 ...
                           max(S_cardinal_1,S_cardinal_2),             ...
                           branch_DIM,                                 ...
                           n_INP);
                       
Trans_Extended_2=zeros(max(S_cardinal_1,S_cardinal_2),                 ...
                           max(S_cardinal_1,S_cardinal_2),             ...
                           branch_DIM,                                 ...
                           n_INP);

for mnrIndex=1:size(TRELLIS,1)
    for mncIndex=1:n_INP
       %------------------------------------------------------------------                   
       Trans_Extended_1(TRELLIS(mnrIndex,1),                           ...
                        TRELLIS(mnrIndex,end),                         ...
                        BRANCH_VEC(mnrIndex),                          ...
                        mncIndex)                                      ...
                                                                       ...
            = Trans_Minimal_1(Trellis_Extended_1(mnrIndex,mncIndex),   ...
                              Trellis_Extended_1(mnrIndex,mncIndex+1));
       %------------------------------------------------------------------                   
       Trans_Extended_2(TRELLIS(mnrIndex,1),                           ...
                        TRELLIS(mnrIndex,end),                         ...
                        BRANCH_VEC(mnrIndex),                          ...
                        mncIndex)                                      ...
                                                                       ...
            = Trans_Minimal_2(Trellis_Extended_2(mnrIndex,mncIndex),   ...
                              Trellis_Extended_2(mnrIndex,mncIndex+1));  
       %------------------------------------------------------------------                   
    end
end

end
