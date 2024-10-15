function [states, branches_indexed, emissions]= JSCgenerate(L,trCube,emCube,Trellis_Index_io)

numBranches = size(trCube,3);
numStates = size(trCube,1);
numEmission = size(emCube,4);

branches_indexed = zeros(1,L);
states = zeros(1,L);
emissions = zeros(1,numEmission*L);


statechange = rand(1,L);

currentstate = 1;

em_branch = zeros(numStates,length(SpecialiZe(Trellis_Index_io(Trellis_Index_io(:,1)==1,2)))*numBranches);
for iStateIndex=1:numStates
    auxBranchIndex = 1;
    for jStateIndex = [Trellis_Index_io(Trellis_Index_io(:,1)==iStateIndex,2)].'
        for branchIndex=1:numBranches
            em_branch(iStateIndex,auxBranchIndex) = trCube(iStateIndex,jStateIndex,branchIndex);
            auxBranchIndex=auxBranchIndex + 1;
        end
    end
end
clear branchIndex;

CUM_em_branch = cumsum(em_branch,2);
emIndex = 1;
for count = 1:L
    post_states = Trellis_Index_io(Trellis_Index_io(:,1)==currentstate,2);
    state = post_states(1);
    branch = 1;
    stateVal = statechange(count);
        for branchIndex = length(SpecialiZe(Trellis_Index_io(Trellis_Index_io(:,1)==1,2)))*numBranches - 1 : -1 : 1
            if stateVal > CUM_em_branch(currentstate,branchIndex)
                branch = branchIndex + 1;
                state  = post_states(ceil(branch/numBranches));
                break;
            end
        end
%     branches(count) = branch;
    branches_indexed(count) = modif_mod(branch,numBranches);
    emissions(emIndex:emIndex+numEmission-1)=(permute(emCube(currentstate,state,branches_indexed(count),:),[4 3 2 1])).';
    states(count) = state;
    currentstate = state;
    emIndex = emIndex + numEmission;
end