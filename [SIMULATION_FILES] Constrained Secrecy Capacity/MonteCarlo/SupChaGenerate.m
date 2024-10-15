function [states, branches_indexed, emissions]= SupChaGenerate(L,trCube,emCube,Trellis_Index_io)

    numBranches = size(trCube,3);
    numStates = size(trCube,1);
    numEmission = size(emCube,4);

    branches_indexed = zeros(1,L);
    states = zeros(1,L);
    emissions = zeros(1,numEmission*L);


    statechange = rand(1,L);

    currentstate = 1;

    em_branch = zeros(numStates,length(Trellis_Index_io(Trellis_Index_io(:,1)==1,2)));
    for iStateIndex=1:numStates
        auxBranchIndex = 1;
        for jStateIndex = SpecialiZe([Trellis_Index_io(Trellis_Index_io(:,1)==iStateIndex,2)].')
            for branchIndex=1:numBranches
                if ~isnan(trCube(iStateIndex,jStateIndex,branchIndex))
                    em_branch(iStateIndex,auxBranchIndex) = trCube(iStateIndex,jStateIndex,branchIndex);
                    auxBranchIndex=auxBranchIndex + 1;
                end
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
            for branchIndex = size(CUM_em_branch,2) - 1 : -1 : 1
                if stateVal > CUM_em_branch(currentstate,branchIndex)
                    branch = branchIndex + 1;
                    state  = post_states(branch);
                    break;
                end
            end
%         branches(count) = branch;
        branches_indexed(count) = 1;
        if length(SpecialiZe(post_states))<length(post_states)
            aux_branch_index=find(post_states==state);            
            branches_indexed(count)=find(branch==aux_branch_index);
        end
        emissions(emIndex:emIndex+numEmission-1)=(permute(emCube(currentstate,state,branches_indexed(count),:),[4 3 2 1])).';
        states(count) = state;
        currentstate = state;
        emIndex = emIndex + numEmission;
    end
end