function [ alpha , norma ] = BCJR_FW( rx_waveform ,VarEyAns , Trans , Trans_Pr , Trellis_Index_io)

 State_Cardinality=size(Trans,1);
 numBranches_par=size(Trans,3);
 Number_of_inputs=size(Trans,4);
 numBranches_tot=size(Trellis_Index_io,1)*numBranches_par;
 rx_waveform_pac=[reshape(rx_waveform,Number_of_inputs,length(rx_waveform)/Number_of_inputs)].';

%% 

gamma=zeros(numBranches_tot,size(rx_waveform_pac,1));
s_map=zeros(numBranches_tot,3);
%% gamma
ScaleFactor=0;
            
for state_i=1:State_Cardinality
    for state_j=(Trellis_Index_io(Trellis_Index_io(:,1)==state_i,2)).'
        for branchIndex=1:numBranches_par
            ScaleFactor=ScaleFactor+1;
            s_map(ScaleFactor,:)=[state_i,state_j,branchIndex];
%             for X=1:size(rx_noisless,1)
%                 if [permute(Trans_EX(state_i,state_j,branchIndex,:),[4 1 2 3])].'==rx_noisless(X,:)
                      ar=( ...
                         prod( ...
                                (1./(2*pi*VarEyAns))...
                                .*exp( -( ...
                                            ( (abs( rx_waveform_pac - (permute(Trans(state_i,state_j,branchIndex,:),[4 1 2 3])).' )).^2 ) ...
                                            /(2*VarEyAns) ...
                                     )  ) ...
                            ,2) ...
                         ).';
                      gamma(ScaleFactor,:)=gamma(ScaleFactor,:)+ar*Trans_Pr(state_i,state_j,branchIndex);
%                 end
%             end
        end
    end
end


%% alpha
 % Initial Value
alpha=zeros(State_Cardinality,size(rx_waveform_pac,1));
% alpha(:,1)=1./size(alpha,1);
alpha(:,1)=eps;
alpha(1,1)=1-eps*(size(alpha,1)-1);

 % alpha Calculation
 norma = zeros(1,size(alpha,2));
for t=2:size(alpha,2)
    for state_i=1:State_Cardinality
        for state_j=[Trellis_Index_io(Trellis_Index_io(:,1)==state_i,2)].'
            alpha(state_j,t)...
                =alpha(state_j,t)...
                 +alpha(state_i,t-1)...
                 .*sum( gamma((s_map(:,1)==state_i)&(s_map(:,2)==state_j),t) );
        end
    end
norma(t)=1/sum(alpha(:,t)); % Normalization factor for alpha
alpha(:,t)=alpha(:,t).*norma(t);
end

end