function [ ProbabilityDomain ] = BCJR_fast_cmplx( Trellis_Index_io , sigma , TRANS , Trans_EX , rx_waveform , Window_Size)
 State_Cardinality=size(Trans_EX,1);
 numBranches_par=size(Trans_EX,3);
 Number_of_inputs=size(Trans_EX,4);
 numBranches_tot=State_Cardinality*2^(Number_of_inputs);
 rx_waveform_pac=[reshape(rx_waveform,Number_of_inputs,length(rx_waveform)/Number_of_inputs)].';

%% 



gamma=zeros(numBranches_tot,size(rx_waveform_pac,1));
s_map=zeros(numBranches_tot,3);
%% gamma
ScaleFactor=0;

% for state_i=1:State_Cardinality
%     for state_j=SpecialiZe(SupCha_CODEBOOK(SupCha_CODEBOOK(:,1)==state_i,2))
%         aux_find_bIndex=(SupCha_CODEBOOK(:,[1,2])==[state_i,state_j]);
%         find_bIndex=SupCha_CODEBOOK((aux_find_bIndex(:,1) & aux_find_bIndex(:,2))==1,3);
%         for branchIndex=find_bIndex.'
            
for state_i=1:State_Cardinality
    for state_j=(Trellis_Index_io(Trellis_Index_io(:,1)==state_i,2)).'
        for branchIndex=1:numBranches_par
            ScaleFactor=ScaleFactor+1;
            s_map(ScaleFactor,:)=[state_i,state_j,branchIndex];
%             for X=1:size(rx_noisless,1)
%                 if [permute(Trans_EX(state_i,state_j,branchIndex,:),[4 1 2 3])].'==rx_noisless(X,:)
                      ar=( ...
                         prod( ...
                                (1./(2*pi*(sigma^2)))...
                                .*exp( -( ...
                                            ( (abs( rx_waveform_pac - (permute(Trans_EX(state_i,state_j,branchIndex,:),[4 1 2 3])).' )).^2 ) ...
                                            /(2*(sigma^2)) ...
                                     )  ) ...
                            ,2) ...
                         ).';
                      gamma(ScaleFactor,:)=gamma(ScaleFactor,:)+ar*TRANS(state_i,state_j,branchIndex);
%                 end
%             end
        end
    end
end



%%
% alpha=zeros(State_Cardinality,size(rx_waveform_pac,1));
% beta=zeros(State_Cardinality,size(rx_waveform_pac,1));
%% alpha (forward recursion)
% alpha(:,1)=1./State_Cardinality;
% alpha(1,1)=1;
% for t=2:size(alpha,2)
%     for state_i=1:State_Cardinality
%         for state_j=SpecialiZe(Trellis_Index_io(Trellis_Index_io(:,1)==state_i,2))
%             alpha(state_j,t)=alpha(state_j,t)+alpha(state_i,t-1).*sum(gamma((s_map(:,1)==state_i)&(s_map(:,2)==state_j),t));
%         end
%     end
% alpha(:,t)=alpha(:,t)/sum(alpha(:,t));
% end
% Atild=log(alpha);
%% beta (backward recursion)
% beta(:,end)=1./State_Cardinality;
% for t=size(beta,2)-1:-1:1
%     for state_i=1:State_Cardinality
%         for state_j=SpecialiZe(Trellis_Index_io(Trellis_Index_io(:,1)==state_i,2))
%             beta(state_j,t)=beta(state_j,t)+beta(state_i,t+1).*sum(gamma((s_map(:,1)==state_i)&(s_map(:,2)==state_j),t+1));
%         end
%     end
% beta(:,t)=beta(:,t)/sum(beta(:,t));
% end
% Btild=log(beta);

%% fast forward backward
post_state_MAR=zeros(numBranches_tot/(2.^Number_of_inputs),2.^Number_of_inputs);
past_state_MAR=zeros(size(post_state_MAR));
for iIndex=1:State_Cardinality
    post_state_MAR(iIndex,:)=[s_map(s_map(:,1)==iIndex,2)].'-1;
    past_state_MAR(iIndex,:)=[s_map(s_map(:,2)==iIndex,1)].'-1;
end
past_s_map=zeros(1,numBranches_tot);
nw_sc_index=1;
for auxRindex=1:size(past_state_MAR,1)
    for past_State=SpecialiZe(past_state_MAR(auxRindex,:))
        aux_past_index=[past_State+1,auxRindex]==s_map(:,[1,2]);
        find_aux_past_index=find(aux_past_index(:,1) & aux_past_index(:,2));
        
        past_s_map(nw_sc_index:nw_sc_index+length(find_aux_past_index)-1)=find_aux_past_index-1;
        nw_sc_index=nw_sc_index+length(find_aux_past_index);
    end
end
alpha_fst=zeros(State_Cardinality,Window_Size+1);
beta_fst=zeros(State_Cardinality,Window_Size+1);
alpha_fst(:,1)=1./State_Cardinality;
beta_fst(:,end)=1./State_Cardinality;
[~,~,APP_sigma]=bcjr_fw_bw(alpha_fst,beta_fst,post_state_MAR,past_state_MAR,gamma,Window_Size,...
                           int32(past_s_map));
%% APP Pr{i,j|Y}
for index=1:size(APP_sigma,2)
        APP_sigma(:,index)=APP_sigma(:,index)./sum(APP_sigma(:,index));
end
APP_sigma_bounded(:,:)=APP_sigma(:,Window_Size+2:end-Window_Size);
%% APP Pr{i|Y}
APP_lambda=zeros(size(APP_sigma,1)/2,size(APP_sigma,2));
for index=1:size(APP_lambda,1)
       APP_lambda(index,:)=sum(APP_sigma(s_map(:,2)==index,:));
end
APP_lambda_bounded(:,:)=APP_lambda(:,Window_Size+2:end-Window_Size);

%%
ProbabilityDomain.APPsigma=APP_sigma_bounded;
ProbabilityDomain.APPlambda=APP_lambda_bounded;
end