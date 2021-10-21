function [ LogDomain , ProbabilityDomain , s_map ] = BCJR( TRANS , rx_noisless , PR_Trellis , rx_waveform , Trellis_Index_io ,sigma)
 State_Cardinality=size(PR_Trellis,1);
 numBranches_par=size(PR_Trellis,3);
 Number_of_inputs=size(PR_Trellis,4);
 numBranches_tot=State_Cardinality*2^(Number_of_inputs);
 rx_waveform_pac=[reshape(rx_waveform,Number_of_inputs,length(rx_waveform)/Number_of_inputs)]';
%% Initialization
% INF=1e6;
 % Size
gamma=zeros(numBranches_tot,size(rx_waveform_pac,1));
s_map=zeros(numBranches_tot,3);
alpha=zeros(size(TRANS,1),size(rx_waveform_pac,1));
beta=zeros(size(alpha));
% Gtild=zeros(size(gamma));
% Atild=zeros(size(alpha));
% Btild=zeros(size(alpha));
%% gamma
ScaleFactor=0;
for state_i=1:State_Cardinality
    for state_j=[Trellis_Index_io(Trellis_Index_io(:,1)==state_i,2)]'
        for branchIndex=1:numBranches_par
                ScaleFactor=ScaleFactor+1;
                s_map(ScaleFactor,:)=[state_i,state_j,branchIndex];
                if TRANS(state_i,state_j,branchIndex)==0
                    gamma(ScaleFactor,:)=0;
                else
    %% gamma Calculation
                    for X=1:size(rx_noisless,1)
                        if [permute(PR_Trellis(state_i,state_j,branchIndex,:),[4 1 2 3])]'==rx_noisless(X,:)
                              ar=[prod((1./sqrt(2*pi*(sigma^2))).*exp(-(((rx_waveform_pac-rx_noisless(X,:)).^2)/(2*(sigma^2)))),2)]';
                              gamma(ScaleFactor,:)=gamma(ScaleFactor,:)+ar*TRANS(state_i,state_j,branchIndex);
                        end
                    end
                end
        end
    end
end
%     gamma(gamma>1)=1;
    Gtild=log(gamma);
%     Gtild(abs(Gtild)>abs(INF))=-INF;
   
% end
%% alpha
 % Initial Value
% alpha(:,1)=1./size(alpha,1);
% Atild=log(alpha);
alpha(1,1)=1;
% Atild(:,:)=-INF;
% Atild(1,1)=0;
 % alpha and Atild Calculation
for t=2:size(alpha,2)
    for state_i=1:State_Cardinality
        for state_j=[Trellis_Index_io(Trellis_Index_io(:,1)==state_i,2)]'
            alpha(state_j,t)=alpha(state_j,t)+alpha(state_i,t-1).*sum(gamma((s_map(:,1)==state_i)&(s_map(:,2)==state_j),t));
            
%             Gtild_aux=log(sum(gamma((s_map(:,1)==state_i)&(s_map(:,2)==state_j),t)));
            
%             Atild(state_j,t)=Max_Star(Atild(state_j,t),Atild(state_i,t-1)+log(sum(gamma((s_map(:,1)==state_i)&(s_map(:,2)==state_j),t))) );
        end
    end
% Atild(:,t)=Atild(:,t)- sum(Atild(:,t));
alpha(:,t)=alpha(:,t)/sum(alpha(:,t));
end
Atild=log(alpha);
% alpha=exp(Atild);
%% beta
 % Initial Value
beta(:,end)=1./size(beta,1);
% Btild=log(beta);
% beta(1,end)=1;
% Btild(:,:)=-INF;
% Btild(1,end)=0;
 % bata and Btild Calculation
for t=size(beta,2)-1:-1:1
    for state_i=1:State_Cardinality
        for state_j=[Trellis_Index_io(Trellis_Index_io(:,1)==state_i,2)]'
            beta(state_j,t)=beta(state_j,t)+beta(state_i,t+1).*sum(gamma((s_map(:,1)==state_i)&(s_map(:,2)==state_j),t+1));
            
%             Gtild_aux=log(sum(gamma((s_map(:,1)==state_i)&(s_map(:,2)==state_j),t+1)));

%             Btild(state_j,t)=Max_Star(Btild(state_j,t),Btild(state_i,t+1)+log(sum(gamma((s_map(:,1)==state_i)&(s_map(:,2)==state_j),t+1))) );
        end
    end
beta(:,t)=beta(:,t)/sum(beta(:,t));
% Btild(:,t)=Btild(:,t)- sum(Btild(:,t));
end
Btild=log(beta);
% beta=exp(Btild);
%% Lambda tild
  Ltild=Atild+Btild;
%% Sigma tild
Stild=zeros(size(Gtild,1),size(Gtild,2));
for state_i=1:State_Cardinality
    for state_j=[Trellis_Index_io(Trellis_Index_io(:,1)==state_i,2)]'
        Stild((s_map(:,1)==state_i)&(s_map(:,2)==state_j),2:end)...
              =Atild(state_i,1:end-1)...
              +Gtild((s_map(:,1)==state_i)&(s_map(:,2)==state_j),2:end)...
              +Btild(state_j,2:end);
    end
end
%% APP Pr{Sr,Ss|Y}
LAPP_sigma=Stild-Ltild(1,end);
APP_sigma=exp(LAPP_sigma);
for index=1:size(APP_sigma,2)
%     if sum(APP_sigma(:,index))>1
        APP_sigma(:,index)=APP_sigma(:,index)./sum(APP_sigma(:,index));
%     end
end
LAPP_sigma=log(APP_sigma);
%% APP Pr{Ss|Y}
LAPP_lambda=Ltild-Ltild(1,end);
APP_lambda=exp(LAPP_lambda);
for index=1:size(APP_lambda,2)
%     if sum(APP_lambda(:,index))>1
        APP_lambda(:,index)=APP_lambda(:,index)./sum(APP_lambda(:,index));
%     end
end


LAPP_lambda=log(APP_lambda);
%% Output Structure
 % Probability Domain
ProbabilityDomain.alpha=alpha;
ProbabilityDomain.beta=beta;
ProbabilityDomain.gamma=gamma;
ProbabilityDomain.APPsigma=APP_sigma;
ProbabilityDomain.APPlambda=APP_lambda;
  % LOG Domain
  LogDomain.Atild=Atild;
  LogDomain.Btild=Btild;
  LogDomain.Gtild=Gtild;
  LogDomain.Lambda=Ltild;
  LogDomain.Sigma=Stild;
  LogDomain.APPsigma=LAPP_sigma;
  LogDomain.APPlambda=LAPP_lambda;
end