function [TRANSTAR_vec,c] = Expectation_Maximization_MET_ii( TRANS , ...
                                               MAIN , WRTP , ...
                                               K , q , stpCrtrn , ...
                                               TRELLIS_io )
%%
% + DUE TO SOME MODIFICATIONS OF THE ALGORITHM ON THE PAPER, THIS METHOD IS 
%   NOT PUBLISHED. BUT IN SOME CASES IT COULD FIND BETTER RESULTS!
%   THIS ALGORITHM IS NOT NUMERICALLY STABLE IN GENERAL
%%
M=size(TRANS,1);
%% Optimization Procedure
Flag=1;
cIndex=1;
TRANSTAR_vec(1,:,:,:)=TRANS;
[SSPr] = SState_Pr(TRANS); % "Mio": (Column vector that returnes Steady state probabilities of Markov model)
Qij=repmat(SSPr,1,length(SSPr)).*TRANS;
while Flag
    
 TRANSTAR_vec(cIndex+1,:,:,:)=TRANS;
 [V,c,nu,omega]=xlnx_Eval(Qij,q,SSPr);
 
 [~,~,MAIN.tx_waveform,WRTP.tx_waveform] = JSWCgenerate(K,TRANS,MAIN.Trans_Extended,WRTP.Trans_Extended,TRELLIS_io);

[~,MAIN.ProbabilityDomain,~]...
    = BCJR(TRANS,MAIN.rx_noisless,MAIN.SNRdB,...
           MAIN.Trans_Extended,MAIN.tx_waveform,TRELLIS_io);
    
[~,WRTP.ProbabilityDomain,~]...
        = BCJR(TRANS,WRTP.rx_noisless,WRTP.SNRdB,...
               WRTP.Trans_Extended,WRTP.tx_waveform,TRELLIS_io);

%% Tij Calculations
T_B=zeros(size(TRANS));
T_E=zeros(size(TRANS));
ScaleFactor=0;
for state_i=1:M
    for state_j=[TRELLIS_io(TRELLIS_io(:,1)==state_i,2)]'
        for branchIndex=1:size(TRANS,3)
        ScaleFactor=ScaleFactor+1;
        %% T1
                T_B(state_i,state_j,branchIndex)...
                    =(1/K)*sum(log(...
                        (MAIN.ProbabilityDomain.APPsigma(ScaleFactor,:).^(MAIN.ProbabilityDomain.APPsigma(ScaleFactor,:)...
                            ./(eps+(SSPr(state_i)*TRANS(state_i,state_j,branchIndex)))))...
                                 ./((MAIN.ProbabilityDomain.APPlambda(state_i,:)).^(MAIN.ProbabilityDomain.APPlambda(state_i,:)...
                                     ./(eps+SSPr(state_i))))+eps));
        %% T2
                T_E(state_i,state_j,branchIndex)...
                    =(1/K)*sum(log(...
                        (WRTP.ProbabilityDomain.APPsigma(ScaleFactor,:).^(WRTP.ProbabilityDomain.APPsigma(ScaleFactor,:)...
                            ./(eps+(SSPr(state_i)*TRANS(state_i,state_j,branchIndex)))))...
                                 ./((WRTP.ProbabilityDomain.APPlambda(state_i,:)).^(WRTP.ProbabilityDomain.APPlambda(state_i,:)...
                                     ./(eps+SSPr(state_i))))+eps));
        end
    end
        %%   
T_B(T_B<-10e6)=-10e6;
T_B(T_B>+10e6)=+10e6;
T_B(isnan(T_B))=0;
T_E(T_E<-10e6)=-10e6;
T_E(T_E>+10e6)=+10e6;
T_E(isnan(T_E))=0;
end
C2=sum(sum((T_B-T_E).*Qij));
        %%
T=T_B-T_E;
T(T<-10e6)=-10e6;
T(T>+10e6)=+10e6;
A=exp(T/q);
A(Qij==0)=0;
if abs(sum(sum(A)))==inf || isnan(sum(sum(A)))
    cIndex=cIndex+1;
    c1(cIndex)=c1(cIndex-1);
    disp(['Iteration*: ',num2str(cIndex),' | Capacity: ',num2str(c1(cIndex))]);
    symbline(length(['Iteration*: ',num2str(cIndex),' | Capacity: ',num2str(c1(cIndex))]),'=');
%     if (abs(mean((sum((reshape(nu.*omega.^2,size(Qij,1),size(Qij,2)).*EqualiZer)')))-eta)>1e-4)
        Flag=0;
%     end
    continue;
end
[EigenVector, EigenValue] = eig(A,'nobalance');
EigenVector=(EigenVector);
EigenValue=(EigenValue);
rho=max(max(EigenValue));
[~,x]=find(EigenValue==rho);clear y;
gamma=EigenVector(:,x);
% gamma_NUM=repmat(gamma' ,M,1);
% gamma_DNUM=repmat(gamma,1,M);
[COEF,SOLVEC] = set_linear_system(gamma,omega,rho,A,V,TrellisIndex,q);
[vX,~,vRREF]=RREF([COEF,SOLVEC]);
stpCnt=1;
for iIndex=1:2*size(TRANS,1)
%     for jIndex=[1,2]
        TRANS(TrellisIndex(iIndex,1),TrellisIndex(iIndex,2))=vX(stpCnt);
        stpCnt=stpCnt+1;
%     end
end
TRANS(TRANS>1)=1;
TRANS(TRANS<0)=0;
TRANS=TRANS./repmat((sum(TRANS'))',1,size(TRANS,2));
SSPr_POST=SSPr;
[SSPr] = SState_Pr(TRANS); % "Mio": (Column vector that returnes Steady state probabilities of Markov model)
omegaVec=1./(repmat(SSPr,length(SSPr),1));
Qij=repmat(SSPr,1,length(SSPr)).*TRANS;
C1=sum(sum((T_B-T_E).*Qij))-xlnx(nu,c,q,Qij,omegaVec);

cIndex=cIndex+1;

fract_ss_pr_Vec(cIndex,:)=SSPr./SSPr_POST;
R_B(cIndex)=sum(sum((T_B).*Qij));
R_E(cIndex)=sum(sum((T_E).*Qij));
c1(cIndex)=C1;
c2(cIndex)=C2;
disp(['Iteration: ',num2str(cIndex),' | Capacity: ',num2str(c2(cIndex))]);
symbline(length(['Iteration: ',num2str(cIndex),' | Capacity: ',num2str(c2(cIndex))]),'=');
if (cIndex>stpCrtrn)
    Flag=0;
end
end
TRANS_all(cIndex+1,:,:)=TRANS;
end