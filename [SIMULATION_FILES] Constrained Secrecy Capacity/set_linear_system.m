function [COEF,SOLVEC] = set_linear_system(gamma,rho,aA,TrellisIndex,kappa,M,SSPr,Qij)
gamma_NUM=repmat(gamma' ,M,1);
gamma_DNUM=repmat(gamma,1,M);
pHat=(gamma_NUM./(gamma_DNUM+eps)).*(aA./(rho+eps));%\hat{p}_{ij}^*
%% old version
% % % A=eye(size(V,1)*2);
% % % auxC=eye(size(V,1)*2);
% % % C=zeros(size(V,1),size(V,1)*2);
% % % for i=0:size(V,1)-1
% % %    C(i+1,:)= auxC(2*i+1,:)+auxC(2*i+2,:);
% % % end
% % % D=-eye(size(V));
% % % B=zeros(size(V,1)/2,size(V,1));
% % % % auxB=repmat(omega,1,length(omega)).*V;
% % % for j=1:size(V,1)
% % %     iN=2*(j-1)+1:2*(j-1)+1+1;
% % %     iM=mod(2*j-1:2*j,length(SSPr));
% % %     iM(iM==0)=length(SSPr);
% % %     B(iN,j)=V(j,iM);
% % % end
% % % COEF=[A,-B;C,D];
% % % 
% % % COEF=[COEF;...
% % %     [ones(1,2*length(SSPr)),zeros(1,length(SSPr))];...
% % %     [zeros(1,2*length(SSPr)),ones(1,length(SSPr))]];
% % % 
% % % 
% % % SOLVEC=zeros((length(SSPr))*2+length(SSPr),1);
% % % 
% % % SOLVEC=[SOLVEC;1;1];
% % % 
% % % coefA=ones(size(aA));
% % % for i=1:length(SSPr)
% % %     for j=1:length(SSPr)
% % %         coefA(i,j)=((1-kappa)/(kappa))*(V(i,j)*SSPr(i)-Qij(i,j));
% % %     end
% % % end
% % % vec_coefA=ones(length(SSPr)*2,1);
% % % for yInd=1:size(TrellisIndex,1)
% % % %     for xInd=1:size(TrellisIndex,2)
% % %         INDEX=TrellisIndex(yInd,:);
% % %         vec_coefA(yInd)=coefA(INDEX(1),INDEX(2));
% % % %     end
% % % end
% % % % coefA=coefA';
% % % % auxcoefA=coefA(:);
% % % SOLVEC(1:(length(SSPr))*2,1)=vec_coefA;
%% Q_{ij} based version
COEF=zeros(2*size(Qij));
for yInd=1:size(TrellisIndex,1)
        INDEX=TrellisIndex(yInd,:);
        COEF(yInd,yInd)=1-pHat(INDEX(1),INDEX(2));
        if rem(yInd,2)==1
            COEF(yInd,yInd+1)=-pHat(INDEX(1),INDEX(2));
        elseif rem(yInd,2)==0
            COEF(yInd,yInd-1)=-pHat(INDEX(1),INDEX(2));
        end
end
COEF=[COEF;ones(1,size(COEF,2))];

auxSStage=eye(size(Qij,1)*2);
SStage=zeros(size(Qij,1),size(pHat,1)*2);
for i=0:size(Qij,1)-1
   SStage(i+1,:)= auxSStage(2*i+1,:)+auxSStage(2*i+2,:);
end

COEF=[COEF;SStage+[-eye(length(SSPr)),-eye(length(SSPr))]];

SOLVEC=zeros((length(SSPr))*2,1);
coefA=ones(size(aA));
for i=1:length(SSPr)
    for j=1:length(SSPr)
        coefA(i,j)=((1-kappa)/(kappa))*(pHat(i,j)*SSPr(i)-Qij(i,j));
    end
end

for yInd=1:size(TrellisIndex,1)
%     for xInd=1:size(TrellisIndex,2)
        INDEX=TrellisIndex(yInd,:);
        SOLVEC(yInd)=coefA(INDEX(1),INDEX(2));
%     end
end
SOLVEC=[SOLVEC;1;zeros(length(SSPr),1)];
end

