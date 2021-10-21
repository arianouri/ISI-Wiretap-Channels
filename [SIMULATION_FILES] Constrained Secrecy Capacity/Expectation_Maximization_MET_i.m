function [TRANSTAR_vec,c]=Expectation_Maximization_MET_i( TRANS , ...
                                               MAIN , WRTP , ...
                                               K , kappa , kappaPrime , stpCrtrn , ...
                                               TRELLIS_io , wind )
M=size(TRANS,1);

%% Optimization Procedure
Flag=1;
cIndex=1;
TRANSTAR_vec(1,:,:,:)=TRANS;
[SSPr] = SState_Pr(TRANS);
% Qij=repmat(SSPr,1,length(SSPr)).*TRANS;
    while Flag
        
        PASTQij=repmat(SSPr,1,size(TRANS,2),size(TRANS,3)).*TRANS;
        PASTSSPr=SSPr;
        
        TRANSTAR_vec(cIndex+1,:,:,:)=TRANS;

        [~,~,MAIN.tx_waveform,WRTP.tx_waveform] = JSWCgenerate(K+2*wind+1,TRANS,MAIN.Trans_Extended,WRTP.Trans_Extended,TRELLIS_io);
        %K+2*wind+1 is for setting out the boundaries 
        
        MAIN.sigma=10.^(-MAIN.SNRdB/20);  
        MAIN.rx_waveform=MAIN.tx_waveform+normrnd(0,MAIN.sigma,1,length(MAIN.tx_waveform));
        WRTP.sigma=10.^(-WRTP.SNRdB/20);  
        WRTP.rx_waveform=WRTP.tx_waveform+normrnd(0,WRTP.sigma,1,length(WRTP.tx_waveform));
        
        %% BCJR

%         [~,MAIN.ProbabilityDomain,~]...
%         = BCJR(TRANS,MAIN.rx_noisless,...
%                MAIN.Trans_Extended,MAIN.rx_waveform,TRELLIS_io,MAIN.sigma);
% 
%         [~,WRTP.ProbabilityDomain,~]...
%         = BCJR(TRANS,WRTP.rx_noisless,...
%                WRTP.Trans_Extended,WRTP.rx_waveform,TRELLIS_io,WRTP.sigma);
           
           %% boosted windowed BCJR
      
        [ MAIN.ProbabilityDomain ]...
        = BCJR_fast( TRELLIS_io , MAIN.sigma , TRANS , MAIN.Trans_Extended , MAIN.rx_waveform , wind);

        [ WRTP.ProbabilityDomain ]...
        = BCJR_fast( TRELLIS_io , WRTP.sigma , TRANS , WRTP.Trans_Extended , WRTP.rx_waveform , wind);


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
                            =(1/K)*sum(log2(...
                                (MAIN.ProbabilityDomain.APPsigma(ScaleFactor,:).^(MAIN.ProbabilityDomain.APPsigma(ScaleFactor,:)...
                                    ./(eps+(SSPr(state_i)*TRANS(state_i,state_j,branchIndex)))))...
                                         ./((MAIN.ProbabilityDomain.APPlambda(state_i,:)).^(MAIN.ProbabilityDomain.APPlambda(state_i,:)...
                                             ./(eps+SSPr(state_i))))+eps));
                %% T2
                        T_E(state_i,state_j,branchIndex)...
                            =(1/K)*sum(log2(...
                                (WRTP.ProbabilityDomain.APPsigma(ScaleFactor,:).^(WRTP.ProbabilityDomain.APPsigma(ScaleFactor,:)...
                                    ./(eps+(SSPr(state_i)*TRANS(state_i,state_j,branchIndex)))))...
                                         ./((WRTP.ProbabilityDomain.APPlambda(state_i,:)).^(WRTP.ProbabilityDomain.APPlambda(state_i,:)...
                                             ./(eps+SSPr(state_i))))+eps));
                end
            end
        end
        q=kappa*kappaPrime;
        T=T_B-T_E;

        unit_mat=zeros(size(TRANS));
        unit_mat(TRELLIS_io(:,1),TRELLIS_io(:,2),:)=1;
        A_Cube=TRANS.*(2.^(T/q));
        A_Cube(unit_mat==0)=0;
        A=sum(A_Cube,3);
        if abs(sum(sum(A)))==inf || isnan(sum(sum(A)))
            cIndex=cIndex+1;
            c(cIndex)=c(cIndex-1);
            disp(['Iteration*: ',num2str(cIndex),' | Capacity: ',num2str(c(cIndex))]);
            symbline(length(['Iteration*: ',num2str(cIndex),' | Capacity: ',num2str(c(cIndex))]),'=');
        %     if (abs(mean((sum((reshape(nu.*omega.^2,size(Qij,1),size(Qij,2)).*EqualiZer)')))-eta)>1e-4)
                Flag=0;
        %     end
            continue;
        end
        [EigenVector, EigenValue] = eig(A,'nobalance');
    %     EigenVector=(EigenVector);
    %     EigenValue=(EigenValue);
        rho=max(max(EigenValue));
        [~,x]=find(EigenValue==rho);clear y;
        gamma=EigenVector(:,x);
    %     [COEF,SOLVEC] = set_linear_system(gamma,rho,A,TRELLIS_io,kappa,M,SSPr,Qij);
    %     [vX,~,vRREF]=RREF([COEF,SOLVEC]);
    %     if (min(vX(1:2*length(SSPr)))<0)
    %         Flag=0;
    %         continue;
    %     end
    %     stpCnt=1;
    %     PASTQij=Qij;
    %     for iIndex=1:2*size(Qij,1)
    %             Qij(TrellisIndex(iIndex,1),TrellisIndex(iIndex,2))=vX(stpCnt);
    %             stpCnt=stpCnt+1;
    %     end

        gamma_NUM=repmat(gamma' ,M,1);
        gamma_DNUM=repmat(gamma,1,M);

        TRANS=(gamma_NUM./(gamma_DNUM+eps)).*((TRANS.*(2.^(T/q)))./(rho+eps));%\hat{p}_{ij}^*
        
        [SSPr] = SState_Pr(TRANS);
        Qij=repmat(SSPr,1,size(TRANS,2),size(TRANS,3)).*TRANS;

    %     [SSPr] = sum(Qij,2);
    %     TRANS=Qij./repmat(SSPr,1,length(SSPr));

        C=sum(sum(sum((T_B-T_E).*Qij)))-cvX_func(Qij,SSPr,PASTQij,PASTSSPr,kappa,kappaPrime);

        cIndex=cIndex+1;
        c(cIndex)=C;
        disp(['Iteration: ',num2str(cIndex),' | Capacity: ',num2str(c(cIndex))]);
        symbline(length(['Iteration: ',num2str(cIndex),' | Capacity: ',num2str(c(cIndex))]),'=');
        if (cIndex>stpCrtrn)
            Flag=0;
        end
    end
end