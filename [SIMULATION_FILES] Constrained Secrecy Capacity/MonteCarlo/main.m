clear all;close all;clc;
load 4.mat
% load trns_pr_ht.mat
collector_samples=0;
Input_Size = 2e5 ;  
MeanIter = 10;
% Trans_Pr=zeros(size(Trans_Extended,1),size(Trans_Extended,2),size(Trans_Extended,3));

for trnsDex = 197:size(TRANSTAR_vec,1)-1
    TRANS(:,:,:) = TRANSTAR_vec(trnsDex+1,:,:,:);

    %% MAin Channel Parameters

    MAIN.VarEyAns=10^(-(MAIN.SNRdB)/10);
    WRTP.VarEyAns=10^(-(WRTP.SNRdB)/10);

    %%
    format long
    C_acc=zeros(1,MeanIter);
    C_B=zeros(1,MeanIter);
    C_E=zeros(1,MeanIter);
        for mIndex=1:MeanIter
                
            C_B(mIndex)=binary_rate( Input_Size , MAIN.VarEyAns , MAIN.Trans_Extended , TRANS , TRELLIS_io );
            C_E(mIndex)=binary_rate( Input_Size , WRTP.VarEyAns , WRTP.Trans_Extended , TRANS , TRELLIS_io );
            
            C_acc(mIndex)=C_B(mIndex)-C_E(mIndex);
            
            disp(['Capacity Sim','(',num2str(mIndex),') :',num2str(C_acc(mIndex))]);
            symbline(length(['Capacity Sim','(',mIndex,') :',num2str(C_acc(mIndex))]),'-');
        
        end

   C(trnsDex,2) = mean(C_acc);
   
   disp(['Eb/N0 BOB: ',num2str(MAIN.SNRdB),'Eb/N0 EVE: ',num2str(WRTP.SNRdB),' | Capacity Avg: ',num2str(mean(C_acc))]);
   symbline(length(['Eb/N0 BOB: ',num2str(MAIN.SNRdB),'Eb/N0 EVE: ',num2str(WRTP.SNRdB),' | Capacity Avg: ',num2str(mean(C_acc))]),'=');

end

C
