function [ Capacity_iid ] = C_iid( Input_Size , VarEyAns , Trans , Trans_Pr )
format long;
sigma = abs(sqrt(VarEyAns));
% ISI_Length=size(FrequencyResponse,2)-1;
% [ outputVec , S_inputVec ] = PR_Dicode( Input_Size , Input_Pr , VarEyAns );
%% Channel Simulation
 [~,OutStates]=hmmgenerate(Input_Size,Trans_Pr,eye(size(Trans_Pr,1)));
OutStates(end)=1;
 tx_waveform=MarkovOutput(OutStates,1,Trans);
rx_waveform=tx_waveform+normrnd(0,sigma,1,length(tx_waveform));
% OutSymb=-2*(rand(1,Input_Size)<0.5)+1; %Generate I.U.D. BPSK sequence
% [outputVec,~] = PR_Channel(OutSymb,FrequencyResponse,sigma);
%%
%alfa = alfa_Calculation( outputVec , Input_Size,SNRdb);
%o_PDF = outputPDF( outputVec , alfa,SNRdb );
[ ~ , o_PDF ] = alfaPDF_Calculation( rx_waveform , VarEyAns , Trans , Trans_Pr );
Capacity_iid = C_iid_R( o_PDF , VarEyAns );
end