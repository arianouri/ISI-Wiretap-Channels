function [ info_rate ] = bit_rate( Input_Size , VarEyAns , PR_Trellis , TRANS , Trellis_Index_io )
format long;
sigma = abs(sqrt(VarEyAns));

%% Channel Simulation
 [~,~,tx_waveform] = JSCgenerate(Input_Size,TRANS,PR_Trellis,Trellis_Index_io);
 rx_waveform=tx_waveform+normrnd(0,sigma,1,length(tx_waveform));

%%
[ ~ , o_PDF ] = alfaPDF_Calculation( rx_waveform , VarEyAns , PR_Trellis , TRANS , Trellis_Index_io );
info_rate = bit_rate_R( o_PDF , VarEyAns , size(PR_Trellis,4) );

end