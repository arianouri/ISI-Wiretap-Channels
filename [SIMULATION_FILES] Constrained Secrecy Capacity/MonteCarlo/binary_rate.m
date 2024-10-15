function [ info_rate ] = binary_rate( Input_Size , VarEyAns , PR_Trellis , TRANS , Trellis_Index_io )
format long;
sigma = abs(sqrt(VarEyAns));

%% Channel Simulation
 [~,~,tx_waveform] = JSCgenerate(Input_Size,TRANS,PR_Trellis,Trellis_Index_io);
 % [~,~,tx_waveform_2] =  SupChaGenerateJSCgenerate(Input_Size,TRANS,PR_Trellis,Trellis_Index_io);

 rx_waveform...
        = tx_waveform...
        + normrnd(0,sigma,1,length(tx_waveform))...
        + 1i*normrnd(0,sigma,1,length(tx_waveform));

%%
% [ ~ , o_PDF ] = alfaPDF_Calculation_up( rx_waveform , VarEyAns , PR_Trellis , TRANS , Trellis_Index_io , states(1));
[ ~ , o_PDF ] = BCJR_FW( rx_waveform , VarEyAns , PR_Trellis , TRANS , Trellis_Index_io);

info_rate = binary_rate_R( o_PDF , VarEyAns , size(PR_Trellis,4) );

end