%% This file calculates CSC in BITs/sec (log2),
%  also the TFP of the channel is considered to be complex.

clear all;clc;

[FIR_rx1,FIR_rx2] = phase_array_FIR(pi/4,pi/6);

%% 
N=4; %number of equidistributed initial points for optimization
Input_Size = 1e5;
optIteration=200;

% MONTE_simLength=5e4;
% MONTE_simAvg=10;

kappa=1;
kappaPrime=6;

wind_length=64;

%%
MAIN.State_Cardinality=4;
WRTP.State_Cardinality=2;
State_Cardinality=max(MAIN.State_Cardinality,WRTP.State_Cardinality);
Number_of_inputs=4;

%% Signal to Noise Ratio
MAIN.SNRdB=-4;
WRTP.SNRdB=-4;

%% Channel Frequency Response
% Main Channel
% MAIN.FrequencyResponse=[0.792,0.610];
% MAIN.FrequencyResponse=[1 1 -1 -1];
% MAIN.FrequencyResponse=MAIN.FrequencyResponse./sqrt(sum(MAIN.FrequencyResponse.^2));
MAIN.FrequencyResponse = FIR_rx1./sqrt(sum( (abs(FIR_rx1)).^2 ));

% Eavesdropper's channel
% WRTP.FrequencyResponse=[0.445516026180429,0.633021994668546,0.633086585454355];
% WRTP.FrequencyResponse=[1 -1];
% WRTP.FrequencyResponse=WRTP.FrequencyResponse./sqrt(sum(WRTP.FrequencyResponse.^2));
WRTP.FrequencyResponse = FIR_rx2./sqrt(sum( (abs(FIR_rx2)).^2 ));

%% Joint Source/Wiretap Channel models trellis (EXTENDED)

[ MAIN.Trans_Minimal    , MAIN.Trans_Extended      , MAIN.Trellis_Minimal ,...
  MAIN.Trellis_Extended , MAIN.Trellis_Extended_io , MAIN.branchIndex_Vec ]...
                                                                           ...
     =PR_Trellis_EX(MAIN.FrequencyResponse,Number_of_inputs,MAIN.State_Cardinality);
 

[ WRTP.Trans_Minimal    , WRTP.Trans_Extended      , WRTP.Trellis_Minimal ,...
  WRTP.Trellis_Extended , WRTP.Trellis_Extended_io , WRTP.branchIndex_Vec ]...
                                                                           ...
     =PR_Trellis_EX(WRTP.FrequencyResponse,Number_of_inputs,WRTP.State_Cardinality);
 
 
[ MAIN.Trans_Extended         , WRTP.Trans_Extended,                       ...
  MAIN.Trellis_Extended       , WRTP.Trellis_Extended, TRELLIS,            ...
  MAIN.Trellis_Extended_io    , WRTP.Trellis_Extended_io, TRELLIS_io,      ...
  MAIN.branchIndex_Vec        , WRTP.branchIndex_Vec ]                     ...
                                                                           ...
= ChannelJoinTrellis( MAIN.Trans_Minimal   , WRTP.Trans_Minimal,           ...
                      MAIN.Trellis_Extended , WRTP.Trellis_Extended,       ...
                      MAIN.branchIndex_Vec  , WRTP.branchIndex_Vec );

%% generating equidistributed initial points

[p_i_vec]=equidistributed_generator(N,size(MAIN.Trellis_Extended+WRTP.Trellis_Extended,1));
TRANS_0=zeros([N,State_Cardinality,State_Cardinality,size(MAIN.Trans_Extended+WRTP.Trans_Extended,3)]);
TRANSVEC=zeros([State_Cardinality,State_Cardinality,size(MAIN.Trans_Extended+WRTP.Trans_Extended,3)]);

for MainTR_index=1:N
    for iIndex=1:length(MAIN.branchIndex_Vec + WRTP.branchIndex_Vec)
            TRANSVEC(iIndex)=p_i_vec(MainTR_index,iIndex);
    end
    
    TRANS_0(MainTR_index,:,:,:)=TRANSVEC./sum(TRANSVEC,[2,3]);
end

%% !!! TEMP !!!
% 
% [p_i_vec]=equidistributed_generator(N,State_Cardinality^2);
% TRANS_0=zeros(N,State_Cardinality,State_Cardinality);
% for gIndex=1:size(p_i_vec,1)
%     sCnt=0;
%     for tiIndex=1:State_Cardinality
%         for tjIndex=1:State_Cardinality
%             sCnt=sCnt+1;
%             TRANS_0(gIndex,tiIndex,tjIndex)=p_i_vec(gIndex,sCnt);
%         end
%     end
%     TRANS_0(gIndex,:,:)=TRANS_0(gIndex,:,:)./sum(TRANS_0(gIndex,:,:),3);
% end

%% Noiseless Vector

MAIN_Trans_Extended=MAIN.Trans_Extended;
WRTP_Trans_Extended=WRTP.Trans_Extended;

vecIndex=1;
MAIN_rx_noisless_aux=zeros(size(MAIN_Trans_Extended,1)*size(MAIN_Trans_Extended,2)*size(MAIN_Trans_Extended,3),Number_of_inputs);
WRTP_rx_noisless_aux=zeros(size(WRTP_Trans_Extended,1)*size(WRTP_Trans_Extended,2)*size(WRTP_Trans_Extended,3),Number_of_inputs);
for iIndex=1:size(MAIN_Trans_Extended,1)
    for jIndex=1:size(MAIN_Trans_Extended,2)
        for bIndex=1:size(MAIN_Trans_Extended,3)
            MAIN_rx_noisless_aux(vecIndex,:)=[permute(MAIN_Trans_Extended(iIndex,jIndex,bIndex,:),[4 1 2 3])].';
            WRTP_rx_noisless_aux(vecIndex,:)=[permute(WRTP_Trans_Extended(iIndex,jIndex,bIndex,:),[4 1 2 3])].';
            vecIndex=vecIndex+1;
        end
    end
end
clear vecIndex;

MAIN.rx_noisless=unique(MAIN_rx_noisless_aux,'rows');
WRTP.rx_noisless=unique(WRTP_rx_noisless_aux,'rows');
%%


for transIndex=1:size(TRANS_0,1)
    [transIndex]
    if length(size(TRANS_0))==3
        Trans_0(:,:)=TRANS_0(transIndex,:,:);
    elseif length(size(TRANS_0))==4
        Trans_0(:,:,:)=TRANS_0(transIndex,:,:,:);
    end

    %% Secrecy Capacity
      disp(['Main Channel SNRdB: ',num2str(MAIN.SNRdB)]);
      disp(['Wiretappers Channel SNRdB: ',num2str(WRTP.SNRdB)]);
              
    % [TRANSTAR_vec,c1]=Expectation_Maximization_MET_i_real( Trans_0 , ... 
     [TRANSTAR_vec,c1]=Expectation_Maximization_MET_i_complex( Trans_0 , ...   
                          MAIN , WRTP , ...
                          Input_Size , kappa , kappaPrime , optIteration , ...
                          TRELLIS_io , wind_length);
C=zeros(length(c1),2);  
C(1:end,1)=c1.';
C(1,1)=NaN;

% TRANS=TRANSTAR;
% PERSENT_AGE=0.92;
    
    SAVEINDEX=transIndex;
    save(num2str(SAVEINDEX));
%     save('MonteCarlo\cube.mat');
%         run('MonteCarlo\Main_Of_Mains.m');
end