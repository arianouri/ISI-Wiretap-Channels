clear all;clc;

%% 
N=100; %number of the initial points
Input_Size = 1e3;
optIteration=80;

MONTE_simLength=1e3;
MONTE_simAvg=10;

kappa=1;
kappaPrime=3;

wind_length=16;

%%
MAIN.State_Cardinality=2;
WRTP.State_Cardinality=4;

State_Cardinality=max(MAIN.State_Cardinality,WRTP.State_Cardinality);
Number_of_inputs=1;

%% Signal to Noise Ratio
MAIN.SNRdB=-6.5;
WRTP.SNRdB=-6.0;

%% Channel Frequency Response
% Main Channel
MAIN.FrequencyResponse=[0.792,0.610];
% MAIN.FrequencyResponse=[1 1 -1 -1];
MAIN.FrequencyResponse=MAIN.FrequencyResponse./sqrt(sum(MAIN.FrequencyResponse.^2));
% Eavesdropper's channel
WRTP.FrequencyResponse=[0.445516026180429,0.633021994668546,0.633086585454355];
% WRTP.FrequencyResponse=[1 -1];
WRTP.FrequencyResponse=WRTP.FrequencyResponse./sqrt(sum(WRTP.FrequencyResponse.^2));

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

[p_i_vec]=equidistributed_generator(N,State_Cardinality);
TRANSVEC=zeros(N,2*State_Cardinality);
for iIndex=1:State_Cardinality
        TRANSVEC(:,2*iIndex-1)=p_i_vec(:,iIndex);
        TRANSVEC(:,2*iIndex)=1-p_i_vec(:,iIndex);
end
stpCnt=1;
for MainTR_index=1:N
    for iIndex=1:2*State_Cardinality
            TRANS_0(MainTR_index,TRELLIS_io(iIndex,1),TRELLIS_io(iIndex,2))=TRANSVEC(MainTR_index,stpCnt);
            stpCnt=stpCnt+1;
    end
    stpCnt=1;
end

%% Noiseless Vector

MAIN_Trans_Extended=MAIN.Trans_Extended;
WRTP_Trans_Extended=WRTP.Trans_Extended;

vecIndex=1;
MAIN_rx_noisless_aux=zeros(size(MAIN_Trans_Extended,1)*size(MAIN_Trans_Extended,2)*size(MAIN_Trans_Extended,3),Number_of_inputs);
WRTP_rx_noisless_aux=zeros(size(WRTP_Trans_Extended,1)*size(WRTP_Trans_Extended,2)*size(WRTP_Trans_Extended,3),Number_of_inputs);
for iIndex=1:size(MAIN_Trans_Extended,1)
    for jIndex=1:size(MAIN_Trans_Extended,2)
        for bIndex=1:size(MAIN_Trans_Extended,3)
            MAIN_rx_noisless_aux(vecIndex,:)=[permute(MAIN_Trans_Extended(iIndex,jIndex,bIndex,:),[4 1 2 3])]';
            WRTP_rx_noisless_aux(vecIndex,:)=[permute(WRTP_Trans_Extended(iIndex,jIndex,bIndex,:),[4 1 2 3])]';
            vecIndex=vecIndex+1;
        end
    end
end
clear vecIndex;

MAIN.rx_noisless=unique(MAIN_rx_noisless_aux,'rows');
WRTP.rx_noisless=unique(WRTP_rx_noisless_aux,'rows');
%%


for transIndex=1:size(TRANS_0,1)
    Trans_0(:,:)=TRANS_0(transIndex,:,:);

    %% Secrecy Capacity
      disp(['Main Channel SNRdB: ',num2str(MAIN.SNRdB)]);
      disp(['Wiretappers Channel SNRdB: ',num2str(WRTP.SNRdB)]);
              
    [TRANSTAR_vec,c1]=Expectation_Maximization_MET_i( Trans_0 , ...
                      MAIN , WRTP , ...
                      Input_Size , kappa , kappaPrime , optIteration , ...
                      TRELLIS_io , wind_length );
C=zeros(optIteration+1,2);  
C(1:end,1)=c1';
mC(1,1)=NaN;

% TRANS=TRANSTAR;
PERSENT_AGE=0.6;
    
    SAVEINDEX=transIndex;
    save('MonteCarlo\cube.mat');
        run('MonteCarlo\Main_Of_Mains.m');
end