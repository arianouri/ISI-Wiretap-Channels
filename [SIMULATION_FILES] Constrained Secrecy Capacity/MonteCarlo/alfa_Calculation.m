function [ alfa ] = alfa_Calculation( outputVec , Input_Size , SNRdb)
SNR = 10.^((SNRdb)/10);
alfa = zeros(Input_Size-1 ,2);
% alfa(State,t) | States -> {1,2};
% first Column -> S=0 ;
% second Column -> S=1 ;
alfa(1,1) = 1;
alfa(1,2) = 0;
% alfa Formula Numerator
% alfa_Cal_numerator = zeros(Input_Size ,2);
% alfa Formula Denumerator
%alfa_Cal_denominator = zeros(Input_Size ,2);
 flag=1;
 ctnTime=2;
 s=1;
    while(flag==1)
        alfa_Cal_numerator  = 0.5.*((alfa(ctnTime-1,1).*outputPDF_R( outputVec(ctnTime) , 1 , s , SNR))+(alfa(ctnTime-1,2).*outputPDF_R( outputVec(ctnTime) , 2 , s , SNR )));
        %alfa_Cal_numerator(ctnTime,s)  = 0.5.*((alfa(ctnTime-1,1).*outputPDF_R( outputVec(ctnTime) , 1 , s , SNRdb))+(alfa(ctnTime-1,2).*outputPDF_R( outputVec(ctnTime) , 2 , s , SNRdb)));
        alfa_Cal_denominator = 0.5*(alfa(ctnTime-1,1)*(outputPDF_R(outputVec(ctnTime),1,1,SNR)+outputPDF_R(outputVec(ctnTime),1,2,SNR))+alfa(ctnTime-1,2)*(outputPDF_R(outputVec(ctnTime),2,1,SNR)+outputPDF_R(outputVec(ctnTime),2,2,SNR)));
        alfa(ctnTime,s) = alfa_Cal_numerator/alfa_Cal_denominator ;
        ctnTime = ctnTime + 1 ;
            if  s == 1 && ctnTime==Input_Size
                s = 2;
                ctnTime=2;
            end
            if s == 2 && ctnTime==Input_Size
                flag = 0;
            end
    end
end

