function [ alfa , o_PDF ] = alfaPDF_Calculation( rx_waveform ,VarEyAns , Trans , Trans_Pr , Trellis_Index_io )
rx_waveform_pac=[reshape(rx_waveform,size(Trans,4),length(rx_waveform)/size(Trans,4))]';
State_Cardinal=size(Trans_Pr,1);
Branch_Cardinal=size(Trans_Pr,3);
alfa = zeros(size(rx_waveform_pac,1),State_Cardinal);
o_PDF = zeros(size(rx_waveform_pac,1),1);
alfa_Cal_numerator = zeros(size(rx_waveform_pac,1),State_Cardinal);
alfa(1,1) = 1;
for ctnTime=2:size(rx_waveform_pac,1)
    for si=1:State_Cardinal
        for sj=[Trellis_Index_io(Trellis_Index_io(:,1)==si,2)]'
            for branchIndex=1:Branch_Cardinal
                o_PDF(ctnTime)=o_PDF(ctnTime)+Trans_Pr(si,sj,branchIndex).*alfa(ctnTime-1,si).*outputPDF_R(rx_waveform_pac(ctnTime,:),si,sj,branchIndex,Trans,VarEyAns);
                alfa_Cal_numerator(ctnTime,sj)=alfa_Cal_numerator(ctnTime,sj)...
                    +Trans_Pr(si,sj,branchIndex).*alfa(ctnTime-1,si).*outputPDF_R(rx_waveform_pac(ctnTime,:),si,sj,branchIndex,Trans,VarEyAns);
            end
        end
    end
    alfa(ctnTime,:) = alfa_Cal_numerator(ctnTime,:) ./ (repmat(o_PDF(ctnTime),1,size(alfa,2))+eps) ;
    alfa(ctnTime,:) = alfa(ctnTime,:)./sum(alfa(ctnTime,:));
end
end