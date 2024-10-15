function SeCap_Curve( Capacity , SNRdB ,StepSize,SNRdB_min,SNRdB_max)
Capacity=smooth(smooth(smooth(Capacity)));
[X,Y]=meshgrid(SNRdB_min:StepSize:SNRdB_max);
% xIndex for EVE's SNRdB
% yindex for BOB's SNRdB
Z=zeros(size(X));
XI=0;
YI=0;
for yIndex=Y(:,1)'
    YI=YI+1;
    for xIndex=X(1,:)
        XI=XI+1;
        if yIndex>xIndex
        % the capacity is zero otherwise
            Z(XI,YI)=Capacity(round(SNRdB*10)/10==round(yIndex*10)/10)...
                    -Capacity(round(SNRdB*10)/10==round(xIndex*10)/10);
        end
    end;XI=0;
end;YI=0;
figure
h=surf(Y,X,Z);
xlabel 'EVE''s SNR (dB)'
ylabel 'BOB''s SNR (dB)'
zlabel 'Secrecy Capacity [bits/channel-use]'
end

