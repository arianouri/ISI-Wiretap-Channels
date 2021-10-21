clear;  clc;

SNR_dB=-6;
transfer_poly_bob=[0.792,0.610];
los_degradation=5; %(db)

%%

Es_N0_dB=SNR_dB-10*log10(2);
N0=10^(-Es_N0_dB/10);

transfer_poly_bob=transfer_poly_bob/sqrt(sum((transfer_poly_bob).^2));
power_bob_dB=20*log10(transfer_poly_bob);

max_los_eve=10^((power_bob_dB(1)-los_degradation)/20);

Es_target=1;

%%
Windex=1;
min_theta=eps;
max_theta=1e4;
tol=1e-4;

%
tp_lvl=100;
theta_quant=1000;
BW_quant=1e4;

%

W=180e3;
T=1./(2.*W);
f=linspace(-W,W,BW_quant);


N=tp_lvl^2;
theta_ans=zeros(1,N);
Es=zeros(1,theta_quant);
capacity=zeros(1,N);
transfer_poly=zeros(N,3);

%%

    for r0=linspace(0,max_los_eve,tp_lvl)
        for r1=linspace(0,sqrt(1-r0^2),tp_lvl)
        transfer_poly(Windex,:)=[r0,r1,sqrt(1-r0^2-r1^2)];
        
        H=0;
        for rIndex=1:length(transfer_poly(Windex,:))
            H=H+transfer_poly(Windex,rIndex)*exp(-1i*rIndex*2*pi*f*T);
        end
        mag_H=(abs(H)).^2;

           lv=min_theta;
           hv=max_theta;
           while(1)
               theta=mean([lv hv]);
%                [lv hv]

                    S_int_arg=theta-((N0/2)./(mag_H));
                    S_int_arg(S_int_arg<0)=0;
                    Es=trapz(f,S_int_arg);


               if Es-Es_target<-tol
                   lv=theta;
               elseif Es-Es_target>+tol
                   hv=theta;
               else
                   theta_ans(Windex)=theta;
                   break
               end
           end
            
                    capacity_int_arg=log2(theta_ans(Windex)./((N0/2)./(mag_H)));
                    capacity_int_arg(capacity_int_arg<0)=0;
                    capacity(Windex)=0.5*trapz(f,capacity_int_arg);

            [Es capacity(Windex)]
% % 
            Windex=Windex+1;
        end
    end
    transfer_poly=real(transfer_poly);