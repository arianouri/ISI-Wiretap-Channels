clear;  clc;

SNR_dB=-5.0;
% transfer_poly=[0.792,0.610];
transfer_poly=[0.445516026180429,0.633021994668546,0.633086585454355];

min_bw=1e-4;
max_bw=500;

%%

Es_N0_dB=SNR_dB-10*log10(2);
N0=10^(-Es_N0_dB/10);

transfer_poly=transfer_poly/sqrt(sum((transfer_poly).^2));
Es_target=1;

%%
Windex=1;
min_theta=eps;
max_theta=1e4;
tol=1e-4;

%

N=500;
theta_quant=1000;
BW_quant=1e4;

%

bandwidth_range=linspace(min_bw,max_bw,N);

W_vec=zeros(1,N);
theta_ans=zeros(1,N);
theta_vec=zeros(1,N);
Es=zeros(1,theta_quant);
capacity=zeros(1,N);

%%

    for W=bandwidth_range
%         W=180e3;
        W_vec(Windex)=W;        
        T=1./(2.*W);
        f=linspace(-W,W,BW_quant);
        H=0;
        for rIndex=1:length(transfer_poly)
            H=H+transfer_poly(rIndex)*exp(-1i*rIndex*2*pi*f*T);
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
        
                capacity_int_arg=log(theta_ans(Windex)./((N0/2)./(mag_H)));
                capacity_int_arg(capacity_int_arg<0)=0;
                capacity(Windex)=0.5*trapz(f,capacity_int_arg);
        
        [Es capacity(Windex)]
        
        Windex=Windex+1;
    end
    plot(W_vec,capacity)