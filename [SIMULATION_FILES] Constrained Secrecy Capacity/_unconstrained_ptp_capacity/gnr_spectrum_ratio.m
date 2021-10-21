clear; close all; clc

%% Bob's str

BOB.SNR_dB=-6.5;
BOB.transfer_poly=[0.445516026180429,0.633021994668546,0.633086585454355];

%% Eve's str

EVE.SNR_dB=-6;
EVE.transfer_poly=[0.792,0.610];

%% 
BOB.Es_N0_dB=BOB.SNR_dB-10*log10(2);
BOB.N0=10^(-BOB.Es_N0_dB/10);
BOB.transfer_poly=BOB.transfer_poly/sqrt(sum((BOB.transfer_poly).^2));
%
EVE.Es_N0_dB=EVE.SNR_dB-10*log10(2);
EVE.N0=10^(-EVE.Es_N0_dB/10);
EVE.transfer_poly=EVE.transfer_poly/sqrt(sum((EVE.transfer_poly).^2));

%%

W=1;
T=1./(2.*W);
f=linspace(-W,W,1e4);

G_BOB=BOB.transfer_poly;
G_EVE=EVE.transfer_poly;

%%
H_BOB=0;
for rIndex=1:length(G_BOB)
    H_BOB=H_BOB+G_BOB(rIndex)*exp(-1i*rIndex*2*pi*f*T);
end
mag_H_BOB=(abs(H_BOB)).^2;
%
H_EVE=0;
for rIndex=1:length(G_EVE)
    H_EVE=H_EVE+G_EVE(rIndex)*exp(-1i*rIndex*2*pi*f*T);
end
mag_H_EVE=(abs(H_EVE)).^2;

%%
plot(f,10*log10((mag_H_BOB)./(BOB.N0/2)),'r'); hold on;
plot(f,10*log10((mag_H_EVE)./(EVE.N0/2)),'black')