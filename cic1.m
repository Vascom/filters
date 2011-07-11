Fs=200;             %Sampling frequency
d=4;                %Divider after CIC filtering
F_pass = 10/(Fs/d/2);
F_stop = 13/(Fs/d/2);

m=[8 8 10 12 12];
m=[1 4 5 6 7];
m=[7 8 10 10 11 12];
m=[4 5 5 6 7 7];
m=[8 8 10 12 12];
%m=[1 9 11 13 14 15]; %20
%m=[3 5 5 7 7 7]; %40
m=[5 6 7 7 9 9]; %30
m=[1 2 4 4 5 6]; %55 - 190
m=[2 9 10 12 14 14]; %20 - 190 MHz
m=[4 5 5 7 9 9]; %33 - 190 MHz 
m=[5 5 5 5 7 7]; %41 - 190 MHz
m=[5 6 7 8 8 8]; %30 - 170 MHz
m=[1 5 5 5 5 7];
m=[2 4 4 5 5 6];
m=[8 10 11 13 15 15];

tt=1:1:512;
g=tt/512/d/2*F_pass;
ach1=1:1:512;
ach1=abs(sin(pi*m(1).*g)./sin(pi.*g)).*abs(sin(pi*m(2).*g)./sin(pi.*g)).*abs(sin(pi*m(3).*g)./sin(pi.*g)).*abs(sin(pi*m(4).*g)./sin(pi.*g)).*abs(sin(pi*m(5).*g)./sin(pi.*g)).*abs(sin(pi*m(6).*g)./sin(pi.*g));

N = 59;  % Order
% Frequency Vector
F = [linspace(0,F_pass,512) F_stop 1];
% Amplitude Vector
A = [1./ach1 0 0];
% Weight Vector
W = [ones(1,255) 0.1 0.1];
% Calculate the coefficients using the FIRLS function.
b  = firls(N, F, A, W);
Hd = dfilt.dffir(b);
set(Hd, 'Arithmetic', 'double');
hn = coeffs(Hd);
[h,w]=freqz(hn.Numerator);

g=tt/512/d/2*1;
ach1=1:1:512;
ach=abs(sin(pi*m(1).*g)./sin(pi.*g)).*abs(sin(pi*m(2).*g)./sin(pi.*g)).*abs(sin(pi*m(3).*g)./sin(pi.*g)).*abs(sin(pi*m(4).*g)./sin(pi.*g)).*abs(sin(pi*m(5).*g)./sin(pi.*g)).*abs(sin(pi*m(6).*g)./sin(pi.*g));
ach=ach.';

plot(g*Fs,20*log10(abs(h)/(abs(h(1))))+20*log10(ach/max(ach))); hold on;
plot(g*Fs,20*log10(abs(h)/(abs(h(1)))),'r');


plot(g*Fs,20*log10(ach/ach(1)),'g');

%hd = floor(hn.Numerator(1:30)/max(abs(hn.Numerator(1:30)))*(2^11) + 0.5);
%hd'