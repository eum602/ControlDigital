clear all; close all; clc;

% Modelo de la Planta
%      b1 s + b2
%   ----------------
%   s^2 + a1 s + a2

b1=64.85;
b2=1.378e5;
a1=417;
a2=1.373e5;
Gp=tf([b1 b2],[1 a1 a2]);

[num,den]=tfdata(Gp,'v');
sd=-1+1j*0.5;
T=0;
PM=0;
figure(1)
[kc,numc,denc,numcl,dencl]=pd(num,den,T,sd,PM);
% ------------------------------------------------------------------------
% Root Locus Tiempo Discreto
% ------------------------------------------------------------------------
T=0.01;
%T=tau/5;
[dnum,dden]=c2dm(num,den,T);
dsd=(exp(sd*T)-1)/T;
figure(2)
[dk,dnumc,ddenc,dnumcl,ddencl]=pd(dnum,dden,T,dsd,PM);
% ------------------------------------------------------------------------
% Bode Tiempo Continuo
% ------------------------------------------------------------------------
T=0;
sd=-2.5+1j*0.4;
sigmad=-real(sd);
wd=imag(sd);
% ------------------------------------------------------------------------
% Bode Tiempo Discreto
% ------------------------------------------------------------------------
teta=atan(-wd/sigmad);
zeta=cos(teta);
wn=sigmad/zeta;
PM=atan(2*zeta/(sqrt(-2*zeta^2+sqrt(1+zeta^4))));
wgc=(2*zeta*wn)/(tan(PM));
sd=1j*wgc;
figure(3)
[kf,numcf,dencf,numclf,denclf]=pd(num,den,T,sd,PM);
T=0.01;
%T=tau/5;
dsd=(exp(1j*wgc*T)-1)/T;
figure(4)
[dkf,dnumcf,ddencf,dnumclf,ddenclf]=pd(dnum,dden,T,dsd,PM);
