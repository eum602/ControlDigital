clear all; close all;clc
% Diseño 1
% -------------------------------------------------------------------------
% Modelo Continuo Motor DC
% -------------------------------------------------------------------------
%Ktb=0.042; Ra=5.3; J=4.3e-6; 
%K=1/Ktb; tau=(Ra*J)/Ktb^2; T=tau/5;
%Pc=tf(K,[tau 1]);

% Modelo de la Planta
%      b1 s + b2
%   ----------------
%   s^2 + a1 s + a2

b1=3.9495e4;
b2=8.4429e7;
a1=408.1;
a2=4.56383e4;

Pc=tf([b1 b2],[1 a1 a2]);
T=0.001;

[num,den]=tfdata(Pc,'v');
[dnum,dden]=c2dm(num,den,T);
% -------------------------------------------------------------------------
% Consideraciones de Diseño
% -------------------------------------------------------------------------
ts=0.00124;
Mp=5/100;
% -------------------------------------------------------------------------
% Formulas de Diseño de Polos Deseados
% -------------------------------------------------------------------------
zeta = sqrt((log(Mp)/pi)^2/(1+(log(Mp)/pi)^2));
wn = 4/(zeta*ts); 
% sd=-2.5+1j*0.4;
% sigmad=-real(sd);
% wd=imag(sd);
% teta=atan(-wd/sigmad);
% zeta=cos(teta);
% wn=sigmad/zeta;
PM=atan(2*zeta/(sqrt(-2*zeta^2+sqrt(1+zeta^4))));
wgc=(2*zeta*wn)/(tan(PM));
dsd=(exp(1j*wgc*T)-1)/T;
%dsd=exp(1j*wgc*T);
% figure
% [dkf,dnumcf,ddencf,dnumclf,ddenclf]=pd(dnum,dden,T,dsd,PM);
[num,den]=tfchk(dnum,dden);
% Sistema ang de G_sys(z0)
X=(polyval(num,dsd,T))/(polyval(den,dsd,T));
% El angulo deseado en el angulo limite
phi=-pi+abs(PM);
% Angulo theta_c
thc=phi-angle(X);
% parte real del punto diseñado
x0=-real(dsd);
% parte imaginaria del punto diseñado
y0=imag(dsd);
% thetap es el angulo del denominador del compensador
thp=atan(T*y0/(1-T*x0));
% thetaz es un angulo del numerador del compensador
th_alpha=thc+thp;
% Evaluacion del cero de compensador
z0=x0+y0/tan(th_alpha);
% El diseño del compensadro PD del sistema unificado
numc=[1 z0]
%denc=[T 1]
denc=[1 0];
% Evaluacion del compensador PD del sistema unificado y diseñado en elpunto
Y=(polyval(numc,dsd,T))/(polyval(denc,dsd,T));
% La ganancia de control
k=1/(abs(X)*abs(Y));
% FT LA del sistema despues del diseño del compensador
numol=conv(numc,num);
denol=conv(denc,den);
% Verificando dimension del sistema en lazo abierto
[numol,denol]=tfchk(numol,denol);
% FTLC del sistema despues del diseño
numcl=k*numol
dencl=denol+numcl
step(numcl,dencl,0.2)
figure()
step(Pc);
%%
%syms z
%numec=poly2sym(numc,'z');
%denoc=poly2sym(denc,'z');
%Gcz=k*numec/denoc;
Gcz=k*tf(numc,denc,T);
Gpz=c2d(Pc,T,'zoh');
% [n,d]=tfdata(Gpz,'v');
% n=poly2sym(n,'z');
% d=poly2sym(d,'z');
% Gpz=n/d;
pcont=feedback(Gpz*Gcz,1);
t=0:T:0.5;
u=ones(1,length(t));
figure
lsim(pcont,u,t)
