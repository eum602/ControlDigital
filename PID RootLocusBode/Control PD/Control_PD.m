clear all;close all;clc;
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
%compesador Proprocional derivativo de la forma (z-alfa)/z
wn=700%sqrt(a2);%900
zeta=a1/(2*wn)%0.85;
sd=-zeta*wn+i*wn*sqrt(1-zeta^2);%polo deseado
zd=exp(sd*T);%x0+y0*i;
x0=real(zd);
y0=imag(zd);
[n,d]=tfdata(c2d(Pc,T,'zoh'),'v');
syms z
n=vpa(poly2sym(n,'z'),4);
d=vpa(poly2sym(d,'z'),4);
G2=n/d;
teta_c=pi-angle(subs(G2,'z',zd))%ang. controlador-crit. angulo
ang_den=angle(zd);
%ang_num-ang_den=teta_c
ang_num=teta_c+ang_den;
%sea 
%ang_num=atan(y0/(x0-alfa))
alfa=x0-y0/tan(ang_num)
%controlador
Gcpdz=tf([1 -alfa],[1 0],T);
Gpz=c2d(Pc,T,'zoh');
G1=Gcpdz*Gpz;
[n,d]=tfdata(G1,'v');
syms z
n=vpa(poly2sym(n,'z'),4);
d=vpa(poly2sym(d,'z'),4);
G2=n/d;
k_c=1/abs(subs(G2,'z',zd));
G_lcerrado=feedback(k_c*G1,1);
step(G_lcerrado)

