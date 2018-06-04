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
wn=1000%207.97;%208%sqrt(a2);
zeta=a1/(2*wn)%0.85;
sd=-zeta*wn+i*wn*sqrt(1-zeta^2);%polo deseado
zd=exp(sd*T);%x0+y0*i;
x0=real(zd);
y0=imag(zd);
[n,d]=tfdata(c2d(Pc,T,'zoh'),'v');
d1=poly2sym(d,'z');
r=solve(d1);
alfa1=real(r(1));
syms z
n=vpa(poly2sym(n,'z'),4);
d=vpa(poly2sym(d,'z'),4);
G2=n/d;
teta_c=pi-angle(subs(G2,'z',zd))%ang. controlador-crit. angulo
ang_den=angle(zd*(zd-1));
%ang_num-ang_den=teta_c
ang_alfa1=angle(zd-0.8138);
alfa1=vpa(alfa1,4);
ang_alfa2=teta_c+ang_den-pi-ang_alfa1;
%sea 
%ang_alfa2=atan(y0/(x0-alfa))
alfa2=x0-y0/tan(ang_alfa2)
%controlador
Gcpdz=tf(conv([1 -0.8138],[1 -0.8138]),conv([1 0],[1 -1]),T);
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

