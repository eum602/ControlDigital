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
wn=300%300%207%sqrt(a2);%207.97;%208%
zeta=a1/(2*wn);%0.85;
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
teta_c=pi-angle(subs(G2,'z',zd));%ang. controlador-crit. angulo
ang_alfa1=angle(zd-0.8138);
angulo_alfa1=ang_alfa1*180/pi
teta_c_prima=teta_c-ang_alfa1+angle(zd-1);
teta_cprima=teta_c_prima*180/pi
alfa2=x0-y0/(tan(teta_c_prima))+pi/40;
ang_alfa2=angle(zd-alfa2);
ang_alfados=ang_alfa2*180/pi
beta=x0-y0/tan(ang_alfa2-teta_c_prima);
ang_denominador=angle(zd-beta)*180/pi
%controlador
Gcpdz=tf(conv([1 -0.8138],[1 -alfa2]),conv([1 -1],[1 -beta]),T);
Gpz=c2d(Pc,T,'zoh');
G1=Gcpdz*Gpz;
[n,d]=tfdata(G1,'v');
syms z
n=vpa(poly2sym(n,'z'),4);
d=vpa(poly2sym(d,'z'),4);
G2=n/d;
k_c=1/abs(subs(G2,'z',zd));
G_lcerrado=feedback(k_c*G1,1);
t=0:T:0.5;
bb=step(Gpz,t)/1000;
plot(t,bb);
hold 
step(G_lcerrado,t)


