clear all;close all;clc;
%diseño de control PID por el método del lugar de las raíces
%planta
%      b1 s + b2
%   ----------------
%   s^2 + a1 s + a2

b1=3.9495e4;
b2=8.4429e7;
a1=408.1;
a2=4.56383e4;
Gp=tf([b1,b2],[1 a1 a2])
T=0.001;
Gpz=c2d(Gp,T,'zoh')
%Variables de diseño: wn, zeta Mp=5%, tr=0.01seg
%Mp=100*exp(pi*zeta/(sqrt(1-zeta^2))
Mp=10;tr=0.01;
a=(log(Mp/100)/pi)^2;
zeta=sqrt(a/(1+a))
wn=1.8/tr
sd=-zeta*wn+i*wn*sqrt(1-zeta^2);
zd=exp(sd*T)
%Grafica Root locus
rlocus(Gpz),title('Lugar de las raíces de la planta en lazo abierto')
%controlador GcPID= kp+ki*T*z/(z-1)+kd/T*(z-1)/z;
%={z^2*(kp+T*ki+kd)+z*(-kp-2*kd/T)+kd/T}/((z-1)*z);
%condición de ángulo: ang_cont+ang_planta=180
[n,d]=tfdata(Gpz,'v');
syms z
n=poly2sym(n,'z');
d=poly2sym(d,'z');
Gplanta=n/d
ang_cont=-pi-angle(subs(Gplanta,'z',zd));
ang_den_c=angle(subs((z-1)*z,'z',zd))
ang_num_c=ang_cont+ang_den_c;
ang_cero_c=ang_num_c/2;
%cero es mag_cero_c=a---> imag(zd)/(real(zd)-a)=tan(ang_cero_c)
a=real(zd)-imag(zd)/tan(ang_cero_c)
Gcpidz=tf(conv([1 a],[1 a]),conv([1 -1],[1 0]),T);
[nc,dc]=tfdata(Gcpidz,'v');
nc=poly2sym(nc,'z');
dc=poly2sym(dc,'z');
Gc=nc/dc;
K=1/(abs(subs(Gplanta,'z',zd))*abs(subs(Gc,'z',zd)));
Gcerrado=feedback(K*Gpz*Gcpidz,1)
figure
rlocus(Gcerrado)
axis([-3 1.5 -3,3])
title('Lugar de las raíces de la plata en lazo cerrado')
figure
t=0:T:0.2;
step(Gcerrado,t),title('Respuesta de la planta controlada ante la entrada step')
title('Respuesta al escalon de la planta controlada')
u=sin(t*2*pi/0.2);
figure
y=lsim(Gcerrado,u,t);
plot(t,u,'b'),hold,stairs(t,y,'g')
title('Respuesta a una señal senoidal de la planta controlada')
