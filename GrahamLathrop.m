clear all; close all; clc
%(39495.3 s+8.44259×?10?^7)/(s^2+408.1 s+45638.3)
 %Ziegler Nichols
 K=8.44259e7;
 %planta: K/a1.s^2+a2.s+a3
 a1=1;
 a2=408.1;
 a3=45638.3;
 %controlador:kp+ki+/s+kd*s;
 %ec. carcterística en lazo cerrado:
 %s^3+s^2.(a2+k.kd)/a1+s.(k.kp+a3)/a1+k.ki/a1;
 %Ecuacion de Graham Lathrop: s^3+1.75w0^2.s^2+2.15w0^2.s+w0^3
Gp=tf(K,[1 408.1 45638.3])%zeta=0.56,wn=370.54-->2pi/wd=0.0205-->Ts=0.00205
syms s 
T=0.001;
tfin=0.3;
[n,d]=tfdata(Gp,'v');
n=poly2sym(n,'s');
d=poly2sym(d,'s');
Ys=n/d*(1/s);
yt=ilaplace(Ys);
y1=subs(yt,'t',0:T:tfin,'b');
%cálculo ante una entrada de 1000 RPM
%ecuación de trensformacion es:(w+585.28)/963.28=voltaje
t=0:T:tfin;
u1=ones(size(t));%entrada de 1 voltio
%entrada a la planta sin controlar
y=lsim(Gp,u1,t);
plot(t,y,'b')
xlabel('\bf t(seg)')
ylabel('\bf y(t)')
legend('Respuesta sin controlador ante entrada escalon de 1 voltio',4)
print -f -dbitmap sincontrol

%% parametros de control
w0=sqrt(a3);
kp=(2.15*w0^2*a1-a3)/K;%0.0059317003679554894251615232103691
ki=(w0^3*a1)/K;%kp/Ti;%2.1969826759301782687311571712613
kd=(1.75*w0*a1-a2)/K;%kp*Td;%0.0000040037945725159527967637526555352
Gc=tf([kd kp ki],[1 0]);
Gtotal=feedback(Gc*Gp,1);
u=1000*ones(size(t));%entrada en RPM
figure
ycont=lsim(Gtotal,u,t);
plot(t,ycont,'g'),legend('controlada');
xlabel('\bf t(seg)')
ylabel('\bf y(t)')
print -f -dbitmap concontrolgraham
%%  discretizando la planta
figure
Gpz=c2d(Gp,T,'zoh');
[ncp,dcp]=tfdata(Gpz,'v');
dato1=[ncp,dcp];
save coef_planta_discreto.lvm dato1 -ascii -tabs
[nc,dc]=tfdata(Gc,'v');
nc=poly2sym(nc);
dc=poly2sym(dc);
syms z
c=nc/dc;
c=subs(c,'x',(2/T)*(z-1)/(z+1));
[nc,dc]=numden(c);
nc=sym2poly(nc);
dc=sym2poly(dc);
Gcz=tf(nc,dc,T); % Función de transferencia del controlador
[num den]=tfdata(Gcz,'v');
dato=[num den];
save coef_PID_discreto.lvm dato -ascii -tabs
Gtotalz=feedback(Gpz*Gcz,1);
yk=lsim(Gtotalz,u,t);
stairs(t,yk,'r'),legend('control discreto',4);
print -f -dbitmap controlGrahamdiscreto
%
plot(t,y,'b',t,ycont,'g'),hold,stairs(t,yk,'r'),legend('Y sin control','Y controlada','Y controlada discreta',4)
print -f -dbitmap comparativaGraham;

