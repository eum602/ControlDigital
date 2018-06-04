clear all; close all; clc
%(39495.3 s+8.44259×?10?^7)/(s^2+408.1 s+45638.3)
 %Ziegler Nichols
 Kp=8.44259e7;
Gp=tf(Kp,[1 408.1 45638.3])%zeta=0.56,wn=370.54-->2pi/wd=0.0205-->Ts=0.00205
syms s 
T=0.0001;
tfin=1;
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
pendiente=diff(yt,'t');
h=max(subs(pendiente,'t',0:T:tfin));
p=diff(diff(yt,'t'));
tcrit=vpa(solve(p,'t'));
%tcrit de punto de inflexion en 0.0032
m=vpa(subs(pendiente,'t',tcrit));
b=vpa((subs(yt,'t',tcrit)-m*tcrit));
L=-b/m;%corte de la curva con el eje t
a=-b;
t=0:T:tfin;
ecu=m*t+b;
%% parametros de control
kp=1.2/a;%0.0059317003679554894251615232103691
Ti=2*L;
Td=0.5*L;
ki=0.5/L*1/a;%kp/Ti;%2.1969826759301782687311571712613
kd=0.5*L*1/a;%kp*Td;%0.0000040037945725159527967637526555352
Gc=tf([0.0000040037945725159527967637526555352 0.0059317003679554894251615232103691 2.1969826759301782687311571712613],[1 0]);
Gtotal=feedback(Gc*Gp,1);
u=1000*ones(size(t));%entrada en RPM
figure
ycont=lsim(Gtotal,u,t);
plot(t,ycont,'g'),legend('controlada');
xlabel('\bf t(seg)')
ylabel('\bf y(t)')
print -f -dbitmap concontrolzieglercontinuo
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
print -f -dbitmap controlzieglerdiscreto
%
plot(t,y,'b',t,ycont,'g'),hold,stairs(t,yk,'r'),legend('Y sin control','Y controlada','Y controlada discreta',4)
print -f -dbitmap comparativa;

