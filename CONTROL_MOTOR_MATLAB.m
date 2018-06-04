clear all; close all; clc;
% Kp=4.531e04;
% Gp=tf(Kp,[1 126.8 4.502e04]);
% Planta Identificada
% Kp=92857.972;
% Gp=tf(Kp,[1 228.2267 92716.7327]);
Kp= 8.44259e7;
Gp=tf(Kp,[1 408.1 45638.3]);

[P,Z]=pzmap(Gp);
figure
t=0:0.001:0.5;
u=ones(size(t));
y=lsim(Gp,u,t);%Y lazo abierto
plot(t,y,'r','LineWidth',2)
hold
plot(t,u,'k','LineWidth',2)
xlabel('\bf t(seg)')
ylabel('\bf y(t)')
legend('Respuesta sin controlador',4)
print -f -dbitmap sincontrol
% Consideraciones de diseño
zeta=0.5066;
ts=0.06;
Mp=exp(-zeta*pi/(sqrt(1-zeta^2)))*100;
wn=5/(ts*0.5);
wd=wn*sqrt(1-zeta^2);
sigma=zeta*wn;
Pd=-sigma + 1j*wd;
Td=2*pi/wd;
% Calculo de magnitudes
r1=sqrt(real(Pd)^2+imag(Pd)^2);
theta1=atan(-wd/sigma);
pr2=-(P(2)-(Pd));
r2=sqrt(real(pr2)^2+imag(pr2)^2);
theta2=atan(imag(pr2)/real(pr2));%59.9731
pr3=-(P(1)-(Pd));
r3=sqrt(real(pr3)^2+imag(pr3)^2);
theta3=atan((imag(pr3)/real(pr3)));%33.9189
% ----------------------------------------------------------------------
theta1=pi+theta1; %atan(-wd/sigma);120.4376°
theta4=pi+theta1+theta2+theta3;
theta4=theta4-2*pi; 
a=wd/tan(theta4)+sigma;%cero del controlador
K=abs((r1*r2*r3)/(Kp*a));
% ----------------------------------------------------------------------
Gc=tf(K*[1 a],[1 0]); % Controlador PI
S=series(Gc,Gp);
F=S/(1+S);
u1=1000*ones(size(t));
yc=lsim(F,u1,t);
figure
plot(t,yc,'k')
hold
plot(t,u1,'b')
xlabel('\bf t(seg)')
ylabel('\bf y(t)')
legend('Respuesta con controlador PI')
print -f -dbitmap concontrolPI
% ----------------------------------------------------------------------
% Rediseño tustin
% ----------------------------------------------------------------------
figure
T=0.001;
%T=Td/10;
% Discretizando la planta
Gpz=c2d(Gp,T,'zoh');
[numpz,denpz]=tfdata(Gpz,'v');
% ----------------------------------------------------------------------
% Discretizando el controlador por Tustin
syscz2=c2d(Gc,T,'tustin');
[nt,dt]=tfdata(syscz2,'v');
dato=[nt dt];
save coef_tustin.lvm dato -ascii -tabs
syspz2=feedback(syscz2*Gpz,1);
tk=0:T:0.5;
u2=ones(size(tk));
yk=1000*step(syspz2,tk);
stairs(tk,yk,'g')
hold
plot(tk,u2,'r')
legend('Control con Tustin')
print -f -dbitmap concontrolportustin
%axis([0 0.5 0 1.2])
figure
plot(t,y,'r',t,yc,'g',t,yk,'b')
legend('Sin control','Control PI','Control PI discreto Tustin',4)
print -f -dbitmap comparativa
%% sea Gp=h/(s^2+a*s+b)=y/u
% Planta discretizada y(k)=m*u(k)+n*y(k-1)+p*y(k-2)
h=8.44259e7;
a=408.1;
b=45638.3;
m=h/(1/T^2+b+a/T);
n=(2/T^2+a/T)/(1/T^2+b+a/T);
p=(1/T^2)/(1/T^2+b+a/T);