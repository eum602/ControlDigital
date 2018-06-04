clc,clear all, close all

% Modelo de la Planta
%           a(1) z^2 + a(2) z + a(3)
% Gp(z) =  -------------------------
%           b(1) z^2 + b(2) z + b(3)

% Modelo del Controlador
%           c(1) z^2 + c(2) z + c(3)
% Gc(z) =  -------------------------
%           d(1) z^2 + d(2) z + d(3)

% Datos Planta
a = [0 69.06 0.001401];
b = [1 -1.628 0.6649];
% Datos Controlador
c = [0.0126 -0.0205 0.0083];
d = [2 -2 0];

Tf = 500;   % numero de muestras
T = 0.001;  % Periodo de muestreo
tk = 0:T:(Tf-1)*T;

%% Simulacion de la Planta Lazo Abierto
% Valores Iniciales
ypa(1)=0;ypa(2)=0;
u(1)=0;u(2)=0;

for i=3:Tf
    u(i)=1;
    ypa(i)=(a(1)*u(i)+a(2)*u(i-1)+a(3)*u(i-2)-b(2)*ypa(i-1)-b(3)*ypa(i-2))/b(1);
end
%figure()
%stairs(ypa);

%% Simulacion de la Planta Lazo Cerrado
% Valores Iniciales
e(1)=0;e(2)=0;
ypc(1)=0;ypc(2)=0;
u(1)=0;u(2)=0;

for i=3:Tf
    u(i)=1;
    e(i)= u(i)-ypc(i-1);
    ypc(i)=(a(1)*e(i)+a(2)*e(i-1)+a(3)*e(i-2)-b(2)*ypc(i-1)-b(3)*ypc(i-2))/b(1);
end
% figure()
% stairs(ypc);

%% Simulacion de la Planta y el Controlador

e(1)=0;e(2)=0;
y(1)=0;y(2)=0;
u(1)=0;u(2)=0;
Ct(1)=0;Ct(2)=0;
for i=3:Tf
    u(i)=1;
    e(i)= u(i)-y(i-1);
    Ct(i)=(c(1)*e(i)+c(2)*e(i-1)+c(3)*e(i-2)-d(2)*Ct(i-1)-d(3)*Ct(i-2))/d(1);
    y(i)=(a(1)*Ct(i)+a(2)*Ct(i-1)+a(3)*Ct(i-2)-b(2)*y(i-1)-b(3)*y(i-2))/b(1);
end
%figure()
%stairs(y);

%% Graficando
figure()
stairs(tk,ypa,'b','LineWidth',1);
xlabel('\bf t(seg)')
ylabel('\bf y(t)')
legend('Planta Lazo Abierto',4);
figure()
stairs(tk,y,'r','LineWidth',1);
xlabel('\bf t(seg)')
ylabel('\bf y(t)')
legend('Planta Controlada',4);
