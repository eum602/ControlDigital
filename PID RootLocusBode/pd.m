function [k,numc,denc,numcl,dencl]=pd(num,den,T,sd,PM)
% Funcion de transferencia en lazo abierto
[num,den]=tfchk(num,den);
% Sistema ang de G_sys(z0)
X=(polyval(num,sd,T))/(polyval(den,sd,T));
% El angulo deseado en el angulo limite
phi=-pi+abs(PM);
% Angulo theta_c
thc=phi-angle(X);
% parte real del punto diseñado
x0=-real(sd);
% parte imaginaria del punto diseñado
y0=imag(sd);
% thetap es el angulo del denominador del compensador
thp=atan(T*y0/(1-T*x0));
% thetaz es un angulo del numerador del compensador
th_alpha=thc+thp;
% Evaluacion del cero de compensador
z=x0+y0/tan(th_alpha);
% El diseño del compensadro PD del sistema unificado
numc=[1 z];
denc=[T 1];
% Evaluacion del compensador PD del sistema unificado y diseñado en elpunto
Y=(polyval(numc,sd,T))/(polyval(denc,sd,T));
% La ganancia de control
k=1/(abs(X)*abs(Y));
% FT LA del sistema despues del diseño del compensador
numol=conv(numc,num);
denol=conv(denc,den);
% Verificando dimension del sistema en lazo abierto
[numol,denol]=tfchk(numol,denol);
% FTLC del sistema despued del diseño
numcl=k*numol;
dencl=denol+numcl;

if T==0
    % FTLA
    Gop=tf(numol,denol);
else
    % Convierte
    [numoz,denoz]=c2dm(numol,denol,T);
    Gop=tf(numoz,denoz,T);
end
Gcl=minreal(k*Gop/(1+k*Gop));
subplot(211)
step(numcl,dencl);
title('Respuesta Step de Compensador PD')
subplot(212)
if PM==0
rlocus(numol,denol);
title('Root Locus del compensadro PD')
hold on
pol=roots(dencl);
zer=roots(numcl);
plot(real(pol),imag(pol),'^',real(zer),imag(zer),'s')
hold off
else
    bode(k*numol,denol)
    title('Respuesta Bode del Compensador PD')
end

