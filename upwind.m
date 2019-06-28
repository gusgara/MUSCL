%% Buckley Leverett upwind terceira ordem (MUSCL k=1/3)
clear all; clc;
%% Dados da malha
ne = 2000;           % Saturacao Irredutivel da Agua
xl = 0;
xr = 1;
dx = (xr - xl)/ne;
x = xl:dx:xr; 
xe=x(1)+dx/2:dx:x(end)-dx/2;
CFL = 0.01;          % Condicao de CFL
dt = CFL*dx;  t =0; tmax = 0.2;
nsteps = round(tmax/dt);

%% Variáveis de entrada
Swi = 0.1;          % Saturacao Irredutivel da Agua 
Sor = 0.1;          % Saturacao Residual do Oleo 
v = 1;              % Velocidade Superficial 
mi_o = 1;           % Viscosidade do Oleo 
mi_w = 1;           % Viscosidade da Agua 

M = mi_w/mi_o;      % Constante - Razao de Viscosidade 
%% Condições de Contorno e Condições iniciais
Sw(1) = 1 - Sor;                 % Condicao de Contorno - Saturacao Prescrita 
Sw(2:ne) = Swi;           % Condicao Inicial 
%% Loop
f = @(Sw,Swi,M) ((Sw-Swi).^2)./(((Sw-Swi).^2) + M.*(1 - Sw - Swi).^2); % função do fluxo fracional da agua em cada no 
%df =@(Sw,Swi,M) 2*M.*(Sw+Swi-1).*(Sw-Swi)./((M.*(Sw+Swi-1).^2)+(Sw-Swi).^2).^2; % derivada do fluxo fracional
df =@(Sw,Swi,M) -2.*((Sw-1).*(((Sw - Swi).^2)+M.*((1 - Sw - Swi).^2))-((Sw - Swi).^2).*((1+M).*Sw-M+(M-1).*Swi))./(((Sw - Swi).^2)+M.*((1 - Sw - Swi).^2)).^2;
 n = 1;
while t < tmax
 krw = (Sw - Swi).^2;      % Permeabilidade Relativa da Agua
 kro = (1 - Sw - Swi).^2;  % Permeabilidade Relativa do Oleo
 DF = df(Sw,Swi,M);
    for i = 1:ne
        dt_aux = (CFL*dx)/(v*DF(i));
        if (dt_aux > 0) && (dt_aux < dt)
            dt = dt_aux;
        end
    end
     t=t+dt;
    [fR,fL] = MUSCL(Sw,f,df,Swi,M,dx);
    Sw = Sw - v*(dt/dx)*(fR-fL);
    n = n +1;
end
ordem = 3;
figure(1);
hold on;
plot(xe,Sw,'k.-');
axis([xl xr 0 1]);
title('MUSCL - upwind de 3ª ordem')
xlabel('Posição (x)');
ylabel('Saturação de Água (Sw)');
grid on;

BL_Semi_Analitic_Solution;

