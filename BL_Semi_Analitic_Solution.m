%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Este programa constroi o perfil de saturacao de um escoamento bifasico,
%  uinidimensional (oleo e agua), utilizando a equacao de Buckley-Leverett 
%  e o metodo semi-analitico de Welge para sua resolucao                                                                     %                                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

%   Dados de entrada
npoin=1000;  % Numero de pontos
Swr=0.1;     % Saturação residual da água
Sor=0.1;     % Saturação residual do óleo
mi_w=1;      % Viscosiade da água
mi_o=1;      % Viscosiade do óleo
A=1;         % Area normal ao escoamento
fi=1.0;      % Porosidade do meio
q=1;         % Taxa de injeçao de água
t=0.2;       % Tempo 
L=1;         % Distância entre injetor e produtor

%   Saturaçao no choque (Método da secante)

%   Valores iniciais para o metodo
Swf=1-Sor;   
Sw0=Swf/2;
    
while abs(Swf-Sw0)>1e-12
    Sw=Swf;
    
%   Calculo dos pontos da funcao de fluxo fracionario
    fw0=1/(1+mi_o/mi_w*((1-Sw0-Swr)/(Sw0-Swr))^2);
    fw=1/(1+mi_o/mi_w*((1-Sw-Swr)/(Sw-Swr))^2);
    
%   Calculo da derivada da funcao de fluxo fracionario
    dfw_dSw0=-2*mi_w/mi_o*((1-Sw0-Swr)/(Sw0-Swr)^3)*(2*Swr-1)/(1+mi_w/mi_o*((1-Sw0-Swr)/(Sw0-Swr))^2)^2;
    dfw_dSw=-2*mi_w/mi_o*((1-Sw-Swr)/(Sw-Swr)^3)*(2*Swr-1)/(1+mi_w/mi_o*((1-Sw-Swr)/(Sw-Swr))^2)^2;
    
%   Equaçao da reta tangente    
    y0=dfw_dSw0*(Sw0-Swr);
    y=dfw_dSw*(Sw-Swr);
    
%   Swf=> Zero da funçao g=y-fw
    g0=y0-fw0;
    g=y-fw;
    Swf=Sw-(Sw-Sw0)*g/(g-g0);
    Sw0=Sw;
end
    
%   Eixo x
delt_x=L/(npoin);
x=0+delt_x/2:delt_x:L-delt_x/2;

% for t=0.2:0.05:0.6
%   Posicao do choque
dfw_dSw_Swf=-2*mi_w/mi_o*((1-Swf-Swr)/(Swf-Swr)^3)*(2*Swr-1)/(1+mi_w/mi_o*((1-Swf-Swr)/(Swf-Swr))^2)^2;
xf=((q*t)/(A*fi))*dfw_dSw_Swf;

Sw(2:length(x))=Swf;
Sw(1)=1-Sor;
ii=1;
while x(ii)<=xf
    ii=ii+1;
    Sw0=.9;
    while abs(Sw(ii)-Sw0)>1e-12
%       Calculo da derivada da funcao de fluxo fracionario
        dfw_dSw0=-2*mi_w/mi_o*((1-Sw0-Swr)/(Sw0-Swr)^3)*(2*Swr-1)/(1+mi_w/mi_o*((1-Sw0-Swr)/(Sw0-Swr))^2)^2;
        dfw_dSw=-2*mi_w/mi_o*((1-Sw(ii)-Swr)/(Sw(ii)-Swr)^3)*(2*Swr-1)/(1+mi_w/mi_o*((1-Sw(ii)-Swr)/(Sw(ii)-Swr))^2)^2;
%       Saturacao em cada x
        Sw1=Sw(ii)-(Sw(ii)-Sw0)*(dfw_dSw-A*fi/(q*t)*x(ii))/(dfw_dSw-dfw_dSw0);
        Sw0=Sw(ii);
        Sw(ii)=Sw1;
    end
end
Sw(ii:length(x))=Swr;

% Grafico
% hold all;
% plot(x,Sw,'r-');
% legend('Solução Semi-analítica de Welge');
% xlabel('x');
% ylabel('Sw');
% grid on;

% end

hold all;
plot(x,Sw,'r-');
legend('Solucao Aproximada: Upwind 1ª Ordem','Solucao Semi-analitica de Welge');


% % Fluxo x Saturação
% figure(2)
% Sw_all = Swr:0.01:1-Sor;
% fw=1./(1+mi_o/mi_w.*((1-Sw_all-Swr)./(Sw_all-Swr)).^2);
% y=dfw_dSw.*(Sw_all-Swr);
% 
% plot(Sw_all, fw, Sw_all, y); axis([0 1 0 1]);
% xlabel('S_w'); ylabel('F_w');
% title('Fluxo Fracionário x Saturação de Água');
% grid on; 