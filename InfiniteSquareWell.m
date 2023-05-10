% SCRIPT PARA SIMULAR LA EVOLUCIÓN TEMPORAL DE LA FUNCIÓN DE ONDA EN UN
% POZO INFINITO DE POTENCIAL

clearvars
close all

% definición de constantes
m = 1; % masa del electrón (9.1e-31 kg)
a = 1; % ancho del pozo (1e-9 m);
A = sqrt(2/a); % constante de normalización de psi
hbar = 1; % constante de Planck (6.626e-34 J s);

% posiciones
x = linspace(0,a,500);

% estados estacionarios del pozo infinito
psi1 = A*sin(pi/a*x);
psi2 = A*sin(2*pi/a*x);
% psi3 = A*sin(3*pi/a*x);
% etc... se pueden agregar los que se quieran

% autovalores de las energías de estados estacionarios
E1 = 1^2*pi^2*hbar^2/2/m/a^2;
E2 = 2^2*pi^2*hbar^2/2/m/a^2;
% E3 = 3^2*pi^2*hbar^2/2/m/a^2;
% etc... se pueden agregar las que se quieran

% coeficientes del desarrollo del estado inicial Psi(x,0) = sum c_i psi_i
c1 = 1/sqrt(2);
c2 = sqrt(1 - c1^2);
% c3 = se pueden agregar más con la condición sum c_i^2 = 1

% estado inicial Psi(x,0)
psi0 = c1*psi1.*exp(-1i*E1/hbar)+c2*psi2.*exp(-1i*E2/hbar);

% gráfico del estado inicial
figure,
h1 = plot(x,real(psi0));                        % parte real
hold on
h2 = plot(x,imag(psi0));                        % parte imaginaria
h3 = plot(x,conj(psi0).*psi0,'LineWidth',2);    % módulo cuadrado
xlim([0 a])                                     % rango de x
xticks([0 a/4 a/2 3*a/4 a])
xticklabels({'0','a/4','a/2','3a/4','a'})       % etiquetas de x
ylim([-2 3.5])                                  % rango de y
legend('Re(\Psi)','Im(\Psi)','|\Psi|^2')                % leyenda
yticks('auto')
% ylabel('Re(Psi), Im(Psi), |Psi|^2')             % etiquetas de y

tfin = 10; % tiempo de simulación
dt = 0.01; % paso de tiempo

for t = 0:dt:tfin
    % estado Psi(x,t)
    psi = c1*psi1.*exp(-1i*E1*t/hbar)+c2*psi2.*exp(-1i*t*E2/hbar);
    h1.YData = real(psi);       % parte real
    h2.YData = imag(psi);       % parte imaginaria
    h3.YData = conj(psi).*psi;  % módulo cuadrado
    
    refreshdata                 % refresca el gráfico
    drawnow                     % ídem
    pause(0.01)                 % pausa breve para visualización
end