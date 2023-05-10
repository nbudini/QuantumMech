% SCRIPT TO SIMULATE THE TIME EVOLUTION OF WAVEFUNCTION
% OF A PARTICLE IN AN INFINITE SQUARE WELL POTENTIAL

clearvars
close all

% constants
m = 1; % masa del electrón (9.1e-31 kg)
a = 1; % ancho del pozo (1e-9 m);
A = sqrt(2/a); % constante de normalización de psi
hbar = 1; % constante de Planck (6.626e-34 J s);

% positions
x = linspace(0,a,500);

% stationary states of infinite square well
psi1 = A*sin(pi/a*x);
psi2 = A*sin(2*pi/a*x);
% and so on... add all you want in the form psin = A*sin(n*pi/a*x)

% energy eigenvalues of stationary states
E1 = 1^2*pi^2*hbar^2/2/m/a^2;
E2 = 2^2*pi^2*hbar^2/2/m/a^2;
% and so on... add all you want in the form En = n^2pi^2*hbar^2/2/m/a^2

% coefficients of the expansion for the initial wavefunction, Psi(x,0) = sum c_i psi_i
c1 = 1/sqrt(2);
c2 = sqrt(1 - c1^2);
% add all you want but keeping in mind normalization, i.e. that sum |c_i|^2 = 1

% initial state, Psi(x,0)
psi0 = c1*psi1.*exp(-1i*E1/hbar)+c2*psi2.*exp(-1i*E2/hbar);

% plot of initial state
figure,
h1 = plot(x,real(psi0));                        % real part plot
hold on
h2 = plot(x,imag(psi0));                        % imaginary part plot
h3 = plot(x,conj(psi0).*psi0,'LineWidth',2);    % square modulus plot
xlim([0 a])                                     % range for x axis
xticks([0 a/4 a/2 3*a/4 a])                     % add some ticks for x axis
xticklabels({'0','a/4','a/2','3a/4','a'})       % add some labels for x axis
ylim([-2 3.5])                                  % range for y axis
legend('Re(\Psi)','Im(\Psi)','|\Psi|^2')        % add legend
yticks('auto')                                  % add auto ticks for y axis
% ylabel('Re(Psi), Im(Psi), |Psi|^2')             % add some labels for y axis

tfin = 10; % simulation time
dt = 0.01; % timestep

for t = 0:dt:tfin
    % Psi(x,t) state
    psi = c1*psi1.*exp(-1i*E1*t/hbar)+c2*psi2.*exp(-1i*t*E2/hbar);
    h1.YData = real(psi);       % update real part plot
    h2.YData = imag(psi);       % update imaginary part plot
    h3.YData = conj(psi).*psi;  % update squared modulus plot
    
    refreshdata                 % refresh plots
    drawnow                     % draw and show
    pause(0.01)                 % pause for visualization
end
