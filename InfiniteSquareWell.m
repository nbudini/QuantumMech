clearvars
close all

m = 1; %9.1e-31;
a = 1; %1e-9;
A = sqrt(2/a);
hbar = 1; %6.626e-34;

%t = linspace(0,1e-10,500);
x = linspace(0,a,500);

psi1 = A*sin(pi/a*x);
psi2 = A*sin(2*pi/a*x);

E1 = 1^2*pi^2*hbar^2/2/m/a^2;
E2 = 2^2*pi^2*hbar^2/2/m/a^2;

c1 = 1/sqrt(2);
c2 = sqrt(1 - c1^2);%1/sqrt(2);

psi = c1*psi1.*exp(-1i*E1/hbar)+c2*psi2.*exp(-1i*E2/hbar);

figure,
h1 = plot(x,real(psi));
hold on
h2 = plot(x,imag(psi));
h3 = plot(x,conj(psi).*psi);
xlim([0 a])
ylim([-2 3.5])
legend('real','imag','|\psi|^2')

for t = 0:0.01:10
%     tesc = t*1e-12;
    psi = c1*psi1.*exp(-1i*E1*t/hbar)+c2*psi2.*exp(-1i*t*E2/hbar);
    h1.YData = real(psi);
    h2.YData = imag(psi);
    h3.YData = conj(psi).*psi;
    refreshdata
    drawnow
    pause(0.01)
end