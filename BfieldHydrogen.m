% Simple script to roughly estimate the magnetic field intensity felt by
% the electron due to the proton in the hydrogen atom for different values
% of n by using Biot-Savart law with r_n and v_n deduced from Bohr model

hbar = h / (2*pi);                   % Reduced Planck constant (J·s)
e = 1.602176634e-19;                 % Elementary charge (C)
epsilon0 = 8.8541878128e-12;         % Vacuum permittivity (F/m)
mu0 = 4*pi*1e-7;                     % Vacuum permeability (H/m)
a0 = 5.29177210903e-11;              % Bohr radius (m)

rn = @(n) a0*n.^2;                   % Radius for n-th level (Bohr)
vn = @(n) e^2/4/pi/epsilon0./n/hbar; % Velocity of electron for n-th level (Bohr)

B = @(n) (mu0*e/4/pi)./rn(n).^2.*vn(n);

n = 1:10;                            % values for n levels
semilogy(n,B(n),'.','MarkerSize',25)
xlabel('n value')
ylabel('Magnetic field (T)')

% Add text labels
for i = 1:length(n)-2
    text(n(i) + 0.2, B(n(i)), sprintf('%.4f %s', B(n(i)), 'T'), ...
        'FontSize', 10, 'VerticalAlignment', 'middle');
end