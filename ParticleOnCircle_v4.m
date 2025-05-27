clearvars
close all

% define constants (setting them to 1 is enough to capture the physics)
m = 1; % mass
hbar = 1; % Planck's constant (6.626e-34 J s)
A = 1/sqrt(2*pi); % eigenstate normalization constant

% define range for polar angle and grid size
npts = 1000; % number of grid points
theta = linspace(0,2*pi,npts);

% define eigenstates as a function of quantum number n and polar angle
psi_n_theta = @(n,theta) A*exp(1i*n*theta);

% define initial wavefunction Psi(theta,0) as a gaussian, but change at
% your taste
sigma = pi/200; % gaussian width (std. dev.)
xcen = pi; % gaussian maximum
Psi0 = @(x) 1/sigma/sqrt(2*pi)*exp(-(x-xcen).^2/2/sigma^2); % gaussian (change at your taste)
normfactor = sqrt(integral(@(theta) conj(Psi0(theta)).*Psi0(theta),0,2*pi)); % calculate norm factor
Psi0n = @(x) 1/sigma/sqrt(2*pi)/normfactor*exp(-(x-xcen).^2/2/sigma^2); % normalize wavefunction

% get normalized coefficients of the eigenstates' expansion of Psi(theta,0)
ncoeff = 20; % number of coefficients in the expansion
coeff = getCoeffs_fun_(ncoeff,Psi0,psi_n_theta); % call external function to get them

% calculate expanded wavefunction at t = 0 for theta in [0,2pi]
Psi0_exp = 0; % initialize vector
for i = 1:numel(coeff)
    Psi0_exp = Psi0_exp + coeff(i)*psi_n_theta(i,theta);
end
Psi0_exp = Psi0_exp/sqrt(sum(abs(Psi0_exp).^2)); % renormalize expansion

% calculate expanded wavefunction as a function of time and plot for each
% time until tf in dt steps

% initialize and configure plots
f1 = figure(1);
subplot(221) % real part
h1 = plot(theta,real(Psi0_exp),'LineWidth',2);
xlim([0 2*pi])
ylim([-0.15 0.17])
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
ylabel('Re(\Psi)')

subplot(222) % imaginary part
h2 = plot(theta,imag(Psi0_exp),'r','LineWidth',2);
xlim([0 2*pi])
ylim([-0.15 0.17])
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
ylabel('Im(\Psi)')

subplot(223) % probability amplitude
h3 = plot(theta,conj(Psi0_exp).*Psi0_exp,'k','LineWidth',2);
xlim([0 2*pi])
ylim([0 0.025])
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
xlabel('\theta coordinate')
ylabel('|\Psi|^2')

r = ones(size(theta)); % radial coordinate 

subplot(224) % probability amplitude in polar coordinates
h4 = polarscatter(theta,r,50,conj(Psi0_exp).*Psi0_exp,'filled');
rlim([0 1.1])
colormap(turbo);
clim([0 max(conj(Psi0_exp).*Psi0_exp)]);

% Get the underlying polar axes
pax = gca; % or explicitly with polaraxes if needed

% Hide r-axis grid lines and ticks
pax.RColor = 'none';          % Removes r-axis lines and labels
pax.ThetaColor = 'k';         % Keep theta labels visible (set color to black)
pax.GridColor = 'none';       % Remove grid lines if you want

% Show only selected angle labels (0, pi/2, pi, 3pi/2, 2pi)
pax.ThetaTick = [0 90 180 270]; 
pax.ThetaTickLabel = {'0','\pi/2','\pi','3\pi/2'};

% begin loop
tf = 12.565;    % end time (system returns to initial state after this time!)
dt = 0.005;     % timestep
k = 0;          % index for image filename only
for t = 0:dt:tf 
    Psi_t = 0; % initialize vector
    for j = 1:numel(coeff)
        % add e^(-i E_n t/hbar) factor to each expansion term
        Psi_t = Psi_t + coeff(j)*psi_n_theta(j,theta).*exp(-1i*(j^2*hbar^2/2/m/r(j)^2)*t/hbar);
    end
    Psi_t = Psi_t/sqrt(sum(abs(Psi_t).^2));

    % refresh plots after each timestep
    h1.YData = real(Psi_t);
    h2.YData = imag(Psi_t);
    h3.YData = conj(Psi_t).*Psi_t;
    h4.CData = conj(Psi_t).*Psi_t;
    
    txt.String = sprintf('%.3f', t);

    refreshdata
    drawnow

    % have a breathe (only if necessary...)
    % pause(0.01)
    
    % Uncomment to save frames as png images
    % fname = sprintf('images/frame_%04d.png', k); % e.g., frame_001.png
    % exportgraphics(f1, fname, 'Resolution', 300); % 300 DPI PNG
    % k = k+1;

    % To create mp4 video from these png images use this command afterwards
    % in a Linux terminal:
    % $ ffmpeg -framerate 20 -i frame_%04d.png -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -pix_fmt yuv420p output.mp4
    % 
    % (Install ffmpeg with: $ sudo apt install ffmpeg)
end