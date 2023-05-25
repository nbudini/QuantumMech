clearvars
close all

% define constants (setting them to 1 is enough to capture the physics)
m = 1; % mass
r = 1; % circle radius
hbar = 1; % Planck's constant (6.626e-34 J s)
A = 1/sqrt(2*pi); % eigenstate normalization constant

% define range for polar angle and grid size
npts = 500; % number of grid points
theta = linspace(0,2*pi,npts);

% define eigenstates as a function of quantum number n and polar angle
psi_n_theta = @(n,theta) A*exp(1i*n*theta);

% define initial wavefunction Psi(theta,0)
sigma = pi/20; % gaussian width (std. dev.)
xcen = 2*pi/2; % gaussian maximum
Psi0 = @(x) 1/sigma/sqrt(2*pi)*exp(-(x-xcen).^2/2/sigma^2); % gaussian (change at your taste)
normfactor = sqrt(integral(@(theta) conj(Psi0(theta)).*Psi0(theta),0,2*pi)); % calculate norm factor
Psi0n = @(x) 1/sigma/sqrt(2*pi)/normfactor*exp(-(x-xcen).^2/2/sigma^2); % normalize wavefunction

% get normalized coefficients of the eigenstates' expansion of Psi(theta,0)
ncoeff = 20; % number of coefficients in the expansion
coeff = getCoeffs_fun(ncoeff,Psi0,psi_n_theta); % call external function to get them

% calculate expanded wavefunction at t = 0 for theta in [0,2pi]
Psi0_exp = 0; % initialize vector
for i = 1:numel(coeff)
    Psi0_exp = Psi0_exp + coeff(i)*psi_n_theta(i,theta);
end
Psi0_exp = Psi0_exp/sqrt(sum(abs(Psi0_exp).^2)); % renormalize expansion

% calculate expanded wavefunction as a function of time and plot for each
% time until tf in dt steps

% initialize and configure plots
figure(1)
subplot(311) % real part
h1 = plot(theta,real(Psi0_exp),'LineWidth',2);
xlim([0 2*pi])
ylim([-1 2])
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
ylabel('Re(\Psi)')

subplot(312) % imaginary part
h2 = plot(theta,imag(Psi0_exp),'r','LineWidth',2);
xlim([0 2*pi])
ylim([-1 2])
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
ylabel('Im(\Psi)')

subplot(313) % probability amplitude
h3 = plot(theta,conj(Psi0_exp).*Psi0_exp,'k','LineWidth',2);
xlim([0 2*pi])
ylim([-0.5 2])
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
xlabel('\theta coordinate')
ylabel('|\Psi|^2')

% initialize gif file to save animation
sflag = 0; % 1 -> save to gif, 0 -> don't save to gif
fname = 'qm.gif';
if sflag == 1 % generate gif with first frame only
    im = frame2im(getframe(gcf)); % capture the current (first) frame
    [imind,cm] = rgb2ind(im,256); % adapt rgb indices and colormap
    imwrite(imind,cm,fname,'gif','delaytime',0.1,'Loopcount',inf); % save
end

% begin loop
tf = 1;
dt = 0.01;
for t = 0:dt:tf 
    Psi_t = 0; % initialize vector
    for j = 1:numel(coeff)
        % add e^(-i E_n t/hbar) factor to each expansion term
        Psi_t = Psi_t + coeff(j)*psi_n_theta(j,theta)*exp(-1i*(j^2*hbar^2/2/m/r^2)*t/hbar);
    end

    % refresh plots after each timestep
    h1.YData = real(Psi_t);
    h2.YData = imag(Psi_t);
    h3.YData = conj(Psi_t).*Psi_t;
    refreshdata
    drawnow

    if sflag == 1 % append frame to gif file
        im = frame2im(getframe(gcf)); % capture the current frame
        [imind,cm] = rgb2ind(im,256); % adapt rgb indices and colormap
        imwrite(imind,cm,fname,'gif','delaytime',0.1,'WriteMode','append'); % save
    end

    % have a breathe (only if necessary...)
    % pause(0.05)
end