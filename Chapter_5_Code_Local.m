%% 1D FOURIER-SPECTRAL SOLVER FOR LOCAL-APPROXIMATION MODEL
% Solves
%
%   u_t = d/dx( alpha*u*u_x - gamma*u*(1-u)*d/dx(M0*u + M2*u_xx) )
%
% on a periodic domain [-L, L).
%
% Uses:
%   - Fourier spectral derivatives
%   - semi-implicit IMEX Euler time stepping
%   - mass tracking
%   - fixed-axis live plotting
%
% Two initial condition options:
%   1) 'cosine'  : first cosine mode
%   2) 'gaussian': small Gaussian perturbation
%
% Both are shifted to have total mass = 1.

clear; close all; clc;

%% ===================== PARAMETERS =====================

% Domain
L  = 10;            % domain is [-L, L)
Nx = 2048;           % number of grid points (even recommended)

% Time stepping
dt        = 0.01;   % time step
tFinal    = 500;    % final time
tolSteady = 1e-8;   % stop if max |u_t| falls below this
plotEvery = 200;    % live plot every this many steps

% Model parameters
a     = 0.1;       % parameter used inside SensoryKernel
alpha = 1.0;
gamma = 1.2;



% Initial condition settings
targetMass = 1.0;       % total mass
icType     = 'cosine';  % 'cosine' or 'gaussian'
amp        = 0.01;      % perturbation amplitude

% Gaussian IC settings (used only if icType = 'gaussian')
gaussCenter = 0.0;
gaussWidth  = 1.0;

% Optional mobility regularization:
% set to 0 for the exact PDE, or small positive number if you want regularization
mobFloor = 0.0;

% Number of quadrature points for M0, M2 computation
Nmom = 2000000;

%% ===================== GRID =====================
tic

dx = (2*L) / Nx;
x  = (-L:dx:L-dx).';     % periodic grid, length Nx

% Fourier wavenumbers for a 2L-periodic domain
k  = (pi/L) * [0:(Nx/2-1), -Nx/2:-1].';
k2 = k.^2;
k4 = k.^4;

% 2/3 dealiasing mask for nonlinear products
kcut    = (2/3) * max(abs(k));
dealias = abs(k) <= kcut;

%% ===================== MOMENTS M0, M2 =====================

% Compute M0 and M2 exactly in the same style as your earlier code
x_moments = linspace(-L, L, Nmom);
w = arrayfun(@(xx) SensoryKernel(xx, a), x_moments);

M0 = trapz(x_moments, w);
M2 = 0.5*trapz(x_moments, x_moments.^2 .* w);

fprintf('Computed moments:\n');
fprintf('  M0 = %.12e\n', M0);
fprintf('  M2 = %.12e\n', M2);

%% ===================== INITIAL CONDITION =====================

baseMean = targetMass / (2*L);   % mean needed for total mass = 1

switch lower(icType)
    case 'cosine'
        % First Fourier cosine mode on [-L, L)
        u = baseMean + amp * cos(pi*x/L);

    case 'gaussian'
        % Small Gaussian perturbation added to the target mean
        u = baseMean + amp * exp(-((x - gaussCenter).^2) / (2*gaussWidth^2));

    otherwise
        error('Unknown icType. Use ''cosine'' or ''gaussian''.');
end

% Shift mean so total mass is exactly targetMass, while preserving shape
u = u + (targetMass/(2*L) - mean(u));

% Initial diagnostics
fprintf('\nInitial condition diagnostics:\n');
fprintf('  mass     = %.12f\n', dx * sum(u));
fprintf('  mean     = %.12f\n', mean(u));
fprintf('  min(u0)  = %.12f\n', min(u));
fprintf('  max(u0)  = %.12f\n', max(u));

%% ===================== STORAGE =====================

NtMax = ceil(tFinal / dt);

tHist    = zeros(NtMax+1, 1);
massHist = zeros(NtMax+1, 1);
uHist    = zeros(NtMax+1, Nx);

tHist(1)    = 0;
massHist(1) = dx * sum(u);
uHist(1,:)  = u.';

%% ===================== LIVE PLOT SETUP =====================

% Fixed y-axis limits so the plot does not rescale every frame
yPad = max(0.05, 0.1 * max(abs(u - mean(u))));
yMin = min(u) - yPad;
yMax = max(u) + yPad;

figure(1);
plotHandle = plot(x, u, 'LineWidth', 1.5);
xlabel('x');
ylabel('u(x,t)');
title(sprintf('t = %.3f', 0));
grid on;
xlim([x(1), x(end)]);
ylim([yMin, yMax]);
drawnow;

%% ===================== TIME STEPPING =====================

t = 0;

for n = 1:NtMax

    % ---------- Current Fourier transform ----------
    u_hat = fft(u);

    % ---------- Spectral derivatives ----------
    ux   = real(ifft(1i * k .* u_hat));
    uxx  = real(ifft(-k2 .* u_hat));
    uxxx = real(ifft(-1i * k.^3 .* u_hat));

    % ---------- Local chemical potential ----------
    % mu = M0*u + M2*u_xx
    % mu_x = M0*u_x + M2*u_xxx
    mux = M0 * ux + M2 * uxxx;

    % ---------- Mobility ----------
    mob = u .* (1 - u) + mobFloor;

    % ---------- Full nonlinear flux ----------
    flux = alpha * u .* ux - gamma * mob .* mux;

    % Dealias before differentiating the nonlinear flux
    flux_hat = fft(flux);
    flux_hat(~dealias) = 0;

    % Full RHS: u_t = d/dx(flux)
    rhsFull = real(ifft(1i * k .* flux_hat));

    % ---------- Semi-implicit linearization ----------
    % Linearize around the current mean ubar:
    %
    %   u_t ≈ alpha*ubar*u_xx - gamma*mbar*(M0*u_xx + M2*u_xxxx)
    %
    ubar = mean(u);
    mbar = ubar * (1 - ubar);

    Lhat = (-alpha * ubar + gamma * mbar * M0) .* k2 ...
           - gamma * mbar * M2 .* k4;

    rhsLin = real(ifft(Lhat .* u_hat));
    rhsNL  = rhsFull - rhsLin;

    % ---------- IMEX Euler update ----------
    u_hat_new = (u_hat + dt * fft(rhsNL)) ./ (1 - dt * Lhat);
    u = real(ifft(u_hat_new));

    % ---------- Advance time ----------
    t = t + dt;

    % ---------- Store ----------
    tHist(n+1)    = t;
    massHist(n+1) = dx * sum(u);
    uHist(n+1,:)  = u.';

    % ---------- Steady-state check on updated solution ----------
    u_hat = fft(u);
    ux    = real(ifft(1i * k .* u_hat));
    uxxx  = real(ifft(-1i * k.^3 .* u_hat));
    mux   = M0 * ux + M2 * uxxx;
    mob   = u .* (1 - u) + mobFloor;

    flux = alpha * u .* ux - gamma * mob .* mux;
    flux_hat = fft(flux);
    flux_hat(~dealias) = 0;
    rhsCheck = real(ifft(1i * k .* flux_hat));

    if max(abs(rhsCheck)) < tolSteady
        fprintf('\nSteady state reached at t = %.6f\n', t);
        tHist    = tHist(1:n+1);
        massHist = massHist(1:n+1);
        uHist    = uHist(1:n+1,:);
        break;
    end

    % ---------- Live plot ----------
    if mod(n, plotEvery) == 0
        set(plotHandle, 'YData', u);
        title(sprintf('t = %.3f', t));
        drawnow;
    end
end

%% ===================== FINAL DIAGNOSTICS =====================

fprintf('\nFinal diagnostics:\n');
fprintf('  final time = %.12f\n', tHist(end));
fprintf('  final mass = %.12f\n', massHist(end));
fprintf('  min(u)     = %.12f\n', min(uHist(end,:)));
fprintf('  max(u)     = %.12f\n', max(uHist(end,:)));
toc

%% ===================== FINAL PLOTS =====================

%figure;
%surf(x, tHist, uHist, 'EdgeColor', 'none');
%xlabel('x');
%ylabel('t');
%zlabel('u(x,t)');
%title('Solution u(x,t)');
%colorbar;
%view(45,30);

figure;
plot(tHist, massHist, 'LineWidth', 2);
xlabel('t');
ylabel('Mass');
title('Mass conservation check');
grid on;

figure;
plot(x, uHist(end,:), 'LineWidth', 2);
xlabel('x');
ylabel('u(x,t_{final})');
title('Final profile');
grid on;

%%

Solution_Energy = zeros(1, length(tHist));

for i = 1:length(tHist)
    u = uHist(i, :).';

    % Compute spatial derivative
    ux = gradient(u, dx);

    % Energy density:
    % (alpha - gamma*M0)/2 * u^2 + (gamma*M2/4) * (u_x)^2
    integrand = ((alpha - gamma*M0)/2) * u.^2 ...
              + (gamma*M2/4) * ux.^2;

    % Integrate
    Solution_Energy(i) = dx * sum(integrand);
end

figure;
plot(tHist, Solution_Energy, 'b-', 'LineWidth', 2)
title('Energy Dissipation for the Local Model')
xlabel('t')
ylabel('Energy')
grid on