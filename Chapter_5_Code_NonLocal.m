%Chapter 5 Code

%%
%Clear Workspace

clear; clc;

%%
%Parameter Definition

L = 10; %Spatial Parameters
N = 200;
x = zeros(N,1); 
dx = (2*L)/N;
for i = 0:N-1
    x(i+1) = (-L) + (i+0.5)*dx; 
end

T = 500; %Temporal Parameters
M = 50000;
t = linspace(0, T, M+1);

alpha = 1; %Model Parameters
a = 0.1;
x_moments = linspace(-L, L, 1000000);
w = zeros(1, 1000000);
for i = 1:length(x_moments)
    x_val = x_moments(i);
    w(i) = SensoryKernel(x_val, a);
end
M0 = trapz(x_moments, w);
M2 = trapz(x_moments, x_moments.^2 .* w);
Gamma = 1;
eps = 0.5;
u0 = 0.05;
m0 = u0*(1-u0);
m1 = 1 - (2*u0);
n_c = 1;
k_c = n_c * pi / L;
w_hat_1 = SensoryKernel_FT(n_c, L, a, 10000000);
gamma_c = (alpha*u0) / (m0*w_hat_1);

%Solver Parameters
nx = 10000; 
maxit = 100;
omega = 0.5;
tol = 1e-10;

%%
%Set Initial Conditions

%Gaussian Perturbation
%sigma_0 = 0.5;
%mu_0 = 0;
%u_0 = 1 + exp( -(x - mu_0).^2 ./ (2*sigma_0^2) ); %Small Gaussian Perturbation
%u_0 = u_0 / (sum(u_0) * dx); %Ensure our starting mass is 1

%Cosine Perturbation
A = 0.01;         
baseMean = 1 / (2*L);             
u_0 = baseMean + A * cos(pi*x/L);

%Homogenous State
%u_0 = ones(N, 1);
%u_0 = u_0 / (sum(u_0) * dx); %Ensure our starting mass is 1

%%
%Non-Local Solve and Plot

%gamma = 1.1;
%Solution = Implicit_Solve(u_0, L, N, T, M, alpha, gamma, maxit, omega, tol, a);

%figure;
%surf(x, t.', Solution.')   
%shading interp;
%xlabel('x');
%ylabel('time');
%zlabel('u(x,t)');
%title('Solution u(x,t)');
%colorbar;
%view(45,30);


%% 
% Nonlocal model just below and just above the critical gamma

gamma_below = 0.8;

% Run the solver below criticality
tic
Solution_below = Implicit_Solve(u_0, L, N, T, M, alpha, gamma_below, maxit, omega, tol, a);
toc
%1823.707316
%2020.790399

%%

gamma_above = 1.2;
% Run the solver above criticality
tic
Solution_above = Implicit_Solve(u_0, L, N, T, M, alpha, gamma_above, maxit, omega, tol, a);
toc
%1765.804893


%%
%Plot 1

figure;

surf(x, t.', Solution_below.');
shading interp;
xlabel('x');
ylabel('time');
zlabel('u(x,t)');
title('Solution u(x,t), \gamma = 0.8');
colorbar;
view(45,30);

%%
%Plot 2

figure;

surf(x, t.', Solution_above.');
shading interp;
xlabel('x');
ylabel('time');
zlabel('u(x,t)');
title('Solution u(x,t), \gamma = 1.2');
colorbar;
view(45,30);

%%

%Tracking Energy
W = zeros(N, 1);
for i = 1:N
    W(i) = SensoryKernel(x(i), a);
end
W = ifftshift(W);
W_hat = fft(W);

Solution_Energy = zeros(1, M+1);
for i = 1:M+1
    conv_u = dx * real(ifft(W_hat .* fft(Solution_above(:,i))));
    integrand = (alpha * 0.5).*Solution_above(:,i).^2 - (gamma_above/2).*Solution_above(:,i).*conv_u;
    F = dx * sum(integrand);
    Solution_Energy(i) = F;
end

figure;
plot(t, Solution_Energy, 'b-', 'LineWidth', 2)
title('Energy Dissipation for the Non-Local Model')
xlabel('t')
ylabel('Energy')
grid on