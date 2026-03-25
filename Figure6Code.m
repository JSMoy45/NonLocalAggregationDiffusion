%Parameters
nx = 1000;              
L = 10;
a = 0.1;
x = linspace(-L, L, nx+1);
x(end) = [];       
dx = 2 * L / nx;
alpha = 1;
%gamma = 1.01;

%Calculate Moments
w = zeros(1, nx);
for i = 1:length(x)
    x_val = x(i);
    w(i) = SensoryKernel(x_val, a);
end
M0 = trapz(x, w);
M2 = trapz(x, x.^2 .* w);
%Predicted = cos(imag(sqrt(2*(alpha - gamma*M0) / (gamma * M2))) * x);
%Pred_Mass = sum(Predicted) * dx;
%Predicted = (Predicted / Pred_Mass);
%Predicted = Predicted(:);


%Preallocate
u = zeros(nx, 1);

%Calculate W_hat
W = zeros(nx, 1);
for i = 1:nx
    W(i) = SensoryKernel(x(i), a);
end
W = ifftshift(W);
W_hat = fft(W);

%Gammas
gamma_nonloc = alpha / SensoryKernel_FT(1, L, a, 1000000);
gamma_loc = alpha / (M0 - ( (M2/2) * (pi /L)^2 ));

function [F, grad] = Energy(u, alpha, gamma, W_hat, dx)
    conv_u = dx * real(ifft(W_hat .* fft(u)));
    integrand = (alpha * 0.5)*u.^2 - (gamma/2)*u.*conv_u;
    F = dx * sum(integrand);
    grad = dx*(alpha*u - gamma*conv_u);
end

function [F, grad] = LocalEnergy(u, alpha, gamma, M0, M2, dx)
    ux = (circshift(u,-1) - circshift(u,1)) / (2*dx);
    uxx = (circshift(u,-1) - 2*u + circshift(u,1)) / dx^2;
    integrand = ((alpha - M0*gamma) * 0.5)*u.^2 + (gamma * M2 * 0.25)*ux.^2;
    F = dx * sum(integrand);
    grad = dx*((alpha-M0*gamma)*u - gamma*0.5*M2.*uxx);
end

%Initial Condition
u0 = 1 + ( 1*exp(-(x).^2 / (2^2)) ); 
Initial_Mass = sum(u0) * dx;
u0 = (u0 / Initial_Mass);
u0 = u0(:);  

%Cosine Perturbation
%A = 0.01;         
%baseMean = 1 / (2*L);             
%u0 = baseMean + A * cos(pi*x/L);
%u0 = u0(:);  


%options = optimoptions('fmincon',...
%    'Algorithm','interior-point',...
%    'SpecifyObjectiveGradient',true,...
%    'Display','iter');

options = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'SpecifyObjectiveGradient',true, ...
    'OptimalityTolerance',1e-12, ...
    'StepTolerance',1e-14, ...
    'ConstraintTolerance',1e-12, ...
    'MaxIterations',1000, ...
    'MaxFunctionEvaluations',1e8, ...
    'Display','iter');

[u_min, E_min] = fmincon(@(u) Energy(u,alpha,gamma_nonloc,W_hat,dx), ...
                u0, ...
                [],[], dx*ones(1,nx), 1, ...   
                zeros(nx, 1), [], ...         
                [], ...            
                options);

[u_min_Local, E_min_Local] = fmincon(@(u) LocalEnergy(u,alpha,gamma_loc,M0,M2,dx), ...
                u0, ...
                [],[], dx*ones(1,nx), 1, ...   
                zeros(nx, 1), [], ...         
                [], ...            
                options);



figure;
plot(x, u0, 'k--', 'LineWidth', 1.2); hold on;
plot(x, u_min, 'b-', 'LineWidth', 2);
plot(x, u_min_Local, 'g--', 'LineWidth', 2);
xlabel('x'); ylabel('u(x)');
legend('Initial State', 'Non-Local Energy Minimiser', 'Local Energy Minimiser', 'Location','Best');
title('Numerically Minimised Energy');
grid on;

mass = sum(u_min) * dx;
fprintf('mass (should equal M) = %.15g\n', mass);