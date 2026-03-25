%Plotting Parameters
n_max = 30;
nx = 10001;

%Model Parameters
alpha = 1;
gamma_1 = 1.0527;
gamma_2 = 1.0528;
gamma_3 = 1.0529;
L = 10;
a = 0.01;
u0 = 0.05;
m0 = 0.05 * 0.95;


%Pre-Allocate for Lambda and k
n_vals = linspace(0, n_max, nx);
k_vals = zeros(1,nx);
lambda_1 = zeros(1, nx);
lambda_2 = zeros(1, nx);
lambda_3 = zeros(1, nx);

%Calculate Lambda
for i = 1:length(n_vals)
    n = n_vals(i);
    k_vals(i) = n * pi / L;
    lambda_1(i) = - (n*pi/L)^2 * ((alpha*u0) - (gamma_1*m0*SensoryKernel_FT(n, L, a, nx)));
    lambda_2(i) = - (n*pi/L)^2 * ((alpha*u0) - (gamma_2*m0*SensoryKernel_FT(n, L, a, nx)));
    lambda_3(i) = - (n*pi/L)^2 * ((alpha*u0) - (gamma_3*m0*SensoryKernel_FT(n, L, a, nx)));
end

%Plot
figure;
plot(k_vals, lambda_1, 'r-', 'LineWidth', 1.5); hold on;
plot(k_vals, lambda_2, 'g-', 'LineWidth', 1.5); 
plot(k_vals, lambda_3, 'b-', 'LineWidth', 1.5); 
xline(pi / L, '-', 'Color', [0.7 0.7 0.7 0.2], 'LineWidth', 1.5)
xline(2 * pi / L, '-', 'Color', [0.7 0.7 0.7 0.2],'LineWidth', 1.5)
xline(3 * pi / L, '-', 'Color', [0.7 0.7 0.7 0.2],'LineWidth', 1.5)
xline(4 * pi / L, '-', 'Color', [0.7 0.7 0.7 0.2],'LineWidth', 1.5)
xline(5 * pi / L, '-', 'Color', [0.7 0.7 0.7 0.2],'LineWidth', 1.5)
xline(6 * pi / L, '-', 'Color', [0.7 0.7 0.7 0.2],'LineWidth', 1.5)
xline(7 * pi / L, '-', 'Color', [0.7 0.7 0.7 0.2],'LineWidth', 1.5)
xline(8 * pi / L, '-', 'Color', [0.7 0.7 0.7 0.2],'LineWidth', 1.5)
xline(9 * pi / L, '-', 'Color', [0.7 0.7 0.7 0.2],'LineWidth', 1.5)
xline(10 * pi / L, '-', 'Color', [0.7 0.7 0.7 0.2],'LineWidth', 1.5)
yline(0, '--k', 'LineWidth', 1.5);   
xlabel('Wave-Number (k_n)');
ylabel('Growth Rate (\lambda_n)');
legend(['\gamma = ', num2str(gamma_1)], ['\gamma = ', num2str(gamma_2)], ['\gamma = ', num2str(gamma_3)], 'Location', 'southwest');
title('Dispersion Relation for the Non-Local Model');
axis([0 2.5 -2e-5 2e-5]);


