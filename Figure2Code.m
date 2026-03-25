%Plotting Parameters
n_max = 1000;
nx = 10001;

%Model Parameters
L = 10;
a = 0.01;


%Pre-Allocate 
n_vals = linspace(0, n_max, nx);
k_vals = zeros(1,nx);
w = zeros(1, n_max);
W_hat = zeros(1, nx);
W_hat_approx = zeros(1, nx);

%Calculate Moments
x = linspace(-L, L, nx);
for i = 1:length(x)
    x_val = x(i);
    w(i) = SensoryKernel(x_val, a);
end
M0 = trapz(x, w);
M2 = trapz(x, x.^2 .* w);

%Calculate W_hat and W_hat_approx
for i = 1:length(n_vals)
    n = n_vals(i);
    k = n* pi / L;
    k_vals(i) = k;
    W_hat(i) = SensoryKernel_FT(n, L, a, nx);
    W_hat_approx(i) = M0 - (0.5 * M2 * k^2);
end

figure;
plot(k_vals, W_hat, 'b-', 'LineWidth', 1.5); hold on;
plot(k_vals, W_hat_approx, 'r--', 'LineWidth', 1.5); 
yline(0, '--k', 'LineWidth', 1.5);   
xlabel('Wave-Number (k_n)');
ylabel('$\widehat{W}\left(k\right)$ vs $M_0 - \frac{M_2}{2}k^2$', Interpreter='latex');
legend('$\widehat{W}\left(k\right)$ ',  '$M_0 - \frac{M_2}{2}k^2$', 'Location', 'northeast', interpreter = 'latex');
title('Fourier Tranform of W vs its Apporximation');
axis([1 250 -2 2]);
