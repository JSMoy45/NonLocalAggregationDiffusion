function [w_hat] = SensoryKernel_FT(n, L, a, dx)
    k = n * pi / L;
    x = linspace(-L, L, dx);
    W = SensoryKernel(x, a);
    integrand = W .* exp(-1i * k * x);
    w_hat = trapz(x, integrand);
    w_hat = real(w_hat);
end