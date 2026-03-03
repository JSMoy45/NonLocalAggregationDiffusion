function [W] = SensoryKernel(x, a)
    Coef = (a * sqrt(2*pi))^(-1);
    Expo = exp( - (x.^2) ./ ( 2* (a^2)) );
    W = Coef * Expo;
end