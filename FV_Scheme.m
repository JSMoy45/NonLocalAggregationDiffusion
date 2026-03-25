%First is to define our implicit map for u^n+1
function [v_1] = FV_Scheme(v_guess, v_0, N, dx, dt, alpha, gamma, a)

x = zeros(N,1); 
for i = 0:N-1
    domain_start = -(N*dx)/2;
    x(i+1) = domain_start + (i+0.5)*dx; %Spatial Cell Centers
end

L_2 = N*dx;

%Define v** as in paper
v_conv = (v_guess + v_0) / 2;

%Define convolution kernal W
sigma_w = a; %Scales the Sensing Range
W = @(r) (1/(sqrt(2*pi)*sigma_w)) * exp(-r^2 / (2*(sigma_w^2))); %Perceptive Function

%Calculate (W*v)_i
K = zeros(N,N);
for i = 1:N
    for j = 1:N
        r = x(i) - x(j);
        r = r - L_2*round(r/L_2);
        K(i,j) = W(r) * dx;  
    end
end
W_star_v = K*v_conv;

%Then calculate alpha*v_i - gamma*(W*v)_i
nu = (alpha.*v_guess);
xi = -(gamma.*W_star_v);


%Next calculate the derivative of xi, aka g at the interfaces. To comply
%with boundary conditions, we set g = 0 manually at the boundaries
% xi and nu are column vectors of length N (as in your code)


xi_right = xi;              
xi_left  = circshift(xi, 1);
g = (-1/dx) .* ( xi_right - xi_left );

nu_right = nu;
nu_left  = circshift(nu, 1);
h = (-1/dx) .* ( nu_right - nu_left );  

%Now we use g to calculate the flux at the interfaces as given in the paper
F = zeros(N,1);
Fg = zeros(N, 1);
Fh = zeros(N, 1);
for q = 1:N
    left  = mod(q-2,N) + 1;    
    right = mod(q-1,N) + 1;

    if g(q) >= 0
        Fg(q) = ( v_guess(left) * (1 - v_guess(left)) *  g(q) );
    elseif g(q) < 0
        Fg(q) = ( v_guess(right) * (1 - v_guess(right)) * g(q) );
    else
        disp("ERROR1")
    end

    if h(q) >= 0
        Fh(q) = ( v_guess(left) * h(q) );
    elseif h(q) < 0
        Fh(q) = ( v_guess(right) * h(q) );
    else
        disp("ERROR2")
    end

    F(q) = Fh(q) + Fg(q);

end

%Finally, use these fluxes to calculate our output
F_right = circshift(F,-1);    
F_diff  = (dt/dx).*(F_right - F); 
v_1 = v_0 - F_diff;
end