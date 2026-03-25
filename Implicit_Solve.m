%Finally, we perform many time steps and store our results to graph
function [sol] = Implicit_Solve(u_0, L, N, T, M, alpha, gamma, maxit, omega, tol, a)
dx = (2*L)/N;
sol = zeros(N,M+1); %Preassign size
sol(:,1) = u_0; %First column is initial conditions

for t = 2:M+1
    sol(:,t) = Picard_Timestep(sol(:,t-1), L, N, T, M, alpha, gamma, maxit, omega, tol, a);
    if 100 ~= 0 && mod(t, 100) == 0
        disp(t)
    end
    mass = sum(sol(:,t))*dx;
    if abs(mass - sum(u_0)*dx) > 1e-6
        disp("Mass Consistency Failure");
    end
end
end