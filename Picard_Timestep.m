%Next, we use a Picard iteration to solve for u^n+1 given u^n
function [u_new] = Picard_Timestep(u, L, N, T, M, alpha, gamma, maxit, omega, tol, a)
  dx = (2*L) / N;
  dt = T / M;
  u_k = u;

  for k = 1:maxit
    u_tilde = FV_Scheme(u_k, u, N, dx, dt, alpha, gamma, a);   % one Picard map
    u_new_candidate = (1-omega)*u_k + omega*u_tilde;     % under relax

    if any(u_new_candidate < 0) %Positivity
        disp("positivity forcibly enforced")
        u_new_candidate(u_new_candidate < 1e-12) = 1e-12;
    end

    if norm(u_new_candidate - u_k, inf) < tol %Check Convergence
      u_k = u_new_candidate; break
    end

    u_k = u_new_candidate; %Update

    if k == maxit
        disp("Max Iterations Reached")
    end
  end
  u_new = u_k;
end