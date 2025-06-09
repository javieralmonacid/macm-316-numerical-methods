function [x, iter, xall] = jacobi(A, b, x0, tol, max_iter)
% JACOBI Solves the system Ax=b using Jacobi iterative method
%
% x = jacobi(A, b, x0, tol, max_iter) computes the solution to the
% system Ax=b using x0 as a starting vector. The convergence criterion is
% set to:
%
%          || x^(k) - x^(k-1) ||_{inf}
%         ------------------------------ < tol.
%                  || x^(k) ||_{inf}
%
% The solver stops iterating once the maximum number of iterations has been
% reached, in case convergence was not detected earlier.
%
% [x, iter] = jacobi(A, b, x, tol, max_iter) works as above but the solver
% now also returns the number of iterations that the solver took
%
% [x, iter, xall] = jacobi(A, b, x, tol, max_iter) works as above but now
% the solver also returns *all* Jacobi iterations.
%
% Author: Javier Almonacid, for SFU's MACM 316 Class (Summer 2025)
% Based on Algorithm 7.1, Burden & Faires (10th Edition).

n = length(b);
iter = 0;
x = zeros(n,1);
xall = x0;

while iter < max_iter
    iter = iter + 1;
    % Main iterates
    for i = 1:n
        x(i) = (- A(i,[1:i-1 i+1:n]) * x0([1:i-1 i+1:n]) + b(i) ) / A(i,i);
    end
    % This is to collect all iterations, if requested
    if nargout == 3
        xall = [xall x];
    end
    % Convergence criterion
    if norm(x - x0, 'inf') / norm(x0) < tol
        break;
    end
    % Store new iterate
    x0 = x;
end

% Display a warning if convergence was not reached
if iter == max_iter
    warning('Jacobi did not converge within the maximum number of iterations (%d).', max_iter);
end