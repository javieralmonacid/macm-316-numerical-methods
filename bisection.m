function [p,summary] = bisection(f,a,b,max_iter,tol)
% BISECTION Attempts to find a root of a function using bisection's method.
%
% p = bisection(f,a,b,max_iter,tol) attempts to find a root of f over the 
% interval [a,b]. The solver stops when the the root has been found (|f(p)| 
% < eps), the relative error between iterations falls below the tolerance
% tol, or when the maximum number of iterations has been reached.
%
% [p,summary] = bisection(f,a,b,max_iter,tol) works as before but this time, the
% solver returns a table with the iteration number, a_k, b_k, and p_k
% values in the algorithm.
%
% Author: Javier Almonacid, for SFU's MACM 316 Class (Summer 2025)

iter = 0;
table = [];

% Check that neither endpoint is a root, and if f(a) and f(b) have the 
% same sign, output an error.
if (f(a) == 0)
    p = a;
    return;
elseif (f(b) == 0)
    p = b;
    return;
elseif (sign(f(a)) * sign(f(b)) > 0)
    error('Function must have opposite signs at a and b');
end

% Perform bisection
p_old = Inf;
while (iter < max_iter)
    iter = iter + 1;
    p = a + (b-a)/2;
    fp = f(p);
    table = [table; iter, a, b, p, fp];
    if (fp == 0 || abs(p-p_old)/abs(p) < tol)
        break; % p is a root
    elseif (sign(fp)*sign(f(a))>0)
        a = p; % Use midpoint as the left endpoint for the next iteration
    else
        b = p; % Use midpoint as the right endpoint for the next iteration
    end
    p_old = p;
end

% Export table
summary = array2table(table, 'VariableNames', ...
            {'Iteration k', 'ak', 'bk', 'pk', 'f(pk)'});
