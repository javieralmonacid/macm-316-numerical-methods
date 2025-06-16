function [p,summary] = secant(f, p0, p1, max_iter, tol)
% SECANT Attempts to find a root using the secant method
%
% [p, summary] = secant(f, p0, p1, max_iter, tol) Attempts to find a 
% root of the equation f(x) = 0 using the secant method with initial
% solution p0. The solver stops if |p_n - p_{n-1}|/|p_n| < tol or the
% maximum number of iterations max_iter has been achieved. Here, f and df
% are function handles for f and its derivative, respectively.
%
% Author: Javier Almonacid, for SFU's MACM 316 Class (Summer 2025)

iter = 0;
table = [iter, p0, f(p0); iter, p1, f(p1)];

% Check that neither neither of the first two points are a root
if (f(p0) == 0)
    p = p0;
    return;
elseif (f(p1) == 0)
    p = p1;
    return;
end

while (iter < max_iter)
    iter = iter + 1;
    slope = (f(p1)-f(p0))/(p1-p0);
    p = p0 - f(p0)/slope;
    fp = f(p);
    table = [table; iter, p, fp];
    if (fp == 0 || abs(p-p0)/abs(p) < tol)
        break;
    end
    p0 = p1;
    p1 = p;
end

% Export table
summary = array2table(table, 'VariableNames', ...
            {'Iteration k', 'pk', 'f(pk)'});
