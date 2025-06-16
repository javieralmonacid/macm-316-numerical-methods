function [p,summary] = newton_raphson(f,df,p0,max_iter,tol)
% NEWTON_RAPHSON Attempts to find a root using Newton-Rhapson's method
%
% [p, summary] = newton_raphson(f,df,p0,max_iter,tol) Attempts to find a 
% root of the equation f(x) = 0 using Newton-Raphson's method with initial
% solution p0. The solver stops if |p_n - p_{n-1}|/|p_n| < tol or the
% maximum number of iterations max_iter has been achieved. Here, f and df
% are function handles for f and its derivative, respectively.
%
% Author: Javier Almonacid, for SFU's MACM 316 Class (Summer 2025)

iter = 0;
table = [iter, p0, f(p0)];
while (iter < max_iter)
    iter = iter + 1;
    p = p0 - f(p0)/df(p0);
    fp = f(p);
    table = [table; iter, p, fp];
    if (fp == 0 || abs(p-p0)/abs(p) < tol)
        break;
    end
    p0 = p;
end

% Export table
summary = array2table(table, 'VariableNames', ...
            {'Iteration k', 'pk', 'f(pk)'});
