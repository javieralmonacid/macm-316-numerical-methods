% MACM 316: Numerical Analysis I
% Demo: Three different problems to illustrate the performance of
% algorithms to solve nonlinear equations.
clear; close all; clc;

%% Example 1

f = @(x) exp(x)-5;
df = @(x) exp(x);
p_exact = log(5);
max_iter = 100;
tol = 1e-16;

% Case 1: p0 is close to p_exact
[~, tbis] = bisection(f, 0, 2, max_iter, tol);
[~, tnew] = newton_raphson(f, df, 1.5, max_iter, tol);
[~, tsec] = secant(f, 1, 2, max_iter, tol);

figure(1)
fig = gcf;
fig.Position(3:4) = [586 438];
fs = 14; % Font Size
semilogy(tbis.("Iteration k"), abs(tbis.pk - p_exact), '-o', ...
         tnew.("Iteration k"), abs(tnew.pk - p_exact), '-s', ...
         tsec.("Iteration k"), abs(tsec.pk - p_exact), '-d', ...
         'LineWidth', 1.5)
legend({'Bisection (a = 0, b = 2)', ...
        'Newton (p_0 = 1.5)', ...
        'Secant (p_0 = 1, p_1 = 2)'}, 'FontSize', fs)
xlabel('Iteration', 'FontSize', fs)
ylabel('Absolute Error', 'FontSize', fs)
ax = gca;
ax.FontSize = fs;
grid on

% Case 2: p0 is far from p_exact
[pbis, tbis] = bisection(f, 0, 2, 100, 1e-16);
[pnew, tnew] = newton_raphson(f, df, -1, 100, 1e-16);
[psec, tsec] = secant(f, -2, 0, 100, 1e-16);

figure(2)
fig = gcf;
fig.Position(3:4) = [586 438];
fs = 14;
semilogy(tbis.("Iteration k"), abs(tbis.pk - p_exact), '-o', ...
         tnew.("Iteration k"), abs(tnew.pk - p_exact), '-s', ...
         tsec.("Iteration k"), abs(tsec.pk - p_exact), '-d', ...
         'LineWidth', 1.5)
legend({'Bisection (a = 0, b = 2)', ...
        'Newton (p_0 = -1)', ...
        'Secant (p_0 = -2, p_1 = 0)'}, 'FontSize', fs)
xlabel('Iteration', 'FontSize', fs)
ylabel('Absolute Error', 'FontSize', fs)
ax = gca;
ax.FontSize = fs;
grid on

%% Example 2

pol = [1 0 -2 2]; % Coefficients of the polynomial
f = @(x) polyval(pol, x);
df = @(x) polyval([3 0 -2], x);
all_roots = roots(pol);
p_exact = all_roots(1);
max_iter = 100;
tol = 1e-16;

% Case 1: p0 is close to p_exact
[~, tbis] = bisection(f, -2, 2, max_iter, tol);
[~, tnew] = newton_raphson(f, df, -2, max_iter, tol);
[~, tsec] = secant(f, 0, 5, max_iter, tol);

figure(3)
fig = gcf;
fig.Position(3:4) = [586 438];
fs = 14; % Font Size
semilogy(tbis.("Iteration k"), abs(tbis.pk - p_exact), '-o', ...
         tnew.("Iteration k"), abs(tnew.pk - p_exact), '-s', ...
         tsec.("Iteration k"), abs(tsec.pk - p_exact), '-d', ...
         'LineWidth', 1.5)
legend({'Bisection (a = -2, b = 2)', ...
        'Newton (p_0 = -2)', ...
        'Secant (p_0 = 0, p_1 = 5)'}, 'FontSize', fs)
xlabel('Iteration', 'FontSize', fs)
ylabel('Absolute Error', 'FontSize', fs)
ax = gca;
ax.FontSize = fs;
grid on

% Case 2: p0 is far from p_exact and Newton diverges
[~, tbis] = bisection(f, -2, 2, max_iter, tol);
[~, tnew] = newton_raphson(f, df, 0, max_iter, tol);
[~, tsec] = secant(f, 0, 5, max_iter, tol);

figure(4)
fig = gcf;
fig.Position(3:4) = [586 438];
fs = 14; % Font Size
semilogy(tbis.("Iteration k"), abs(tbis.pk - p_exact), '-o', ...
         tnew.("Iteration k"), abs(tnew.pk - p_exact), '-s', ...
         tsec.("Iteration k"), abs(tsec.pk - p_exact), '-d', ...
         'LineWidth', 1.5)
legend({'Bisection (a = -2, b = 2)', ...
        'Newton (p_0 = 0)', ...
        'Secant (p_0 = 0, p_1 = 5)'}, 'FontSize', fs)
xlabel('Iteration', 'FontSize', fs)
ylabel('Absolute Error', 'FontSize', fs)
ax = gca;
ax.FontSize = fs;
grid on

%% Example 3

f = @(x) x./(1+x.^2);
df = @(x)  (1-x.^2)./(x.^2+1).^2;
p_exact = 0;
max_iter = 100;
tol = 1e-16;

% Case 1: p0 is close to p_exact
[~, tbis] = bisection(f, -1, 2, max_iter, tol);
[~, tnew] = newton_raphson(f, df, 0.56, max_iter, tol);
[~, tsec] = secant(f, 0.5, 1.5, max_iter, tol);

figure(5)
fig = gcf;
fig.Position(3:4) = [586 438];
fs = 14; % Font Size
semilogy(tbis.("Iteration k"), abs(tbis.pk - p_exact), '-o', ...
         tnew.("Iteration k"), abs(tnew.pk - p_exact), '-s', ...
         tsec.("Iteration k"), abs(tsec.pk - p_exact), '-d', ...
         'LineWidth', 1.5)
legend({'Bisection (a = -1, b = 2)', ...
        'Newton (p_0 = 0.56)', ...
        'Secant (p_0 = 0.5, p_1 = 1.5)'}, 'FontSize', fs)
xlabel('Iteration', 'FontSize', fs)
ylabel('Absolute Error', 'FontSize', fs)
ax = gca;
ax.FontSize = fs;
grid on

% Case 2: p0 is slighlty less close to p_exact
[~, tbis] = bisection(f, -1, 2, max_iter, tol);
[~, tnew] = newton_raphson(f, df, 0.58, max_iter, tol);
[~, tsec] = secant(f, 0.5, 1.5, max_iter, tol);

figure(6)
fig = gcf;
fig.Position(3:4) = [586 438];
fs = 14; % Font Size
semilogy(tbis.("Iteration k"), abs(tbis.pk - p_exact), '-o', ...
         tnew.("Iteration k"), abs(tnew.pk - p_exact), '-s', ...
         tsec.("Iteration k"), abs(tsec.pk - p_exact), '-d', ...
         'LineWidth', 1.5)
legend({'Bisection (a = -1, b = 2)', ...
        'Newton (p_0 = 0.58)', ...
        'Secant (p_0 = 0.5, p_1 = 1.5)'}, 'FontSize', fs)
xlabel('Iteration', 'FontSize', fs)
ylabel('Absolute Error', 'FontSize', fs)
ax = gca;
ax.FontSize = fs;
grid on
