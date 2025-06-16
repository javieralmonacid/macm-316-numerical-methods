% MACM 316: Numerical Analysis I
% Demo: Performance of Newton's method for the case of roots with
% multiplicity >= 1
clear; close all; clc;

f1 = @(x) exp(x)-1;
df1 = @(x) exp(x);
f2 = @(x) exp(x) - 1 - x;
df2 = @(x) exp(x) - 1;
f3 = @(x) exp(x) - 1 - x - x^2/2;
df3 = @(x) exp(x) - 1 - x;

[p1,t1] = newton_raphson(f1, df1, 1, 25, 1e-16);
[p2,t2] = newton_raphson(f2, df2, 1, 25, 1e-16);
[p3,t3] = newton_raphson(f3, df3, 1, 25, 1e-16);

figure(1)
fig = gcf;
fig.Position(3:4) = [586 438];
fs = 14;
semilogy(t1.("Iteration k"), t1.pk, '-o', ...
         t2.("Iteration k"), t2.pk, '-s', ...
         t3.("Iteration k"), t3.pk, '-d', ...
         'LineWidth', 1.5)
legend({'$\exp(x)-1$','$\exp(x)-1-x$','$\exp(x)-1-x-x^2/2$'}, 'FontSize', fs, 'Interpreter','latex', 'Location', 'SouthEast')
xlabel('Iteration', 'FontSize', fs)
ylabel('Absolute Error', 'FontSize', fs)
ax = gca;
ax.FontSize = fs;
grid on

