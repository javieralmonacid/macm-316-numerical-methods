% MACM 316: Numerical Analysis I
% Demo: An application of bisection method, including an error plot
clear;

% Define f(x) and the root p
f = @(x) sin(x);
pexact = pi;

% Parameters
a = 1; % Left endpoint a
b = 5; % Right endpoint b
N = 100; % Maximum number of iterations
tol = 1e-16; % Tolerance

% Call the function bisection.m
[p,table] = bisection(f,a,b,N,tol);

% Plot the error versus iteration number
figure(1);
semilogy(abs(table.('pk') - pexact*ones(height(table),1)), '-*')
xlabel('Iteration Number k','FontSize',16)
ylabel('Absolute Error |p-p_k|','FontSize',16)
title('Bisection Method Error','FontSize',16)
ax = gca;
ax.FontSize = 16;
grid on
grid minor