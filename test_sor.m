function tests = test_sor
    % Unit tests for the custom SOR solver by Javier Almonacid
    tests = functiontests(localfunctions);
end

function testKnownSystem(testCase)
    A = [4 -1 0; -1 4 -1; 0 -1 4];
    b = [15; 10; 10];
    x0 = zeros(3,1);
    omega = 1.1;  % Relaxation factor
    tol = 1e-6;
    max_iter = 100;

    [x, iter, xall] = sor(A, b, x0, omega, tol, max_iter);
    expected = A \ b;

    verifyLessThanOrEqual(testCase, norm(x - expected, inf), tol*5);
    verifyLessThanOrEqual(testCase, iter, max_iter);
    verifyEqual(testCase, size(xall,1), length(b));
    verifyGreaterThanOrEqual(testCase, size(xall,2), 1);
end

function testConvergenceWithOmega1(testCase)
    A = eye(2);
    b = [3; 5];
    x0 = [0; 0];
    omega = 1.0;  % Equivalent to Gauss-Seidel
    tol = 1e-12;
    max_iter = 10;

    [x, iter, ~] = sor(A, b, x0, omega, tol, max_iter);
    verifyEqual(testCase, x, b);
    verifyLessThanOrEqual(testCase, iter, 2);
end

function testNonConvergenceWithBadOmega(testCase)
    A = [2 -1; -1 2]; % Diagonally dominant
    b = [1; 0];
    x0 = zeros(2,1);
    omega = 1.95;  % Too aggressive (may fail to converge)
    tol = 1e-12;
    max_iter = 5;

    [~, iter, ~] = sor(A, b, x0, omega, tol, max_iter);
    verifyEqual(testCase, iter, max_iter);
end

function testFinalXMatchesXall(testCase)
    A = diag([4 4 4]) + diag([-1 -1],1) + diag([-1 -1],-1);
    b = [10; 20; 30];
    x0 = zeros(3,1);
    omega = 1.2;
    tol = 1e-8;
    max_iter = 100;

    [x, ~, xall] = sor(A, b, x0, omega, tol, max_iter);
    verifyEqual(testCase, x, xall(:,end), 'AbsTol', 1e-10);
end
