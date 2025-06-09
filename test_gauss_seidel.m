function tests = test_gauss_seidel
    % Unit tests for the custom Gauss-Seidel solver by Javier Almonacid
    tests = functiontests(localfunctions);
end

function testKnownSolution(testCase)
    A = [4 -1 0; -1 4 -1; 0 -1 4];
    b = [15; 10; 10];
    x0 = zeros(3,1);
    tol = 1e-6;
    max_iter = 100;

    [x, iter, xall] = gauss_seidel(A, b, x0, tol, max_iter);

    expected = A \ b;
    verifyLessThanOrEqual(testCase, norm(x - expected, inf), tol*5);
    verifyLessThanOrEqual(testCase, iter, max_iter);
    verifyEqual(testCase, size(xall,1), 3);
    verifyGreaterThanOrEqual(testCase, size(xall,2), 1);
end

function testConvergenceInOneStep(testCase)
    A = eye(3);
    b = [1; 2; 3];
    x0 = zeros(3,1);
    tol = 1e-12;
    max_iter = 5;

    [x, iter, ~] = gauss_seidel(A, b, x0, tol, max_iter);

    verifyEqual(testCase, x, b);
    verifyLessThanOrEqual(testCase, iter, 2);
end

function testNonConvergence(testCase)
    A = [1 3; 2 1]; % Not diagonally dominant
    b = [5; 5];
    x0 = [0; 0];
    tol = 1e-10;
    max_iter = 5;

    [~, iter, ~] = gauss_seidel(A, b, x0, tol, max_iter);

    % Expect no convergence within 5 iterations
    verifyEqual(testCase, iter, max_iter);
end

function testFinalXMatchesXall(testCase)
    A = diag([5 5 5]) + diag([-1 -1], 1) + diag([-1 -1], -1);
    b = [5; 10; 15];
    x0 = zeros(3,1);
    tol = 1e-8;
    max_iter = 100;

    [x, ~, xall] = gauss_seidel(A, b, x0, tol, max_iter);
    verifyEqual(testCase, x, xall(:,end), 'AbsTol', 1e-10);
end
