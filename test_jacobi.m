function tests = test_jacobi
    % Unit tests for the custom Jacobi function by Javier Almonacid
    tests = functiontests(localfunctions);
end

function testSimpleSystem(testCase)
    A = [4 -1 0; -1 4 -1; 0 -1 4];
    b = [15; 10; 10];
    x0 = zeros(3, 1);
    tol = 1e-6;
    max_iter = 100;

    [x, iter, xall] = jacobi(A, b, x0, tol, max_iter); 

    % Reference solution
    expected = A \ b;

    % Check solution accuracy.
    verifyLessThanOrEqual(testCase, norm(x - expected, inf), tol*5);

    % Check iteration count within bounds
    verifyLessThanOrEqual(testCase, iter, max_iter);

    % Check that xall has correct shape
    verifyEqual(testCase, size(xall, 1), length(b));
    verifyGreaterThanOrEqual(testCase, size(xall, 2), 1);
end

function testConvergenceCriterion(testCase)
    A = eye(2);
    b = [1; 2];
    x0 = zeros(2,1);
    tol = 1e-10;
    max_iter = 10;

    [x, iter, ~] = jacobi(A, b, x0, tol, max_iter);

    verifyEqual(testCase, x, b);  % should converge in one iteration
    verifyLessThanOrEqual(testCase, iter, 2);
end

function testNoConvergence(testCase)
    % Use a matrix that is not diagonally dominant (may not converge)
    A = [1 2; 3 1];
    b = [3; 4];
    x0 = [0; 0];
    tol = 1e-12;
    max_iter = 5;

    [~, iter, ~] = jacobi(A, b, x0, tol, max_iter);

    % It should stop due to max_iter
    verifyEqual(testCase, iter, max_iter);
end

function testXallHasCorrectNumberOfColumns(testCase)
    A = [10 -1 0; -1 10 -1; 0 -1 10];
    b = [9; 7; 8];
    x0 = zeros(3,1);
    tol = 1e-8;
    max_iter = 100;

    [~, iter, xall] = jacobi(A, b, x0, tol, max_iter);
        
    verifyEqual(testCase, size(xall, 2), iter+1, 'AbsTol', 1e-10);
end

