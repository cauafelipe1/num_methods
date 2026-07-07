function [x, iter] = gauss_seidel(A, b, x0, tol, max_iter)
% gauss_seidel
%  solves the linear system Ax = b using the Gauss-Seidel iterative method.
%   x = gauss_seidel(A, b) returns the solution vector x using default parameters.
%   [x, iter] = gauss_seidel(A, b, x0, tol, max_iter) allows you to specify the 
%   initial guess, tolerance, and maximum number of iterations.
%
%   hint: Gauss-Seidel requires the matrix to have non-zero elements on the 
%   main diagonal. It is guaranteed to converge for strictly diagonally dominant 
%   or symmetric positive-definite matrices.

    [n, m] = size(A);
    if n != m
        error("gauss_seidel wraps only square matrices.");
    endif

    % default arguments handling
    if nargin < 3 || isempty(x0)
        x0 = zeros(n, 1);
    endif
    if nargin < 4 || isempty(tol)
        tol = 1e-6;
    endif
    if nargin < 5 || isempty(max_iter)
        max_iter = 1000;
    endif

    % check for zeros on the main diagonal to prevent division by zero
    for i = 1:n
        if abs(A(i, i)) < 1e-15
            error("zeros on the main diagonal are not allowed for gauss-seidel. consider row swapping.");
        endif
    endfor

    x = x0;
    x_prev = x0;

    for iter = 1:max_iter
        for i = 1:n
            % sum_knowns uses updated values for j < i and old values for j > i
            sum_knowns = ( 
            A(i, 1:i-1) * x(1:i-1)
            +
            A(i, i+1:n) * x_prev(i+1:n));
            
            x(i) = (b(i) - sum_knowns) / A(i, i);
        endfor
        
        % check for convergence using the infinity norm (max absolute difference)
        err = norm(x - x_prev, inf);
        if err < tol
            return;
        endif
        
        x_prev = x;
    endfor

    warning("gauss_seidel did not converge within the maximum number of iterations.");
endfunction