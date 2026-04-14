function [P, L, U] = palu_decomposition(A)
% palu_decomposition
%  computes the PA=LU decomposition of matrix A using partial pivoting.
%   [P, L, U] = palu_decomposition(A) returns permutation matrix P, 
%   lower triangular matrix L and upper triangular matrix U such that P * A = L * U.
%   
%   hint: the permutation matrix P handles row swaps to ensure numerical stability 
%   and prevent division by zero (especially when a pivot on the main diagonal is zero). 
%   The following derivation to solve the system Ax = b preceeds:
%    Ax = b => 
%    PAx = Pb =>
%    LUx = Pb =>                      (substitute PA for LU)
%    Ly = Pb =>     (solve for y using forward substitution)
%    Ux = y            (solve for x using back substitution)

    [n, ~] = size(A);
    P = eye(n);
    L = eye(n);
    U = A;

    for j = 1:n-1
        % partial pivoting
        [max_val, max_idx] = max(abs(U(j:n, j)));
        pivot_row = max_idx + j - 1;

        if max_val < 1e-15
            error("can't apply PA=LU method on a singular or near-singular matrix");
        endif

        # row swap
        if pivot_row != j
            U([j, pivot_row], :) = U([pivot_row, j], :);
            P([j, pivot_row], :) = P([pivot_row, j], :);

            if j > 1
                L([j, pivot_row], 1:j-1) = L([pivot_row, j], 1:j-1);
            endif
        endif
        
        % elimination
        for i = j+1:n
            factor = U(i, j) / U(j, j);
            L(i, j) = factor;
            U(i, j:end) = U(i, j:end) - factor * U(j, j:end);
        endfor
    endfor
endfunction