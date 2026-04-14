function [x, Ab] = gauss_jordan(A, b)
% gauss_jordan
%  solves the system Ax = b by reducing A to the identity matrix.
%   x = gauss_jordan(A, b) returns the solution vector x.
%   [x, Ab] = gauss_jordan(A, b) also returns the final augmented matrix.
%   
%   hint: you can also apply the Gauss-Jordan method to inverse a matrix by passing eye(m) as b
%   where m is the number of columns of A.

    [n, m] = size(A);
    
    % augmented matrix [A | b]
    Ab = [A, b];

    for j = 1:n
        % partial pivoting
        [max_val, max_idx] = max(abs(Ab(j:n, j)));
        pivot_row = max_idx + j - 1;

        % swap current row with pivot row
        if pivot_row != j
            Ab([j, pivot_row], :) = Ab([pivot_row, j], :);
        endif

        if max_val < 1e-15
            error("can't apply gauss-jordan method on a singular or near-singular matrix");
        endif

        % normalize the pivot row to have 1 on the diagonal
        Ab(j, :) = Ab(j, :) / Ab(j, j);

        % full elimination
        for i = 1:n
            if i != j
                factor = Ab(i, j);
                Ab(i, j:end) = Ab(i, j:end) - factor * Ab(j, j:end);
            endif
        endfor
    endfor
    x = Ab(:, n+1:end);
endfunction