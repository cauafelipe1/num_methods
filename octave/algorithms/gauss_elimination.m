function [output_1, output_2] = gauss_elimination(A, b)
% gauss_elimination
%  performs Gaussian Elimination with partial pivoting.
%   Ab = gauss_elimination(A, b) returns the augmented matrix [A|b] in upper triangular form.
%   [A, b] = gauss_elimination(A, b) returns the upper triangular matrix A and modified vector b separately. 
% 
%   hint: if you chose to work with the augmented matrix Ab, call Ab(:, 1:m) to get A and Ab(:, end) to get b,
%   where m is the number of columns of A.

    [n, m] = size(A);

    % augmented matrix [A | b]
    Ab = [A, b];

    for j = 1:n-1
        % partial pivoting
        [max_val, max_idx] = max(abs(Ab(j:n, j)));
        pivot_row = max_idx + j - 1;

        % swap current row with pivot row
        if pivot_row != j
            Ab([j, pivot_row], :) = Ab([pivot_row, j], :);
        endif

        if max_val < 1e-15
            error("can't apply gauss elimination method on a singular or near-singular matrix");
        endif

        % elimination
        for i = j+1:n
            factor = (Ab(i, j) / Ab(j, j));
            Ab(i, j:end) = Ab(i, j:end) - factor * Ab(j, j:end);
        endfor
    endfor
    
    % check output arguments
    if nargout <= 1
        % Ab = gauss_elimination(...)
        output_1 = Ab;
    else
        % [A, b] = gauss_elimination(...)
        output_1 = Ab(:, 1:m);      % triangular matrix A
        output_2 = Ab(:, end);      % modified vector b
    endif
endfunction