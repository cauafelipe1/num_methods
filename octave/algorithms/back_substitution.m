function x = back_substitution(R, b)
% back_substitution 
%  solves the system Rx = b for an upper triangular matrix R.
%   x = back_substitution(R, b) returns the solution vector x.

    [n, m] = size(R);
    if n != m
        error("back_substitution wraps only square matrices.");
    endif

    x = zeros(n, 1);

    % start from the back (last row)
    for i = n:-1:1
        % makes use of already calculated x values
        sum_prev = R(i, i+1:n) * x(i+1:n);
        x(i) = (b(i) - sum_prev) / R(i, i);
    
    endfor
endfunction