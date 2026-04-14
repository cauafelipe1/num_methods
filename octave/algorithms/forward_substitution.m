function x = forward_substitution(L, b)
% back_substitution 
%  solves the system Lx = b for a lower triangular matrix L.
%   x = forward_substitution(L, b) returns the solution vector x.

    [n, m] = size(L);
    if n != m
        error("forward_substitution wraps only square matrices.");
    endif

    x = zeros(n, 1);

    % start from the first row and goes forward
    for i = 1:n
        sum_prev = L(i, 1:i-1) * x(1:i-1);   
        x(i) = (b(i) - sum_prev) / L(i, i);
    endfor
endfunction