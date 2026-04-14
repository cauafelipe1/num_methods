function L = cholesky(A)
% cholesky decomposition
%  decomposes a positive-definite matrix into the product of a lower triangular matrix and its conjugate transpose.
%    L = cholesky(A) returns the lower triangular matrix L which is equal to A when multiplied by its on tranpose U (upper triangular)
%    U = L'
%    L * U = L * L' = A
    [n, m] = size(A);
    if n != m
        error("cholesky algorith wraps only square matrices");
    endif

    L = zeros(n, n);

    for i = 1:n
        for j = 1:i
            s = L(i, 1:j-1) * L(j, 1:j-1)';
            if (i == j)
                L(i, j) = sqrt(A(i, i) - s);
		if (L(i, j) < 0)
			error("not a positive-definite matrix");
		endif
            else
                L(i, j) = (A(i, j) - s) / L(j, j);
            endif
        endfor
    endfor
endfunction
