function G = gram_schmidt(A)
% gram_schimidt 
% computes an orthogonal basis from the columns of matrix A.
%   G = gram_schmidt(A) returns matrix G with orthogonal columns using 
%   the classical gram-schmidt process.
	A = [1 2 3 ; 4 5 5 ; 7 8 9 ];

    [rows, cols] = size(A);
    G = zeros(rows, cols);

    for k = 1:cols
        current_v = A(:, k);
        g = current_v;

        for j = 1:k-1
            prev = G(:, j);
            g = g - proj_v(current_v, prev);
        endfor

        G(:, k) = g;
    endfor
endfunction