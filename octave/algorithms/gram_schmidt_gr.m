function [G, R] = gram_schmidt_gr(A)
% gram_schmidt_gr
%  gramm schmidt decomposition
%  computes an orthogonal basis from the columns of matrix A
%   [G, R] = gram_schmidt_gr(A) returns orthogonal matrix G and 
%   upper triangular matrix R such that A = G * R.

    [rows, cols] = size(A);

    G = zeros(rows, cols);    
    R = eye(cols); % upper triangular matrix with 1s on the diagonal

    for k = 1:cols
        current_v = A(:, k);
        g = current_v;

        % orthogonalization step
        for j = 1:k-1
            prev = G(:, j);
            
            R(j, k) = dot(current_v, prev) / dot(prev, prev);        
            g = g - proj_v(current_v, prev);
        endfor

        G(:, k) = g;
    endfor
endfunction