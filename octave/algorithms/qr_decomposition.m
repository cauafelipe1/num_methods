function [Q, R] = qr_decomposition(A)
% qr_decomposition
%  computes the QR decomposition of matrix A.
%   [Q, R] = qr_decomposition(A) returns orthonormal matrix Q and 
%   upper triangular matrix R such that A = Q * R.
%   
%   hint: since Q is orthonormal, Q' = inv(Q) => Q'Q = I, which allows the following derivation:
%    Ax = b => 
%    QRx = b =>
%    Q'QRx = Q'b =>
%    Rx = Q'b

    [rows, cols] = size(A);

    Q = zeros(rows, cols);
    R = zeros(cols, cols);

    for k = 1:cols
        current_v = A(:, k);
        q = current_v;
        
        % orthogonalization step
        for j = 1:k-1
            R(j, k) = dot(current_v, Q(:, j));
            q = q - R(j, k) * Q(:, j);
        endfor

        % normalization step
        R(k, k) = norm(q);
        Q(:, k) = q / R(k, k);
    endfor
endfunction