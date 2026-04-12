clear; clc;
addpath('algorithms');

% tests

disp("[Cholesky decomposition]");
disp("input:")
% symetric and defined positive matrix
A = [4, 12, -16; 12, 37, -43; -16, -43, 98]

disp("resulting matrix lower triangular");
L = cholesky(A)

disp("resulting matrix upper triangular (U = L')");
U = L'

disp("verification (L * U) = (L * L') = A:")
L * U

disp("[Gram-Schmidt decomposition]");
disp("input:"); A
[G, R] = gram_schmidt_gr(A)

disp("verification (G * R = A):")
GR = G * R


disp("[QR decomposition]");
disp("input:"); A
[Q, R] = qr_decomposition(A)

disp("Q verification (Q' = inv(Q)):")
I = inv(Q) * Q
I = Q' * Q

disp("QR verification (A = Q * R):")
QR = Q * R

% applies QR and backsubstitution to solve Ax = b
disp("[Ax = b system using QR decomposition]")
disp("input:"); A
b = [1; 1; 1]

disp("QR decomposition:")
[Q, R] = qr_decomposition(A)

disp("since Q is orthonormal, the following derivation preceeds:")
disp("Ax = b =>");
disp("QRx = b =>");
disp("Q'QRx = Q'b =>");
disp("Rx = Q'b");

y = Q' * b;

disp("solving x using back_substitution:");
x = back_substitution(R, y)

disp("octave solution (A \\ b):"), disp(A \ b);