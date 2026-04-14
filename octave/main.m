clear; clc;
addpath('algorithms');

% tests
% symetric and defined positive matrix
A = [4, 12, -16; 12, 37, -43; -16, -43, 98]
b = [1; 1; 1]

disp("MAIN TARGET: solve the system Ax = b, which can be described by augmented matrix [A|b]:")
[A, b]

disp("[Gauss Elimination with partial pivoting]");
disp("input:")
A
b

disp("augmented matrix after gauss elimination:")
Ab = gauss_elimination(A, b)

disp("solving the system Ax = b using back substitution:")
x = back_substitution(Ab(:, 1:end-1), Ab(:, end))

disp("[Gauss-Jordan elimination]");
disp("input:")
A
b

disp("solving the system Ax = b using gauss-jordan method (it goes directly!)")
[x, Ab] = gauss_jordan(A, b); x
disp("augmented matrix after gauss-jordan elimination:"); Ab

disp("you can also find the inverse of a matrix by using gauss-jordan method")
[A_inv, GJ] = gauss_jordan(A, eye(3)); A_inv
disp("augmented matrix after inverse calculation using gauss-jordan method"); GJ

disp("[Cholesky decomposition]");
disp("input:"); A

disp("resulting matrix lower triangular");
L = cholesky(A)

disp("resulting matrix upper triangular (U = L')");
U = L'

disp("verification (L * U) = (L * L') = A:")
L * U

disp("solving LUx = b using forward and backward substitution:")
x = back_substitution(U, forward_substitution(L, b))

disp("here we basically derivated:")
disp("Ax = b =>                   p.")
disp("LUx = b =>        cholesky(A).")
disp("Lx = y =>                  fs.")
disp("Uy = b =>                  bs.")

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
disp("input:"); 
A 
b

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

disp("[PA=LU decomposition using partial pivoting]")
disp("input:"); A
[P, L, U] = palu_decomposition(A)

disp("PA = LU derivation");
disp("Ax = b => ");
disp("PAx = Pb =>");
disp("LUx = Pb =>                      (substitute PA for LU)");
disp("Ly = Pb =>     (solve for y using forward substitution)");
disp("Ux = y            (solve for x using back substitution)");
disp("solving x");

x = back_substitution(U, forward_substitution(L, P * b))

disp("octave solution (A \\ b):"), disp(A \ b);

