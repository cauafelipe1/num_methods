clear; clc;
addpath('../algorithms');

disp("=============================================");
disp("       TESTE DO METODO DE GAUSS-SEIDEL       ");
disp("=============================================");

% Matriz estritamente diagonal dominante para garantir convergência
A = [4, 1, -1; 
     2, 7,  1; 
     1, -3, 12];

% Vetor b
b = [3; 19; 31];

disp("Matriz de entrada A:");
disp(A);

disp("Vetor de entrada b:");
disp(b);
disp("Matriz aumentada [A b]");
disp([A b]);

disp("---------------------------------------------");
disp("Teste 1: Parametros padrao (default)");
disp("---------------------------------------------");
[x1, iter1] = gauss_seidel(A, b);

disp("Solucao x:");
disp(x1);
disp(["Convergencia atingida em iteracoes: ", num2str(iter1)]);

disp("---------------------------------------------");
disp("Teste 2: Parametros customizados");
disp("---------------------------------------------");
% Definindo chute inicial, tolerancia rigorosa e limite de iterações
x0 = [0; 0; 0];
tol = 1e-10;
max_iter = 100;

[x2, iter2] = gauss_seidel(A, b, x0, tol, max_iter);

disp("Solucao x (com tol = 1e-10):");
disp(x2);
disp(["Convergencia atingida em iteracoes: ", num2str(iter2)]);

disp("---------------------------------------------");
disp("Verificacao dos Resultados");
disp("---------------------------------------------");
disp("Calculando A * x (deve ser igual a b):");
disp(A * x1);

disp("Solucao embutida do Octave (A \\ b):");
disp(A \ b);