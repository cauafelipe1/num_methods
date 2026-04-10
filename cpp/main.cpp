#include "matrix.hpp"
#include <cassert>
#include <iostream>
using namespace std;

int main() {
    // identity matrix
    auto I = Matrix<3, 3, double>::identity();
    cout << "Identity I:\n" << I;

    Matrix<3, 3, double> A;
    A(0, 0) = 4;   A(0, 1) = 12;  A(0, 2) = -16;
    A(1, 0) = 12;  A(1, 1) = 37;  A(1, 2) = -43;
    A(2, 0) = -16; A(2, 1) = -43; A(2, 2) = 98;
    cout << "\nMatrix A:\n" << A << "\n";
    
    // ELEMENTARY OPERATIONS
    // row swapping
    A.row_swap(0, 1);
    cout << "\nMatrix A row_swap(0, 1):\n" << A << "\n";
    A.row_swap(0, 1);

    // row multiplication
    A.row_multiply(0, 2);
    cout << "\nMatrix A row_multiply(0, 2):\n" << A << "\n";
    A.row_multiply(0, 0.5);
    cout << "\nMatrix A row_multiply(0, 0.5):\n" << A << "\n";
    
    // row addition
    A.row_add(0, 1);
    cout << "\nMatrix A row_add(0, 1):\n" << A << "\n";
    A.row_add(0, 1, -1);
    cout << "\nMatrix A row_add(0, 1, -1):\n" << A << "\n";

    // cholesky decomposition
    auto L = A.cholesky();
    cout << "\nMatrix L (lower triangular) after cholesky decomposition:\n" << L;
    auto U = L.transpose();
    cout << "\nMatrix U (upper triangular) after cholesky decomposition:\n" << U;
    auto LU = L * U;
    cout << "\nMatrix A = LU = L * L':\n";
    I *= LU;
    cout << "\n" << LU;
    assert(A == LU);

    // inverse test using gauss jordan
    Matrix<3, 3, double> M;
    M(0,0) = 2; M(0,1) = 1; M(0,2) = 1;
    M(1,0) = 1; M(1,1) = 2; M(1,2) = 1;
    M(2,0) = 1; M(2,1) = 1; M(2,2) = 2;
    cout << "\nMatrix M:\n" << M << "\n";
    cout << "\nMatrix M.inverse():\n" << M.inverse() << "\n";
    cout << "\nMatrix M * M.inverse():\n" << M * M.inverse() << "\n";
    assert(M * M.inverse() == M.inverse() * M);
    return 0;
}