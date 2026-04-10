#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <cmath>
#include <algorithm>

template <int R, int C, typename T>
struct Matrix {
public:    
    T data[R * C];

    Matrix(const T &v = T{0});
    ~Matrix();
    
    // @brief access matrix element at (i, j)
    T &operator()(int i, int j);

    // @brief access matrix element at (i, j)
    T operator()(int i, int j) const;

    // @brief sum of matrices
    Matrix<R, C, T> operator+(const Matrix& other) const;

    // @brief subtraction of matrices
    Matrix<R, C, T> operator-(const Matrix& other) const;
    
    // @brief multiplication of matrices
    template <int K>
    Matrix<R, K, T> operator*(const Matrix<C, K, T>& other) const;
    
    // @brief in-place sum of matrices
    Matrix<R, C, T> &operator+=(const Matrix& other);
    
    // @brief in-place subtraction of matrices
    Matrix<R, C, T> &operator-=(const Matrix& other);
    
    // @brief multiplication by scalar
    // @param scalar value
    Matrix<R, C, T> &operator*=(const T& scalar);
    
    // @brief division by scalar
    Matrix<R, C, T> &operator/=(const T& scalar);
    
    // @brief unary minus / negative matrix
    Matrix<R, C, T> operator-() const;
    
    // @brief in-place multiplication of matrices
    Matrix<R, C, T> &operator*=(const Matrix<R, C, T> &other);
    
    // @brief equality operator
    bool operator==(const Matrix& other) const;
    
    // @brief inequality operator
    bool operator!=(const Matrix& other) const;
        
    // @brief elementary row multiplication operation
    // @concept r -> k * r1
    // @param  r multiplied row
    // @param k number which the row is multiplied 
    void row_multiply(int r, T k);
    
    // @brief elementary row addition operation
    // @concept r1 -> r1 + k * r2
    // @param r1 operated row
    // @param r2 added row
    // @param k scalar multipling added to r2, 1 by default
    void row_add(int r1, int r2, T k = T{1});
    
    // @brief elementary row swap operation
    // @concept r1 -> r2 & r2 -> r1
    // @param r1 first row
    // @param r2 second row
    void row_swap(int r1, int r2);
    
    // @brief identity matrix
    // @returns static identity matrix
    static Matrix<R, C, T> identity();

    // @brief transposed matrix
    // @returns a new matrix that is the transpose of the current matrix
    Matrix<C, R, T> transpose();

    // @brief gauss elimination without pivot
    // @concept transforms the matrix into upper triangular form (U) using elementary row operations 
    // @returns a reference of the matrix after elimination
    Matrix<R, C, T> &unpivoted_gauss_elimination();

    // @brief gauss elimination using partial pivoting
    // @concept still transforms the matrix into upper triangular form (U) using elementary row operations
    //  but instead finds the row with the largest absolute value in the current column
    // and swap it with the current row.
    // @returns a reference of the matrix after elimination
    Matrix<R, C, T> &gauss_elimination();
    
    // @brief inverse matrix using gauss-jordan elimination
    // @concept augments the matrix with the identity matrix and applies row operations
    // until the original matrix becomes the identity.
    // @returns new matrix which is the inverse of the current matrixs
    Matrix<R, C, T> gauss_jordan_inverse();

    // @brief inverse matrix currently using gauss-jordan elimination
    // @returns new matrix which is the inverse of the current matrix
    Matrix<R, C, T> inverse();
    
    // @brief cholesky decomposition
    // @concept decomposes a positive-definite matrix into the product of a lower triangular matrix
    //  and its conjugate transpose.
    // @example: 
    //  Matrix<3, 3, double> A;
    //  ...
    //  auto L = A.cholesky();
    //  auto U = L.transpose();
    //  assert(A = L * U); true!
    // @returns lower triangular matrix L which is equal to A when multiplied by its on tranpose U (upper triangular)
    Matrix<R, C, T> cholesky() const;
    
    // @brief output stream operator for matrix printing
    // @returns a reference to the output stream
    template <int R_orig, int C_orig, typename T_orig>
    friend std::ostream& operator<<(std::ostream& os, const Matrix<R_orig, C_orig, T_orig>& mat);
};

template <int R, int C, typename T>
Matrix<R, C, T>::Matrix(const T &v) {
    std::fill(data, data + (R * C), v);
}
template <int R, int C, typename T>
Matrix<R, C, T>::~Matrix() {}

template <int R, int C, typename T>
T& Matrix<R, C, T>::operator()(int i, int j) {
    return data[i * C + j];
}
template <int R, int C, typename T>
T Matrix<R, C, T>::operator()(int i, int j) const {
    return data[i * C + j];
}
template <int R, int C, typename T>
Matrix<R, C, T> Matrix<R, C, T>::operator+(const Matrix<R, C, T> &other) const {
    Matrix<R, C, T> res;
    for (int i = 0; i < R * C; i++) {
        res.data[i] = this->data[i] + other.data[i];
    }
    return res;
}
template <int R, int C, typename T>
Matrix<R, C, T> Matrix<R, C, T>::operator-(const Matrix<R, C, T> &other) const {
    Matrix<R, C, T> res;
    for (int i = 0; i < R * C; i++) {
        res.data[i] = this->data[i] - other.data[i];
    }
    return res;
}
template <int R, int C, typename T>
template <int K>
Matrix<R, K, T> Matrix<R, C, T>::operator*(const Matrix<C, K, T> &other) const {
    Matrix<R, K, T> result;    
    for (int i = 0; i < R; i++) {
        for (int j = 0; j < K; j++) {
            T sum = 0;
            for (int k = 0; k < C; k++) {
                sum += (*this)(i, k) * other(k, j);
            }
            result(i, j) = sum;
        }
    }
    return result;
}
template <int R, int C, typename T>
Matrix<R, C, T> &Matrix<R, C, T>::operator+=(const Matrix<R, C, T> &other) {
    for (int i = 0; i < R * C; i++) {
        this->data[i] = this->data[i] + other.data[i];
    }
    return *this;
}
template <int R, int C, typename T>
Matrix<R, C, T> &Matrix<R, C, T>::operator-=(const Matrix<R, C, T> &other) {
    for (int i = 0; i < R * C; i++) {
        this->data[i] = this->data[i] - other.data[i];
    }
    return *this;
}
template <int R, int C, typename T>
Matrix<R, C, T> &Matrix<R, C, T>::operator*=(const T &scalar) {
    for (int i = 0; i < R * C; i++) {
        this->data[i] *= scalar;
    }
    return *this;
}
template <int R, int C, typename T>
Matrix<R, C, T> &Matrix<R, C, T>::operator/=(const T &scalar) {
    if (scalar == T{0}) throw std::runtime_error("division by zero");
    for (int i = 0; i < R * C; i++) {
        this->data[i] /= scalar;
    }
    return *this;
}
template <int R, int C, typename T>
bool Matrix<R, C, T>::operator==(const Matrix<R, C, T> &other) const {
    for (int i = 0; i < R * C; i++) {
        if (data[i] != other.data[i]) return false;
    }
    return true;
}

template <int R, int C, typename T>
bool Matrix<R, C, T>::operator!=(const Matrix<R, C, T> &other) const {
    return !(*this == other);
}
template <int R, int C, typename T>
Matrix<R, C, T> Matrix<R, C, T>::operator-() const {
    Matrix<R, C, T> res;
    for (int i = 0; i < R * C; i++) res.data[i] = -data[i];
    return res;
}
template <int R, int C, typename T>
Matrix<R, C, T> &Matrix<R, C, T>::operator*=(const Matrix<R, C, T> &other) {
    static_assert(R == C, "In-place multiplication (*=) only allowed for square matrices");
    
    *this = (*this) * other; 
    return *this;
}

template <int R, int C, typename T>
void Matrix<R, C, T>::row_multiply(int r, T k) {
    if (r < 0 || r >= R) {
        throw std::out_of_range("Row index out of range");
    }
    T* rowPtr = &data[r * C];
    for (int j = 0; j < C; ++j) {
        rowPtr[j] *= k;
    }
}

template <int R, int C, typename T>
void Matrix<R, C, T>::row_add(int r1, int r2, T k) {
    if (r1 < 0 || r1 >= R || r2 < 0 || r2 >= R) {
        throw std::out_of_range("Row indexes out of range");
    }
    T* row1 = &data[r1 * C]; // destination row
    T* row2  = &data[r2 * C]; // added row

    for (int j = 0; j < C; ++j) {
        // r1 -> k * r2
        row1[j] += k * row2[j];
    }
}

template <int R, int C, typename T>
void Matrix<R, C, T>::row_swap(int r1, int r2) {
    if (r1 == r2) return;
    
    if (r1 < 0 || r1 >= R || r2 < 0 || r2 >= R) {
        throw std::out_of_range("Row indexes out of range");
    }

    T* row1 = &data[r1 * C];
    T* row2 = &data[r2 * C];

    for (int j = 0; j < C; ++j) {
        T temp = row1[j];
        row1[j] = row2[j];
        row2[j] = temp;
    }
}

template <int R, int C, typename T>
Matrix<R, C, T> &Matrix<R, C, T>::unpivoted_gauss_elimination() {
    int steps = (R < C) ? R : C;

    for (int j = 0; j < steps; ++j) {
        T pivot = (*this)(j, j);

        if (std::abs(pivot) <= 1e-15) throw std::runtime_error("Near zero pivot found... Consider using gauss_elimination()");

        for (int i = j + 1; i < R; ++i) {
            T factor = -(*this)(i, j) / pivot;
            this->row_add(i, j, factor);
        }
    }
    return *this;
}

template <int R, int C, typename T>
Matrix<R, C, T> &Matrix<R, C, T>::gauss_elimination() {
    int steps = (R < C) ? R : C;

    for (int j = 0; j < steps; ++j) {

        // partial pivoting
        int max_row = j;
        T max_value = std::abs((*this)(j, j));

        for (int k = j + 1; k < R; ++k) {
            T current = std::abs((*this)(k, j));
            if (current > max_value) {
                max_value = current;
                max_row = k;
            }
        }

        if (max_row != j) {
            this->row_swap(j, max_row);
        }
        
        // proceeds to standard elimination
        T pivot = (*this)(j, j);

        if (std::abs(pivot) <= 1e-15) throw std::runtime_error("Can't apply gauss elimination on a singular or near-singular Matrix");

        for (int i = j + 1; i < R; ++i) {
            T factor = -(*this)(i, j) / pivot;
            this->row_add(i, j, factor);
        }
    }
    return *this;
}

template <int R, int C, typename T>
Matrix<R, C, T> Matrix<R, C, T>::gauss_jordan_inverse() {
    Matrix<R, C, T> A = *this;
    Matrix<R, C, T> inv = Matrix<R, C, T>::identity();

    for (int j = 0; j < R; ++j) {
        // partially pivoting
        int max_row = j;
        T max_val = std::abs(A(j, j));

        for (int k = j + 1; k < R; ++k) {
            if (std::abs(A(k, j)) > max_val) {
                max_val = std::abs(A(k, j));
                max_row = k;
            }
        }

        if (max_val < 1e-15) {
            throw std::runtime_error("Can't inverse using gauss-jordan method on a singular matrix");
        }

        // elementary operations
        // has to be applied on both matrices
        if (max_row != j) {
            A.row_swap(j, max_row);
            inv.row_swap(j, max_row);
        }

        T pivot_val = A(j, j);
        A.row_multiply(j, T{1} / pivot_val);
        inv.row_multiply(j, T{1} / pivot_val);

        for (int i = 0; i < R; ++i) {
            if (i != j) {
                T factor = -A(i, j);
                A.row_add(i, j, factor);
                inv.row_add(i, j, factor);
            }
        }
    }
    return inv;
}

template <int R, int C, typename T>
Matrix<R, C, T> Matrix<R, C, T>::inverse() {
    return gauss_jordan_inverse();
}

template <int R, int C, typename T>
Matrix<R, C, T> Matrix<R, C, T>::identity() {
    static_assert(R == C, "Identity matrix must be square");
    Matrix<R, C, T> res;
    for (int i = 0; i < R; i++) {
        res(i, i) = T{1};
    }
    return res;
}

template<int R, int C, typename T>
Matrix<C, R, T> Matrix<R, C, T>::transpose() {
    Matrix<C, R, T> res;
    for (int i = 0; i < R; i++) {
        for (int j = 0; j < C; j++) {
            res(j, i) = (*this)(i, j);
        }
    }
    return res;
}

template<int R, int C, typename T>
Matrix<R, C, T> Matrix<R, C, T>::cholesky() const {
        static_assert(R == C, "cholesky algorith wraps only square matrices");
        
        Matrix<R, C, T> L;
        for (int row = 0; row < R; row++) {
            for (int col = 0; col <= row; col++) {
                T sum = 0;
                for (int i = 0; i < col; i++) {
                    sum += L(row, i) * L(col, i);
                }
                if (row == col) {
                    T val = (*this)(row, row) - sum;
                    if (val < 0) throw std::runtime_error("not a positive-definite matrix");
                    L(row, col) = std::sqrt(val);
                } else {
                    L(row, col) = ((*this)(row, col) - sum) / L(col, col);
                }
            }
        }
        return L;
    }
template <int R, int C, typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<R, C, T>& mat) {
    for (int i = 0; i < R; i++) {
        os << "[ ";
        for (int j = 0; j < C; j++) {
            os << mat(i, j) << (j == C - 1 ? "" : " ");
        }
        os << " ]\n";
    }
    return os;
}
#endif