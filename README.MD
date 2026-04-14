# num_methods

C++ and octave based repository for academic studies related to Numerical Methods and it's algorithms

## C++

### Features
- **Linear Systems**
  - Gaussian Elimination
  - Cholesky decomposition

### Requirements
- C++17 compiler (GCC/Clang)
- CMake 3.10+

### Build

```bash
# on cpp dir
mkdir build && cd build
cmake ..
make
```

## Octave

### Features
- Back substitution
- Cholesky decomposition
- Forward susbtitution
- Gauss elimination
- Gauss-Jordan elimination
- Gram Schmidt process
- Gram Schmidt GR decomposition
- Projections:
  - Scalar projection
  - Vector projection
  - Orthogonal complement
- QR decomposition

### Requirements
- GNU Octave, version 10.2.0

### Build

```bash
# on octave dir
octave
octave:1> run main.m
```