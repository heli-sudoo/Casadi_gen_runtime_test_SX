The function `Read_Dyn_Info_Deprecated` has been deprecated. This function is computationally expensive since it fills in the sparse matrix every time it is called, resulting in complexity O((nnz+1)nnz/2) where nnz is the total number of nonzeros. To fill in the sparse matrix in a more efficient way. Use, for example,
```cpp
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Matrix<double, 14,1> StateVec;
typedef Eigen::Matrix<double, 4, 1> ContrlVec;

SpMat H(xsize, xsize);
SpMat tau(qsize, 1);
StateVec x;
ContrlVec u;
Initiailize(H, tau);
Read_Dyn_Info(x, u, H, tau);
```

`Initialize(H,tau)` looks for sparsity info for H and tau and fill in with zeros. `Read_Dyn_Info()` updates with runtime values by using e.g., H.valPtr().