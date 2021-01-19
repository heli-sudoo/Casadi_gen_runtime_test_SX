
// #include <iostream>
// #include <vector>
// #include <string>
// #include <cstring>
// #include <chrono> // system clock
// #include <memory> // smart pointer 

// #include "Dyn_FL_cas.h"
// #include "Dyn_FL_par_cas.h"
// #include "Jacob_F_cas.h"
// #include "Jacob_B_cas.h"
// #include "Dyn_FL_sub.h"
// #include "Dyn_FL_sub_par.h"
// #include "dynamics.h"



// int main()
// {   
//     cout << "This code benchmark runtime of casadi_gen func and matlab_gen fun." << endl;    

//     x.setZero();
//     u.setZero();
//     y.setZero();
//     q.setZero();


//     auto start = high_resolution_clock::now();
//     Forward_dyn(x,u,x_next,y);
//     auto stop = high_resolution_clock::now();
//     auto duration = duration_cast<microseconds>(stop - start);
//     cout << "Runtime : " << duration.count() << " microseconds " <<endl;
    
//     // Initialize_Mat(H, tau);
//     Initialize_Mat_all(); 
//     auto start2 = high_resolution_clock::now();
//     Eigen::BiCGSTAB<SpMat> solver;
//     Read_Dyn_Info(x, u, H, tau);
//     solver.compute(H); 
//     auto q_next = solver.solve(tau);
//     auto stop2 = high_resolution_clock::now();
//     auto duration2 = duration_cast<microseconds>(stop2 - start2);
//     cout << "Runtime : " << duration2.count() << " microseconds " <<endl; 
//     Eigen::saveMarket(H, "H.mtx");
//     Eigen::saveMarketVector(CoorVec(tau), "H_b.mtx");
//     cout << "J = " << "\n" << MatrixXd(JF) << endl;

// }

// void Initialize_Mat_all()
// {
//     vector<SpMat*>Mat{&H, &tau};
//     Initialize_Mat_new(Dyn_FL_sub_sparsity_out, Mat);
//     cout << "H = " << "\n" << H << endl;

//     Mat = {&JF, &JFd};
//     Initialize_Mat_new(Jacob_F_sparsity_out, Mat);
//     cout << "JF = " << "\n" << JFd << endl;

// }

// void Initialize_Mat_new(const int_T *f_sparse_out(int_T), vector<SpMat*>&mat)
// {
//     const int n_out = mat.size();
//     size_t ou, i_out, col, col_temp;
//     size_t row, i_nz = 0;

//     const casadi_int* ou_sp_info;
//     const casadi_int* ou_index[n_out];
//     const casadi_int* ou_colptr_info[n_out];
//     const casadi_int* ou_row_info[n_out];
//     casadi_int nnz_ou[n_out];
//     for (ou = 0; ou < n_out; ou++)
//     {
//         ou_sp_info = f_sparse_out(ou);
//         ou_index[ou] = ou_sp_info;
//         ou_colptr_info[ou] = ou_sp_info + 2;
//         ou_row_info[ou] = ou_colptr_info[ou] + ou_index[ou][1] + 1;
//         nnz_ou[ou] = *(ou_row_info[ou] - 1);
//         cout << "NNZ = " << nnz_ou[ou] << endl;
//     }

//     for (i_out = 0; i_out < n_out; i_out++)
//     {
//         i_nz = 0;
//         mat[i_out]->reserve(nnz_ou[i_out]);
//         for (col = 0; col < ou_index[i_out][1]; col++)
//         {
//             col_temp = col;
//             while ((col_temp+1< ou_index[i_out][1]) && (ou_colptr_info[i_out][col_temp+1] == 0))
//             {
//                 col_temp = col_temp + 1;
//             }
            
            
//             for (row = 0;i_nz < ou_colptr_info[i_out][col_temp+1] ; ++i_nz, ++row)
//             {
//                 mat[i_out]->insert(ou_row_info[i_out][i_nz], col) = 0;                         
//             }
//         }
//         mat[i_out]->makeCompressed();
//     }
//     // cout << "H = " << "\n" << mat[0] << endl;
// }

// void Initialize_Mat(SpMat& H, SpMat& tau)
// {
//     const int n_out = 2;
//     size_t ou, i_out, col, col_temp;
//     size_t row, i_nz = 0;

//     const casadi_int* ou_sp_info;
//     const casadi_int* ou_index[n_out];
//     const casadi_int* ou_colptr_info[n_out];
//     const casadi_int* ou_row_info[n_out];
//     casadi_int nnz_ou[n_out];

//     for (ou = 0; ou < n_out; ou++)
//     {
//         ou_sp_info = Dyn_FL_sub_sparsity_out(ou);
//         ou_index[ou] = ou_sp_info;
//         ou_colptr_info[ou] = ou_sp_info + 2;
//         ou_row_info[ou] = ou_colptr_info[ou] + ou_index[ou][1] + 1;
//         nnz_ou[ou] = *(ou_row_info[ou] - 1);
//     }

//     H.reserve(nnz_ou[0]);
//     tau.reserve(nnz_ou[1]);
//     for (i_out = 0; i_out < n_out; i_out++)
//     {
//         i_nz = 0;
//         for (col = 0; col < ou_index[i_out][1]; col++)
//         {
//             col_temp = col;
//             while ((col_temp+1< ou_index[i_out][1]) && (ou_colptr_info[i_out][col_temp+1] == 0))
//             {
//                 col_temp = col_temp + 1;
//             }
            
            
//             for (row = 0;i_nz < ou_colptr_info[i_out][col_temp+1] ; ++i_nz, ++row)
//             {
//                 if (i_out == 0)
//                 {
//                     H.insert(ou_row_info[i_out][i_nz], col) = 0;
//                 }
//                 else
//                 {
//                     tau.insert(ou_row_info[i_out][i_nz], col) = 0;
//                 }                
                
//             }
//         }
//     }
//     H.makeCompressed();
//     tau.makeCompressed();
// }

// void Read_Dyn_Info(StateVec&x, ContrlVec&u, SpMat&H, SpMat&tau)
// {
//     double *arg[2], *res[2], *w;
//     casadi_int *iw;
//     casadi_int sz_arg, sz_res, sz_iw, sz_w;
//     Dyn_FL_sub_work(&sz_arg, &sz_res, &sz_iw, &sz_w);

//     arg[0] = x.data();
//     arg[1] = u.data();
//     res[0] = H.valuePtr();
//     res[1] = tau.valuePtr();
//     iw = new casadi_int[sz_iw];
//     w = new double[sz_w];

        
//     Dyn_FL_sub((const double**)arg, res, iw, w, 1);

//     delete[] iw;
//     delete[] w;

// }

// void test()
// {
//     const int a[5] = {1,2,3,4,5};
//     const char* b = "Hello";
//     cout << "a = " << a << endl;
//     cout << "b = " << b << endl;
// }

// void Read_Dyn_Info_Deprecated(StateVec&x, ContrlVec&u, SpMat&H, SpMat&tau)
// {
//     /*This function is deprecated */
//     double **arg;
//     double **res,  *w;
//     int_T *iw;
//     int_T sz_arg, sz_res, sz_iw, sz_w;
//     const size_t n_in=2, n_out=2;  
//     size_t ou;

//     Dyn_FL_sub_work(&sz_arg, &sz_res, &sz_iw, &sz_w);
   
//     const casadi_int* ou_sp_info;
//     const casadi_int* ou_index[n_out];
//     const casadi_int* ou_colptr_info[n_out];
//     const casadi_int* ou_row_info[n_out];
//     casadi_int nnz_ou[n_out];
   
//     arg = new double *[n_in];
//     res = new double *[n_out];
//     iw = new int_T[sz_iw];
//     w = new double[sz_w];

//     arg[0] = x.data();
//     arg[1] = u.data();


//     for (ou = 0; ou < n_out; ou++)
//     {
//         ou_sp_info = Dyn_FL_sub_sparsity_out(ou);
//         ou_index[ou] = ou_sp_info;
//         ou_colptr_info[ou] = ou_sp_info + 2;
//         ou_row_info[ou] = ou_colptr_info[ou] + ou_index[ou][1] + 1;
//         nnz_ou[ou] = *(ou_row_info[ou] - 1);

//         res[ou] = new double[nnz_ou[ou]];
//     }

//     Dyn_FL_sub((const double**)arg, res, iw, w, 1);
   
//     size_t col, col_temp, row, i_nz = 0, i_out;
//     H.reserve(nnz_ou[0]);
//     tau.reserve(nnz_ou[1]);
//     for (i_out = 0; i_out < n_out; i_out++)
//     {
//         i_nz = 0;
//         for (col = 0; col < ou_index[i_out][1]; col++)
//         {
//             col_temp = col;
//             while ((col_temp+1< ou_index[i_out][1]) && (ou_colptr_info[i_out][col_temp+1] == 0))
//             {
//                 col_temp = col_temp + 1;
//             }
            
            
//             for (row = 0;i_nz < ou_colptr_info[i_out][col_temp+1] ; ++i_nz, ++row)
//             {
//                 if (i_out == 0)
//                 {
//                     H.insert(ou_row_info[i_out][i_nz], col) = res[i_out][ou_colptr_info[i_out][col]+ row];
//                 }
//                 else
//                 {
//                     tau.insert(ou_row_info[i_out][i_nz], col) = res[i_out][ou_colptr_info[i_out][col]+ row];
//                 }                
                
//             }
//         }
//     }

   
// }

// void Forward_dyn(StateVec &x, ContrlVec &u, StateVec &x_next, OutputVec &y)
// {
//     double **arg;
//     double **res,  *w;
//     int_T *iw;
//     int_T sz_arg, sz_res, sz_iw, sz_w;
//     size_t n_in, n_out;  

//     n_in = Dyn_FL_n_in();
//     n_out = Dyn_FL_n_out();
//     Dyn_FL_work(&sz_arg, &sz_res, &sz_iw, &sz_w);

//     arg = new double *[n_in];
//     res = new double *[n_out];
//     iw = new int_T[sz_iw];
//     w = new double[sz_w];

//     arg[0] = x.data();
//     arg[1] = u.data();
//     res[0] = x_next.data();
//     res[1] = y.data();
    

//     Dyn_FL((const double**) arg, res, iw, w, 0);
//     delete[] arg;
//     delete[] res;
// }

// void Jacobian(StateVec&x, SpMat&J, SpMat&Jd)
// {
//     double **arg, **res, *w;
//     int_T *iw;
//     int_T sz_arg, sz_res, sz_iw, sz_w;
//     size_t n_in, n_out, ou;  

//     n_in = Jacob_F_n_in();
//     n_out = Jacob_F_n_out();
//     Jacob_F_work(&sz_arg, &sz_res, &sz_iw, &sz_w);
//     const casadi_int* ou_sp_info;
//     const casadi_int** ou_index = new const casadi_int*[n_out];
//     const casadi_int** ou_colptr_info = new const casadi_int*[n_out];
//     const casadi_int** ou_row_info = new const casadi_int*[n_out];
//     casadi_int *nnz_ou = new casadi_int [n_out]; // nnz of output variables

//     arg = new double *[n_in];
//     res = new double *[n_out];
//     iw = new int_T[sz_iw];
//     w = new double[sz_w];

//     arg[0] = x.data();
    
//     for (ou = 0; ou < n_out; ou++)
//     {
//         ou_sp_info = Jacob_F_sparsity_out(ou);
//         ou_index[ou] = ou_sp_info;
//         ou_colptr_info[ou] = ou_sp_info + 2;
//         ou_row_info[ou] = ou_colptr_info[ou] + ou_index[ou][1] + 1;
//         nnz_ou[ou] = *(ou_row_info[ou] - 1);

//         res[ou] = new double[nnz_ou[ou]];
//     }
//     Jacob_F((const double**)arg, res, iw, w, 0);

//     vector<T> *triplist = new vector<T>[n_out];
    
//     size_t col, col_temp, row, i_nz, i_out;
    
//     for (i_out = 0; i_out < n_out; i_out++)
//     {
//         i_nz = 0;
//         for (col = 0; col < ou_index[i_out][1]; col++)
//         {
//             col_temp = col;
//             while ((col_temp+1< ou_index[i_out][1]) && (ou_colptr_info[i_out][col_temp+1] == 0))
//             {
//                 col_temp = col_temp + 1;
//             }
            
            
//             for (row = 0;i_nz < ou_colptr_info[i_out][col_temp+1] ; ++i_nz, ++row)
//             {
//                 triplist[i_out].push_back(T(ou_row_info[i_out][i_nz], col, res[i_out][ou_colptr_info[i_out][col]+ row]));
//             }
//         }
//     }
//     J.setFromTriplets(triplist[0].begin(), triplist[0].end());
//     Jd.setFromTriplets(triplist[1].begin(), triplist[1].end());
// }


// void Forward_dyn_par(StateVec &x, ContrlVec &u, SpMat&fx, SpMat&fu, SpMat&gx, SpMat&gu)
// {
//     double **arg, **res, *w;
//     int_T *iw;
//     size_t n_in, n_out, ou;  
//     int_T sz_arg, sz_res, sz_iw, sz_w;

//     n_in = Dyn_FL_par_n_in();
//     n_out = Dyn_FL_par_n_out();
//     Dyn_FL_par_work(&sz_arg, &sz_res, &sz_iw, &sz_w);    

//     const casadi_int* ou_sp_info;
//     const casadi_int** ou_index = new const casadi_int*[n_out];
//     const casadi_int** ou_colptr_info = new const casadi_int*[n_out];
//     const casadi_int** ou_row_info = new const casadi_int*[n_out];
//     casadi_int *nnz_ou = new casadi_int [n_out]; // nnz of output variables

//     arg = new double *[n_in];
//     res = new double *[n_out];
//     iw = new int_T[sz_iw];
//     w = new double[sz_w];
//     arg[0] = x.data();
//     arg[1] = u.data();

//     std::cout << "Number of output: " << n_out<<endl;
//     std::cout << "Number of input: " << n_in<<endl;

//     for (ou = 0; ou < n_out; ou++)
//     {
//         ou_sp_info = Dyn_FL_par_sparsity_out(ou);
//         ou_index[ou] = ou_sp_info;
//         ou_colptr_info[ou] = ou_sp_info + 2;
//         ou_row_info[ou] = ou_colptr_info[ou] + ou_index[ou][1] + 1;
//         nnz_ou[ou] = *(ou_row_info[ou] - 1);

//         res[ou] = new double[nnz_ou[ou]];
//     }
       
//     Dyn_FL_par((const double**)arg, res, iw, w, 0);

//     vector<T> *triplist = new vector<T>[n_out];
    
//     size_t col, col_temp, row, i_nz = 0, i_out;
    
//     for (i_out = 0; i_out < n_out; i_out++)
//     {
//         i_nz = 0;
//         for (col = 0; col < ou_index[i_out][1]; col++)
//         {
//             col_temp = col;
//             while ((col_temp+1< ou_index[i_out][1]) && (ou_colptr_info[i_out][col_temp+1] == 0))
//             {
//                 col_temp = col_temp + 1;
//             }
            
            
//             for (row = 0;i_nz < ou_colptr_info[i_out][col_temp+1] ; ++i_nz, ++row)
//             {
//                 triplist[i_out].push_back(T(ou_row_info[i_out][i_nz], col, res[i_out][ou_colptr_info[i_out][col]+ row]));
//             }
//         }
//     }
    
//     fx.setFromTriplets(triplist[0].begin(), triplist[0].end());
//     fu.setFromTriplets(triplist[1].begin(), triplist[1].end());
//     gx.setFromTriplets(triplist[2].begin(), triplist[2].end());
//     gu.setFromTriplets(triplist[3].begin(), triplist[3].end());
    
//     for (ou = 0; ou < n_out; ou++)
//     {
//         delete[] res[ou];
//     }
//     delete[] res;
//     delete[] arg;
//     delete[] iw;
//     delete[] w;
//     delete[] nnz_ou;
//     delete[] ou_index;
//     delete[] ou_row_info;
//     delete[] ou_colptr_info;
//     delete[] triplist;
// }

