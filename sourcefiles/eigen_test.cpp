// #include <iostream>
// #include <Eigen/Dense>
// #include <Eigen/Sparse>
 
// using Eigen::MatrixXd;
// typedef Eigen::SparseMatrix<double> SpMat;
// typedef Eigen::Triplet<double> T;
 
// int main()
// {
//   MatrixXd DA(2,2);
//   DA(0,0) = 3;
//   DA(1,0) = 2.5;
//   DA(0,1) = -1;
//   DA(1,1) = DA(1,0) + DA(0,1);
//   std::cout<<"Dense matrix DA = " << std::endl;
//   std::cout <<DA << std::endl;

//   int m = 2, n = 2;
//   std::vector<T> tripletlist;
//   SpMat SA(m, n);
//   tripletlist.reserve(2);
//   for (int i = 0; i < m; i++)
//   {
//     for (int j = 0; j < n; j++)
//     {
//       tripletlist.push_back(T(i,j,i+j));
//     }
//   }
//   SA.setFromTriplets(tripletlist.begin(), tripletlist.end());
//   std::cout<<"Sparse matrix SA = "<<std::endl;
//   std::cout<<SA<<std::endl;
//   std::cout<<"DA*SA = "<<std::endl;
//   std::cout<< DA*SA <<std::endl;  
// }