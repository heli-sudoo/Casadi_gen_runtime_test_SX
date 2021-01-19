#ifndef DYNAMICS_DDP
#define DYNAMICS_DDP

#include <vector>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>
#include <eigen3/unsupported/Eigen/SparseExtra>
#include <eigen3/Eigen/IterativeLinearSolvers>

#define int_T long long int

using namespace std;
using namespace std::chrono;

using Eigen::MatrixXd;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;
typedef Eigen::Matrix<double, 14,14> StateTransMat;
typedef Eigen::Matrix<double, 14, 4> ContrlMat;
typedef Eigen::Matrix<double, 2, 14> ObservMat;
typedef Eigen::Matrix<double, 2, 4> DirectMat; 
typedef Eigen::Matrix<double, 14,1> StateVec;
typedef Eigen::Matrix<double, 2, 7> JacobMat;
typedef Eigen::Matrix4d ContrlVec;
typedef Eigen::Vector2d OutputVec;
typedef Eigen::Matrix<double, 7, 1> CoorVec;

void test();
void Initialize_Mat(const int_T *f_sparse_out(int_T), vector<SpMat*>&);
void Initialization();
void Read_Dyn_Info(StateVec&, ContrlVec&, SpMat&, SpMat&);
void Forward_dyn(StateVec&, ContrlVec&, StateVec&, OutputVec&);
void Forward_dyn_par(StateVec&, ContrlVec&, SpMat&, SpMat&, SpMat&, SpMat&);


void Forward_dyn_old(StateVec&, ContrlVec&, StateVec&, OutputVec&);
void Forward_dyn_par_old(StateVec&, ContrlVec&, SpMat&, SpMat&, SpMat&, SpMat&);
void Jacobian_old(StateVec&, SpMat&, SpMat&);
void Initialize_Mat_old(SpMat&, SpMat&);

const int qsize = 7, xsize = 2*qsize, usize = 4, ysize = 2;

StateVec x, x_next;
ContrlVec u;
OutputVec y;
CoorVec q, qd_next_f;

SpMat fx(xsize, xsize);
SpMat fu(xsize, usize);
SpMat gx(ysize, xsize);
SpMat gu(ysize, usize);
SpMat JF(ysize, qsize);
SpMat JFd(ysize, qsize);
SpMat JB(ysize, qsize);
SpMat JBd(ysize, qsize);
SpMat H(qsize, qsize);
SpMat tau(qsize, 1);
Eigen::BiCGSTAB<SpMat> solver;


#endif