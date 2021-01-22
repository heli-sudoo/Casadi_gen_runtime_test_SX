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
void initialize_mat(const int_T *f_sparse_out(int_T), vector<SpMat*>&);
void initialization();
void read_dyn_info(StateVec&, ContrlVec&, SpMat&, SpMat&);
void dynamics(StateVec&, ContrlVec&, StateVec&, OutputVec&, int);
void dynamics_par(StateVec&, ContrlVec&, SpMat&, SpMat&, SpMat&, SpMat&, int);
void resetmap(StateVec&, StateVec&, OutputVec&, int);
void resetmap_par(StateVec&, SpMat&, SpMat&, SpMat&, SpMat&, int);
void footJacobian(StateVec&, SpMat&, SpMat&, int);
void sparseIdentityMat(SpMat&);

void dynamics_old(StateVec&, ContrlVec&, StateVec&, OutputVec&);
void dynamics_par_old(StateVec&, ContrlVec&, SpMat&, SpMat&, SpMat&, SpMat&);
void Jacobian_old(StateVec&, SpMat&, SpMat&);
void initialize_mat_old(SpMat&, SpMat&);

const int qsize = 7, xsize = 2*qsize, usize = 4, ysize = 2;
float dt;

StateVec x, x_next, fc;
ContrlVec u;
OutputVec y;
CoorVec q, qd, qd_next_f;

// StateTransMat fx;
// ContrlMat fu;

SpMat fx(xsize, xsize);  // dynamics partial w.r.t. x in discrete time
SpMat fu(xsize, usize);  // dynamics partial w.r.t. u in discrete time
SpMat fcx(xsize, xsize); // dynamics partial w.r.t. x in continuous time
SpMat fcu(xsize, usize); // dynamics partial w.r.t. u in continuous time
SpMat Px(xsize, xsize);  // resetmap partial
SpMat Pu(xsize, usize);
SpMat gx(ysize, xsize);
SpMat gu(ysize, usize);
SpMat J(ysize, qsize);
SpMat Jd(ysize, qsize);
SpMat H(qsize, qsize);
SpMat tau(qsize, 1);
SpMat SpI14(xsize, xsize);
SpMat SpI7(qsize, qsize);
SpMat SpI4(usize, usize);
Eigen::BiCGSTAB<SpMat> solver;


#endif