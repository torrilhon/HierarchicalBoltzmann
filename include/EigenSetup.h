/*************************************************************************\
 EigenSetup.h  - declares namespace, some routines and defines relevant 
                 data types. 
                 Note: All double vectors in the code are represented 
                 as one column matrices.
 This is based on the Eigen template library, see
 http://eigen.tuxfamily.org/
\*************************************************************************/
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include "defs.h"

using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXi;
using Eigen::SparseMatrix;
typedef SparseMatrix<double,Eigen::RowMajor> SpMatrix;
typedef SparseMatrix<double,Eigen::RowMajor,INT> SpMatrixL;

void SpInit( int n, int m, const VectorXi& Ai, const VectorXi& Aj, const MatrixXd& Aval, SpMatrix& A );
void SpBlock( int idx, const MatrixXd& P, SpMatrix& T );
