/*
 * Linear algebra definitions for the copl interior point solver 
 * Santiago Akle 
 */
//#define DEBUG_TO_FILE
#ifndef COPL_LINALG_H
#define COPL_LINALG_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <random>
#include <iostream>
#include <fstream>
#include <time.h>
#include <copl_debug.h>
using namespace std;

namespace copl_ip 
{

typedef Eigen::SparseMatrix<double> copl_matrix;
typedef Eigen::VectorXd copl_vector;

}
#endif
