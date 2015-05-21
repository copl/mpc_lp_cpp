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
#include <copl_debug.h>
using namespace std;

namespace copl_ip 
{

// #ifdef DEBUG_TO_FILE
// extern ofstream LOG_FILE_VARIABLE;
// #define OUTPUT LOG_FILE_VARIABLE
// #else
// #define OUTPUT cout
// #endif
typedef Eigen::SparseMatrix<double> copl_matrix;
typedef Eigen::VectorXd copl_vector;

}
#endif
