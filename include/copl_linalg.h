/*
 * Linear algebra definitions for the copl interior point solver 
 * Santiago Akle 
 */
#define DEBUG_TO_FILE
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
typedef Eigen::SparseMatrix<double> EigenSpMat_t;
typedef std::vector<double> copl_vector;

//Wrapper for our matrix classes
class copl_matrix 
{
public: 
	EigenSpMat_t* eigenMat;	

	copl_matrix(int m, int n);
	
    //Random sparse matrix with dist p
	copl_matrix(int m, int n, double p);

    //Insert element to sparse matrix, this is slow and we should not use it
	void insert_at(int m, int n, double val);

    //set value of element in sparse matrix, this is slow and we should not use it
	//void set_value(int m, int n, double val);
	// number of rows
	int num_rows();
	// number of columns
	int num_cols();
	// value at
	double value_at(int i, int j);
	// var dump
	void var_dump();
	
	//Return the number of non-zero entires in the matrix	
	int nnz();
	//Destructor 
	~copl_matrix();

};


//These functions do not belong to either vector or matrices

//y<- alpha Ax + beta y with sparse A
void sp_dgemv(double alpha, double beta, copl_matrix &copl_A, copl_vector &copl_x, copl_vector &copl_y);

//Matrix vector multiply and accumulate in y
//y<- alpha A^Tx + beta y with sparse A (The input matrix is A not A^T)
void sp_dgemtv(double alpha, double beta, copl_matrix &copl_A, copl_vector &copl_x, copl_vector &copl_y);

//Scale the vector
//x<-alpha *x
void scal(copl_vector &copl_x, double alpha );

//y<- alpha*x + y
void axpy(double alpha, copl_vector &copl_x, copl_vector &copl_y);

// x^Ty
double dot(copl_vector &copl_y, copl_vector &copl_x);

//Zero out 
void zeros(copl_vector &y);

//Two norm 
double norm2(copl_vector &y);

//Infinity norm 	
double normInf(copl_vector &y);

}
#endif
