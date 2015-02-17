/*
 * Linear algebra definitions for the copl interior point solver 
 * Santiago Akle 
 */
#ifndef COPL_LINALG_H
#define COPL_LINALG_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <random>

namespace copl_ip 
{
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
	insert_at(int m, int n, double val)
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

//k_newton_copl_matrix
class k_newton_copl_matrix {
    private: 
    	//The assembled eigen matrix 
    	EigenSpMat_t* eigenKMat;	
		//Indices in the permuted matrix that correspond
	    //to the entries of the Hessian as read columnwise
		vector<int> permuted_indices;
		//Called from the constructor to assemble the first
		//version of the matrix
		void assemble(lp_input problem_definition);

		//Calls the symbolic analysis function
		void permute();
			
	public:
		k_newton_copl_matrix(lp_input problem_data);
		~k_newton_copl_matrix();
		void update(lp_variables variables);
		void factor();
	  	void solve(vector<double> &solution, vector<double> rhs);

};

}
#endif
