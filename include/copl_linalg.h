#ifndef MATRIX_POINT_DATA_STRUCTURES
#define MATRIX_DATA_STRUCTURES

#include <Eigen/Core>
#include <Eigen/Sparse>

typedef Eigen::SparseMatrix<double> EigenSpMat_t;
typedef Eigen:Eigen::VectorXd EigenVector_t;
//These functions do not belong to either vector or matrices

//y<- alpha Ax + beta y with sparse A
void sp_dgemv(double alpha, double beta, copl_matrix A, copl_vector x, copl_vector y);

//Matrix vector multiply and accumulate in y
//y<- alpha A^Tx + beta y with sparse A (The input matrix is A not A^T)
void sp_dgemtv(double alpha, double beta, copl_matrix A, copl_vector x, copl_vector y);

//Scale the vector
//x<-alpha *x
void scal(copl_vector x, double alpha );

//y<- a*x + y
void axpy(double alpha, copl_vector x, copl_vector y);

// x^Ty
double dot(copl_vector y, copl_vector x);

//Zero out 
double zeros(copl_vector y);

//Wrapper for our matrix classes
class copl_matrix 
{
public: 
	EigenSpMat_t eigenMat;	

	//Constructs a matrix from the CSR format vectors
	copl_matrix(vector<int> &col_counts, vector<int> &row_ix, vector<double> &vals);

	//This returns a copl_matrix from 3 vectors or rows columns and values
	static copl_matrix from_coo_format(vector<int> &row_ix, vector<int> &col_ix, vector<double> &vals);
};

class copl_vector {
	EigenVector_t *vec;	

	//Construct from an STL vector	
	copl_vector(vector<double> vec);

	//Construct form a raw memory vector
	copl_vector(int n, double*);

	~copl_vector()
	{
		delete(vec;)			
	}
};

#endif
