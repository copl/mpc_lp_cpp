#include <copl_core.h>
#include <copl_linalg.h>
#include <iostream>
#include <gtest/gtest_prod.h>

using namespace std;

namespace copl_ip {

//This class is used by the optimization code 
//to solve the linear systems of the iterations
class k_newton_copl_matrix {
        
	public:        
        //TODO: Move this to settings
        double DELTA = 1.e-4;
        
        bool isFactored = false;	
      
        //Hessian entries 
        std::vector<int> *hessianIx;	
     
        //Set up the eigen solver object
    	Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::AMDOrdering<int> >   solver;
               //Private empty constructor for testing

        //The assembled eigen matrix 
        EigenSpMat_t* eigenKMat;    
  
        //TODO: move these two to private and  
        //This function assembles the K newton matrix with identities in the diagonals
        //K = [-I A' G']
        //    [A -I    ]
        //    [G     I ]
        void assemble_matrix(copl_matrix &A, copl_matrix &G);          
        k_newton_copl_matrix(int m, int n);

        int nnz();
		k_newton_copl_matrix(copl_matrix& A, copl_matrix& G);
	  	void solve(copl_vector &solution, copl_vector &rhs);
	  	void update(lp_variables &variables);
	  	~k_newton_copl_matrix();
};

}
