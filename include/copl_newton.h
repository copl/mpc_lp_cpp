#include <copl_core.h>
#include <copl_linalg.h>
#include <iostream>
#include <gtest/gtest_prod.h>

using namespace std;

namespace copl_ip {

//This class is used by the optimization code 
//to solve the linear systems of the iterations
class k_newton_copl_matrix {
    private: 
        //TODO: Move this to settings
        double DELTA = 1.e-4;
        
        bool isFactored = false;	
        //The assembled eigen matrix 
    	EigenSpMat_t* eigenKMat;	
        
        //Hessian entries 
        std::vector<int> *hessianIx;	

        //We keep a reference for generating the matices 
    	//during the updates.
    	lp_input* pData; 	
       
        //Set up the eigen solver object
    	Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::AMDOrdering<int> >   solver;
        
        //This function assembles the K newton matrix with identities in the diagonals
        //K = [-I A' G']
		//    [A -I    ]
        //    [G     I ]
        void assemble_matrix(copl_matrix &A, copl_matrix &G);          

        //Private empty constructor for testing
        k_newton_copl_matrix();

        //Testing classes
        FRIEND_TEST(KNEWTON,Assemble);
	public:
		k_newton_copl_matrix( lp_input &problem_data);
		void factor();
	  	void solve(copl_vector &solution, copl_vector &rhs);
	  	void update(lp_variables &variables);
	  	~k_newton_copl_matrix();
};

}
