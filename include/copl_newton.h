#include <copl_core.h>
#include <copl_linalg.h>
#include <iostream>
#include <gtest/gtest_prod.h>

using namespace std;

namespace copl_ip {

//This class is used by the optimization code 
//to solve the linear systems of the iterations
class k_newton_copl_matrix {
        
	protected:        
        //TODO: Move this to settings
        double DELTA = 1.e-4; 
        bool isFactored = false;	

        int m,n,p;      
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

        //Friend tests
        FRIEND_TEST(KNEWTON,Assemble);
        FRIEND_TEST(KNEWTON,NonZeroPerCols);
        FRIEND_TEST(KNEWTON,nnz);
        FRIEND_TEST(KNEWTON,Constructor);
        FRIEND_TEST(KNEWTON,Update);
        FRIEND_TEST(KNEWTON,solve);

    public: 
        int nnz();
		k_newton_copl_matrix(copl_matrix& A, copl_matrix& G);
	  	void solve(copl_vector &solution, copl_vector &rhs);
	  	void update(lp_variables &variables);
	  	~k_newton_copl_matrix();
};

//Extends the k_newton matrix and implements the methods to solve homogeneous systems.
class homogeneous_solver : protected k_newton_copl_matrix {
    
    copl_vector &_c, &_h, &_b;
    copl_vector rhs_2;
    copl_vector sol_1;
    copl_vector sol_2;

    //The homogeneous system is 
    //[0   A' G' c]
    //[A        -b]
    //[G        -h]
    //[-c' b' h   ]
    homogeneous_solver(copl_matrix &A, 
                       copl_matrix &G,
                       copl_vector &c,
                       copl_vector &b,
                       copl_vector &h); 

    void update(lp_variables &variables);

    void solve(lp_direction &dir, linear_system_rhs& rhs);

};

}
