#include <copl_core.h>
#include <copl_linalg.h>
#include <iostream>
#include <gtest/gtest_prod.h>

using namespace std;

namespace copl_ip {

/*  
   Defines the methods used to solve 
   for the search directions.
*/

/*
  Used to solve with the symmetric  
  quasi definte system defined by 
  [  A' G'][dx]   
  [A      ][dy]
  [G   -H ][dz].
  This class regularizes the above system to generate 
  a strongly quasidefinite matrix that can be factored without 
  pivoting. 
  TODO: Use the factored matrix as a preconditioner for MINRES
  The present implementation is inefficient in two fronts, 
  the factorization calculates the permutation at each iteration and 
  uses an algorithm that will potentially pivot (super lu)
*/
class k_newton_copl_matrix {
        
	protected:        
        //TODO: Move this to settings
        double DELTA = 1.e-7; 
        bool isFactored = false;	

        //Problem size
        int m,n,k;      

        //Hessian entries 
        std::vector<int> *hessianIx;	
     
        //Set up the eigen solver object with AMD ordering
    	Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::AMDOrdering<int> >   solver;

        //The assembled eigen matrix 
        copl_matrix* eigenKMat;    
  
        //This function assembles the K newton matrix with identities in the diagonals
        //K = [-I A' G']
        //    [A -I    ]
        //    [G     I ]
        void assemble_matrix(copl_matrix &A, copl_matrix &G);          
       	//This constructor is used in some testing routines 
	k_newton_copl_matrix(int m, int n);
        
	//Tests
        FRIEND_TEST(KNEWTON,Assemble);
        FRIEND_TEST(KNEWTON,NonZeroPerCols);
        FRIEND_TEST(KNEWTON,nnz);
        FRIEND_TEST(KNEWTON,Constructor);
        FRIEND_TEST(KNEWTON,Update);
        FRIEND_TEST(KNEWTON,solve);

    public: 
	//Number of nonzeros
        int nnz();
	//Constructor
	k_newton_copl_matrix(copl_matrix& A, copl_matrix& G);
	//Solve with the symmetric system
	void solve(copl_vector &solution, copl_vector &rhs);
	//Update the variables and recalculate the hessian block
	void update(lp_variables &variables);
        //Destructor	
	~k_newton_copl_matrix();
};

/* Extends the k_newton matrix and implements methods to solve 
 * the system of the simplfied homogeneous self dual search directions (homogeneous system)
 
 * [0   A'   G' c]dz         q1
 * [-A          b]dy        =q2
 * [-G          h]dz  -ds    q3
 * [-c' -b' -h   ]dt  -dk    q4
 *  H dz + ds                q5 
 *  kappa/tau dt + dk        q6
*/
class homogeneous_solver : protected k_newton_copl_matrix {
    
    protected:
        copl_vector &_c, &_h, &_b;
        copl_vector rhs_1;
        copl_vector sol_1;
        copl_vector sol_2;
        double tau, kappa, dtau_denom;
        void reduce_rhs(linear_system_rhs  &rhs);

	//Solves with the system
        //[0   A'   G'  c]dz      q1
        //[A           -b]dy    = q2
        //[G       -H  -h]dz      q3
        //[-c' -b' -h k/t]dt      q4 
        void solve_reduced(lp_direction &dir, linear_system_rhs &rhs);

	//Given dx,dy,dz,dt calculate ds,dk
        void back_substitute(lp_direction &dir, linear_system_rhs  &rhs, lp_variables &var);

        FRIEND_TEST(KNEWTON,reduce_rhs_test);
        FRIEND_TEST(KNEWTON,back_substitute_test);
	FRIEND_TEST(HOMOGENEOUS_SOLVER,solve);
	FRIEND_TEST(HOMOGENEOUS_SOLVER,solve_reduced);
	FRIEND_TEST(HOMOGENEOUS_SOLVER,solve_fixed_rhs);

    public: 
    	homogeneous_solver(lp_input &prob);
   	
    	//Updates the value of the Hessian block and the tau and kappa variables 
    	//in the matrix. Factors a matrix derived from the homogenos system
    	// and solves one system with the factored matrix.
    	void update(lp_variables &variables);
    	
    	//Solves the homogeneous system
    	void solve(lp_direction &dir, linear_system_rhs& rhs, lp_variables &var);

};

}
