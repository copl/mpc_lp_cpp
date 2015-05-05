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

        //Problem size
        int m,n,k;      

        //Hessian entries 
        std::vector<int> *hessianIx;	
     
        //Set up the eigen solver object
    	Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::AMDOrdering<int> >   solver;

        //The assembled eigen matrix 
        copl_matrix* eigenKMat;    
  
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
