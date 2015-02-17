#ifndef COPL_NEWTON_H
#define COPL_NEWTON_H

#include <copl_core.h>
#include <copl_linalg.h>

//k_newton_copl_matrix
class k_newton_copl_matrix {
	
    private: 
    	//The assembled eigen matrix 
    	EigenSpMat_t* eigenKMat;	
		//Indices in the permuted matrix that correspond
	    //to the entries of the Hessian as read columnwise
		std::vector<int> permuted_indices;
		//Called from the constructor to assemble the first
		//version of the matrix
		void assemble(copl_matrix A, copl_matrix G, int m, int n, int p);
		//Calls the symbolic analysis function
		void permute();
			
	public:
		k_newton_copl_matrix(lp_input problem_data);
		~k_newton_copl_matrix();
		void update(lp_variables variables);
		void factor();
	  	void solve(vector<double> &solution, vector<double> rhs);

};

#endif