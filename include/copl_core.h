#ifndef INTERIOR_POINT_DATA_STRUCTURES
#define PREVENT_COPY_CONSTRUCTOR
#define INTERIOR_POINT_DATA_STRUCTURES

#include <copl_linalg.h>
#include <iostream>
namespace copl_ip{

// preliminary definitions //
class lp_direction;
class k_newton_copl_matrix;

//lp_input
class lp_input {
public:
	int m, n, k_var;
	
	copl_matrix A, G;
	
	copl_vector c, h, b;

	lp_input(int m, int n, int k_var);
	lp_input(const lp_input &obj);
	~lp_input();
	void var_dump()  ;
};
//--------End lp_input--------

//lp_variables
class lp_variables {
public:
	copl_vector x, s, z, y;

	double tau, kappa;
	
	lp_variables(int n, int m, int k_var);
	lp_variables(const lp_variables &obj);
	void take_step(lp_direction &direction);
	~lp_variables();
};
//--------End lp_variables--------


//algorithm_state
class algorithm_state {
public:
	double mu;
	double sigma;
	double gap;


	algorithm_state();
	void update_mu (lp_variables &variables, lp_input &problem_data);
	void update_gap (lp_variables &variables,  lp_input &problem_data);
};


//--------End algorithm_state--------
	
//settings
class lp_settings {
	int max_iter;
	double linear_feas_tol;	//This is a relative tolerance w.r.t. 
                            //some normalizing norms
	double comp_tol; // How small must s^Tz must be when we stop
	
	//ant length of the 
    //maximum combined step to the boundary to use
    double bkscale;
public:
	//Configuration for solver
	//linear_solver_lp_settings 

	lp_settings(
		int max_iter,
		double linear_feas_tol,	
		double comp_tol,	
		double bkscale
		);
	lp_settings(const lp_settings &obj); //copy constructor - 
	int get_max_iter();
	double get_linear_feas_tol();
	double get_comp_tol();
	double get_bkscale();
};

//--------End settings--------


//lp_residuals

class lp_residuals {
public:
	copl_vector r1;
	copl_vector r2;
	copl_vector r3;
	double r4;
	
	lp_residuals( lp_input &problem_data);
	void compute_residuals( lp_input &problem_data, lp_variables &variables);

	double get_r1_norm();
	double get_r2_norm();
	double get_r3_norm();
	double get_norm_squared();
	
	void var_dump();
};
//--------End lp_residuals--------


//linear_system_rhs
class linear_system_rhs {
	copl_vector q1;
	copl_vector q2;
	copl_vector q3;
	double q4;
	copl_vector q5;
	double q6;
public:
	linear_system_rhs( lp_input &problem_data);
	linear_system_rhs(const linear_system_rhs &obj);
	~linear_system_rhs();
	
	void compute_affine_rhs(lp_residuals &residuals, lp_variables &variables);
	void compute_corrector_rhs(lp_residuals &residuals, lp_variables &variables, algorithm_state &state, lp_direction &direction,  lp_input &problem_data);

	void var_dump();
};

//--------End linear_system_rhs--------


//lp_direction
class lp_direction {
public:
	copl_vector dx;
	copl_vector dy;
	copl_vector dz;
	copl_vector ds;
	double dtau;
	double dkappa;
	double alpha;
	
	lp_direction(lp_variables &variables);
	lp_direction(const lp_direction &obj);

	void compute_affine_direction(	
		linear_system_rhs &affine_rhs,
		 lp_input &problem_data,
		lp_variables &variables,
		k_newton_copl_matrix &K_matrix
		);

	void compute_corrector_direction(
		linear_system_rhs &corrector_rhs,
		lp_input &problem_data,
		lp_variables &variables,
		algorithm_state &state,
		lp_settings &settings,
		k_newton_copl_matrix &K_matrix
		);
	
	void compute_min_ratio_alpha(copl_vector &var, copl_vector &dvar, double& alpha_val);
	void compute_min_ratio_alpha(double var, double dvar, double& alpha_val);
	
	double get_alpha();
	double get_dtau(); 
	double get_dkappa(); 
	copl_vector get_dx(); 
	copl_vector get_dy(); 
	copl_vector get_dz(); 
	copl_vector get_ds(); 

};

//k_newton_copl_matrix
class k_newton_copl_matrix {
    private: 
    	//The assembled eigen matrix 
    	EigenSpMat_t* eigenKMat;	
		//Indices in the permuted matrix that correspond
	    //to the entries of the Hessian as read columnwise
		std::vector<int> permuted_indices;
		//Called from the ructor to assemble the first
		//version of the matrix
		//void assemble(copl_matrix A, copl_matrix G, int m, int n, int p);
		//Calls the symbolic analysis function
		//void permute();
			
	public:
		k_newton_copl_matrix( lp_input &problem_data);
		//~k_newton_copl_matrix();
		//void update(lp_variables variables);
		void factor();
	  	void solve(copl_vector &solution, copl_vector &rhs);
	  	void update(lp_variables &variables);
};

//--------End lp_direction--------

//lp_result
class lp_result {

};


class copl_utility {
      public:
          // Load problem from "The University of Florida Sparse Matrix Collection"
          // the object would be created and stored in "problem_data"
		  static lp_input* Trivial_Test1();
          static lp_input* loadFromUF(string UF_group, string UF_name);
		  
};

//--------End lp_result--------
}

#endif
