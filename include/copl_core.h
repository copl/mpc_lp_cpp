#ifndef INTERIOR_POINT_DATA_STRUCTURES
#define PREVENT_COPY_CONSTRUCTOR
#define INTERIOR_POINT_DATA_STRUCTURES


#include <copl_linalg.h>
#include <iostream>
#include <gtest/gtest_prod.h>


namespace copl_ip{

// preliminary definitions //
class lp_direction;
class k_newton_copl_matrix;

//lp_input
class lp_input {
public:
	int m, n, k_var;
	copl_matrix &A, &G;
	copl_external_vector &c, &h, &b;

public:
	lp_input(copl_matrix &_A, copl_external_vector & _b, copl_external_vector  &c, copl_matrix &G, copl_external_vector &h);	
	void var_dump()  ;
};
//--------End lp_input--------

//lp_variables
class lp_variables {
public:
	copl_external_vector &x;
	copl_external_vector &s; 
	copl_external_vector &z; 
	copl_external_vector &y;
	double &tau;
	double &kappa;
	lp_variables(copl_external_vector &_x,
		     copl_external_vector &_y,
		     copl_external_vector &_s,
		     copl_external_vector &_z, double &_tau, double &_kappa);
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
	public:
		int max_iter;
		double linear_feas_tol;	//This is a relative tolerance w.r.t. 
								//some normalizing norms
		double comp_tol; // How small must s^Tz must be when we stop
		
		//ant length of the 
		//maximum combined step to the boundary to use
		double bkscale;
		//Regularization for the SQD matrix 
		double regularization;

	lp_settings(
		int max_iter,
		double linear_feas_tol,	
		double comp_tol,	
		double bkscale,
		double regularization
		);
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

    //Variables for the norms 
    double hn1, hn2, hn3; //Homogeneous residual norms A'y+G'z, Ax, Gx+s
    double n1,  n2,  n3, n4;  //Residual norms A'y+G'z + tc , Ax-tb , Gx+s-th

	void var_dump();
};
//--------End lp_residuals--------

//linear_system_rhs
class linear_system_rhs {
    //Use one vector for [q1, q2, q3] and q4
	copl_vector q123;
	double q4;
	copl_vector q5;
	double q6;
    int m,n,k;
    friend class homogeneous_solver;
    FRIEND_TEST(KNEWTON,reduce_rhs_test);
    FRIEND_TEST(KNEWTON,back_substitute_test);
    FRIEND_TEST(HOMOGENEOUS_SOLVER,solve_reduced);
    FRIEND_TEST(HOMOGENEOUS_SOLVER,solve_fixed_rhs);
    FRIEND_TEST(HOMOGENEOUS_SOLVER,solve_full);
    FRIEND_TEST(CORE,affine_rhs);
    FRIEND_TEST(CORE, residual_test);

public:
	linear_system_rhs( lp_input &problem_data);
	linear_system_rhs(const linear_system_rhs &obj);
	~linear_system_rhs();
	
	void compute_affine_rhs(lp_residuals &residuals, lp_variables &variables);
	void compute_combined_rhs(lp_residuals &residuals, lp_variables &variables, 
                                                       lp_direction &dir_affine,
                                                       double sigma,
                                                       double mu);
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

	void compute_direction(
		linear_system_rhs &affine_rhs,
		lp_input &problem_data,
		lp_variables &variables,
		algorithm_state &state,
		lp_settings &settings,
		k_newton_copl_matrix &K_matrix
		);
		
	void solve_linear_system_for_new_direction(linear_system_rhs& rhs, k_newton_copl_matrix& K_matrix);	
	void compute_step_size(lp_variables& variables, lp_settings& settings);
};

//--------End lp_direction--------

//lp_result
class lp_result {

};


class copl_utility {
      public:
          // Load problem from "The University of Florida Sparse Matrix Collection"
          // the object would be created and stored in "problem_data"
          static lp_input* loadFromUF(string UF_group, string UF_name);
		  
};

//--------End lp_result--------
}

#endif
