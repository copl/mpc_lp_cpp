#ifndef INTERIOR_POINT_DATA_STRUCTURES
#define INTERIOR_POINT_DATA_STRUCTURES

#include "copl_linalg.h"


namespace copl_ip{

// preliminary definitions //
class lp_direction;

//lp_input
class lp_input {
public:
	int m, n, k;
	
	copl_matrix A, G;
	
	copl_vector c, h, b;

	lp_input(int m, int n, int k);
};
//--------End lp_input--------

//lp_variables
class lp_variables {
public:
	copl_vector x;
	copl_vector s;
	copl_vector z;
	copl_vector y;

	float tau;
	float kappa;
	
	lp_variables(lp_input problem_data);
	void take_step(lp_direction dir);
};
//--------End lp_variables--------

//k_newton_copl_matrix
class k_newton_copl_matrix {
	public:
		k_newton_copl_matrix(lp_input problem_data);
		void update(lp_variables variables);

};

//--------End k_newton_copl_matrix--------



//algorithm_state
class algorithm_state {

	float mu;
	float sigma;
	float gap;

public:
	algorithm_state();
	void update_mu (lp_variables variables, lp_input progblem_data); //TODO
	void update_gap (lp_variables variables, lp_input problem_data); //TODO
};


//--------End algorithm_state--------
	
//settings
class lp_settings {
	int max_iter;
	float linear_feas_tol;	//This is a relative tolerance w.r.t. 
                            //some normalizing norms
	float comp_tol; // How small must s^Tz must be when we stop
	
	//Constant length of the 
    //maximum combined step to the boundary to use
    float bkscale;
	
	//Configuration for solver
	//linear_solver_lp_settings 
public:
	lp_settings(
		int max_iter,
		float linear_feas_tol,	
		float comp_tol,	
		float bkscale
		);
	int get_max_iter();
	float get_linear_feas_tol();
	float get_comp_tol();
	float get_bkscale();
};

//--------End settings--------

//residuals

class residuals {
public:
	copl_vector r1;
	copl_vector r2;
	copl_vector r3;
	copl_vector r4;

	float r1_norm;
	float r2_norm;
	float r3_norm;
	float normed_squared;
	
	residuals();
	void update_values();
	void compute_residuals();
};
//--------End residuals--------

//lp_residuals

class lp_residuals {
public:
	copl_vector r1;
	copl_vector r2;
	copl_vector r3;
	float r4;

	float r1_norm;
	float r2_norm;
	float r3_norm;
	float normed_squared;
	
	lp_residuals(lp_input problem_data);
	void compute_residuals(lp_input in_problem_data, lp_variables variables);

	float get_r1_norm();
	float get_r2_norm();
	float get_r3_norm();
};
//--------End lp_residuals--------


//linear_system_rhs
class linear_system_rhs {
	copl_vector q1;
	copl_vector q2;
	copl_vector q3;
	copl_vector q4;
	copl_vector q5;
	copl_vector q6;
	void update_values(
		copl_vector q1, 
		copl_vector q2, 
		copl_vector q3, 
		copl_vector q4, 
		copl_vector q5, 
		copl_vector q6
	); //Not called in algorithm
public:
	linear_system_rhs(lp_input problem_data);
	
	void compute_affine_rhs(lp_residuals residuals, lp_variables variables);
	void compute_corrector_rhs(lp_residuals residuals, lp_variables variables, algorithm_state state, lp_direction direction, lp_input problem_data);

};

//--------End linear_system_rhs--------


//lp_direction
class lp_direction {
public:
	copl_vector dx;
	copl_vector dy;
	copl_vector dz;
	copl_vector ds;
	float dtau;
	float dkappa;
	float alpha;
	
	lp_direction(lp_variables variables);
	void update_values(
		copl_vector dx,
		copl_vector dy,
		copl_vector dz,
		copl_vector ds, 
		float dtau,
		float dkappa,
		float alpha
		     );
	void compute_affine_direction(	
		linear_system_rhs affine_rhs,
		lp_input problem_data,
		lp_variables variables,
		k_newton_copl_matrix K_newton_copl_matrix
		);
	void compute_corrector_direction(
		linear_system_rhs corrector_rhs,
		lp_input problem_data,
		lp_variables variables,
		algorithm_state state,
		lp_settings settings,
		k_newton_copl_matrix K_newton_copl_matrix
		);
	void compute_alpha(
		algorithm_state state,
		lp_settings settings
		);
	void compute_min_ratio_alpha (
		copl_vector var,
		copl_vector dvar
		);
	float get_alpha();
	float get_dtau(); 
	float get_dkappa(); 
	copl_vector get_dx(); 
	copl_vector get_dy(); 
	copl_vector get_dz(); 
	copl_vector get_ds(); 

};
//--------End lp_direction--------

//lp_result
class lp_result {

};
//--------End lp_result--------
}


#endif
