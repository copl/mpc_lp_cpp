#ifndef INTERIOR_POINT_DATA_STRUCTURES
#define INTERIOR_POINT_DATA_STRUCTURES

#include "copl_linalg.h"


namespace copl_ip{
//lp_settings
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
//--------End lp_settings--------

//lp_input
class lp_input {
	copl_ip::matrix A;
	copl_ip::matrix G;
	
	copl_ip::vector c;
	copl_ip::vector h;
	copl_ip::vector b;
	
	int m;
	int n;
	int k;
public:
	lp_input();
};
//--------End lp_input--------

//lp_residuals

class lp_residuals {
	copl_ip::vector r1;
	copl_ip::vector r2;
	copl_ip::vector r3;
	copl_ip::vector r4;

	float r1_norm;
	float r2_norm;
	float r3_norm;
	float normed_squared;
public:
	lp_residuals(lp_input problem_data);
	void update_values(
		copl_ip::vector r1,
		copl_ip::vector r2, 
		copl_ip::vector r3,
		copl_ip::vector r4
		);
	void compute_residuals(lp_input in_problem_data, lp_variables variables);

	float get_r1_norm();
	float get_r2_norm();
	float get_r3_norm();
};
//--------End lp_residuals--------


//linear_system_rhs

class linear_system_rhs {
	copl_ip::vector q1;
	copl_ip::vector q2;
	copl_ip::vector q3;
	copl_ip::vector q4;
	copl_ip::vector q5;
	copl_ip::vector q6;
	void update_values(
		copl_ip::vector q1, 
		copl_ip::vector q2, 
		copl_ip::vector q3, 
		copl_ip::vector q4, 
		copl_ip::vector q5, 
		copl_ip::vector q6
	); //Not called in algorithm
public:
	linear_system_rhs(lp_input problem_data);
	
	void compute_affine_rhs(lp_residuals residuals, lp_variables variables);
	void compute_corrector_rhs(lp_residuals residuals, lp_variables variables);

};

//--------End linear_system_rhs--------

//lp_direction
class lp_direction {
	copl_ip::vector dx;
	copl_ip::vector dy;
	copl_ip::vector dz;
	copl_ip::vector ds;
	float dtau;
	float dkappa;
	float alpha;
public:
	lp_direction();
	void update_values(
		copl_ip::vector dx,
		copl_ip::vector dy,
		copl_ip::vector dz,
		copl_ip::vector ds, 
		float dtau,
		float dkappa,
		float alpha
		     );
	void compute_affine_direction(	
		linear_system_rhs affine_rhs,
		lp_input problem_data,
		lp_variables variables,
		k_newton_matrix K_newton_matrix
		);
	void compute_corrector_direction(
		linear_system_rhs corrector_rhs,
		lp_input problem_data,
		lp_variables variables,
		algorithm_state state,
		lp_settings settings,
		k_newton_matrix K_newton_matrix
		);
	void compute_alpha(
		algorithm_state state,
		lp_settings settings
		);
	void compute_min_ratio_alpha (
		copl_ip::vector var,
		copl_ip::vector dvar
		);
	float get_alpha();
	float get_dtau(); 
	float get_dkappa(); 
	copl_ip::vector get_dx(); 
	copl_ip::vector get_dy(); 
	copl_ip::vector get_dz(); 
	copl_ip::vector get_ds(); 

};
//--------End lp_direction--------

//lp_variables
class lp_variables {
	copl_ip::vector x;
	copl_ip::vector s;
	copl_ip::vector z;
	copl_ip::vector y;

	float tau;
	float kappa;
public:
	lp_variables(lp_input problem_data);
	void take_step(lp_direction dir);
};
//--------End lp_variables--------

//k_newton_matrix
class k_newton_matrix {
	public:
		k_newton_matrix(lp_input input);
		void update(lp_variables variables);

};

//--------End k_newton_matrix--------

//lp_result
class lp_result {

};
//--------End lp_result--------



//algorithm_state
class algorithm_state {

	float mu;
	float sigma;
	float gap;

public:
	algorithm_state();
	void update_mu (lp_variables variables, lp_input problem_data); //TODO
	void update_gap (lp_variables variables, lp_input problem_data); //TODO
};


//--------End algorithm_state--------


}


#endif
