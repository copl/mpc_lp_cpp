#ifndef INTERIOR_POINT_DATA_STRUCTURES
#define INTERIOR_POINT_DATA_STRUCTURES

#include "copl_linalg.h"


namespace copl_ip{

// preliminary definitions //
class direction;

//lp_input
class lp_input {
	copl_ip::copl_matrix A;
	copl_ip::copl_matrix G;
	
	copl_ip::copl_vector c;
	copl_ip::copl_vector h;
	copl_ip::copl_vector b;
	
	int m;
	int n;
	int k;
public:
	lp_input();
};
//--------End lp_input--------

//lp_variables
class lp_variables {
	copl_ip::copl_vector x;
	copl_ip::copl_vector s;
	copl_ip::copl_vector z;
	copl_ip::copl_vector y;

	float tau;
	float kappa;
public:
	lp_variables(lp_input problem_data);
	void take_step(direction dir);
};
//--------End lp_variables--------

//k_newton_matrix
class k_newton_matrix {
	public:
		k_newton_matrix();
		void update(lp_variables variables);

};

//--------End k_newton_matrix--------



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
	copl_ip::copl_vector r1;
	copl_ip::copl_vector r2;
	copl_ip::copl_vector r3;
	copl_ip::copl_vector r4;

	float r1_norm;
	float r2_norm;
	float r3_norm;
	float normed_squared;
public:
	residuals();
	void update_values();
	void compute_residuals();
};
//--------End residuals--------

//lp_residuals

class lp_residuals {
	copl_ip::copl_vector r1;
	copl_ip::copl_vector r2;
	copl_ip::copl_vector r3;
	copl_ip::copl_vector r4;

	float r1_norm;
	float r2_norm;
	float r3_norm;
	float normed_squared;
public:
	lp_residuals(lp_input problem_data);
	void update_values(
		copl_ip::copl_vector r1,
		copl_ip::copl_vector r2, 
		copl_ip::copl_vector r3,
		copl_ip::copl_vector r4
		);
	void compute_residuals(lp_input in_problem_data, lp_variables variables);

	float get_r1_norm();
	float get_r2_norm();
	float get_r3_norm();
};
//--------End lp_residuals--------


//linear_system_rhs
class linear_system_rhs {
	copl_ip::copl_vector q1;
	copl_ip::copl_vector q2;
	copl_ip::copl_vector q3;
	copl_ip::copl_vector q4;
	copl_ip::copl_vector q5;
	copl_ip::copl_vector q6;
	void update_values(
		copl_ip::copl_vector q1, 
		copl_ip::copl_vector q2, 
		copl_ip::copl_vector q3, 
		copl_ip::copl_vector q4, 
		copl_ip::copl_vector q5, 
		copl_ip::copl_vector q6
	); //Not called in algorithm
public:
	linear_system_rhs(lp_input problem_data);
	
	void compute_affine_rhs(lp_residuals residuals, lp_variables variables);
	void compute_corrector_rhs(lp_residuals residuals, lp_variables variables);

};

//--------End linear_system_rhs--------


//lp_direction
class lp_direction {
	copl_ip::copl_vector dx;
	copl_ip::copl_vector dy;
	copl_ip::copl_vector dz;
	copl_ip::copl_vector ds;
	float dtau;
	float dkappa;
	float alpha;
public:
	lp_direction();
	void update_values(
		copl_ip::copl_vector dx,
		copl_ip::copl_vector dy,
		copl_ip::copl_vector dz,
		copl_ip::copl_vector ds, 
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
		copl_ip::copl_vector var,
		copl_ip::copl_vector dvar
		);
	float get_alpha();
	float get_dtau(); 
	float get_dkappa(); 
	copl_ip::copl_vector get_dx(); 
	copl_ip::copl_vector get_dy(); 
	copl_ip::copl_vector get_dz(); 
	copl_ip::copl_vector get_ds(); 

};
//--------End lp_direction--------

//lp_result
class lp_result {

};
//--------End lp_result--------
}


#endif
