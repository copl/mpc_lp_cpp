#include "../include/copl_core.h"
#include "../include/copl_linalg.h"

// lp_settings
copl_ip::lp_settings::lp_settings (int input_max_iter, float input_linear_feas_tol, float input_comp_tol, float input_bkscale) {
	max_iter 			= input_max_iter;
	linear_feas_tol 	= input_linear_feas_tol;
	comp_tol 			= input_comp_tol;
	bkscale 			=input_bkscale;
}
int copl_ip::lp_settings::get_max_iter()			{return max_iter;}
float copl_ip::lp_settings::get_linear_feas_tol() 	{return linear_feas_tol;}
float copl_ip::lp_settings::get_comp_tol() 			{return comp_tol;}
float copl_ip::lp_settings::get_bkscale() 			{return bkscale;}

//--------End lp_settings--------

// lp_residuals
copl_ip::lp_residuals::lp_residuals(){ }

void copl_ip::lp_residuals::update_values(
	copl_ip::vector r1,
	copl_ip::vector r2,
	copl_ip::vectpr r3,
	copl_ip::vector r4
	){
	this->r1 = r1;
	this->r2 = r2;
	this->r3 = r3;
	this->r4 = r4;

	// TODO
}

void copl_ip::lp_residuals::compute_residuals(lp_input in_problem_data, lp_variables variables){
	// TODO
}

float copl_ip::lp_residuals::get_r1_norm(){return r1_norm;}
float copl_ip::lp_residuals::get_r2_norm(){return r2_norm;}
float copl_ip::lp_residuals::get_r3_norm(){return r3_norm;}



//--------End lp_residuals--------



// lp_variables
void copl_ip::lp_variables::lp_variables(lp_input problem_data){

	tau = 1;
	kappa = 1;
}

void lp_variables::take_step(lp_direction direction){
	float alpha = direction.get_alpha();
	//TODO: Implement copl_ip::vector.multiply and copl_ip::vector.add
	// x = x.add(direction.get_dx.multiply(alpha));
	// s = s.add(direction.get_ds.multiply(alpha));
	// z = z.add(direction.get_dz.multiply(alpha));
	// y = y.add(direction.get_dy.multiply(alpha));
	tau = tau + alpha * direction.get_dtau();
	kappa = kappa + alpha * direction.get_dkappa();

}
//--------End lp_variables--------

// algorithm_state

algorithm_state::algorithm_state() {
	
}
void algorithm_state::update_gap(lp_variables variables, lp_input problem_data){

}

void algorithm_state::update_mu(lp_variables variables, lp_input problem_data){
	
}
//--------End algorithm_state--------


// linear_system_rhs

linear_system_rhs::linear_system_rhs(lp_input problem_data){

}

void linear_system_rhs::update_values(
	copl_ip::vector q1, 
	copl_ip::vector q2, 
	copl_ip::vector q3, 
	copl_ip::vector q4, 
	copl_ip::vector q5, 
	copl_ip::vector q6
	){

}

void linear_system_rhs::compute_affine_rhs(lp_residuals residuals, lp_variables variables){
	// TOOD
}
void linear_system_rhs::compute_correctpr_rhs(lp_residual residuals, lp_variables variables){
	// TOOD
}
//--------End linear_system_rhs--------


