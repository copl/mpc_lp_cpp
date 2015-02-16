#include "../include/copl_core.h"
#include "../include/copl_linalg.h"

namespace copl_ip {

//lp_input
lp_input::lp_input(int _m, int _n, int _k) 
	: A(_n,_m), G(_n,_m), c(_n,0.0), h(_n,0.0), b(_n,0.0) { // ************
	m = _m;
	n = _n;
	k = _k;
};
//--------End lp_input--------

// lp_settings
lp_settings::lp_settings (int input_max_iter, float input_linear_feas_tol, float input_comp_tol, float input_bkscale) {
	max_iter 			= input_max_iter;
	linear_feas_tol 	= input_linear_feas_tol;
	comp_tol 			= input_comp_tol;
	bkscale 			=input_bkscale;
}
int lp_settings::get_max_iter()			{return max_iter;}
float lp_settings::get_linear_feas_tol() 	{return linear_feas_tol;}
float lp_settings::get_comp_tol() 			{return comp_tol;}
float lp_settings::get_bkscale() 			{return bkscale;}

//--------End lp_settings--------

// lp_residuals
lp_residuals::lp_residuals(lp_input problem_data){ }

void lp_residuals::compute_residuals(lp_input problem_data, lp_variables variables){
	// r1 = -pd.A'*variables.y - pd.G'*variables.z - pd.c*variables.tau;
	
	zeros(r1);
	sp_dgemtv(-1.0, 1.0, problem_data.A, variables.y, r1);
	sp_dgemtv(-1.0, 1.0, problem_data.G, variables.z, r1);
	axpy(-variables.tau, problem_data.c, r1);
	
	// r2 = pd.A*variables.x - pd.b*variables.tau;
	zeros(r2);
	sp_dgemv(1.0, 1.0, problem_data.A, variables.x, r2);
	axpy(-variables.tau, problem_data.c, r2);
	
	// r3 = variables.s + pd.G*variables.x - variables.tau*pd.h;
	zeros(r3);
	axpy(1.0, variables.s, r3);
	sp_dgemv(1.0, 1.0, problem_data.G, variables.x, r3);
	axpy(-variables.tau, problem_data.h, r3);
	
	//r4 = variables.kappa + pd.c'*variables.x + pd.b'*variables.y + + pd.h'*variables.z;
	r4 = variables.kappa;
	r4 += dot(problem_data.c, variables.x);
	r4 += dot(problem_data.b, variables.y);
	r4 += dot(problem_data.h, variables.z);
}

float lp_residuals::get_r1_norm(){return r1_norm;}
float lp_residuals::get_r2_norm(){return r2_norm;}
float lp_residuals::get_r3_norm(){return r3_norm;}



//--------End lp_residuals--------

// lp direction
lp_direction::lp_direction(lp_variables variables) {
	// allocate memory
	copl_vector dx(variables.x.size(),0);
	copl_vector dy(variables.y.size(),0);
	copl_vector dz(variables.z.size(),0);
	copl_vector ds(variables.s.size(),0);
}

float lp_direction::get_alpha() { return alpha; }
float lp_direction::get_dtau() { return dtau; }
float lp_direction::get_dkappa() { return dkappa; } 
copl_vector lp_direction::get_dx() { return dx; }
copl_vector lp_direction::get_dy() { return dy; } 
copl_vector lp_direction::get_dz() { return dz; }
copl_vector lp_direction::get_ds() { return ds; }

void  lp_direction::compute_affine_direction(linear_system_rhs affine_rhs,
		lp_input problem_data,
		lp_variables variables,
		k_newton_copl_matrix K_newton_copl_matrix) {
	
}

void lp_direction::compute_corrector_direction(
		linear_system_rhs corrector_rhs,
		lp_input problem_data,
		lp_variables variables,
		algorithm_state state,
		lp_settings settings,
		k_newton_copl_matrix K_newton_copl_matrix
		) {
			
}

//-----------End lp direction



// lp_variables
lp_variables::lp_variables(lp_input problem_data){

	tau = 1;
	kappa = 1;
}

void lp_variables::take_step(lp_direction direction){
	float alpha = direction.get_alpha();
	//TODO: Implement copl_vector.multiply and copl_vector.add
	// x = x.add(direction.get_dx.multiply(alpha));
	// s = s.add(direction.get_ds.multiply(alpha));
	// z = z.add(direction.get_dz.multiply(alpha));
	// y = y.add(direction.get_dy.multiply(alpha));
	tau = tau + alpha * direction.get_dtau();
	kappa = kappa + alpha * direction.get_dkappa();

}
//--------End lp_variables--------

//k_newton_copl_matrix
k_newton_copl_matrix::k_newton_copl_matrix(lp_input problem_data) {
	
}
void k_newton_copl_matrix::update(lp_variables variables) {
	
}

//--------End k_newton_copl_matrix--------

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
	copl_vector q1, 
	copl_vector q2, 
	copl_vector q3, 
	copl_vector q4, 
	copl_vector q5, 
	copl_vector q6
	){

}

void linear_system_rhs::compute_affine_rhs(lp_residuals residuals, lp_variables variables){
	// TOOD
}
void linear_system_rhs::compute_corrector_rhs(lp_residuals residuals, lp_variables variables, algorithm_state state, lp_direction direction, lp_input problem_data){
	// TOOD
}
//--------End linear_system_rhs--------
}

