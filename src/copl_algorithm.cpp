#include "../include/copl_linalg.h"
#include "../include/copl_core.h"
#include "../include/copl_algorithm.h"

bool termination_criterion_met(lp_settings settings, algorithm_state state, lp_residuals residuals){
	// TODO
	retrn false;
}

void interior_point_algorithm(copl_ip::lp_input problem_data, copl_ip::lp_settings settings){
	copl_ip::lp_variables variables (problem_data);	
	copl_ip::algorithm_state state();
	state.update_mu(variables, problem_data);

	copl_ip::k_newton_matrix K_newton_matrix(problem_data);
	copl_ip::linear_system_rhs rhs(problem_data);
	copl_ip::lp_direction direction(problem_data);
	copl_ip::lp_residuals residuals(problem_data);

	// Begin iteration
	for (int itr = 1; itr <=settings.get_max_iter; i++){
		// To be sent to Tiago's Linear Solver
		K_newton_matrix.update(variables);

		// compute residuals
		residuals.compute_residuals(problem_data, variables);

		if (termination_criterion_met(settings, state, residuals)){
			break;
		}
		// compute affine rhs
		rhs.compute_affine_rhs(residuals, variables);

		// compute affine direction using new affine rhs
		direction.compute_affine_direction(rhs,problem_data,variables,K_newton_matrix); //Incomplete??
		
		// update corrector rhs using new affine direction
		rhs.compute_corrector_rhs(residuals,variables,state,direction,problem_data)

		// update corrector direction using new corrector rhs
	    direction.compute_corrector_direction(
	    	rhs,
			problem_data,
			variables,
			state,
			settings,
			K_newton_matrix
			);
		
		
		// take step in the corrector direction
		variables.take_step(direction)
	    
		// compute the gap after the step
		state.update_mu(variables,problem_data)
		state.update_gap(variables,problem_data)


	}
}
