
#include <copl_algorithm.h>



using namespace std;
namespace copl_ip {
	void interior_point_algorithm(lp_input &problem_data, lp_settings &settings){
	
		problem_data.var_dump();
		
		// create data structures
		lp_variables variables (problem_data.n,problem_data.m,problem_data.k_var);	
		algorithm_state state;
		k_newton_copl_matrix K_matrix(problem_data);
		lp_direction direction(variables);
		linear_system_rhs rhs(problem_data);
		lp_residuals residuals(problem_data);
		
		state.update_mu(variables, problem_data);
		
		// Begin iteration
		int MAX_IT = settings.get_max_iter();
		for (int itr = 1; itr <= MAX_IT; itr++){
			// To be sent to Tiago's Linear Solver
			K_matrix.update(variables);
			// compute residuals
			residuals.compute_residuals(problem_data, variables);
			residuals.var_dump();
			
			if (termination_criteria_met(settings, state, residuals)){
				break;
			}
			
			
			// compute affine rhs
			rhs.compute_affine_rhs(residuals, variables);
			rhs.var_dump();
			// --- BEGIN NOT WORKING --- //
			// compute affine direction using new affine rhs
			direction.compute_affine_direction(rhs,problem_data,variables,K_matrix); //Incomplete??
			
			// update corrector rhs using new affine direction
			rhs.compute_corrector_rhs(residuals,variables,state,direction,problem_data);

			// update corrector direction using new corrector rhs
			direction.compute_corrector_direction(
				rhs,
				problem_data,
				variables,
				state,
				settings,
				K_matrix
				);
			
			// --- END NOT WORKING --- //
			
			// take step in the corrector direction
			variables.take_step(direction);
			
			// compute the gap after the step
			state.update_mu(variables,problem_data);
			state.update_gap(variables,problem_data);
			
			print_status(state, direction, variables, residuals, itr);
		}
		cout << "IP algorithm finished" << endl;
	}

	bool termination_criteria_met(lp_settings &settings, algorithm_state &state, lp_residuals &residuals){
		/*
		# TO DO
		# store a bunch of norms
		
		#Evaluate termination criteria	
		#TODO: This part will have to be a more sofisticated test to detect 
		#unbounded and infeasible problems.
		*/
		
		if (residuals.get_r1_norm() < settings.get_linear_feas_tol() && 
			residuals.get_r2_norm() < settings.get_linear_feas_tol() &&
			residuals.get_r3_norm() < settings.get_linear_feas_tol() && 
			state.mu < settings.get_comp_tol())
			 return true;
		else
			return false;
		
		return false;
	}
	
	void print_status(algorithm_state &state, lp_direction &direction, lp_variables &variables, lp_residuals &residuals, int itr) {
		cout << "it:" << itr << " gap:" << state.gap << " mu:" << state.mu << " alpha:" << direction.alpha << " tau:" << variables.tau << " residuals:" << residuals.get_norm_squared() << endl;
		//@printf("%3i\t%3.3e\t%3.3e\t%3.3e\t%3.3e\t%3.3e\n", itr, state.gap ,state.mu, direction.alpha, variables.tau, residuals.normed_squared)
	}
}