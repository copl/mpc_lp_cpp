
#include <copl_algorithm.h>
#include <copl_core.h>

namespace copl_ip {
	
	void interior_point_algorithm(lp_input problem_data, lp_settings settings){
		// create data structures
		lp_variables variables (problem_data.n,problem_data.m,problem_data.k_var);	
		algorithm_state state;
		
		cout << "here1" << endl;
		
		k_newton_copl_matrix K_matrix(problem_data);
		cout << "here2" << endl;
		lp_direction direction(variables);
		
		cout << "here3" << endl;
		
		linear_system_rhs rhs(problem_data);
		
		cout << "here4" << endl;
		
		lp_residuals residuals(problem_data);
		
		state.update_mu(variables, problem_data);

		cout << "here" << endl;

		// Begin iteration
		
		int MAX_IT = settings.get_max_iter();
		for (int itr = 1; itr <= MAX_IT; itr++){
			// To be sent to Tiago's Linear Solver
			K_matrix.update(variables);

			// compute residuals
			residuals.compute_residuals(problem_data, variables);

			if (termination_criteria_met(settings, state, residuals)){
				break;
			}
			// compute affine rhs
			rhs.compute_affine_rhs(residuals, variables);

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
			
			
			// take step in the corrector direction
			variables.take_step(direction);
			
			// compute the gap after the step
			state.update_mu(variables,problem_data);
			state.update_gap(variables,problem_data);
		}
	}

	bool termination_criteria_met(lp_settings settings, algorithm_state state, lp_residuals residuals){
		/*
		# TO DO
		# store a bunch of norms
		
		#Evaluate termination criteria	
		#TODO: This part will have to be a more sofisticated test to detect 
		#unbounded and infeasible problems.
		
		if (residuals.r1_norm < settings.linear_feas_tol && 
			residuals.r2_norm < settings.linear_feas_tol &&
			residuals.r3_norm < settings.linear_feas_tol && 
			state.mu < settings.comp_tol)
			 println("Ended");
			 return true
		else
			return false
		end
		*/
		return false;
	}
}