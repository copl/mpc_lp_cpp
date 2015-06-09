#include <copl_algorithm.h>

using namespace std;
namespace copl_ip {
	void interior_point_algorithm(lp_input &problem_data, lp_settings &settings){
			OUTPUT << "start" << endl;
			
			lp_timer interior_point_timer;

			interior_point_timer.start();

        	problem_data.var_dump();
        		
        	//Allocate the variables
        	lp_variables variables (problem_data.m,problem_data.n,problem_data.k_var);	
        	algorithm_state state;
                
			//Create and analyze the newton matrix 
			homogeneous_solver K_matrix(problem_data,settings);
			
			//This stores the search directions
			lp_direction direction(variables);
			
			//Helps generate the rhs for the linear solves	
			linear_system_rhs rhs(problem_data);
                
            //Contains the linear residuals
        	lp_residuals residuals(problem_data);
        	    
            //Compute initial mu
        	state.update_mu(variables, problem_data);
		
		// Begin iteration
		print_status(state, direction, variables, residuals, 0);

		lp_timer iteration_timer;
		for (int itr = 1; itr <= settings.max_iter; itr++){
			OUTPUT << "Beginning Iteration " << itr << endl;
			iteration_timer.start();
			//Update the linear system with the present value
			K_matrix.update(variables);

			// compute residuals
			residuals.compute_residuals(problem_data, variables);

			//Debug call
            		residuals.var_dump();
			
			if (termination_criteria_met(settings, state, residuals)){
				break;
			}
			
			// compute affine rhs
			rhs.compute_affine_rhs(residuals, variables);
			rhs.var_dump();

			//Solve the system and compute the affine direction	
			K_matrix.solve(direction,rhs,variables);
     			
			//Calculate sigma 
			direction.compute_step_size(variables,settings);		
			
			state.sigma =  (1-direction.alpha);
			state.sigma *= (1-direction.alpha);
			state.sigma *= (1-direction.alpha);
			
			// update corrector rhs using new affine direction
			rhs.compute_combined_rhs(residuals,variables,direction,state.sigma,state.mu);

			//Solve the system and compute the combined direction	
			K_matrix.solve(direction,rhs,variables);
		
			direction.compute_step_size(variables,settings);				
			
			// take step in the corrector direction
			variables.take_step(direction);
			
	    	        // compute the gap after the step
			state.update_mu(variables,problem_data);
			state.update_gap(variables,problem_data);
			iteration_timer.end();
			print_status(state, direction, variables, residuals, itr);
			OUTPUT << "Iteration " << itr << " elapsed time: " << iteration_timer.get_total_time() << " seconds" << endl;
			OUTPUT << "------------------------------" << endl;
		}
		interior_point_timer.end();
		OUTPUT << "IP algorithm finished" << endl;
		OUTPUT << "Total time elapsed: " << interior_point_timer.get_total_time() << " seconds" << endl;
	}

	bool termination_criteria_met(lp_settings &settings, algorithm_state &state, lp_residuals &residuals){
		/*
		# TO DO
		# store a bunch of norms
		
		#Evaluate termination criteria	
		#TODO: This part will have to be a more sofisticated test to detect 
		#unbounded and infeasible problems.
		*/
		
		if (residuals.n1 < settings.linear_feas_tol && 
			residuals.n2 < settings.linear_feas_tol &&
			residuals.n3 < settings.linear_feas_tol && 
			state.mu < settings.comp_tol)
			 return true;
		else
			return false;
		
		return false;
	}
	
	void print_status(algorithm_state &state, lp_direction &direction, lp_variables &variables, lp_residuals &residuals, int itr) {
		cout << "it:" << itr << " gap:" << state.gap << " mu:" << state.mu << " sigma:" << state.sigma << " alpha:" << direction.alpha << " tau:" << variables.tau << " residuals:" << residuals.n1<<","<<residuals.n2<<","<<residuals.n3 << ","<<residuals.n4 << endl;
		//@printf("%3i\t%3.3e\t%3.3e\t%3.3e\t%3.3e\t%3.3e\n", itr, state.gap ,state.mu, direction.alpha, variables.tau, residuals.normed_squared)
	}
}
