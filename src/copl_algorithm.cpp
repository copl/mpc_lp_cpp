#include <copl_algorithm.h>
#include <vector>
#include <stdio.h>
using namespace std;
namespace copl_ip {

	void interior_point_algorithm_no_answer(lp_input &problem_data, lp_settings &settings)
	{
		std::vector<double> x(problem_data.n);
		std::vector<double> y(problem_data.k_var);
		std::vector<double> s(problem_data.m);
		std::vector<double> z(problem_data.m);
		double tau, kappa;
		copl_external_vector cx(&x[0],problem_data.n);	
		copl_external_vector cy(&y[0],problem_data.k_var);
		copl_external_vector cs(&s[0],problem_data.m);
		copl_external_vector cz(&z[0],problem_data.m);

        	lp_variables variables (cx,cy,cs,cz,tau,kappa);	
 		interior_point_algorithm(problem_data,settings,variables);
	}
	
	double c_callable_interior_point_algorithm(int m, int n, int k,\
						 double* Gp, int* Gi, int* Gj, double* h,\
						 double* c,\
						 double* Ap, int* Ai, int* Aj, double* b,\
						 double* x, double* y, double* s, double* z, double* tau, double* kappa,\
						 int max_iter, double linear_feas_tol, double comp_tol){
		//Construct the problem data structure
		copl_external_vector cc(c,n); 
		copl_external_vector cb(b,k); 
		copl_external_vector ch(h,m); 
		copl_matrix cA(k,n,Aj[n],Aj,Ai,Ap);
        	copl_matrix cG(m,n,Gj[n],Gj,Gi,Gp);
		lp_input lp_problem(cA,cb,cc,cG,ch);
		//Wrap the variables
		copl_external_vector vars_cx(x,n), vars_cy(y,k), vars_cs(s,m), vars_cz(z,m);
		lp_variables vars(vars_cx,vars_cy,vars_cs,vars_cz,*tau,*kappa);
		lp_settings settings(max_iter,linear_feas_tol,comp_tol,0.95,1e-7);	
		interior_point_algorithm(lp_problem, settings, vars);
		return 0.0;

	}	


	void interior_point_algorithm(lp_input &problem_data, lp_settings &settings, lp_variables &variables){
		
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
		print_header();
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
            		//residuals.var_dump();
			
			if (termination_criteria_met(settings, state, residuals)){
				break;
			}
			
			// compute affine rhs
			rhs.compute_affine_rhs(residuals, variables);
			//rhs.var_dump();


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

	bool termination_criteria_met(lp_settings &settings, algorithm_state &state, lp_residuals &residuals) {
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
	
	void print_header() {
		printf("%2s  %6s  %5s  %8s  %8s  %6s  %8s  %8s  %8s  %8s \n","it","gap","mu","sigma","alpha","tau","r1","r2","r3","r4");
		printf("--------------------------------------------------\n");
	}	

	void print_status(algorithm_state &state, lp_direction &direction, lp_variables &variables, lp_residuals &residuals, int itr) {
		printf("%2i %3.3e %2.3e %5.3e %5.3e %3.3e %5.3e %5.3e %5.3e %5.3e \n",itr\
								   ,state.gap
								   ,state.mu
								   ,state.sigma
								   ,direction.alpha
								   ,variables.tau
								   ,residuals.n1
								   ,residuals.n2
								   ,residuals.n3
								   ,residuals.n4);
	}
}
