#ifndef INTERIOR_POINT_ALGORITHM
#define INTERIOR_POINT_ALGORITHM

#include <copl_core.h>
#include <copl_linalg.h>
#include <copl_newton.h>
#include <copl_debug.h>

namespace copl_ip{

	bool termination_criteria_met(lp_settings &settings, algorithm_state &state, lp_residuals &residuals);

	void interior_point_algorithm(lp_input &problem_data, lp_settings &settings, lp_variables &variables);

	void interior_point_algorithm_no_answer(lp_input &problem_data, lp_settings &settings);

	void print_status(algorithm_state &state, lp_direction &direction, lp_variables &variables, lp_residuals &residuals, int itr);

	void print_header();
	
	extern "C" {
			

	double c_callable_interior_point_algorithm(int m, int n, int k,\
						 double* Gp, int* Gi, int* Gj, double* h,\
						 double* c,\
						 double* Ap, int* Ai, int* Aj, double* b,\
						 double* x, double* y, double* s, double* z, double* tau, double* kappa,\
						 int max_iter, double linear_feas_tol, double comp_tol);

	double c_callable_test(int n, double* A);


	}

}

#endif
