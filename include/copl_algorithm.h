#ifndef INTERIOR_POINT_ALGORITHM
#define INTERIOR_POINT_ALGORITHM

#include <copl_core.h>
#include <copl_linalg.h>
#include <copl_newton.h>
#include <copl_debug.h>

namespace copl_ip{

	bool termination_criteria_met(lp_settings &settings, algorithm_state &state, lp_residuals &residuals);

	void interior_point_algorithm(lp_input &problem_data, lp_settings &settings, lp_variables &variables);

	void print_status(algorithm_state &state, lp_direction &direction, lp_variables &variables, lp_residuals &residuals, int itr);

}

#endif
