#ifndef INTERIOR_POINT_ALGORITHM
#define INTERIOR_POINT_ALGORITHM

#include <copl_core.h>
#include <copl_linalg.h>

namespace copl_ip{

bool termination_criteria_met(lp_settings settings, algorithm_state state, lp_residuals residuals);

// return type subject to change
void interior_point_algorithm(lp_input problem_data, lp_settings settings);

}

#endif
