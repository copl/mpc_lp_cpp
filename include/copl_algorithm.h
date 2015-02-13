#ifndef INTERIOR_POINT_ALGORITHM
#define INTERIOR_POINT_ALGORITHM

#include "copl_core.h"
<<<<<<< HEAD

void interior_point_algorithm(copl_ip::lp_input problem_data, copl_ip::settings settings);
=======
#include "copl_linalg.h"

bool termination_criteria_met(lp_settings settings, algorithm_state state, lp_residuals residuals);

// return type subject to change
void interior_point_algorithm(lp_input problem_data, lp_settings settings);
>>>>>>> e56f22bdfdbaaba154a7681f7b79d59c7dba3341


#endif
