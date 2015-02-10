copl_ip::class_settings::class_settings (int input_max_iter, float input_linear_feas_tol, float input_comp_tol, float input_bkscale) {
	max_iter 			= input_max_iter;
	linear_feas_tol 	= input_linear_feas_tol;
	comp_tol 			= input_comp_tol;
	bkscale 			=input_bkscale;
}
//--------End class_settings--------
void copl_ip::class_linear_program_variables::class_linear_program_variables(class_linear_program_input problem_data){

	tau = 1;
	kappa = 1;
}

void class_linear_program_variables::take_step(class_direction direction){
	float alpha = direction.get_alpha();
	//TODO: Implement copl_ip::vector.multiply and copl_ip::vector.add
	// x = x.add(direction.get_dx.multiply(alpha));
	// s = s.add(direction.get_ds.multiply(alpha));
	// z = z.add(direction.get_dz.multiply(alpha));
	// y = y.add(direction.get_dy.multiply(alpha));
	tau = tau + alpha * direction.get_dtau();
	kappa = kappa + alpha * direction.get_dkappa();

}
//--------End class_lienar_program_variables--------

class_algorithm_state::class_algorithm_state() {
	
}
void class_algorithm_state::update_gap(class_linear_program_variables variables, class_linear_program_input problem_data){

}

void class_algorithm_state::update_mu(class_linear_program_variables variables, class_linear_program_input problem_data){
	
}
//--------End class_algorithm_state--------

class_linear_system_rhs::class_linear_system_rhs(class_linear_program_input problem_data){

}

void class_linear_system_rhs::update_values(copl_ip::vector q1, copl_ip::vector q2, copl_ip::vector q3, copl_ip::vector q4, copl_ip::vector q5, copl_ip::vector q6){

}

//--------End class_lienar_system_rhs--------


