/*
 * main.cpp
 *
 *  Created on: 28/01/2015
 *      Author: Oliver
 */
#include <iostream>
#include <copl_algorithm.h>
#include <copl_linalg.h>


using namespace std;
using namespace copl_ip;



int main()
{
	#ifdef DEBUG_TO_FILE
	LOG_FILE_VARIABLE.open("./bin/log.txt"); //Global Variable - defined in copl_debug.h
	#endif
	const int max_iter				= 5;
	const double linear_feas_tol 	= 1e-8; //Assuming possible Integer Overflow
	const double comp_tol			= 1e-8; //Assuming possible Integer Overflow
	const double bkscale			= 0.95;

	//copl_vector test(10,10);

	OUTPUT << "COPL 2015" << endl;
	OUTPUT << "Interior point algorithm coming" << endl;
	
	// Initialize configuration variable
	lp_settings settings(max_iter,linear_feas_tol,comp_tol,bkscale);

	// We are creating an instance of LP. In practice, we should read the problem data from input stream.
	/*lp_input problem_data(2,3,4); //= construct_instance1();
	
	*/
	lp_input * problem_data = copl_utility::Trivial_Test1();
	//copl_utility::loadFromUF("UF_group", "name", &p_problem_data);
	//lp_input problem_data = *p_problem_data;
	
	problem_data->var_dump();
		
	// The main function that run interior point algorithm.
	interior_point_algorithm(*problem_data,settings);
	
	delete problem_data;
};



/*

# Constructing a random LP
function construct_instance1()
	n = 5;
	k = 10;
	m = 10;
	
	problem_data = class_linear_program_input()
	x0 = rand(k,1)

	A = rand(n, k);
	G = -diagm(ones(m));
	c = rand(k,1)
	h = zeros(m);
	b = A*x0;


	problem_data.A = A
	problem_data.G = G
	problem_data.c = c
	problem_data.h = h
	problem_data.b = b
	problem_data.m = m
	problem_data.k = k
	problem_data.n = n
	
	return(problem_data)
end

######################
#  run the program   #
######################

*/
