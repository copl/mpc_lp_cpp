/*
 * main.cpp
 *
 *  Created on: 28/01/2015
 *      Author: Oliver
 */
#include <iostream>
#include "copl_algorithm.h"
#include "copl_linalg.h"

#include "copl_core.cpp"
#include "copl_algorithm.cpp"
#include "copl_linalg.cpp"

using namespace std;
using namespace copl_ip;

lp_input construct_instance1() {
	lp_input problem_data(2,2,2);
	
	problem_data.var_dump();
	
	return problem_data;
};

int main()
{
	const int max_iter			= 20;
	const double linear_feas_tol 	= 1e-8; //Assuming possible Integer Overflow
	const double comp_tol			= 1e-8; //Assuming possible Integer Overflow
	const double bkscale			= 0.95;

	copl_vector test(10,10);

	cout << "COPL 2015" << endl;
	cout << "Interior point algorithm coming" << endl;

	// Initialize configuration variable
	lp_settings settings(max_iter,linear_feas_tol,comp_tol,bkscale);

	// We are creating an instance of LP. In practice, we should read the problem data from input stream.
	lp_input problem_data = construct_instance1();

	// The main function that run interior point algorithm.
	interior_point_algorithm(problem_data,settings);
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
