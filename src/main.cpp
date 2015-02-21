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
	const int max_iter			= 5;
	const double linear_feas_tol 	= 1e-8; //Assuming possible Integer Overflow
	const double comp_tol			= 1e-8; //Assuming possible Integer Overflow
	const double bkscale			= 0.95;

	//copl_vector test(10,10);

	cout << "COPL 2015" << endl;
	cout << "Interior point algorithm coming" << endl;
	
	// Initialize configuration variable
	lp_settings settings(max_iter,linear_feas_tol,comp_tol,bkscale);

	// We are creating an instance of LP. In practice, we should read the problem data from input stream.
	lp_input problem_data(2,3,4); //= construct_instance1();
	{
		problem_data.A.insert_at(0,0,1.0);
		problem_data.A.insert_at(1,1,1.0);
		problem_data.A.insert_at(2,2,1.0);
		problem_data.A.insert_at(2,3,1.0);
		
		problem_data.G.insert_at(0,0,1.0);
		problem_data.G.insert_at(1,1,1.0);
		problem_data.G.insert_at(0,2,1.0);
		problem_data.G.insert_at(1,3,1.0);
		
		problem_data.c[0] = 1;
		problem_data.c[1] = 2;
		problem_data.c[2] = 3;
		problem_data.c[3] = 4;
		
		copl_vector temp_x(4,1.0);
		
		zeros(problem_data.b);
		sp_dgemv(1.0, 1.0, problem_data.A, temp_x, problem_data.b);
		
		zeros(problem_data.h);
		sp_dgemv(1.0, 1.0, problem_data.G, temp_x, problem_data.h);
		
	}
        
	problem_data.var_dump();
		
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
