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


void make_trivial_problem(copl_matrix &A, copl_matrix &G, copl_vector &c, copl_vector &b, copl_vector &h)
{

    c  << 1,2,3,4;
    b  << 5,6;
    h  << 7,8,9;

    A.insert(0,0) = 1.0;
    A.insert(0,1) = 2.0;
    A.insert(1,0) = 3.0;
    A.insert(1,2) = 4.0;    
    A.insert(1,3) = 5.0;
         
    G.insert(0,3) = 6.0;    
    G.insert(1,2) = 7.0;
    G.insert(2,0) = 8.0;
    G.insert(2,1) = 9.0;
    G.insert(2,2) = 10.0;
    G.insert(2,3) = 11.0;
 
}
 

int main()
{
	const int max_iter			    = 5;
	const double linear_feas_tol 	= 1e-8; //Assuming possible Integer Overflow
	const double comp_tol			= 1e-8; //Assuming possible Integer Overflow
	const double bkscale			= 0.95;

	//copl_vector test(10,10);

	OUTPUT << "COPL 2015" << endl;
	OUTPUT << "Interior point algorithm coming" << endl;
	
	// Initialize configuration variable
	lp_settings settings(max_iter,linear_feas_tol,comp_tol,bkscale);	
    
    copl_matrix A(2,4);
    copl_matrix G(3,4);
    copl_vector c(4),b(2),h(3);
    make_trivial_problem(A,G,c,b,h);
    
    lp_input problem_data(A,b,c,G,h);
    
     problem_data.var_dump();
     
	// The main function that run interior point algorithm.
	interior_point_algorithm(problem_data,settings);
	
};

