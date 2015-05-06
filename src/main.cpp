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
	srand (0);
	for (int i = 0; i  < 4; i++)
		for (int j = 0; j  < 5; j++)
			A.insert_at(i,j,rand() % 10 + 1);
	for (int i = 0; i  < 6; i++)
		for (int j = 0; j  < 5; j++)
			G.insert_at(i,j,rand() % 10 + 1);

	for (int i = 0; i  < 5; i++)
		c.at(i) = rand() % 10 + 1;

	for (int i = 0; i  < 4; i++)
		b.at(i) = rand() % 10 + 1;

	for (int i = 0; i  < 6; i++)
		h.at(i) = rand() % 10 + 1;	
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
    
    copl_matrix A(4,5);
    copl_matrix G(6,5);
    copl_vector c(5),b(4),h(6);
    make_trivial_problem(A,G,c,b,h);
    
    lp_input problem_data(&A,&b,&c,&G,&h);
	//copl_utility::loadFromUF("UF_group", "name", &p_problem_data);
	//lp_input problem_data = *p_problem_data;
	
	problem_data.var_dump();
		
	// The main function that run interior point algorithm.
	interior_point_algorithm(problem_data,settings);
	
};

