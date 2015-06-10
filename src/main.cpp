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


void make_trivial_problem(Eigen::SparseMatrix<double> &A, Eigen::SparseMatrix<double> &G, copl_vector &c, copl_vector &b, copl_vector &h)
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
	const int max_iter			    = 50;
	const double linear_feas_tol 	= 1e-8; //Assuming possible Integer Overflow
	const double comp_tol			= 1e-8; //Assuming possible Integer Overflow
	const double bkscale			= 0.95;
	const double regularization     = 1e-7;

	//copl_vector test(10,10);

	OUTPUT << "COPL 2015" << endl;
	OUTPUT << "Interior point algorithm coming" << endl;
	
	// Initialize configuration variable
	lp_settings settings(max_iter,linear_feas_tol,comp_tol,bkscale,regularization);	
    
        Eigen::SparseMatrix<double> A(2,4);
        Eigen::SparseMatrix<double> G(3,4);
        
        copl_vector c(4),b(2),h(3);
        make_trivial_problem(A,G,c,b,h);
        
	copl_external_vector cc(&c[0],c.size()); 
	copl_external_vector cb(&b[0],b.size()); 
	copl_external_vector ch(&h[0],h.size()); 	
	A.makeCompressed();
	G.makeCompressed(); 
	copl_matrix cA(A.rows(),A.cols(),A.nonZeros(),A.innerIndexPtr(),A.outerIndexPtr(),A.valuePtr());
        copl_matrix cG(G.rows(),G.cols(),G.nonZeros(),G.innerIndexPtr(),G.outerIndexPtr(),G.valuePtr());

        lp_input problem_data(cA,cb,cc,cG,ch);
    
        problem_data.var_dump();
     
	// The main function that run interior point algorithm.
	interior_point_algorithm_no_answer(problem_data,settings);
	
};

