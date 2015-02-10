//TODO: Use a test bench for these 

//Using copl linalg 
#include <copl_linalg.h>
#include <random>
#include <iostream>

//Define a random matrix
//and a random rhs

//Compute a residual 
//Compute the residual norm
using namespace copl_ip;
using namespace std;

//Utility function to generate random vectors
void fill_rand(copl_vector &b)
{
    normal_distribution<double> rnorm(0,1);
	default_random_engine generator;
	for(size_t i=0;i<b.size();i++)
		b[i] = rnorm(generator);
}

int main(void)
{

	int m = 10;
	int n = 10;
	double p = 0.5;
	//Make a copl matrix
	matrix A(m,n,p);

	//make a random vector x 
	std::vector<double> x(n);
	fill_rand(x);

	//Make a random vector b	
	std::vector<double> b(m);
	fill_rand(b);

	//Calculate r = b-Ax
	std::vector<double> r(m);	
	cout << "r[0]: " << r[0] << "\n";
    sp_dgemv(-1.0, 0.0, A, x, r);	  	 
   	//Now add b to r and store in r
   	axpy(1.0,b,r);
   	//At this point r = b-Ax
	cout << "r = b-Ax, r[0]: " << r[0] << "\n";

	//Return the norm of r 
	cout << "norm2(r): " << norm2(r) << "\n";		
	cout << "normInf(r): " << normInf(r) << "\n";		
	//Zero out and try again 
	zeros(r);
	cout << "norm2(r): " << norm2(r) << "\n";		
	cout << "normInf(r): " << normInf(r) << "\n";	

}