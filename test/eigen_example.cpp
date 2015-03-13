//For those who are curious about EIGEN
//We should not use Eigen directly, instead we should use the copl_linalg 

//Load a matrix A and vector b 
//Calculate r = b-Ax in a predetermined storage r
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <random>
#include <iostream>

using namespace std;
//Shorthand typedef
typedef Eigen::SparseMatrix<double> SpMat;
//Shorthand for a vector of triplets
typedef std::vector<Eigen::Triplet<double>> triplet_vector_t;

//Generates a random A and returns nnz, sets the values in triplet vector
void  generate_random_A(int m, int n, triplet_vector_t &vals)
{
	default_random_engine generator;
	bernoulli_distribution random(0.3);
	normal_distribution<double> random_normal(0,1.0);
	//Sample each entry with highish probability so we get a full rank matrix
	int nnz_count = 0;
	double val;

	//Sample all entries and select iid	
	for(int in=0; in<n;in++)
	for(int im=0;im<m;im++)
	{
		if(random(generator)) //Sample and entry
		{ 	
			//Put the entry in the list 
			val = random_normal(generator);
			vals.push_back(Eigen::Triplet<double>(im,in,val));
			nnz_count++;
		}
	}
}

void generate_random_b(int m, Eigen::VectorXd &b)
{
    normal_distribution<double> rnorm(0,1);
	default_random_engine generator;
	for(int i=0;i<m;i++)
		b[i] = rnorm(generator);
}

void concatenate_test()
{
	//Build a sprase Eigien matrix A and G 
	//Form IA
 	//     A0
	 Eigen::SparseMatrix<double> A(10,100), B(100,100), C(110,110), Z(100,100);
	 //Call the concatenation operator 	
	 Z.setIdentity();	 
	 cout << " Z ("<<Z.rows() << "," << Z.cols() << ")\n";
	 C.block(0,0,101,101) = Z;
	 //Test the size of C
	 cout << " C ("<<C.rows() << "," << C.cols() << ")\n";
}

int main(void)
{
	int m = 10;
	int n = 20;
	triplet_vector_t entries; 
	generate_random_A(m,n,entries);
	SpMat A(m,n);
	A.setFromTriplets(entries.begin(),entries.end());
	SpMat AAt = A*A.transpose();
	Eigen::SimplicialCholesky<SpMat> solver(AAt);
	//Put b in eigen vector VectorXd 
	Eigen::VectorXd b(m);
	generate_random_b(m,b);

	std::vector<double> xstl(m);
	xstl[0] = 1;
	Eigen::Map<Eigen::VectorXd> x(&xstl[0],m);

	std::vector<double> ystl(m);
	ystl[0] = 1;
	Eigen::Map<Eigen::VectorXd> y(&ystl[0],m);
	
	x = solver.solve(b);
	cout << "x(0): " << x[0] << "\n";
	y = AAt*x;
	cout << "y(0): " << y[0] << "\n";

	concatenate_test();
}