//Load a matrix A and vector b 
//Calculate r = b-Ax in a predetermined storage r
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <random>

using namespace std;
//Shorthand typedef
typedef Eigen::SparseMatrix<double> SpMat;
//Shorthand for a vector of triplets
typedef std::vector<Eigen::Triplet<double>> triplet_vector_t;

//Generates a random A and returns nnz, sets the values in triplet vector
int  generate_random_A(int m, int n, triplet_vector_t &vals)
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
	return nnz_count;
}

void generate_random_b(int m, Eigen::VectorXd &b)
{
    normal_distribution<double> rnorm(0,1);
	default_random_engine generator;
	for(int i=0;i<m;i++)
		b(i) = rnorm(generator);
}

int main(void)
{
	int m = 10;
	int n = 20;
	triplet_vector_t entries; 
	int nnz = generate_random_A(m,n,entries);
	SpMat A(m,n);
	A.setFromTriplets(entries.begin(),entries.end());
	SpMat AtA = A*A.transpose();
	Eigen::SimplicialCholesky<SpMat> solver(AtA);
	//Put b in eigen vector VectorXd 
	Eigen::VectorXd b(m);
	generate_random_b(m,b);
	//solve 
	//Check result 

}