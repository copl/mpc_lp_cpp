/*
 * Linear algebra definitions for the copl interior point solver 
 */

#include <copl_linalg.h>

namespace copl_ip
{

//y<- alpha Ax + beta y with sparse A
void sp_dgemv(double alpha, double beta, copl_matrix &copl_A, copl_vector &copl_x, copl_vector &copl_y)
{
	//Make a map to be able to use y and x in Eigen
	Eigen::Map<Eigen::VectorXd> x(&copl_x[0],copl_x.size());
	Eigen::Map<Eigen::VectorXd> y(&copl_y[0],copl_y.size());	
	EigenSpMat_t A = *(copl_A.eigenMat);
	y = alpha*(A*x)+beta*y;
}

//copl_matrix vector multiply and accumulate in y
//y<- alpha A^Tx + beta y with sparse A (The input copl_matrix is A not A^T)
void sp_dgemtv(double alpha, double beta, copl_matrix &copl_A, copl_vector &copl_x, copl_vector &copl_y)
{
	//Make a map to be able to use y and x in Eigen
	Eigen::Map<Eigen::VectorXd> x(&copl_x[0],copl_x.size());
	Eigen::Map<Eigen::VectorXd> y(&copl_y[0],copl_y.size());
	EigenSpMat_t A = *(copl_A.eigenMat);
	OUTPUT << A.innerSize() << "\t" << A.outerSize() << "\t" << x.size()<<"---\n";
	y = alpha*(A.transpose()*x)+beta*y;	

}

//Scale the vector
//x<-alpha *x
void scal(copl_vector &copl_x, double alpha )
{
	Eigen::Map<Eigen::VectorXd> x(&copl_x[0],copl_x.size());
	x *= alpha;
}

//y<- a*x + y
void axpy(double alpha, copl_vector &copl_x, copl_vector &copl_y)
{
	Eigen::Map<Eigen::VectorXd> x(&copl_x[0],copl_x.size());
	Eigen::Map<Eigen::VectorXd> y(&copl_y[0],copl_y.size());
	y = alpha*x+y;
}

// x^Ty
double dot(copl_vector &copl_y, copl_vector &copl_x)
{
	Eigen::Map<Eigen::VectorXd> x(&copl_x[0],copl_x.size());
	Eigen::Map<Eigen::VectorXd> y(&copl_y[0],copl_y.size());
	return y.dot(x);
}

//Zero out 
void zeros(copl_vector &y)
{
	size_t j = 0;
	for(j=0;j<y.size();j++)
	{
		y[j] = 0.0;	
	}
}

double norm2(copl_vector &y)
{

	double norm = 0.0;
	for(auto iter = y.begin(); iter< y.end(); iter++)
	{
		norm+= (*iter)*(*iter);
	}
	return sqrt(norm);
}

double normInf(copl_vector &y)
{
	double max_val;
	double abs_val;
	auto iter = y.begin();
	max_val = *iter;
	iter++;
	for(; iter< y.end(); iter++)
	{
		abs_val = fabs(*iter);
		if(max_val< abs_val)
			max_val = abs_val;
	}
	return max_val;
}

typedef std::vector<Eigen::Triplet<double>> triplet_vector_t;
//Generates a random sets the values in triplet vector
void  generate_random_A(int m, int n, triplet_vector_t &vals, double p)
{
	std::default_random_engine generator;
	std::bernoulli_distribution random(p);
	std::normal_distribution<double> random_normal(0,1.0);
	//Sample each entry with highish probability so we get a full rank copl_matrix
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
		}
	}
}

//The copl copl_matrix 
	
	copl_matrix::copl_matrix(int m, int n)
	{
		eigenMat = new EigenSpMat_t(m,n); 
	}

	//Generate a random copl_matrix selecting the entries iid with prob p
	copl_matrix::copl_matrix(int m, int n, double p)
	{
		triplet_vector_t entries;
		generate_random_A(m,n,entries,p);
		eigenMat = new EigenSpMat_t(m,n);
		eigenMat->setFromTriplets(entries.begin(),entries.end());
	}

	void copl_matrix::insert_at(int m, int n, double val)
	{
		(eigenMat->insert(m,n)) = val;
	}
	
	//void copl_matrix::set_value(int m, int n, double val)
	//{
		//eigenMat(m,n) = val;
	//}
	
	int copl_matrix::num_rows()
	{
		return eigenMat->rows();
	}
	
	int copl_matrix::num_cols()
	{
		return eigenMat->cols();
	}
	
	double copl_matrix::value_at(int i, int j)
	{
		return eigenMat->coeffRef(i,j);
	}

	void copl_matrix::var_dump()
	{
		OUTPUT << *eigenMat << endl;
	}
	
	
	//Destructor 
	copl_matrix::~copl_matrix()
	{
		OUTPUT << "deleting eigen matrix" << endl;
		if(eigenMat!=NULL)
			delete(eigenMat);
	}

	//Return the number of non zeros
	int copl_matrix::nnz()
	{
		return eigenMat->nonZeros();	
	}
}
