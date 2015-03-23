/*
 * Linear algebra definitions for the copl interior point solver 
 */

#include <copl_linalg.h>

namespace copl_ip
{

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

}
