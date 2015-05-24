#include <copl_newton.h>

namespace copl_ip{

void copl_vector_dump1(copl_vector &vec) {
	
	OUTPUT << "{";
	for(int i = 0; i < vec.size() - 1; i++) {
		OUTPUT << vec.at(i) << ",";
	}
	OUTPUT << vec[vec.size() - 1];
	OUTPUT << "}";
	OUTPUT << endl;
};

int k_newton_copl_matrix::nnz()
{
	return eigenKMat->nonZeros();
}

k_newton_copl_matrix::k_newton_copl_matrix(int m, int n)
{
	cout << "Testing constructor\n";
	hessianIx = new std::vector<int>(m);
	eigenKMat = new EigenSpMat_t(n,n);
}

k_newton_copl_matrix::k_newton_copl_matrix(copl_matrix &A, copl_matrix &G)
{
	
	p = A.num_rows();	
	n = G.num_cols();		
    m = G.num_rows();	

	//Now find the indices of the hessian
	hessianIx = new std::vector<int>(m);
	eigenKMat = new EigenSpMat_t(n+m+p,n+m+p);
	
    //Assemble 	
	assemble_matrix(A, G);
    
    //Compress 
	eigenKMat->makeCompressed();	

    //The indices are the last nonzero per column for the cols n+p to n+p+m-1
	int* outerIx = eigenKMat->outerIndexPtr();		
	for(int i = 0; i < m; i++ )
		 (*hessianIx)[i] = outerIx[i+n+p+1]-1;	
    
    solver.analyzePattern(*eigenKMat);

}

void k_newton_copl_matrix::assemble_matrix(copl_matrix &A, copl_matrix &G)
{
	//Assemble 
	//Count the number of non zeros per row and col of A and G
	std::vector<int> nzRowA(p);		
	std::vector<int> nzRowG(m);	
	std::vector<int> nzColA(n);		
	std::vector<int> nzColG(n);	
	std::vector<int> nnzColsK(m+n+p);	
	//TODO this depends on EIGENs iterators and will need 
	//eventual refactoring	
	//Get the internal matrices 

   	EigenSpMat_t &eA = *(A.eigenMat);	
	EigenSpMat_t &eG = *(G.eigenMat);
	//Count the non zeros per column and row
	for(int i = 0; i < eA.cols(); i++ )
		for(Eigen::SparseMatrix<double>::InnerIterator it(eA,i); it; ++it)	
		{
			nzRowA[it.row()]++;
			nzColA[it.col()]++;
		}
	for(int i = 0; i < eG.cols(); i++ )
		for(Eigen::SparseMatrix<double>::InnerIterator it(eG,i); it; ++it)	
		{
			nzRowG[it.row()]++;
			nzColG[it.col()]++;
		}		
	//Calculate the number of non zeros per column of K
	// d A'  G'
    // A -d  0
	// G 0  -d - H
	//The d represent dI delta times identity 
	int j = 0;
    for(int i = 0; i < n; i++)
    	nnzColsK[j++] += 1 + nzColA[i] + nzColG[i];
    for(int i = 0; i < p; i++)    
    	nnzColsK[j++] += 1 + nzRowA[i];
	for(int i = 0; i < m; i++)    
    	nnzColsK[j++] += 1 + nzRowG[i];
	//Reserve space for the nonzeros
	eigenKMat->reserve(nnzColsK);
    
	//Now copy the values into K
	//Start with the diagonal	
	for(int i = 0; i < n; i++)
		eigenKMat->insert(i,i) = DELTA;

	for(int i = 0; i < n; i++ )
		for(Eigen::SparseMatrix<double>::InnerIterator it(eA,i); it; ++it)	
		{
			int c,r;
			c = it.col();
			r = it.row();
			eigenKMat->insert(r+n,c) = -it.value();
		}
	for(int i = 0; i < n; i++ )
		for(Eigen::SparseMatrix<double>::InnerIterator it(eG,i); it; ++it)	
		{
			int c,r;
			c = it.col();
			r = it.row();
			eigenKMat->insert(r+n+p,c) = -it.value();	
		}

	for(int i = 0; i < n; i++ )
		for(Eigen::SparseMatrix<double>::InnerIterator it(eA,i); it; ++it)	
		{
			int c,r;
			r = it.col();
			c = it.row();
			eigenKMat->insert(r,n+c) = it.value();
		}	

	for(int i = 0; i < p; i++ )		
			eigenKMat->insert(n+i,n+i) = -DELTA;


	for(int i = 0; i < n; i++ )
		for(Eigen::SparseMatrix<double>::InnerIterator it(eG,i); it; ++it)	
		{
			int c,r;
			r = it.col();
			c = it.row();
			eigenKMat->insert(r,n+p+c) = it.value();
		}	

	//These correspond to the diagonals of the hessian
	for(int i = 0; i < m; i++ )		
		eigenKMat->insert(n+p+i,n+p+i) = -1.0 - DELTA;
     
}

k_newton_copl_matrix::~k_newton_copl_matrix()
{
	delete(hessianIx);
	delete(eigenKMat);
}


void k_newton_copl_matrix::solve(copl_vector &solution, copl_vector &rhs) {
	if(!isFactored) {
		std::cout << "TRIED TO USE SOLVE BEFORE UPDATE";
		throw new std::exception();
	}

	int m = rhs.size();
	Eigen::Map<Eigen::VectorXd> rhsEigen(&rhs[0],m), solEigen(&solution[0],m);	
	//TODO: log this
	solEigen = solver.solve(rhsEigen);
	//TODO: iterative refinement with MINRES
}

void k_newton_copl_matrix::solve(lp_direction &dir, linear_system_rhs& rhs, lp_variables &variables) {
	OUTPUT << "IN k_newton_copl_matrix::solve now..." << endl;

}


void k_newton_copl_matrix::update(lp_variables &variables)
{

	//Update the entries in the diagonal starting 
	//use the hessianIx to find the entries that need updating cheaply
	if(!eigenKMat->isCompressed()) 
	{		
		//TODO: log error
		std::cout << "CALLED UPDATE ON UNCOMPRESSED MATRIX \n";
    	throw new std::exception();
    }

	double* Kvals = eigenKMat->valuePtr();
    int m = variables.s.size();
	for(int j = 0; j < m; j++)
	{
	    Kvals[ (*hessianIx)[j] ] = variables.s[j]/variables.z[j]-DELTA;
	}
	//Real deal
	solver.factorize(*eigenKMat);
	//Solve the fixed rhs
	isFactored = true;
}

//Homogeneous solver implementation 
 homogeneous_solver::homogeneous_solver(copl_matrix &A, 
                       				    copl_matrix &G,
                                        copl_vector &c,
                                        copl_vector &b,
                                        copl_vector &h):k_newton_copl_matrix(A,G),
                                        				_c(c),
 														_h(h),
 														_b(b),
 														rhs_1(m+n+p),
 														rhs_2(m+n+p),
 														sol_1(m+n+p),
 														sol_2(m+n+p)
 														 {

//Build the RHS for the first system -[c; b; h]'

 	int j = 0;
 	for(int i = 0; i < n; i++)
 		rhs_1[j++] = -c[i];
 	for(int i = 0; i < p; i++)
 		rhs_1[j++] = -b[i]; 	
 	for(int i = 0; i < m; i++) 		
 		rhs_1[j++] = -h[i]; 	

}

void homogeneous_solver::update(lp_variables &variables) {
	k_newton_copl_matrix::update(variables);
	//Solve the first rhs system
	//k_newton_copl_matrix::solve(sol_1,rhs_2);
}

void homogeneous_solver::solve(lp_direction &dir, linear_system_rhs& rhs, lp_variables &variables) {
	OUTPUT << "IN homogeneous_solver::solve now..." << endl;
	k_newton_copl_matrix::solve(sol_1,rhs_1);//can be moved to homogeneous_solver::update
	
	this->build_rhs2(rhs, variables);
	k_newton_copl_matrix::solve(sol_2,rhs_2);
	
	double q8 = rhs.q4-rhs.q6/variables.tau;
	dir.dtau = (-q8 + dotat(_c, sol_2, 0) + dotat(_b, sol_2, n) + dotat(_h, sol_2, n+p))/
			(variables.kappa/variables.tau - dotat(_c, sol_1, 0) - dotat(_b, sol_1, n) - dotat(_h, sol_1, n+p) );


	addat(dir.dtau, sol_1, sol_2, dir.dx, 0, 0);
	addat(dir.dtau, sol_1, sol_2, dir.dy, n, n);
	addat(dir.dtau, sol_1, sol_2, dir.dz, n+p, n+p);
	dir.dkappa = (rhs.q6-dir.dtau*variables.kappa)/variables.tau;
	
	for(int i = 0; i < m ; i++)
	 		dir.ds[i] = (rhs.q5[i] - dir.dz[i]*variables.s[i])/variables.z[i]; 	
	
}

void homogeneous_solver::build_rhs2(linear_system_rhs& rhs, lp_variables &variables) {
	
	//Build the RHS for the second system -[c; b; h]'

	 	int j = 0;
	 	for(int i = 0; i < n; i++)
	 		rhs_2[j++] = -rhs.q1[i];
	 	for(int i = 0; i < p; i++)
	 		rhs_2[j++] = -rhs.q2[i]; 	
	 	for(int i = 0; i < m; i++) 		
	 		rhs_2[j++] = -(rhs.q3[i]-rhs.q5[i]/variables.z[i]); 	

	
}

}
