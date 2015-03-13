#include <copl_newton.h>

namespace copl_ip{

k_newton_copl_matrix::k_newton_copl_matrix( lp_input &problem_data)
{
	
	int n = problem_data.A.num_cols();		
	int p = problem_data.A.num_rows();	
	int m = problem_data.G.num_rows();	

	//Assemble 	
	assemble_matrix(problem_data.A, problem_data.G);

	//Now find the indices of the hessian
	hessianIx = new std::vector<int>(m);
	eigenKMat = new EigenSpMat_t(n+m+p,n+m+p);
	//Compress 
	eigenKMat->makeCompressed();	
	//The indices are the last nonzero per column for the cols n+p to n+p+m-1
	int* outerIx = eigenKMat->outerIndexPtr();
	for(int i = n+p; i < n+p+m; i++ )
		(*hessianIx)[i] = outerIx[i+1]-1;

	//Now do the symbolic analysis of the matrix 
	cout << "Will analyze pattern\n";
	cout.flush();
	solver.analyzePattern(*eigenKMat);
}

void k_newton_copl_matrix::assemble_matrix(copl_matrix &A, copl_matrix &G)
{

	//Assemble 
	int n = A.num_cols();		
	int p = A.num_rows();	
	int m = G.num_rows();	

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
			eigenKMat->insert(r+n,c) = it.value();
		}
	for(int i = 0; i < n; i++ )
		for(Eigen::SparseMatrix<double>::InnerIterator it(eG,i); it; ++it)	
		{
			int c,r;
			c = it.col();
			r = it.row();
			eigenKMat->insert(r+n+p,c) = it.value();	
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

	int m = variables.s.size();
	double* Kvals = eigenKMat->valuePtr();
	for(int j = 0; j < m; j++)
	{
		Kvals[ (*hessianIx)[j] ] = -variables.s[j]/variables.z[j]-DELTA;
	}
	//Real deal
	solver.factorize(*eigenKMat);
	//Solve the fixed rhs
	isFactored = true;
}

}
