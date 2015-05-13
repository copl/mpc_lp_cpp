#include <copl_newton.h>

namespace copl_ip{

int k_newton_copl_matrix::nnz()
{
	return eigenKMat->nonZeros();
}

k_newton_copl_matrix::k_newton_copl_matrix(int m, int n)
{
	cout << "Testing constructor (m,n) " <<m <<","<<n <<"\n";
	hessianIx = new std::vector<int>(m);
	eigenKMat = new copl_matrix(n,n);
    cout << "This constructor should only be used in testing";
    cout.flush();
}

k_newton_copl_matrix::k_newton_copl_matrix(copl_matrix &A, copl_matrix &G)
{
	
	k = A.rows();	
	n = G.cols();		
    m = G.rows();	

	//Now find the indices of the hessian
	hessianIx = new std::vector<int>(m);
	eigenKMat = new copl_matrix(n+m+k,n+m+k);
	
    //Assemble 	
	assemble_matrix(A, G);
    
    //Compress 
	eigenKMat->makeCompressed();	

    //The indices are the last nonzero per column for the cols n+p to n+p+m-1
	int* outerIx = eigenKMat->outerIndexPtr();		
	for(int i = 0; i < m; i++ )
		 (*hessianIx)[i] = outerIx[i+n+k+1]-1;	
    
    solver.analyzePattern(*eigenKMat);

}

void k_newton_copl_matrix::assemble_matrix(copl_matrix &A, copl_matrix &G)
{
	//Assemble 
	//Count the number of non zeros per row and col of A and G
	std::vector<int> nzRowA(k);		
	std::vector<int> nzRowG(m);	
	std::vector<int> nzColA(n);		
	std::vector<int> nzColG(n);	
	std::vector<int> nnzColsK(m+n+k);	
	//TODO this depends on EIGENs iterators and will need 
	//eventual refactoring	
	//Get the internal matrices 
	
    //Count the non zeros per column and row
	for(int i = 0; i < A.cols(); i++ )
		for(Eigen::SparseMatrix<double>::InnerIterator it(A,i); it; ++it)	
		{
			nzRowA[it.row()]++;
			nzColA[it.col()]++;
		}
	for(int i = 0; i < G.cols(); i++ )
		for(Eigen::SparseMatrix<double>::InnerIterator it(G,i); it; ++it)	
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
    for(int i = 0; i < k; i++)    
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
		for(Eigen::SparseMatrix<double>::InnerIterator it(A,i); it; ++it)	
		{
			int c,r;
			c = it.col();
			r = it.row();
			eigenKMat->insert(r+n,c) = it.value();
		}
	for(int i = 0; i < n; i++ )
		for(Eigen::SparseMatrix<double>::InnerIterator it(G,i); it; ++it)	
		{
			int c,r;
			c = it.col();
			r = it.row();
			eigenKMat->insert(r+n+k,c) = it.value();	
		}

	for(int i = 0; i < n; i++ )
		for(Eigen::SparseMatrix<double>::InnerIterator it(A,i); it; ++it)	
		{
			int c,r;
			r = it.col();
			c = it.row();
			eigenKMat->insert(r,n+c) = it.value();
		}	

	for(int i = 0; i < k; i++ )		
			eigenKMat->insert(n+i,n+i) = -DELTA;


	for(int i = 0; i < n; i++ )
		for(Eigen::SparseMatrix<double>::InnerIterator it(G,i); it; ++it)	
		{
			int c,r;
			r = it.col();
			c = it.row();
			eigenKMat->insert(r,n+k+c) = it.value();
		}	

	//These correspond to the diagonals of the hessian
	for(int i = 0; i < m; i++ )		
			eigenKMat->insert(n+k+i,n+k+i) = -1.0 - DELTA;
     
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
    
	solution = solver.solve(rhs);
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

	double* Kvals = eigenKMat->valuePtr();
    int m = variables.s.size();
	for(int j = 0; j < m; j++)
	{
	    Kvals[ (*hessianIx)[j] ] = -variables.s[j]/variables.z[j]-DELTA;
	}
	//Real deal
	solver.factorize(*eigenKMat);
	isFactored = true;
}

//Homogeneous solver implementation 
 homogeneous_solver::homogeneous_solver(lp_input &prob):_c(prob.c),
							_h(prob.h),
							_b(prob.b),
							rhs_1(m+n+k),
							k_newton_copl_matrix(prob.A,prob.G) {
    //Build the RHS for the first system [-c; b; h]'
    rhs_1.segment(0,n)   = -prob.c;
    rhs_1.segment(n,k)   = prob.b;
    rhs_1.segment(n+k,m) = prob.h;

}

/* Updates the values of the hessian block of the augmented system and 
 * solves the equation 
 * [  A' G'][dx]   [-c]
 * [A      ][dy] = [b]
 * [G   -H ][dz]   [h]
 * And solves the solution in sol_1 */

void homogeneous_solver::update(lp_variables &variables) {
	k_newton_copl_matrix::update(variables);
        tau = variables.tau;
    	kappa = variables.kappa;

	//Solve the first rhs system
	k_newton_copl_matrix::solve(sol_1,rhs_1);
 
        // dtau_denom = kap/tau + (c'*x1 + by1 + h'*z1); 
	dtau_denom = kappa/tau 
        + _c.dot(sol_1.segment(0,n)) +_b.dot(sol_1.segment(n,k)) +_h.dot(sol_1.segment(n+k,m));
}

/* Warning mutates rhs. 
 * This method solves the homogeneous equation defined by var and rhs and saves the result on dir.
 * It assumes that sol_1 is populated with the solution of the system 
 * [  A' G'][dx]   [-c]
 * [A      ][dy] = [b]
 * [G   -H ][dz]   [h]
 * Given the solution */

void homogeneous_solver::solve(lp_direction &dir, linear_system_rhs& rhs, lp_variables &var) {
    
     //Solve for the right hand side
     k_newton_copl_matrix::solve(sol_2,rhs.q123); 
     
     // dtau = (-q4+q6 + c'*x2 + by2 + h'*z2)/dtau_denom
     dir.dtau = -rhs.q4 + rhs.q6 + _c.dot(sol_2.segment(0,n)) + _b.dot(sol_2.segment(n,k)) + _h.dot(sol_2.segment(n+k,m));
     dir.dtau /= dtau_denom;

     //d2+dtau d1
     sol_2+=dir.dtau*sol_1;
     //sol_2 contains dx,dy,dz, copy them to the solution. 
     dir.dx = sol_2.segment(0,n);
     dir.dy = sol_2.segment(n,k);
     dir.dz = sol_2.segment(n+k,m);

 }


//Solves with the system
//[0   A'   G'  c]dz      q1
//[A           -b]dy    = q2
//[G       -H  -h]dz      q3
//[-c' -b' -h k/t]dt      q4 

void homogeneous_solver::solve_reduced(lp_direction &dir, linear_system_rhs &rhs) {
     //Solve with the right hand side formed with the top of the reduced rhs
     k_newton_copl_matrix::solve(sol_2,rhs.q123); 
     
     // dtau  = (q4+q6 + c'*x2 + by2 + h'*z2)/dtau_denom
     dir.dtau = rhs.q4 + _c.dot(sol_2.segment(0,n)) + _b.dot(sol_2.segment(n,k)) + _h.dot(sol_2.segment(n+k,m));
     dir.dtau /= dtau_denom;

     //d2+dtau d1
     sol_2+=dir.dtau*sol_1;
     //sol_2 contains dx,dy,dz, copy them to the solution. 
     dir.dx = sol_2.segment(0,n);
     dir.dy = sol_2.segment(n,k);
     dir.dz = sol_2.segment(n+k,m);
}



void homogeneous_solver::reduce_rhs(linear_system_rhs  &rhs)
{
   //Now calculate -r3+s and place it in the last entry of q123
    rhs.q123.segment(n+k,m) -= rhs.q5;
    //And the same for the scalar entry of the last row
    rhs.q4 -= rhs.q6;
}

void homogeneous_solver::back_substitute(lp_direction &dir, linear_system_rhs  &rhs, lp_variables &var)
{
    //Compute s from Hdz + ds = q5
     dir.ds = rhs.q5-(dir.dz.array()*var.s.array()/var.z.array()).matrix();
     //Compute dk from k/t dt + dk = q6
     dir.dkappa = rhs.q6 - kappa/tau*dir.dtau;
}
       

}
