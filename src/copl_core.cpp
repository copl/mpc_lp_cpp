#include <copl_core.h>
#include <copl_linalg.h>
#include <iostream>
#include <fstream>
using namespace std;

namespace copl_ip {

void copl_vector_dump(copl_vector &vec) {
	cout << "{";
	for(int i = 0; i < vec.size() - 1; i++) {
		cout << vec[i] << ",";
	}
	cout << vec[vec.size() - 1];
	cout << "}";
	cout << endl;
};

//lp_input
lp_input::lp_input(int _m, int _n, int _k_var) // generates things at random this must change!!!
	: A(_n,_k_var), G(_m,_k_var), c(_k_var,2.0), h(_m,2.0), b(_n,2.0) { // ************
	m = _m;
	n = _n;
	k_var = _k_var;
};


lp_input::~lp_input() {
	cout << "deleting lp input" << endl;
	cout << "constructor not complete" << endl;
};

void lp_input::var_dump()  {
	
	cout << "=== BEGIN LP INPUT VAR DUMP ====" << endl;
	cout << "c = ";
	copl_vector_dump(c);
	cout << "b = ";
	copl_vector_dump(b);
	cout << "h = ";
	copl_vector_dump(h);
	
	cout << "A = ";
	A.var_dump();
	cout << "G = ";
	G.var_dump();
	
	// ********* complete for matricies
	cout << "=== END VAR DUMP ====" << endl;
};

//--------End lp_input--------

// lp_settings
lp_settings::lp_settings (int input_max_iter, double input_linear_feas_tol, double input_comp_tol, double input_bkscale) {
	max_iter 			= input_max_iter;
	linear_feas_tol 	= input_linear_feas_tol;
	comp_tol 			= input_comp_tol;
	bkscale 			= input_bkscale;
}
int lp_settings::get_max_iter()			{return max_iter;}
double lp_settings::get_linear_feas_tol() 	{return linear_feas_tol;}
double lp_settings::get_comp_tol() 			{return comp_tol;}
double lp_settings::get_bkscale() 			{return bkscale;}

//--------End lp_settings--------

// lp_residuals
lp_residuals::lp_residuals( lp_input &problem_data) : r1(problem_data.k_var,0.0), r2(problem_data.n,0.0), r3(problem_data.m,0.0) {

}

/*void lp_residuals::update_norms() {
	r1_norm = norm2(r1);
	r2_norm = norm2(r2);
	r3_norm = norm2(r3);
}*/

void lp_residuals::compute_residuals( lp_input &problem_data, lp_variables &variables){
	/*
	cout << "A " << problem_data.A.num_rows() << "," << problem_data.A.num_cols() << endl;
	cout << "G " << problem_data.G.num_rows() << "," << problem_data.G.num_cols() << endl;
	cout << "r1 " << r1.size() << endl;
	cout << "y " << variables.y.size() << endl;
	cout << "z " << variables.z.size() << endl;
	cout << "x " << variables.x.size() << endl;
	cout << "c " << problem_data.c.size() << endl;
	*/
	
	// r1 = -pd.A'*variables.y - pd.G'*variables.z - pd.c*variables.tau;
	zeros(r1);
	sp_dgemtv(-1.0, 1.0, problem_data.A, variables.y, r1);
	sp_dgemtv(-1.0, 1.0, problem_data.G, variables.z, r1);
	axpy(-variables.tau, problem_data.c, r1);
	
	cout << "r1:";
	copl_vector_dump(r1);
	
	// r2 = pd.A*variables.x - pd.b*variables.tau;
	zeros(r2);
	sp_dgemv(1.0, 1.0, problem_data.A, variables.x, r2);
	axpy(-variables.tau, problem_data.b, r2);
	
	cout << "r2:";
	copl_vector_dump(r2);
	
	// r3 = variables.s + pd.G*variables.x - variables.tau*pd.h;
	zeros(r3);
	axpy(1.0, variables.s, r3);
	sp_dgemv(1.0, 1.0, problem_data.G, variables.x, r3);
	axpy(-variables.tau, problem_data.h, r3);
	
	cout << "r3:";
	copl_vector_dump(r3);
	
	//r4 = variables.kappa + pd.c'*variables.x + pd.b'*variables.y + + pd.h'*variables.z;
	r4 = variables.kappa;
	r4 += dot(problem_data.c, variables.x);
	r4 += dot(problem_data.b, variables.y);
	r4 += dot(problem_data.h, variables.z);
	
	cout << "r4:" << r4 << endl;
}

double lp_residuals::get_r1_norm(){return norm2(r1);}
double lp_residuals::get_r2_norm(){return norm2(r2);}
double lp_residuals::get_r3_norm(){return norm2(r3);}
double lp_residuals::get_norm_squared() {return -1.0;} // TO DO



//--------End lp_residuals--------

// lp direction
lp_direction::lp_direction(lp_variables &variables) 
	: dx(variables.x.size(),0.0), ds(variables.s.size(),1.0), dz(variables.z.size(),1.0), dy(variables.y.size(),0.0) 
{

}

double lp_direction::get_alpha() { return alpha; }
double lp_direction::get_dtau() { return dtau; }
double lp_direction::get_dkappa() { return dkappa; } 
copl_vector lp_direction::get_dx() { return dx; }
copl_vector lp_direction::get_dy() { return dy; } 
copl_vector lp_direction::get_dz() { return dz; }
copl_vector lp_direction::get_ds() { return ds; }

void  lp_direction::compute_affine_direction(linear_system_rhs &affine_rhs,
		 lp_input &problem_data,
		lp_variables &variables,
		k_newton_copl_matrix &K_matrix) {
}

void lp_direction::compute_corrector_direction(
		linear_system_rhs &corrector_rhs,
		 lp_input &problem_data,
		lp_variables &variables,
		algorithm_state &state,
		lp_settings &settings,
		k_newton_copl_matrix &K_matrix
		) {
			
}

//-----------End lp direction



// lp_variables
lp_variables::lp_variables(int n, int m, int k_var) :
	x(k_var,0.0), s(m,1.0), z(m,1.0), y(n,0.0) {
	tau = 1;
	kappa = 1;
}

lp_variables::~lp_variables() {
	cout << "lp variables being deleted" << endl;
	cout << "destructor not yet complete" << endl;
}

void lp_variables::take_step(lp_direction &direction){
	double alpha = direction.get_alpha();
	//TODO: Implement copl_vector.multiply and copl_vector.add
	// x = x.add(direction.get_dx.multiply(alpha));
	
	axpy(alpha,direction.dx,x);
	axpy(alpha,direction.ds,s);
	axpy(alpha,direction.dz,z);
	axpy(alpha,direction.dy,y);
	
	tau = tau + alpha * direction.get_dtau();
	kappa = kappa + alpha * direction.get_dkappa();

}

// algorithm_state

algorithm_state::algorithm_state() {
	mu = -1;
	sigma = -1;
	gap = -1;
}
void algorithm_state::update_gap(lp_variables &variables,  lp_input &problem_data){
	// ((problem_data.c)'*(variables.x) + (problem_data.h)'*(variables.z) + (problem_data.b)'*(variables.y))
	gap = dot(problem_data.c,variables.x) + dot(problem_data.h,variables.z) + dot(problem_data.b,variables.y);
}

void algorithm_state::update_mu(lp_variables &variables,  lp_input &problem_data){
	// ( ((variables.s)'*(variables.z) + (variables.tau)*(variables.kappa))/(problem_data.m+1))
	mu = ( dot(variables.s,variables.z) + variables.tau*variables.kappa )/( problem_data.m + 1.0 );
}
//--------End algorithm_state--------


// linear_system_rhs

linear_system_rhs::linear_system_rhs( lp_input &problem_data) :
	q1(1,0), q2(1,0), q3(1,0), q4(1,0), q5(1,0), q6(1,0)
{
	cout << "xx" << endl;
}

void linear_system_rhs::compute_affine_rhs(lp_residuals &residuals, lp_variables &variables){
	// TOOD
}
void linear_system_rhs::compute_corrector_rhs(lp_residuals &residuals, lp_variables &variables, algorithm_state &state, lp_direction &direction,  lp_input &problem_data){
	// TOOD
}
//--------End linear_system_rhs--------


k_newton_copl_matrix::k_newton_copl_matrix( lp_input &problem_data)
{

}
void k_newton_copl_matrix::factor(){

}
void k_newton_copl_matrix::solve(copl_vector &solution, copl_vector &rhs) {

}

void k_newton_copl_matrix::update(lp_variables &variables)
{
	
}



void copl_utility::loadFromUF(string UF_group, string UF_name, lp_input &problem_data){
	// Inserted an script that download file from UF repository
        ifstream A_File("../example_problems/ex3sta1/ex3sta1.mtx");
  	if (!A_File) {
  		cerr << "Error Loading from UF dataset. File not found." << endl;
    		return ;
	}

	string line;
 	while (getline(A_File, line)) {
  		if (line.empty()) continue;
                if (line.at(0) == '%') continue;

		std::istringstream ss(line);
		std::string token;

		while(std::getline(ss, token, ' ')) {
			std::cout << token << '\n';
		}	
	}

	A_File.close();
}


}
