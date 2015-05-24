#include <copl_core.h>
#include <copl_linalg.h>
#include <copl_newton.h>
#include <iostream>


using namespace std;

namespace copl_ip {

ofstream LOG_FILE_VARIABLE;

lp_input::lp_input(copl_matrix *_A, copl_vector* _b, copl_vector* _c, copl_matrix* _G, copl_vector* _h):A(_A),G(_G),c(_c),b(_b),h(_h) {
  if(A!=NULL)
  	k_var = A->num_rows();
  m = G->num_rows();
  n = G->num_cols();
}

void copl_vector_dump(copl_vector &vec) {
	
	OUTPUT << "{";
	for(int i = 0; i < vec.size() - 1; i++) {
		OUTPUT << vec.at(i) << ",";
	}
	OUTPUT << vec[vec.size() - 1];
	OUTPUT << "}";
	OUTPUT << endl;
};

void lp_input::var_dump()  {
	
	OUTPUT << "=== BEGIN LP INPUT VAR DUMP ====" << endl;
	OUTPUT << "c = ";
	copl_vector_dump(*c);
	OUTPUT << "b = ";
	copl_vector_dump(*b);
	OUTPUT << "h = ";
	copl_vector_dump(*h);
	
	OUTPUT << "A = ";
	A->var_dump();
	OUTPUT << "G = ";
	G->var_dump();
	
	// ********* complete for matricies
	OUTPUT << "=== END VAR DUMP ====" << endl;
};

//--------End lp_input--------

// lp_settings
lp_settings::lp_settings (int input_max_iter, double input_linear_feas_tol, double input_comp_tol, double input_bkscale) {
	max_iter 			= input_max_iter;
	linear_feas_tol 	= input_linear_feas_tol;
	comp_tol 			= input_comp_tol;
	bkscale 			= input_bkscale;
}
int lp_settings::get_max_iter()				{return max_iter;}
double lp_settings::get_linear_feas_tol() 	{return linear_feas_tol;}
double lp_settings::get_comp_tol() 			{return comp_tol;}
double lp_settings::get_bkscale() 			{return bkscale;}

//--------End lp_settings--------

// lp_residuals
lp_residuals::lp_residuals( lp_input &problem_data) : r1(problem_data.n,0.0), r2(problem_data.k_var,0.0), r3(problem_data.m,0.0) {

}

/*void lp_residuals::update_norms() {
	r1_norm = norm2(r1);
	r2_norm = norm2(r2);
	r3_norm = norm2(r3);
}*/

void lp_residuals::compute_residuals( lp_input &problem_data, lp_variables &variables){
	// r1 = -pd.A'*variables.y - pd.G'*variables.z - pd.c*variables.tau;
	zeros(r1);
	sp_dgemtv(-1.0, 1.0, *problem_data.A, variables.y, r1);
	sp_dgemtv(-1.0, 1.0, *problem_data.G, variables.z, r1);
	axpy(-variables.tau, *problem_data.c, r1);
	// r2 = pd.A*variables.x - pd.b*variables.tau;
	zeros(r2);
	sp_dgemv(1.0, 1.0, *problem_data.A, variables.x, r2);
	axpy(-variables.tau, *problem_data.b, r2);
	//copl_vector_dump(r2);
	//copl_vector_dump(problem_data.b);
	
	// r3 = variables.s + pd.G*variables.x - variables.tau*pd.h;
	zeros(r3);
	axpy(1.0, variables.s, r3);
	sp_dgemv(1.0, 1.0, *problem_data.G, variables.x, r3);
	axpy(-variables.tau, *problem_data.h, r3);

	//r4 = variables.kappa + pd.c'*variables.x + pd.b'*variables.y + + pd.h'*variables.z;
	r4 = variables.kappa;
	r4 += dot(*problem_data.c, variables.x);
	r4 += dot(*problem_data.b, variables.y);
	r4 += dot(*problem_data.h, variables.z);
}

void lp_residuals::var_dump() {
	OUTPUT << "DUMP RESIDUALS OBJECT" << endl;
	OUTPUT << "r1:";
	copl_vector_dump(r1);
	
	OUTPUT << "r2:";
	copl_vector_dump(r2);
	
	OUTPUT << "r3:";
	copl_vector_dump(r3);
	
	OUTPUT << "r4:" << r4 << endl;
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
lp_direction::lp_direction(const lp_direction &obj){
	#ifdef PREVENT_COPY_CONSTRUCTOR
	OUTPUT << "WARNING: Direction copy constructor called." << "\n";
	#endif
}
double lp_direction::get_alpha() { return alpha; }
double lp_direction::get_dtau() { return dtau; }
double lp_direction::get_dkappa() { return dkappa; } 
copl_vector lp_direction::get_dx() { return dx; }
copl_vector lp_direction::get_dy() { return dy; } 
copl_vector lp_direction::get_dz() { return dz; }
copl_vector lp_direction::get_ds() { return ds; }


void lp_direction::solve_linear_system_for_new_direction(linear_system_rhs& rhs, k_newton_copl_matrix& K_matrix, lp_variables &variables) {
	K_matrix.solve(*this, rhs, variables);
}

void  lp_direction::compute_direction(
		linear_system_rhs &affine_rhs,
		lp_input &problem_data,
		lp_variables &variables,
		algorithm_state &state,
		lp_settings &settings,
		k_newton_copl_matrix &K_matrix
		) {
		
	this->solve_linear_system_for_new_direction(affine_rhs, K_matrix, variables);
	this->compute_step_size(variables,settings);
	/*
	this.compute_affine_direction = function(affine_rhs::class_linear_system_rhs,
													problem_data::class_linear_program_input,	
													variables::class_linear_program_variables,	
													K_newton_matrix::class_K_newton_matrix)
			
		dir = solveLinearEquation(problem_data, variables, affine_rhs, K_newton_matrix)
		
		m = problem_data.m
		n = problem_data.n
		k = problem_data.k
		
		x = variables.x
		z = variables.z
		s = variables.s
		y = variables.y
		tau = variables.tau
		kappa = variables.kappa
		
		dx_a = dir[1:k];
		dy_a = dir[(k+1):(k+n)];
		dz_a = dir[(k+n+1):(k+n+m)];
		dtau_a = dir[(k+n+m+1)];
		ds_a = ( -z.*s - dz_a.*s)./z;
		dkappa_a = (-(tau)*(kappa) - dtau_a*(kappa))/(tau)\
		
		this.update_values(dx_a,dy_a,dz_a,dtau_a,ds_a,dkappa_a,alpha)
	end*/	
	
	// figure out alpha (line search value)
	
}


void lp_direction::compute_step_size(lp_variables& variables, lp_settings& settings) {
	this->alpha = 1;
	compute_min_ratio_alpha(variables.s,this->ds,this->alpha); // TO DO change to dz
	compute_min_ratio_alpha(variables.z,this->dz,this->alpha);
	compute_min_ratio_alpha(variables.kappa,this->dkappa,this->alpha);
	compute_min_ratio_alpha(variables.tau,this->dtau,this->alpha);
	this->alpha = this->alpha*settings.bkscale;
}

void lp_direction::compute_min_ratio_alpha(double var, double dvar, double& alpha_val) {
	if (dvar != 0) { // GREATER THAN TOL ??????
		double candidate_alpha = -var/dvar;
		if (candidate_alpha > 0) {
			alpha_val = min(alpha_val, candidate_alpha);
		}
	}
}

void lp_direction::compute_min_ratio_alpha(copl_vector &var, copl_vector &dvar, double& alpha_val) {
	assert(var.size() == dvar.size());
	
	for (int i = 0; i < var.size(); i++) {
		double var_double = var[i];
		double dvar_double = var[i];
		compute_min_ratio_alpha(var_double, dvar_double, alpha_val);
	}
}

void lp_direction::var_dump() {
	
	OUTPUT << "DUMP DIRECTION OBJECT" << endl;
	OUTPUT << "dx ";
	copl_vector_dump(dx);
	OUTPUT << "dy ";
	copl_vector_dump(dy);
	OUTPUT << "dz ";
	copl_vector_dump(dz);
	OUTPUT << "ds ";
	copl_vector_dump(ds);
	OUTPUT << "dtau " << dtau << endl;
	OUTPUT << "dkappa " << dkappa << endl;
	OUTPUT << "alpha " << alpha << endl;
	
}

//-----------End lp direction


// lp_variables
lp_variables::lp_variables(int n, int m, int k_var) :
	x(n,0.0), s(m,1.0), z(m,1.0), y(k_var,0.0) {
	tau = 1;
	kappa = 1;
}

lp_variables::lp_variables(const lp_variables &obj){
	#ifdef PREVENT_COPY_CONSTRUCTOR
	OUTPUT << "WARNING: Variables copy constructor called." << "\n";
	#endif
}

lp_variables::~lp_variables() {
	OUTPUT << "lp variables being deleted" << endl;
	OUTPUT << "destructor not yet complete" << endl;
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
	gap = dot(*problem_data.c,variables.x) + dot(*problem_data.h,variables.z) + dot(*problem_data.b,variables.y);
}

void algorithm_state::update_mu(lp_variables &variables,  lp_input &problem_data){
	// ( ((variables.s)'*(variables.z) + (variables.tau)*(variables.kappa))/(problem_data.m+1))
	mu = ( dot(variables.s,variables.z) + variables.tau*variables.kappa )/( problem_data.m + 1.0 );
}
//--------End algorithm_state--------


// linear_system_rhs

linear_system_rhs::linear_system_rhs( lp_input &problem_data) :
	q1(problem_data.n,0), q2(problem_data.k_var,0), q3(problem_data.m,0), q5(problem_data.m,0)
{
	
}

linear_system_rhs::linear_system_rhs(const linear_system_rhs &obj){
	#ifdef PREVENT_COPY_CONSTRUCTOR
	OUTPUT << "WARNING: Linear system rhs copy constructor called." << "\n";
	#endif
}
linear_system_rhs::~linear_system_rhs() {
	OUTPUT << "destructor for linear system rhs called" << endl;
}

void linear_system_rhs::compute_affine_rhs(lp_residuals &residuals, lp_variables &variables){
	// TOOD	
	
	// q1 = -r1
	zeros(q1);
	axpy(-1.0, residuals.r1, q1);
	
	// q2 = -r2
	zeros(q2);
	axpy(-1.0, residuals.r2, q2);
	
	// q3 = -r3
	zeros(q3);
	axpy(-1.0, residuals.r3, q3);
	
	// q4 = -r4
	q4 = -residuals.r4;
	
	// q5 = -(variables.z).*(variables.s)
	zeros(q5);
	for (int i = 0; i < q5.size(); i++)
		q5[i] = -variables.z[i]*variables.s[i];
	
	// q6 = -(variables.tau)*(variables.kappa)
	q6 = -variables.tau*variables.kappa;
}
void linear_system_rhs::compute_corrector_rhs(lp_residuals &residuals, lp_variables &variables, algorithm_state &state, lp_direction &direction,  lp_input &problem_data){
	// TODO
	//
	
	// mu_a = ((s+alpha*ds_a)'*(z+alpha*dz_a) + (tau + alpha*dtau_a)*((kappa) + alpha*dkappa_a))/(m+1.0);
	// mu_a = ((s+alpha*ds_a)'*(z+alpha*dz_a) + (tau + alpha*dtau_a)*((kappa) + alpha*dkappa_a))/(m+1.0);
	
	
	/*
	this.compute_corrector_rhs = function(residuals::class_residuals,variables::class_linear_program_variables,state::class_algorithm_state,affine_direction::class_direction,problem_data::class_linear_program_input)
		m = problem_data.m
		
		z = variables.z
		s = variables.s
		tau = variables.tau
		kappa = variables.kappa
		
		dx_a = direction.dx
		dy_a = direction.dy
		ds_a = direction.ds
		dz_a = direction.dz
		dtau_a = direction.dtau
		dkappa_a = direction.dkappa
		alpha = direction.alpha
		
		mu = state.mu
		
		mu_a = ((s+alpha*ds_a)'*(z+alpha*dz_a) + (tau + alpha*dtau_a)*((kappa) + alpha*dkappa_a))/(m+1.0);
		sigma = ((mu_a/(mu))^3)[1]
		
		state.sigma = sigma
		
		this.update_values(-(1-sigma)*residuals.r1, -(1-sigma)*residuals.r2, -(1-sigma)*residuals.r3, -(1-sigma)*residuals.r4, -z.*s -ds_a.*dz_a + sigma*mu,  -tau*kappa-dtau_a*dkappa_a + sigma*mu)
	end
	*/
}

void linear_system_rhs::var_dump() {
	OUTPUT << "DUMP RHS OBJECT" << endl;
	OUTPUT << "q1";
	copl_vector_dump(q1);
	OUTPUT << "q2";
	copl_vector_dump(q2);
	OUTPUT << "q3";
	copl_vector_dump(q3);
	OUTPUT << "q4" << q4 << endl;
	OUTPUT << "q5";
	copl_vector_dump(q5);
	OUTPUT << "q6" << q6 << endl;
}
//--------End linear_system_rhs--------


#ifdef random
lp_input* copl_utility::loadFromUF(string UF_group, string UF_name){
	lp_input* problem_data;
	
	// Inserted an script that download file from UF repository
    ifstream A_File("../example_problems/ex3sta1/ex3sta1.mtx");
	
  	if (!A_File) {
  		cerr << "Error Loading from UF dataset. File not found (A Matrix)" << endl;
    		return ;
	}
	
	int idx = 0;
	string line;
    bool first_line = true;
	
 	while (getline(A_File, line)) {
  		if (line.empty()) continue;
                if (line.at(0) == '%') continue;

		std::istringstream ss(line);
		std::string token;
		string indata[3] = {"", "", ""};

		idx = 0;
		while(std::getline(ss, token, ' ')) {
			indata[idx++] = token;
		}	
		if (first_line){
                	first_line = false;
			int m = stoi(indata[0]);
			int n = m;
			int k_var = stoi(indata[1]);
                        problem_data = new lp_input(m,n,k_var);
                }else{
			int row = stoi(indata[0]) - 1;
			int col = stoi(indata[1]) - 1;
			double value = stod(indata[2]); 
			problem_data->A.insert_at(row,col,value);
		}
		
	}
	A_File.close();


	ifstream b_File("../example_problems/ex3sta1/ex3sta1_b.mtx");
	if (!b_File) {
			cerr << "Error Loading from UF dataset. File not found (b vector)" << endl;
			return ;
	}

	first_line = true;
	int col_idx = 0;
	while (getline(b_File, line)) {
			if (line.empty()) continue;
			if (line.at(0) == '%') continue;

			std::istringstream ss(line);
			std::string token;

			string indata[1] = {""};
			idx = 0;
			while(std::getline(ss, token, ' ')) {
					indata[idx++] = token;
			}
			if (first_line){
					first_line = false;
			}else{
					double value = std::stod(indata[0]);
					problem_data->b[col_idx++] = value;
			}

	}
	b_File.close();


	ifstream c_File("../example_problems/ex3sta1/ex3sta1_c.mtx");
	if (!c_File) {
			cerr << "Error Loading from UF dataset. File not found (c vector)" << endl;
			return ;
	}

	first_line = true;
	col_idx = 0;
	while (getline(c_File, line)) {
			if (line.empty()) continue;
			if (line.at(0) == '%') continue;

			std::istringstream ss(line);
			std::string token;                

			string indata[1] = {""};
			idx = 0;
			while(std::getline(ss, token, ' ')) {
					indata[idx++] = token;
			}
			if (first_line){
					first_line = false;
			}else{
					double value = std::stod(indata[0]);
					problem_data->c[col_idx++] = value;
			}

	}
	c_File.close();

	return problem_data;
}
#endif  


}
