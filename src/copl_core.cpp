#include <copl_core.h>
#include <copl_linalg.h>
#include <iostream>


using namespace std;

namespace copl_ip {

ofstream LOG_FILE_VARIABLE;

lp_input::lp_input(copl_matrix &_A, 
		   copl_external_vector &_b,
		   copl_external_vector &_c,
		   copl_matrix &_G,
		   copl_external_vector &_h):A(_A),G(_G),c(_c),b(_b),h(_h) {
  k_var = A.rows();
  m = G.rows();
  n = G.cols();
  //TODO: This metho should test that the sizes are consistent!
  //There should be a unittest for that
}

void lp_input::var_dump()  {
	
	OUTPUT << "=== BEGIN LP INPUT VAR DUMP ====" << endl;
	OUTPUT << "c = " << c;
	OUTPUT << "b = " << b;
	OUTPUT << "h = " << h;
	
	OUTPUT << "A = " << A;
	OUTPUT << "G = " << G;
	
	// ********* complete for matricies
	OUTPUT << "=== END VAR DUMP ====" << endl;
};

// lp_variables
	lp_variables::lp_variables(copl_external_vector &_x,
				     copl_external_vector &_y,
				     copl_external_vector &_s,
				     copl_external_vector &_z, double &_tau, double &_kappa):
	x(_x),y(_y),s(_s),z(_z),tau(_tau),kappa(_kappa){
 
	x.setZero();
    	y.setZero();
    	s.setConstant(1.0);
    	z.setConstant(1.0);
	tau   = 1.0;
	kappa = 1.0;
}

lp_variables::lp_variables(const lp_variables &obj):
	x(obj.x),y(obj.y),s(obj.s),z(obj.z),tau(obj.tau),kappa(obj.kappa){
	#ifdef PREVENT_COPY_CONSTRUCTOR
	OUTPUT << "WARNING: Variables copy constructor called." << "\n";
	#endif
}

lp_variables::~lp_variables() {
	OUTPUT << "lp variables being deleted" << endl;
	OUTPUT << "destructor not yet complete" << endl;
}

void lp_variables::take_step(lp_direction &direction){
	double alpha = direction.alpha;
    x = x+alpha*direction.dx;	
    y = y+alpha*direction.dy;
    z = z+alpha*direction.dz;
    s = s+alpha*direction.ds;
	
	tau = tau + alpha * direction.dtau;
	kappa = kappa + alpha * direction.dkappa;

}

// lp direction
lp_direction::lp_direction(lp_variables &variables) 
	: dx(variables.x.size()), ds(variables.s.size()), dz(variables.z.size()), dy(variables.y.size()) 
{

}

lp_direction::lp_direction(const lp_direction &obj){
	#ifdef PREVENT_COPY_CONSTRUCTOR
	OUTPUT << "WARNING: Direction copy constructor called." << "\n";
	#endif
}

void lp_direction::solve_linear_system_for_new_direction(linear_system_rhs& rhs, k_newton_copl_matrix& K_matrix) {

}

//void  lp_direction::compute_direction(
//		linear_system_rhs &affine_rhs,
//		lp_input &problem_data,
//		lp_variables &variables,
//		algorithm_state &state,
//		lp_settings &settings,
//		k_newton_copl_matrix &K_matrix
//		) {
//		
//	this->solve_linear_system_for_new_direction(affine_rhs, K_matrix);
//	this->compute_step_size(variables,settings);
//	/*
//	this.compute_affine_direction = function(affine_rhs::class_linear_system_rhs,
//													problem_data::class_linear_program_input,	
//													variables::class_linear_program_variables,	
//													K_newton_matrix::class_K_newton_matrix)
//			
//		dir = solveLinearEquation(problem_data, variables, affine_rhs, K_newton_matrix)
//		
//		m = problem_data.m
//		n = problem_data.n
//		k = problem_data.k
//		
//		x = variables.x
//		z = variables.z
//		s = variables.s
//		y = variables.y
//		tau = variables.tau
//		kappa = variables.kappa
//		
//		dx_a = dir[1:k];
//		dy_a = dir[(k+1):(k+n)];
//		dz_a = dir[(k+n+1):(k+n+m)];
//		dtau_a = dir[(k+n+m+1)];
//		ds_a = ( -z.*s - dz_a.*s)./z;
//		dkappa_a = (-(tau)*(kappa) - dtau_a*(kappa))/(tau)\
//		
//		this.update_values(dx_a,dy_a,dz_a,dtau_a,ds_a,dkappa_a,alpha)
//	end*/	
//	
//	// figure out alpha (line search value)
//	
//}
//--------End lp_input--------

// lp_settings
lp_settings::lp_settings (int _max_iter,
                          double _linear_feas_tol,
                          double _comp_tol,
                          double _bkscale,
			  double _regularization):max_iter(_max_iter),
						  linear_feas_tol(_linear_feas_tol),
						  comp_tol(_comp_tol),
						  bkscale(_bkscale),
						  regularization(_regularization){}

//--------End lp_settings--------

// lp_residuals
lp_residuals::lp_residuals( lp_input &problem_data) : r1(problem_data.n),
                                                      r2(problem_data.k_var),
                                                      r3(problem_data.m) {


}

void lp_residuals::compute_residuals( lp_input &problem_data, lp_variables &variables){
    r1 = problem_data.A.transpose()*variables.y +
         problem_data.G.transpose()*variables.z;

    hn1 = r1.norm();
    r1 += problem_data.c*variables.tau;
    n1  = r1.norm();

    r2 = -problem_data.A*variables.x;	
    hn2 = r2.norm();
    r2 += problem_data.b*variables.tau;
    n2 = r2.norm();

    r3 = -variables.s - problem_data.G*variables.x;
    hn3 = r3.norm();
    r3 += variables.tau*problem_data.h;
    n3 = r3.norm();

    r4 = -variables.kappa - problem_data.c.dot(variables.x) 
                          - problem_data.b.dot(variables.y) 
                          - problem_data.h.dot(variables.z);
    n4 = fabs(r4);
    
}

void lp_residuals::var_dump() {
	OUTPUT << "DUMP RESIDUALS OBJECT" << endl;
	OUTPUT << "r1:" << r1 << endl;
	OUTPUT << "r2:" << r2 << endl;
	OUTPUT << "r3:" << r3 << endl;
	OUTPUT << "r4:" << r4 << endl;
}

////--------End lp_residuals--------
//TODO: Testing routines for these three methods
void lp_direction::compute_step_size(lp_variables &variables, 
                                     lp_settings  &settings) {
	double var, dvar, ratio;
	copl_vector s = variables.s;
	copl_vector z = variables.z;

	this->alpha = 1.0;
	for(int i=0; i<variables.s.size(); i++)
	{
		ratio = -s[i]/ds[i];
		this->alpha = ratio > 0 ? min(this->alpha,ratio) : this->alpha;	
		ratio = -z[i]/dz[i];
		this->alpha = ratio > 0 ? min(this->alpha,ratio) : this->alpha;	

	}
	if(dtau < 0 ) this->alpha = min(-variables.tau/dtau,this->alpha);
	if(dkappa < 0 ) this->alpha = min(-variables.kappa/dkappa,this->alpha);
	this->alpha = this->alpha*settings.bkscale;
}

////-----------End lp direction
//

// algorithm_state
algorithm_state::algorithm_state() {
	mu = -1;
	sigma = -1;
	gap = -1;
}

void algorithm_state::update_gap(lp_variables &variables,  lp_input &problem_data){
	// ((problem_data.c)'*(variables.x) + (problem_data.h)'*(variables.z) + (problem_data.b)'*(variables.y))
    gap = problem_data.c.dot(variables.x) + problem_data.h.dot(variables.z) + problem_data.b.dot(variables.y);
}

void algorithm_state::update_mu(lp_variables &variables,  lp_input &problem_data){
    mu = (variables.s.dot(variables.z) + variables.tau*variables.kappa)/(problem_data.m+1.0);
}

//Right hand side
linear_system_rhs::linear_system_rhs( lp_input &problem_data) :
	q123(problem_data.k_var + problem_data.n + problem_data.m), q5(problem_data.m),
    m(problem_data.m), n(problem_data.n), k(problem_data.k_var)
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
	// TODO: Test	
    	
	// q1 = -r1
	// q2 = -r2
	// q3 = -r3
    q123.segment(0,n)   = -residuals.r1;
    q123.segment(n,k)   = -residuals.r2;
    q123.segment(n+k,m) = -residuals.r3;
	
	// q4 = -r4
	q4 = -residuals.r4;
	
	// q5 = -s
    q5 = -variables.s;
	
	// q6 = -k
	q6 = -variables.kappa;

}

void linear_system_rhs::compute_combined_rhs(lp_residuals &residuals, 
                                             lp_variables &variables,
                                             lp_direction &affine_dir,
                                             double sigma, 
                                             double mu){
    
    q123.segment(0,n)   = -(1-sigma)*residuals.r1;
    q123.segment(n,k)   = -(1-sigma)*residuals.r2;
    q123.segment(n+k,m) = -(1-sigma)*residuals.r3;
	
	q4 = -(1-sigma)*residuals.r4;
	
	// q5 = -s+sigma mu/z + (ds.*dz)./z
        q5 = (-variables.s.array()+sigma*mu/variables.z.array() - affine_dir.ds.array()*affine_dir.dz.array()/variables.z.array()).matrix();
	
	// q6 = -k+sigma*mu/tau
	q6 = -variables.kappa+sigma*mu/variables.tau - affine_dir.dtau*affine_dir.dkappa/variables.tau;
}

void linear_system_rhs::var_dump() {
	OUTPUT << "DUMP RHS OBJECT" << endl;
	OUTPUT << "q123" << q123 << endl;
	OUTPUT << "q4"   << q4 << endl;
	OUTPUT << "q5"   << q5 << endl;
	OUTPUT << "q6" << q6 << endl;
}

//Utility functions
//#ifdef random
//lp_input* copl_utility::loadFromUF(string UF_group, string UF_name){
//	
//	lp_input* problem_data;
//	// Inserted an script that download file from UF repository
//    ifstream A_File("../example_problems/ex3sta1/ex3sta1.mtx");
//	
//  	if (!A_File) {
//  		cerr << "Error Loading from UF dataset. File not found (A Matrix)" << endl;
//    		return ;
//	}
//	
//	int idx = 0;
//	string line;
//    bool first_line = true;
//	
// 	while (getline(A_File, line)) {
//  		if (line.empty()) continue;
//                if (line.at(0) == '%') continue;
//
//		std::istringstream ss(line);
//		std::string token;
//		string indata[3] = {"", "", ""};
//
//		idx = 0;
//		while(std::getline(ss, token, ' ')) {
//			indata[idx++] = token;
//		}	
//		if (first_line){
//                	first_line = false;
//			int m = stoi(indata[0]);
//			int n = m;
//			int k_var = stoi(indata[1]);
//                        problem_data = new lp_input(m,n,k_var);
//                }else{
//			int row = stoi(indata[0]) - 1;
//			int col = stoi(indata[1]) - 1;
//			double value = stod(indata[2]); 
//			problem_data->A.insert_at(row,col,value);
//		}
//		
//	}
//	A_File.close();
//
//
//	ifstream b_File("../example_problems/ex3sta1/ex3sta1_b.mtx");
//	if (!b_File) {
//			cerr << "Error Loading from UF dataset. File not found (b vector)" << endl;
//			return ;
//	}
//
//	first_line = true;
//	int col_idx = 0;
//	while (getline(b_File, line)) {
//			if (line.empty()) continue;
//			if (line.at(0) == '%') continue;
//
//			std::istringstream ss(line);
//			std::string token;
//
//			string indata[1] = {""};
//			idx = 0;
//			while(std::getline(ss, token, ' ')) {
//					indata[idx++] = token;
//			}
//			if (first_line){
//					first_line = false;
//			}else{
//					double value = std::stod(indata[0]);
//					problem_data->b[col_idx++] = value;
//			}
//
//	}
//	b_File.close();
//
//
//	ifstream c_File("../example_problems/ex3sta1/ex3sta1_c.mtx");
//	if (!c_File) {
//			cerr << "Error Loading from UF dataset. File not found (c vector)" << endl;
//			return ;
//	}
//
//	first_line = true;
//	col_idx = 0;
//	while (getline(c_File, line)) {
//			if (line.empty()) continue;
//			if (line.at(0) == '%') continue;
//
//			std::istringstream ss(line);
//			std::string token;                
//
//			string indata[1] = {""};
//			idx = 0;
//			while(std::getline(ss, token, ' ')) {
//					indata[idx++] = token;
//			}
//			if (first_line){
//					first_line = false;
//			}else{
//					double value = std::stod(indata[0]);
//					problem_data->c[col_idx++] = value;
//			}
//
//	}
//	c_File.close();
//
//	return problem_data;
//}
//#endif  


}
