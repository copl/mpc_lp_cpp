/* 
 * Defines variables vars_x, vars_y .... vars_tau, vars_kappa. And 
 * maps them into Eigen vectors vars_cx,... etc
 * defines the vars lp_variables object 
 * 
 */


#define VARS(m,n,k) \
	std::vector<double> vars_x(n);\
	std::vector<double> vars_y(k);\
	std::vector<double> vars_s(m);\
	std::vector<double> vars_z(m);\
		    double vars_tau;\
		    double vars_kappa;\
	copl_external_vector vars_cx(&vars_x[0],n), vars_cy(&vars_y[0],k), vars_cs(&vars_s[0],m), vars_cz(&vars_z[0],m);\
	lp_variables vars(vars_cx,vars_cy,vars_cs,vars_cz,vars_tau,vars_kappa);
/* Defines the variables cA,cG from Eigen Sparse Matrices
 *
 *
 */
#define PROB(A,b,c,G,h) \
	copl_external_vector cc(&c[0],c.size()); \
	copl_external_vector cb(&b[0],b.size()); \
	copl_external_vector ch(&h[0],h.size()); \
	A.makeCompressed();\
	G.makeCompressed();\
	copl_matrix cA(A.rows(),A.cols(),A.nonZeros(),A.outerIndexPtr(),A.innerIndexPtr(),A.valuePtr());\
        copl_matrix cG(G.rows(),G.cols(),G.nonZeros(),G.outerIndexPtr(),G.innerIndexPtr(),G.valuePtr());\
	lp_input lp_problem(cA,cb,cc,cG,ch);


#define MATS(A,G) \
	A.makeCompressed();\
	G.makeCompressed();\
	copl_matrix cA(A.rows(),A.cols(),A.nonZeros(),A.outerIndexPtr(),A.innerIndexPtr(),A.valuePtr());\
        copl_matrix cG(G.rows(),G.cols(),G.nonZeros(),G.outerIndexPtr(),G.innerIndexPtr(),G.valuePtr());\	
