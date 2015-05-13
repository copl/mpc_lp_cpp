#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "copl_linalg.h"
#include "copl_newton.h"
#include <vector>

//We need this for the gmock predicates like Each
using namespace testing;

//The tests have to be in the same namespace as the classes they are friends to
namespace copl_ip {
//Test member of EigenTest named Assemble
TEST(KNEWTON,Assemble)
{
    //Make a copl_matrix A and G
    //With patter 
    //    x x 0 0
    //A = x 0 x x
    //
    //G = 0 0 0 x
    //    0 0 x 0
    //    x x x x
    // 

    copl_matrix A(2,4);
    copl_matrix G(3,4);
    A.insert(0,0) = 1.0;
    A.insert(0,1) = 2.0;
    A.insert(1,0) = 3.0;
    A.insert(1,2) = 4.0;    
    A.insert(1,3) = 5.0;
    
     
    G.insert(0,3) = 6.0;    
    G.insert(1,2) = 7.0;
    G.insert(2,0) = 8.0;
    G.insert(2,1) = 9.0;
    G.insert(2,2) = 10.0;
    G.insert(2,3) = 11.0;
	
    //Assemble a matrix K
    k_newton_copl_matrix k_mat(4,9); 
    k_mat.m = G.rows();
    k_mat.n = G.cols();
    k_mat.k = A.rows();
    k_mat.assemble_matrix(A,G);
    //The new matrix must have 10 + 2*(nnz(A)+nnz(G)
    ASSERT_EQ(k_mat.nnz(),9+2*A.nonZeros()+2*G.nonZeros());
    //Test the nonzero pattern the columns should have 
    //4344 34225
}

TEST(KNEWTON,NonZeroPerCols)
{

    //Make a copl_matrix A and G
    //With patter 
    //    x x 0 0
    //A = x 0 x x
    //
    //G = 0 0 0 x
    //    0 0 x 0
    //    x x x x
    // 

    copl_matrix A(2,4);
    copl_matrix G(3,4);
    A.insert(0,0) = 1.0;
    A.insert(0,1) = 2.0;
    A.insert(1,0) = 3.0;
    A.insert(1,2) = 4.0;    
    A.insert(1,3) = 5.0;
    
    G.insert(0,3) = 6.0;    
    G.insert(1,2) = 7.0;
    G.insert(2,0) = 8.0;
    G.insert(2,1) = 9.0;
    G.insert(2,2) = 10.0;
    G.insert(2,3) = 11.0;
	
    //Assemble a matrix K
    k_newton_copl_matrix k_mat(4,9); 
    k_mat.m = G.rows();
    k_mat.n = G.cols();
    k_mat.k = A.rows();

    k_mat.assemble_matrix(A,G);

    //Test the nonzero pattern the columns should have 
    //4344 34 225
    k_mat.eigenKMat->makeCompressed();
    int *ii = k_mat.eigenKMat->outerIndexPtr();
    EXPECT_EQ(ii[1]-ii[0],4);
    EXPECT_EQ(ii[2]-ii[1],3);
    EXPECT_EQ(ii[3]-ii[2],4);
    EXPECT_EQ(ii[4]-ii[3],4);
    EXPECT_EQ(ii[5]-ii[4],3);
    EXPECT_EQ(ii[6]-ii[5],4);
    EXPECT_EQ(ii[7]-ii[6],2);
    EXPECT_EQ(ii[8]-ii[7],2);
    EXPECT_EQ(ii[9]-ii[8],5);
}

TEST(KNEWTON,nnz){

    //Make a copl_matrix A and G
    //With patter 
    //    x x 0 0
    //A = x 0 x x
    //
    //G = 0 0 0 x
    //    0 0 x 0
    //    x x x x
    // 

    copl_matrix A(2,4);
    copl_matrix G(3,4);
    A.insert(0,0) = 1.0;
    A.insert(0,1) = 2.0;
    A.insert(1,0) = 3.0;
    A.insert(1,2) = 4.0;    
    A.insert(1,3) = 5.0;
    
    G.insert(0,3) = 6.0;    
    G.insert(1,2) = 7.0;
    G.insert(2,0) = 8.0;
    G.insert(2,1) = 9.0;
    G.insert(2,2) = 10.0;
    G.insert(2,3) = 11.0;
	 
    //Assemble a matrix K
    k_newton_copl_matrix k_mat(4,9); 
    k_mat.m = G.rows();
    k_mat.n = G.cols();
    k_mat.k = A.rows();
    
    k_mat.assemble_matrix(A,G);
    ASSERT_EQ(k_mat.nnz(),22+2+3+4);
}

//Test the public constructor
TEST(KNEWTON,Constructor)
{
    //Make a copl_matrix A and G
    //With patter 
    //    x x 0 0
    //A = x 0 x x
    //
    //G = 0 0 0 x
    //    0 0 x 0
    //    x x x x
    // 

    copl_matrix A(2,4);
    copl_matrix G(3,4);
    A.insert(0,0) = 1.0;
    A.insert(0,1) = 2.0;
    A.insert(1,0) = 3.0;
    A.insert(1,2) = 4.0;    
    A.insert(1,3) = 5.0;
    
    G.insert(0,3) = 6.0;    
    G.insert(1,2) = 7.0;
    G.insert(2,0) = 8.0;
    G.insert(2,1) = 9.0;
    G.insert(2,2) = 10.0;
    G.insert(2,3) = 11.0;
	
	
    copl_vector c,b,h;
    
    //Assemble a matrix K
    k_newton_copl_matrix k_mat(A,G); 
   
    //Just make sure it does not crash when constructing
    EXPECT_EQ( (*k_mat.hessianIx)[0], 23);
    EXPECT_EQ( (*k_mat.hessianIx)[1], 25);
    ASSERT_EQ( (*k_mat.hessianIx)[2], 30);
}

//Test the update routine
TEST(KNEWTON,Update)
{
    //Make a copl_matrix A and G
    //With patter 
    //    x x 0 0
    //A = x 0 x x
    //
    //G = 0 0 0 x
    //    0 0 x 0
    //    x x x x
    // 
    copl_matrix A(2,4);
    copl_matrix G(3,4);
    A.insert(0,0) = 1.0;
    A.insert(0,1) = 2.0;
    A.insert(1,0) = 3.0;
    A.insert(1,2) = 4.0;    
    A.insert(1,3) = 5.0;
    
    G.insert(0,3) = 6.0;    
    G.insert(1,2) = 7.0;
    G.insert(2,0) = 8.0;
    G.insert(2,1) = 9.0;
    G.insert(2,2) = 10.0;
    G.insert(2,3) = 11.0;
	 	
    copl_vector c,b,h;
    
    //Assemble a matrix K
    k_newton_copl_matrix k_mat(A,G); 
    int m = 3;
    int n = 4;
    int k = 2;
    
    lp_variables vars(m,n,k);
    
    for(int j = 0; j < m ;j++) {
        vars.s[j] = 3.0;
        vars.z[j] = 2.0;
    }
    
    k_mat.update(vars);
    //After the update the matrix should be 
    //[dI A' G'   ] 
    //[A dI       ]
    //[G    -s/z-d]

    copl_vector x(m+n+k);
    copl_vector y(m+n+k);
    copl_vector w(m);

    //Set x to all ones 
    for(int j = 0; j < m+n+k; j++)
        x[j] = 1.0;

    for(int j = 0; j < m; j++)
        w[j] = vars.s[j]/vars.z[j]*x[j+n+k];

    copl_matrix K = *(k_mat.eigenKMat);
    y = K*x;

    y.segment(0,n)   = y.segment(0,n) - k_mat.DELTA*x.segment(0,n)
                      - A.transpose()*x.segment(n,k) 
                      - G.transpose()*x.segment(n+k,m);
    y.segment(n,k)   = y.segment(n,k) - A*x.segment(0,n) \
                       + k_mat.DELTA*x.segment(n,k);
    y.segment(n+k,m) = y.segment(n+k,m) - G*x.segment(0,n)\
                       + k_mat.DELTA*x.segment(n+k,m) + w; 
    
    cout << "Residual norm: " << y.norm() << "\n";
    ASSERT_LT(y.norm(),1.e-14);

}

TEST(KNEWTON,solve)
{
    //Make a copl_matrix A and G
    //With patter 
    //    x x 0 0
    //A = x 0 x x
    //
    //G = 0 0 0 x
    //    0 0 x 0
    //    x x x x
    // 
    copl_matrix A(2,4);
    copl_matrix G(3,4);
    A.insert(0,0) = 1.0;
    A.insert(0,1) = 2.0;
    A.insert(1,0) = 3.0;
    A.insert(1,2) = 4.0;    
    A.insert(1,3) = 5.0;
    
    G.insert(0,3) = 6.0;    
    G.insert(1,2) = 7.0;
    G.insert(2,0) = 8.0;
    G.insert(2,1) = 9.0;
    G.insert(2,2) = 10.0;
    G.insert(2,3) = 11.0;
	
    copl_vector c,b,h;
    
    //Assemble a matrix K
    k_newton_copl_matrix k_mat(A,G); 
    int m = 3;
    int n = 4;
    int k = 2;
    
    lp_variables vars(m,n,k);
    
    for(int j = 0; j < m ;j++) {
        vars.s[j] = 3.0;
        vars.z[j] = 2.0;
    }
    
    k_mat.update(vars);
    //After the update the matrix should be 
    //[dI A' G'   ] 
    //[A dI       ]
    //[G    -s/z-d]

    copl_vector x(m+n+k);
    copl_vector y(m+n+k);
    copl_vector w(m);

    //Set x to all ones 
    for(int j = 0; j < m+n+k; j++)
        x[j] = 1.0;

    for(int j = 0; j < m; j++)
        w[j] = vars.s[j]/vars.z[j]*x[j+n+k];

    copl_matrix K = *(k_mat.eigenKMat);
    
    k_mat.solve(y,x);
    //y <- K\x
    //r = K*y-x
    y = K*y-x;
                   
    cout << "Residual norm of solution: " << y.norm() << "\n";
    ASSERT_LT(y.norm(),1.e-14);

}

TEST(KNEWTON,solveException)
{
    //Make a copl_matrix A and G
    //With patter 
    //    x x 0 0
    //A = x 0 x x
    //
    //G = 0 0 0 x
    //    0 0 x 0
    //    x x x x
    // 
    copl_matrix A(2,4);
    copl_matrix G(3,4);
    A.insert(0,0) = 1.0;
    A.insert(0,1) = 2.0;
    A.insert(1,0) = 3.0;
    A.insert(1,2) = 4.0;    
    A.insert(1,3) = 5.0;
    
    G.insert(0,3) = 6.0;    
    G.insert(1,2) = 7.0;
    G.insert(2,0) = 8.0;
    G.insert(2,1) = 9.0;
    G.insert(2,2) = 10.0;
    G.insert(2,3) = 11.0;
	
    copl_vector c,b,h;
    
    //Assemble a matrix K
    k_newton_copl_matrix k_mat(A,G); 
    copl_vector x(9);
    copl_vector y(9);
    ASSERT_ANY_THROW(k_mat.solve(y,x));
}

TEST(KNEWTON,reduce_rhs_test)
{
//This method should form the rhs for 
//[    A   G    c]dx    r1
//[A           -b]dy    r2
//[G      -H   -h]dz  = r3 - r5
//[-c' b' h' -k/t]dt    r4 - r6


    copl_matrix A(2,4);
    copl_matrix G(3,4);
    copl_vector c(4),b(2),h(3);
    c  << 1,2,3,4;
    b  << 5,6;
    h  << 7,8,9;

    A.insert(0,0) = 1.0;
    A.insert(0,1) = 2.0;
    A.insert(1,0) = 3.0;
    A.insert(1,2) = 4.0;    
    A.insert(1,3) = 5.0;
         
    G.insert(0,3) = 6.0;    
    G.insert(1,2) = 7.0;
    G.insert(2,0) = 8.0;
    G.insert(2,1) = 9.0;
    G.insert(2,2) = 10.0;
    G.insert(2,3) = 11.0;
            
    lp_input lp_problem(A,b,c,G,h);
    //Form the rhs 
    linear_system_rhs rhs(lp_problem);
    //Fill the rhs with stuff
    rhs.q123.setConstant(2.0);
    rhs.q4 = 2.0;
    rhs.q5.setConstant(10.0);
    rhs.q6 = 10.0;

    lp_variables vars(3,4,2);
    //Set the variables to some state
    vars.x << 0.5,0.6,0.7,0.8;
    vars.y << 0.9,1.1;
    vars.s << 0.5,0.6,0.9;
    vars.z << 0.5,1.e-3,1.7;
    vars.tau = 0.4;
    vars.kappa = 1.e-4;
 
    homogeneous_solver K_solver(lp_problem);
    K_solver.update(vars);
    K_solver.reduce_rhs(rhs);
    //check the last entires of q123, it should now read  -8
    
    //Copy to stl to use gmock 
    std::vector<double> qs(6);
    std::vector<double> qe(3);
    std::vector<double> q5(3);
    std::copy(&rhs.q123[0],&rhs.q123[0]+6,qs.begin());
    std::copy(&rhs.q123[6],&rhs.q123[6]+3,qe.begin());
    std::copy(&rhs.q5[0],&rhs.q5[0]+3,q5.begin());
     
    //Check that the entries are correct
    ASSERT_THAT(qs,Each(2.0));
    ASSERT_THAT(qe,Each(-8.0));
    ASSERT_THAT(rhs.q4,-8.0);
    ASSERT_THAT(q5,Each(10.0));
}

TEST(KNEWTON,back_substitute_test)
{

    copl_matrix A(2,4);
    copl_matrix G(3,4);
    copl_vector c(4),b(2),h(3);
    c  << 1,2,3,4;
    b  << 5,6;
    h  << 7,8,9;

    A.insert(0,0) = 1.0;
    A.insert(0,1) = 2.0;
    A.insert(1,0) = 3.0;
    A.insert(1,2) = 4.0;    
    A.insert(1,3) = 5.0;
         
    G.insert(0,3) = 6.0;    
    G.insert(1,2) = 7.0;
    G.insert(2,0) = 8.0;
    G.insert(2,1) = 9.0;
    G.insert(2,2) = 10.0;
    G.insert(2,3) = 11.0;

    lp_variables vars(3,4,2);
    vars.s << 1.0,1.0,1.0;
    vars.z << 0.5,0.5,0.5;
    vars.tau = 0.5;
    vars.kappa = 1.0;
    
    //Input
    lp_input lp_problem(A,b,c,G,h);

    //Form the rhs 
    linear_system_rhs rhs(lp_problem); 
    lp_direction dir(vars);
    
    //Set the directions
    dir.dz << 1.0,2.0,3.0; 
    dir.dtau = 0.5;

    //Set the rhs
    rhs.q5 << 0.5,1.5,2.5;
    rhs.q6 = 1.5;

    homogeneous_solver K_solver(lp_problem);
    K_solver.update(vars);
    K_solver.back_substitute(dir, rhs, vars);
   
    copl_vector error(3);
    //Hdz + ds = rhs
    error = vars.s.array()/vars.z.array()*dir.dz.array() - rhs.q5.array() + dir.ds.array();
    ASSERT_LT(error.norm(),1.e-15);
}

//Validate that after update sol1 has the solution to 
//[0   A'   G' ]dz      c
//[A           ]dy    = -b
//[G       -H  ]dz      -h

TEST(HOMOGENEOUS_SOLVER,solve_fixed_rhs)
{
    //Populate A,G,c,b,h
    copl_matrix A(2,4);
    copl_matrix G(3,4);
    copl_vector c(4),b(2),h(3);
    c  << 1,2,3,4;
    b  << 5,6;
    h  << 7,8,9;
    A.insert(0,0) = 1.0;
    A.insert(0,1) = 2.0;
    A.insert(1,0) = 3.0;
    A.insert(1,2) = 4.0;    
    A.insert(1,3) = 5.0; 
    G.insert(0,3) = 6.0;    
    G.insert(1,2) = 7.0;
    G.insert(2,0) = 8.0;
    G.insert(2,1) = 9.0;
    G.insert(2,2) = 10.0;
    G.insert(2,3) = 11.0;
    
    lp_input lp_problem(A,b,c,G,h);
    linear_system_rhs rhs(lp_problem);

    //Fill the rhs with stuff
    rhs.q123.setConstant(2.0);
    rhs.q4 = 2.0;
    rhs.q5.setConstant(10.0);
    rhs.q6 = 10.0;

    lp_variables vars(3,4,2);
    //Set the variables to some state
    vars.x << 0.5,0.6,0.7,0.8;
    vars.y << 0.9,1.1;
    vars.s << 0.5,0.6,0.9;
    vars.z << 0.5,1e-3,1.7;
    vars.tau = 0.4;
    vars.kappa = 1.e-4;
 
    //Allocate space for the solution 
    lp_direction dir(vars);
 
    homogeneous_solver K_solver(lp_problem);
    K_solver.update(vars);//Solves and populates sol_1
    //Vector to compute the error 
    copl_vector errorx(4);
    copl_vector errory(2);
    copl_vector errorz(3);
    errorx = A.transpose()*K_solver.sol_1.segment(4,2)  + 
	     K_solver.DELTA*K_solver.sol_1.segment(0,4) +
             G.transpose()*K_solver.sol_1.segment(6,3)  + c;
    errory = A*K_solver.sol_1.segment(0,4) -
             K_solver.DELTA*K_solver.sol_1.segment(4,2) - b;
    errorz = G*K_solver.sol_1.segment(0,4) - 
             K_solver.DELTA*K_solver.sol_1.segment(6,3) -
             (vars.s.array()/vars.z.array()*K_solver.sol_1.segment(6,3).array()).matrix() - h;
    EXPECT_LT(errorx.norm(),1e-14);  
    EXPECT_LT(errory.norm(),1e-14); 
    EXPECT_LT(errorz.norm(),1e-14); 

}


TEST(HOMOGENEOUS_SOLVER,solve_reduced)
{
    //Populate A,G,c,b,h
    copl_matrix A(2,4);
    copl_matrix G(3,4);
    copl_vector c(4),b(2),h(3);
    c  << 1,2,3,4;
    b  << 5,6;
    h  << 7,8,9;
    A.insert(0,0) = 1.0;
    A.insert(0,1) = 2.0;
    A.insert(1,0) = 3.0;
    A.insert(1,2) = 4.0;    
    A.insert(1,3) = 5.0; 
    G.insert(0,3) = 6.0;    
    G.insert(1,2) = 7.0;
    G.insert(2,0) = 8.0;
    G.insert(2,1) = 9.0;
    G.insert(2,2) = 10.0;
    G.insert(2,3) = 11.0;
    
    lp_input lp_problem(A,b,c,G,h);
    linear_system_rhs rhs(lp_problem);

    //Fill the rhs with stuff
    rhs.q123.setConstant(2.0);
    rhs.q4 = 2.0;
    rhs.q5.setConstant(10.0);
    rhs.q6 = 10.0;

    lp_variables vars(3,4,2);
    //Set the variables to some state
    vars.x << 0.5,0.6,0.7,0.8;
    vars.y << 0.9,1.1;
    vars.s << 0.5,0.6,0.9;
    vars.z << 0.5,1e-3,1.7;
    vars.tau = 0.4;
    vars.kappa = 1.e-4;
 
    //Allocate space for the solution 
    lp_direction dir(vars);
 
    homogeneous_solver K_solver(lp_problem);
    K_solver.update(vars);//Solves and populates sol_1
    K_solver.solve_reduced(dir,rhs);

    //Vector to compute the error 
    copl_vector errorx(4);
    copl_vector errory(2);
    copl_vector errorz(3);
    double errort;

    errorx = A.transpose()*dir.dy+G.transpose()*dir.dz+dir.dtau*c +
	     K_solver.DELTA*dir.dx - rhs.q123.segment(0,4);
    errory = A*dir.dx-dir.dtau*b   -
	     K_solver.DELTA*dir.dy -
	     rhs.q123.segment(4,2);
    errorz = G*dir.dx - dir.dtau*h - 
	     K_solver.DELTA*dir.dz -
             rhs.q123.segment(6,3)-(vars.s.array()/vars.z.array()*dir.dz.array()).matrix();
    errort = -c.dot(dir.dx) - b.dot(dir.dy) - h.dot(dir.dz) + vars.kappa/vars.tau*dir.dtau - rhs.q4;
    EXPECT_LT(errorx.norm(),1e-14);  
    EXPECT_LT(errory.norm(),1e-14); 
    EXPECT_LT(errorz.norm(),1e-14); 
    ASSERT_LT(errort,1e-15); 

}

}
