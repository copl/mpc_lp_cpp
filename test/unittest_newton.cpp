#include "gtest/gtest.h"
#include "copl_linalg.h"
#include "copl_newton.h"
#include <vector>

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
    A.insert_at(0,0,1.0);
    A.insert_at(0,1,2.0);
    A.insert_at(1,0,3.0);
    A.insert_at(1,2,4.0);    
    A.insert_at(1,3,5.0);
    
     
    G.insert_at(0,3,6.0);    
    G.insert_at(1,2,7.0);
    G.insert_at(2,0,8.0);
    G.insert_at(2,1,9.0);
    G.insert_at(2,2,10.0);
    G.insert_at(2,3,11.0);
	
    //Assemble a matrix K
    k_newton_copl_matrix k_mat(4,10); 
    k_mat.m = G.num_rows();
    k_mat.n = G.num_cols();
    k_mat.p = A.num_rows();
    k_mat.assemble_matrix(A,G);
    //The new matrix must have 10 + 2*(nnz(A)+nnz(G)
    ASSERT_EQ(k_mat.nnz(),9+2*A.nnz()+2*G.nnz());
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
    A.insert_at(0,0,1.0);
    A.insert_at(0,1,2.0);
    A.insert_at(1,0,3.0);
    A.insert_at(1,2,4.0);    
    A.insert_at(1,3,5.0);
    
     
    G.insert_at(0,3,6.0);    
    G.insert_at(1,2,7.0);
    G.insert_at(2,0,8.0);
    G.insert_at(2,1,9.0);
    G.insert_at(2,2,10.0);
    G.insert_at(2,3,11.0);
    
    //Assemble a matrix K
    k_newton_copl_matrix k_mat(4,10); 
    k_mat.m = G.num_rows();
    k_mat.n = G.num_cols();
    k_mat.p = A.num_rows();

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
    A.insert_at(0,0,1.0);
    A.insert_at(0,1,2.0);
    A.insert_at(1,0,3.0);
    A.insert_at(1,2,4.0);    
    A.insert_at(1,3,5.0);
    
     
    G.insert_at(0,3,6.0);    
    G.insert_at(1,2,7.0);
    G.insert_at(2,0,8.0);
    G.insert_at(2,1,9.0);
    G.insert_at(2,2,10.0);
    G.insert_at(2,3,11.0);
    
    //Assemble a matrix K
    k_newton_copl_matrix k_mat(4,10); 
    k_mat.m = G.num_rows();
    k_mat.n = G.num_cols();
    k_mat.p = A.num_rows();
    
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
    
    A.insert_at(0,0,1.0);
    A.insert_at(0,1,2.0);
    A.insert_at(1,0,3.0);
    A.insert_at(1,2,4.0);    
    A.insert_at(1,3,5.0);
    
     
    G.insert_at(0,3,6.0);    
    G.insert_at(1,2,7.0);
    G.insert_at(2,0,8.0);
    G.insert_at(2,1,9.0);
    G.insert_at(2,2,10.0);
    G.insert_at(2,3,11.0);
    	
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
    
    A.insert_at(0,0,1.0);
    A.insert_at(0,1,2.0);
    A.insert_at(1,0,3.0);
    A.insert_at(1,2,4.0);    
    A.insert_at(1,3,5.0);
    
     
    G.insert_at(0,3,6.0);    
    G.insert_at(1,2,7.0);
    G.insert_at(2,0,8.0);
    G.insert_at(2,1,9.0);
    G.insert_at(2,2,10.0);
    G.insert_at(2,3,11.0);
    	
    copl_vector c,b,h;
    
    //Assemble a matrix K
    k_newton_copl_matrix k_mat(A,G); 
    int m = 3;
    int n = 4;
    int k = 2;
    
    lp_variables vars(n,m,k);
    
    for(int j = 0; j < m ;j++) {
        vars.s[j] = 3.0;
        vars.z[j] = 2.0;
    }
    
    k_mat.update(vars);
    //After the update the matrix should be 
    //[dI A' G'   ] 
    //[A dI       ]
    //[G    -s/z-d]

    copl_vector copl_x(m+n+k);
    copl_vector copl_y(m+n+k);
    copl_vector copl_w(m);

    //Set x to all ones 
    for(int j = 0; j < m+n+k; j++)
        copl_x[j] = 1.0;

    for(int j = 0; j < m; j++)
        copl_w[j] = vars.s[j]/vars.z[j]*copl_x[j+n+k];

    //Make a map to be able to use y and x in Eigen
    Eigen::Map<Eigen::VectorXd> x(&copl_x[0],copl_x.size());
    Eigen::Map<Eigen::VectorXd> y(&copl_y[0],copl_y.size());        
    Eigen::Map<Eigen::VectorXd> w(&copl_w[0],copl_w.size());    
    EigenSpMat_t K = *(k_mat.eigenKMat);
    y = K*x;

    y.segment(0,n)   = y.segment(0,n) - k_mat.DELTA*x.segment(0,n)
                      - (A.eigenMat->transpose())*x.segment(n,k) 
                      - G.eigenMat->transpose()*x.segment(n+k,m);
    y.segment(n,k)   = y.segment(n,k) - *(A.eigenMat)*x.segment(0,n) \
                       + k_mat.DELTA*x.segment(n,k);
    y.segment(n+k,m) = y.segment(n+k,m) - *(G.eigenMat)*x.segment(0,n)\
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
    
    A.insert_at(0,0,1.0);
    A.insert_at(0,1,2.0);
    A.insert_at(1,0,3.0);
    A.insert_at(1,2,4.0);    
    A.insert_at(1,3,5.0);
    
     
    G.insert_at(0,3,6.0);    
    G.insert_at(1,2,7.0);
    G.insert_at(2,0,8.0);
    G.insert_at(2,1,9.0);
    G.insert_at(2,2,10.0);
    G.insert_at(2,3,11.0);
        
    copl_vector c,b,h;
    
    //Assemble a matrix K
    k_newton_copl_matrix k_mat(A,G); 
    int m = 3;
    int n = 4;
    int k = 2;
    
    lp_variables vars(n,m,k);
    
    for(int j = 0; j < m ;j++) {
        vars.s[j] = 3.0;
        vars.z[j] = 2.0;
    }
    
    k_mat.update(vars);
    //After the update the matrix should be 
    //[dI A' G'   ] 
    //[A dI       ]
    //[G    -s/z-d]

    copl_vector copl_x(m+n+k);
    copl_vector copl_y(m+n+k);
    copl_vector copl_w(m);

    //Set x to all ones 
    for(int j = 0; j < m+n+k; j++)
        copl_x[j] = 1.0;

    for(int j = 0; j < m; j++)
        copl_w[j] = vars.s[j]/vars.z[j]*copl_x[j+n+k];

    //Make a map to be able to use y and x in Eigen
    Eigen::Map<Eigen::VectorXd> x(&copl_x[0],copl_x.size());
    Eigen::Map<Eigen::VectorXd> y(&copl_y[0],copl_y.size());        
    Eigen::Map<Eigen::VectorXd> w(&copl_w[0],copl_w.size());    
    EigenSpMat_t K = *(k_mat.eigenKMat);
    
    k_mat.solve(copl_y,copl_x);
    //y <- K\x
    //r = K*y-x
    y = K*y-x;
                   
    cout << "Residual norm of solution: " << y.norm() << "\n";
    ASSERT_LT(y.norm(),1.e-14);

}
}