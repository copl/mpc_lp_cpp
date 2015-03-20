#include "gtest/gtest.h"
#include "copl_linalg.h"
#include "copl_newton.h"
#include <vector>

using namespace copl_ip;
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

TEST(KNEWTON,Construct)
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
        vars.s[j] = 1.0;
        vars.z[j] = 1.0;
    }
    
    k_mat.update(vars);
}

