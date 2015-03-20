#include "gtest/gtest.h"
#include <copl_linalg.h>
#include <copl_core.h>

using namespace copl_ip;
TEST(linalg_copl_vector,access)
{
    copl_vector vec(10);
    EXPECT_EQ(vec.size(),10) << "Vector reports incorrect length";
    copl_vector vec2(0);
    EXPECT_EQ(vec2.size(),0) << "Vector reports incorrect length";
    vec[1] = 100;  
    EXPECT_EQ(vec[1],100) << "Vector[1] should be = 100";
    
}

//Test the copl matrix 
TEST(linalg_copl_matrix,create)
{
    copl_matrix A(10,100);
    EXPECT_EQ(A.num_cols(),100);
    EXPECT_EQ(A.num_rows(),10);
    EXPECT_EQ(A.nnz(),0);
    A.insert_at(5,50,100.0);
    EXPECT_EQ(A.value_at(5,50),100.0);
    EXPECT_EQ(A.value_at(5,51),0.0);
}

TEST(linalg_copl_marix,dgemv)
{
  //y<- alpha Ax + beta y with sparse A
  copl_vector y(10); 

  copl_vector x(10);
  copl_matrix A(10,10);
  //All but the last diag
  for(int j = 0; j < 9; j++)
  {
    A.insert_at(j,j,0.3);
    x[j] = 1.0;
  }
  x[10] = 1.0; 
    
  //sp_dgemv(2.0,0.0,A,x,y);
  //All entries of y should be 2*0.3 except for the last which should be 0

}

////y<- alpha Ax + beta y with sparse A
//void sp_dgemv(double alpha, double beta, copl_matrix &copl_A, copl_vector &copl_x, copl_vector &copl_y);
//
////Matrix vector multiply and accumulate in y
////y<- alpha A^Tx + beta y with sparse A (The input matrix is A not A^T)
//void sp_dgemtv(double alpha, double beta, copl_matrix &copl_A, copl_vector &copl_x, copl_vector &copl_y);
//
////Scale the vector
////x<-alpha *x
//void scal(copl_vector &copl_x, double alpha );
//
////y<- alpha*x + y
//void axpy(double alpha, copl_vector &copl_x, copl_vector &copl_y);
//
//// x^Ty
//double dot(copl_vector &copl_y, copl_vector &copl_x);
//
////Zero out 
//void zeros(copl_vector &y);
//
////Two norm 
//double norm2(copl_vector &y);
//
////Infinity norm 	
//double normInf(copl_vector &y);
//
//
