#include "gtest/gtest.h"
#include "copl_linalg.h"
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>

//Test member of EigenTest named Assemble
TEST(EigenTests,Assemble)
{
	//Build a sprase Eigien matrix A and G 
	//Form IA
 	//     A0
	 Eigen::SparseMatrix<double> A(10,100), B(100,100), C(110,110), Z(100,100);

     //This vector will count the nnz in the rows of A
     std::vector<int> nnzRowsA(A.rows());
     std::vector<int> reserve;

     for(int j=0;j<A.cols();j++)
     {
        for(SparseMatrix<double>::InnerIterator it(A,j); it; it++)
        {
           nnzRows[it.row()]++;
        }
     }
}

