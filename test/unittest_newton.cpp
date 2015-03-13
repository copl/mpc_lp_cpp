#include "gtest/gtest.h"
#include "copl_linalg.h"
#include <vector>

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
    
}


