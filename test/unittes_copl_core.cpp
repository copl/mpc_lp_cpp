#include "gtest/gtest.h"
#include "copl_linalg.h"
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>

//Test member of EigenTest named Assemble
TEST(KNEWTON,Assemble)
{
    //Make a copl_matrix A and G
    //Assemble a matrix K
    lp_input(A,G,c,b,h);
    k_newton_copl_matrix(lp_input);
}

