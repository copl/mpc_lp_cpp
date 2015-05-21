#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <vector>
#include <copl_linalg.h>
#include <copl_core.h>
#include <math.h>

using namespace copl_ip;
using namespace testing;
TEST(Core,input)
{
//TODO:
//
}

//Test the vars constructor
TEST(CORE,lp_vars_constructor)
{
    int m = 10;
    int n = 3;
    int k = 5;
    lp_variables vars(m,n,k);
    //Copy to stl to use gmock 
    std::vector<double> x(n),y(k),z(m),s(m);
    
    std::copy(&vars.x[0],&vars.x[0]+n,x.begin());
    std::copy(&vars.y[0],&vars.y[0]+k,y.begin());
    std::copy(&vars.s[0],&vars.s[0]+m,s.begin());
    std::copy(&vars.z[0],&vars.z[0]+m,z.begin());
    
    //Check that the entries are correct
    ASSERT_THAT(x,Each(0.0));
    ASSERT_THAT(y,Each(0.0));
    ASSERT_THAT(s,Each(1.0));
    ASSERT_THAT(z,Each(1.0));
}

//Test the direction constructor
TEST(CORE,dir_constructor)
{
    int m = 10;
    int n = 3;
    int k = 5;
    lp_variables vars(m,n,k);
    lp_direction dir(vars);
    //Check that the initialized vectors are the same length
    ASSERT_EQ(dir.dx.size(),n);
    ASSERT_EQ(dir.dy.size(),k);
    ASSERT_EQ(dir.ds.size(),m);
    ASSERT_EQ(dir.dz.size(),m);
}

//Test the code that updates the step 
TEST(CORE,test_step)
{
    int m = 10;
    int n = 3;
    int k = 5;
    lp_variables vars(m,n,k);
    lp_direction dir(vars); 
    
    vars.x.setConstant(1.0);
    vars.y.setConstant(2.0);
    vars.s.setConstant(3.0);
    vars.z.setConstant(4.0);
    vars.tau = 5.0;
    vars.kappa = 6.0;
    //Set the direction 
    dir.dx.setConstant(1.0);
    dir.dy.setConstant(2.0);
    dir.ds.setConstant(3.0);
    dir.dz.setConstant(4.0);
    dir.dtau = 5.0;
    dir.dkappa = 6.0; 
    dir.alpha = 0.5;
    //Take a step along the direction  
    vars.take_step(dir);

    //Copy to stl to use gmock 
    std::vector<double> x(n),y(k),z(m),s(m);
    std::copy(&vars.x[0],&vars.x[0]+n,x.begin());
    std::copy(&vars.y[0],&vars.y[0]+k,y.begin());
    std::copy(&vars.s[0],&vars.s[0]+m,s.begin());
    std::copy(&vars.z[0],&vars.z[0]+m,z.begin());
    
    //Check that the entries are correct
    ASSERT_THAT(x,Each(1.5));
    ASSERT_THAT(y,Each(3.0));
    ASSERT_THAT(s,Each(4.5));
    ASSERT_THAT(z,Each(6.0));
    ASSERT_THAT(vars.tau,7.5);
    ASSERT_THAT(vars.kappa,9.0);
}

//Test the code that first calculatest the step size and then updates the step 
TEST(CORE,test_stepsize_andstep)
{
    int m = 10;
    int n = 3;
    int k = 5;
    lp_variables vars(m,n,k);
    lp_direction dir(vars); 
    
    vars.x.setConstant(1.0);
    vars.y.setConstant(2.0);
    vars.s.setConstant(3.0);
    vars.z.setConstant(4.0);
    vars.tau = 5.0;
    vars.kappa = 6.0;
    //Set the direction 
    dir.dx.setConstant(-1.0);
    dir.dy.setConstant(-2.0);
    dir.ds.setConstant(-3.0);
    dir.dz.setConstant(-4.0);
    dir.dtau = -5.0;
    dir.dkappa = -6.0; 
     
    lp_settings settings(10,1.e-10,1.e-10,1.0,1.e-7);
    dir.compute_step_size(vars,settings);
    
    //Take a step along the direction  
    vars.take_step(dir);

    //Copy to stl to use gmock 
    std::vector<double> x(n),y(k),z(m),s(m);
    std::copy(&vars.x[0],&vars.x[0]+n,x.begin());
    std::copy(&vars.y[0],&vars.y[0]+k,y.begin());
    std::copy(&vars.s[0],&vars.s[0]+m,s.begin());
    std::copy(&vars.z[0],&vars.z[0]+m,z.begin());
    
    //Check that the entries are correct
    ASSERT_THAT(x,Each(0.0));
    ASSERT_THAT(y,Each(0.0));
    ASSERT_THAT(s,Each(0.0));
    ASSERT_THAT(z,Each(0.0));
    ASSERT_THAT(vars.tau,0.0);
    ASSERT_THAT(vars.kappa,0.0);
}


TEST(CORE,residuals_constructor)
{
  
    int m = 3;
    int n = 4;
    int k = 2;
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
    lp_residuals res(lp_problem);

    ASSERT_EQ(res.r1.size(),n);
    ASSERT_EQ(res.r2.size(),k);
    ASSERT_EQ(res.r3.size(),m);
}

TEST(CORE,residuals)
{
    int m = 3;
    int n = 4;
    int k = 2;
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
    lp_variables vars(m,n,k);
    //Set the variables to some state
    vars.x << 0.5,0.6,0.7,0.8;
    vars.y << 0.9,1.1;
    vars.s << 0.5,0.6,0.9;
    vars.z << 0.5,1.e-3,1.7;
    vars.tau = 0.4;
    vars.kappa = 1.e-4;
    lp_residuals res(lp_problem);
  
    //Call the residual calculation method 
    res.compute_residuals(lp_problem, vars);
    //the residuals should be 
    copl_vector error1(n), error2(k), error3(m); 
    double error4;

    error1 = res.r1 + (A.transpose()*vars.y + G.transpose()*vars.z + c*vars.tau);
    error2 = res.r2 + (-A*vars.x + vars.tau*b);
    error3 = res.r3 + (-G*vars.x - vars.s + vars.tau*h);
    error4 = res.r4 + (- c.dot(vars.x) - b.dot(vars.y) - h.dot(vars.z) - vars.kappa); 
    double err_norm = error1.norm()+error2.norm()+error3.norm()+fabs(error4);
    ASSERT_LT(err_norm, 1.e-15);
  }

 TEST(CORE,step_size)
 {
    int m = 3;
    int n = 4;
    int k = 2;
    lp_variables vars(m,n,k);
    //Set the variables to some state
    vars.x << 0.5,0.6,0.7,0.8;
    vars.y << 0.9,1.1;
    vars.s << 0.5,0.6,0.9;
    vars.z << 0.5,1.e-3,1.7;
    vars.tau = 0.4;
    vars.kappa = 1.e-4;
    lp_direction dir(vars); 
    dir.ds << -0.1,-0.1,-0.1;
    dir.dz << -0.1,-0.1,-0.1;
    dir.dtau = 0.1;
    dir.dkappa = 0.1;
    lp_settings settings(10,1.e-10,1.e-10,1.0,1.e-7);
    dir.compute_step_size(vars,settings);
    
    //The step size should be 10 times the smallest of s and z
    ASSERT_EQ(1.e-2,dir.alpha);
 }

 TEST(CORE,residual_norms)
 {
 //TODO: Test that the norms are correctly calculated, in compute_residuals
 cout << "Missing test\n";
 }
 TEST(CORE,gap)
 {
 //TODO: Build test for gap 
 cout << "Missing test\n";
 }

 TEST(CORE,mu)
 {
 //TODO test mu update   
 cout << "Missing test\n";
 }
