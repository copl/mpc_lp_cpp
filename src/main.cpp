/*
 * main.cpp
 *
 *  Created on: 28/01/2015
 *      Author: Oliver
 */
#include <iostream>
#include <copl_algorithm.h>
#include <copl_linalg.h>

using namespace std;
using namespace copl_ip;


void make_trivial_problem(copl_matrix &A, copl_matrix &G, copl_vector &c, copl_vector &b, copl_vector &h)
{

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
 
}
 
void getDimensionUF(string problem_name, int* k_var, int* n)
{
    ifstream A_File("../example_problems/lp_afiro/lp_afiro.mtx");
  if (!A_File) {
                cerr << "Error Loading from UF dataset. File not found (A Matrix)" << endl;
                return;
        }

        int idx = 0;
        string line;
        bool first_line = true;

        while (getline(A_File, line)) {
                if (line.empty()) continue;
                if (line.at(0) == '%') continue;

                std::istringstream ss(line);
                std::string token;
                string indata[3] = {"", "", ""};

                idx = 0;
                while(std::getline(ss, token, ' ')) {
                        indata[idx++] = token;
                }
                if (first_line){
                        first_line = false;
                        *k_var = stoi(indata[0]);
                        *n = stoi(indata[1]);

                }

        }
        A_File.close();

} 
void loadFromUF(string UF_group, string UF_name, copl_matrix &A, copl_matrix &G, copl_vector &c, copl_vector &b, copl_vector &h){
		
	int n = 0;
	int m = 0;
	int k_var = 0;
	// Inserted an script that download file from UF repository
        ifstream A_File("../example_problems/lp_afiro/lp_afiro.mtx");
	
  	if (!A_File) {
  		cerr << "Error Loading from UF dataset. File not found (A Matrix)" << endl;
    		return;
	}
	
	int idx = 0;
	string line;
    bool first_line = true;
	
 	while (getline(A_File, line)) {
  		if (line.empty()) continue;
                if (line.at(0) == '%') continue;

		std::istringstream ss(line);
		std::string token;
		string indata[3] = {"", "", ""};

		idx = 0;
		while(std::getline(ss, token, ' ')) {
			indata[idx++] = token;
		}	
		if (first_line){
                	first_line = false;
			n = stoi(indata[1]);
			m = 2*n;
			k_var = stoi(indata[0]);

		}else{
			int row = stoi(indata[0]) - 1;
			int col = stoi(indata[1]) - 1;
			double value = stod(indata[2]); 
			A.insert(row,col) = value;
		}
		
	}
	A_File.close();


	ifstream b_File("../example_problems/lp_afiro/lp_afiro_b.mtx");
	if (!b_File) {
			cerr << "Error Loading from UF dataset. File not found (b vector)" << endl;
			return;
	}

	first_line = true;
	int col_idx = 0;
	while (getline(b_File, line)) {
			if (line.empty()) continue;
			if (line.at(0) == '%') continue;

			std::istringstream ss(line);
			std::string token;

			if (first_line){
					first_line = false;
			}else{
				string indata[1] = {""};
				idx = 0;
				while(std::getline(ss, token, ' ')) {
						indata[idx++] = token;
				}
				double value = std::stod(indata[0]);
				b[col_idx++] = value;
			}

	}
	b_File.close();


	ifstream c_File("../example_problems/lp_afiro/lp_afiro_c.mtx");
	if (!c_File) {
			cerr << "Error Loading from UF dataset. File not found (c vector)" << endl;
			return;
	}

	first_line = true;
	col_idx = 0;
	while (getline(c_File, line)) {
			if (line.empty()) continue;
			if (line.at(0) == '%') continue;

			std::istringstream ss(line);
			std::string token;                

			if (first_line){
					first_line = false;
			}else{
				string indata[1] = {""};
				idx = 0;
				while(std::getline(ss, token, ' ')) {
						indata[idx++] = token;
				}
				double value = std::stod(indata[0]);
					c[col_idx++] = value;
			}

	}
	c_File.close();

        // initialize G
        for (int i = 0; i < n; i++)
        {
           G.insert(i, i) =  1;
           G.insert(n+i, i) =  -1;
        }

        ifstream hi_File("../example_problems/lp_afiro/lp_afiro_hi.mtx");
        if (!hi_File) {
                        cerr << "Error Loading from UF dataset. File not found (high vector)" << endl;
                        return;
        }

        first_line = true;
        col_idx = 0;
        while (getline(hi_File, line)) {
                        if (line.empty()) continue;
                        if (line.at(0) == '%') continue;

                        std::istringstream ss(line);
                        std::string token;

                        if (first_line){
                                        first_line = false;
                        }else{
                            string indata[1] = {""};
                            idx = 0;
                            while(std::getline(ss, token, ' ')) {
                                            indata[idx++] = token;
                            }
                                        double value = std::stod(indata[0]);
                                        h[col_idx++] = value;
                        }
        }
        hi_File.close();


        ifstream lo_File("../example_problems/lp_afiro/lp_afiro_lo.mtx");
        if (!lo_File) {
                        cerr << "Error Loading from UF dataset. File not found (low vector)" << endl;
                        return;
        }

        first_line = true;
        col_idx = 0;
        while (getline(lo_File, line)) {
                        if (line.empty()) continue;
                        if (line.at(0) == '%') continue;

                        std::istringstream ss(line);
                        std::string token;
                        if (first_line){
                               first_line = false;

                        }else{
                        		
                                    string indata[1] = {""};
                                    idx = 0;
                                    while(std::getline(ss, token, ' ')) {

                                                    indata[idx++] = token;
                                    }


                                        double value = std::stod(indata[0]);

                                        h[n + col_idx++] = value;
                        }
        }
        lo_File.close();


}


int main()
{
	const int max_iter			    = 50;
	const double linear_feas_tol 	= 1e-8; //Assuming possible Integer Overflow
	const double comp_tol			= 1e-8; //Assuming possible Integer Overflow
	const double bkscale			= 0.95;
	const double regularization     = 1e-7;

	//copl_vector test(10,10);

	OUTPUT << "COPL 2015" << endl;
	OUTPUT << "Interior point algorithm coming" << endl;
	
	// Initialize configuration variable
	lp_settings settings(max_iter,linear_feas_tol,comp_tol,bkscale,regularization);	
    
    
    string problem_name = "lp_afiro";
    int k_var,n;
    getDimensionUF(problem_name, &k_var, &n);
    OUTPUT << k_var << ":" << n << endl;
    copl_matrix A(k_var,n);
    copl_matrix G(2*n,n);
    copl_vector c(n),b(k_var),h(2*n);
    //make_trivial_problem(A,G,c,b,h);
    loadFromUF("", "", A,G,c,b,h);
    lp_input problem_data(A,b,c,G,h);
    
     //problem_data.var_dump();
     
	// The main function that run interior point algorithm.
	interior_point_algorithm(problem_data,settings);
	
};

