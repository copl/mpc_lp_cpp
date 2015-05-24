% For Oli:
/cygdrive/c/Users/Oliver/Documents/GitHub/mpc_lp_cpp

Davood:
To Run, I had to
comment out line 131 of ./include/Eigen/src/Core/CwiseBinaryOp.h

Kvals[ (*hessianIx)[j] ] = -variables.s[j]/variables.z[j]-DELTA; -->
Kvals[ (*hessianIx)[j] ] = variables.s[j]/variables.z[j]-DELTA;

			eigenKMat->insert(r+n+p,c) = -it.value();	
			eigenKMat->insert(r+n,c) = -it.value();
	 		rhs_2[j++] = -rhs.q1[i]; 	

