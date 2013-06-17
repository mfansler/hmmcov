void IAL(double* y, double** X, double** Xnew, double* b0,   
         double* b, double* v, double delta1, double tau1,    
         int* dims, double *Xj2, double conv, double* resid, 
         double* BIC, int* n2kp, int* w2kp, double* b2kp, 
         int initIAL, int trace, int nReEstimate, 
	 double *weights, double *prop);
