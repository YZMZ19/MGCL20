/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include <cmath>

// BVUNIT computes the vector VUNIT of the input vector V. 
// *** INPUT *** 
//   V(NCD),NCD..... IS VECTOR OF SPACE DIMENSION NCD. 
// *** OUTPUT *** 
//   VUNIT.......... is output unit vector. 
// *** NOTE *** 
// V must not be zero vector, this causes zero division abort. 
void bvunit_(const double *v, int ncd, double *vunit){
    // Builtin functions 

    // Local variables 
    int i;
    double vv;

    // Function Body 
    vv = 0.f;
    for(i=0; i<ncd; ++i)
		vv += v[i]*v[i];
    vv = sqrt(vv);
    for(i=0; i<ncd; ++i)
		vunit[i] = v[i]/vv;
}
