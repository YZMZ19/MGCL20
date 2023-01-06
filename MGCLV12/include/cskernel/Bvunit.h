/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
// BVUNIT computes the vector VUNIT of the input vector V. 
// *** INPUT *** 
//   V(NCD),NCD..... IS VECTOR OF SPACE DIMENSION NCD. 
// *** OUTPUT *** 
//   VUNIT.......... is output unit vector. 
// *** NOTE *** 
// V must not be zero vector, this causes zero division abort. 
void bvunit_(const double *v, int ncd, double *vunit);
