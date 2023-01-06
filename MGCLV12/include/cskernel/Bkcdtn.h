/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
// bkcdtn_ WILL GENERATE DATA POINTS TAU(J) J=0,...,M-1 FROM 
// T(I) I=0,..N-1 SO AS TO BE PROPORTIONALLY DISTRIBUTED IN T. 
// *** INPUT *** 
// T(N) ..... OLD DATA POINT SEQUENCE (MAY BE KNOTVECTOR). 
// N ........ LENGTH OF T. 
// M ........ NUMBER OF NEW DATA POINTS TAU. 
// *** OUTPUT *** 
// TAU(M).... AREA TO STORE NEW DATA POINTS. 
void bkcdtn_(int n,const double *t, int m, double *tau);
