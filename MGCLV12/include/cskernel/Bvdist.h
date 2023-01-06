/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
// BVDIST computes distance of two points P1 and P2. 
// INPUT *** NCD is the space dimension of the two points P1 
//           and P2. 
//           P1(NCD), P2(NCD)....two points. 
// OUTPUT *** BVDIST : DISTANCE OF P1 AND P2. 
double bvdist_(int ncd,const double *p1,const double *p2);
