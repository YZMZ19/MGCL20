/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
//   REAL FUNCTION TO EVALUATE JDERIV-TH DERIVATIVE AT THE PARAMETER 
//   VALUE X OF THE B-REP (N,T,RCOEF). 
// *** INPUT *
//     K,N,T(N+K),RCOEF(N)......B-REP TO EVALUATE,    ORDER, B-REP DIM- 
//     X        VALUE AT WHICH THE DERIVATIVE IS EVALUATED. 
//     JDERIV   ORDER OF THE DERIVATIVE,MAY BE ZERO. 
// *** OUTPUT *
//     BLER    THE VALUE OF THE B-SPLINE AT THE PARMETER X. 
double bler_(int k, int n, const double *t, const double *rcoef, double x, int jderiv);
