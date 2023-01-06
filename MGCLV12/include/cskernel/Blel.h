/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
// REAL FUNCTION TO EVALUATE LEFT-CONTINUOUS JDERIV-TH DERIVATIVE 
// AT THE PARAMETER VALUE X OF THE B-REP. (T,RCOEF). 
// *** INPUT * 
//     K,N,T(N+K),RCOEF(N).......B-REP TO EVALUATE,  ORDER, B-REP DIM- 
//              EMNSION, KNOT VECTOR, AND B-COEFFICIENTS EACH. 
//     X        VALUE AT WHICH THE DERIVATIVE IS EVALUATED. 
//     JDERIV   ORDER OF THE DERIVATIVE,MAY BE ZERO. 
// *** OUTPUT * 
//     BLEL   THE VALUE OF THE B-SPLINE AT THE PARMETER X. 
// *** NOTE * 
//     FUNCTION BLER EVALUATES RIGHT-CONTINUOUS DERIVATIVE WHILE 
//     BLEL DOES LEFT-CONTINUOUS. 
double blel_(int k, int n, const double *t, const double *rcoef, double x, int jderiv);