/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
// BVCANG EVALUATES COSINE ANGLE OF TWO VECTORS. LET THETA IS 
// THE ANGLE OF TWO VECTORS, THEN BVCANG=COS(THETA). 
// INPUT *** 
//     NCD...... is the space dimension of the two vector V1 
//               and V2. 
//     IV1,IV2....ROW DIMENSION OF VECTOR V1 AND V2. 
//     V1(IV1,NCD),V2(IV2,NCD)........ are two vectors, input as 
//               VECTOR V1(1,.) and VECTOR V2(1,.). 
// OUTPUT *** BVCANG : COSINE ANGLE OF THE TWO VECTOR. 
double bvcang_(int ncd,const double *v1,const double *v2);
