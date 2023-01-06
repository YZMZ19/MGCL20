/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
// BVINPD computes inner product of two vectors V1 and V2. 
// *** INPUT *** 
//   NCD,V1(NCD),V2(NCD)..... are two input vectors. 
//         NCD IS THE SPACE DIENSION. 
// *** OUTPUT *** 
//   BVINPD........is scalar data of the inner product. 
double bvinpd_(int ncd,const double *v1,const double *v2);
