/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
// BLGCS2 is a dedicated subroutine of BLGCS, generates data point 
// sequence TAU(.) from data (VAL(j,NCD),j=1...N). VAL(j,.) may 
// include circle inf that is declared by KVAL(j). 
// ******Input****** 
//   NCD,N,KVAL(N),VAL(IV,NCD).....Input data of Space Dimension NCD, and 
//         of length N. KVAL(j) is a knuckle inf of VAL(j,.). 
// ******Output***** 
//   TAU(N).........Data point abssisa obtained of length N. 
void blgcs2_(int ncd, int n, const int *kval, 
	const double *val, int iv, double *tau
);
