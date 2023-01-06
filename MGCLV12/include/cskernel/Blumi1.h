/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
// PRIVATE SUB OF BLUMIX, ADD THE FOLLOWING DATA IN THE SEQ , 
//    (KSNEW(NADD),VALNEW(3,NADD),TAUNEW(NADD))  AT VAL(IPOS,3) 
//   ****DATA INSERTED AFTER IPOS**** 
// IFLAG......=0:SUCCESSFULL RETURN, =2:AREA OF VAL EXAUSTED. 
void blumi1_(int *n, int *kseq, double *val, 
	int iv, double *tau, int ipos, int nadd,const int *ksnew,
	const double *valnew, const double *taunew, int *iflag);

// PRIVATE SUB OF BLUMIX, DELETE NDEL DATA FROM THE SEQ 
void blumi3_(int *n, int *kseq, double *val, 
	int iv, double *tau, int ipos, int ndel);
