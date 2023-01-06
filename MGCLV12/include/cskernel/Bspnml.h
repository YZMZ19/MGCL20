/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
// BSPNML GETS NORMAL VECTOR OF SURFACE B-REP. 
// *** INPUT  * 
//     KU,LUD,UKT(LUD+KU),KV,LVD,VKT(LVD+KV),SURF(ISR1,ISR2,3) 
//            ...SURFACE B-REP 
//     U,V.....PARAMETER VALUES OF THE ABOVE SURFACE B-REP. 
// *** OUTPUT * 
//     VNML(3) : NORMAL VECTOR 
void bspnml_(int ku, int lud,const double *ukt, 
	int kv, int lvd,const double *vkt,const double *surf, int isr1, int isr2,
	double u, double v, double *vnml);
