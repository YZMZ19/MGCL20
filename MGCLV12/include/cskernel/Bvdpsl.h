/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
// FUNCTION TO OBTAIN DISTANCE BETWEEN a LINE(G) AND a POINT(P). 
// INPUT *** 
//    NCD......SPACE DIMENSION OF G AND P, MUST BE 2 OR 3. 
//    G(NCD,2) : PARAMETER OF S.L. AS BELOW 
//            G(.,1):A POINT ON THE S.L., G(.,2):DIRECTIONAL COSINE 
//    P(NCD) : COORDINATE OF A POINT 
// OUTPUT *** BVDPSL : DISTANCE 
double bvdpsl_(int ncd,const double *g,const double *p);
