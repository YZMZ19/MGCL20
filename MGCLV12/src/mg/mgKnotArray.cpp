/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/KnotArray.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
// Implemetation of MGKnotArray Class.

//Friend Function

//Constructor
MGKnotArray::MGKnotArray(const MGKnot& knot)		//From Knot.
:m_vec(1,knot){;}

MGKnotArray::MGKnotArray(double knot, int mult)
//From knot and the multiplicity.
:m_vec(1,MGKnot(knot,mult))
{assert(mult>0);}
