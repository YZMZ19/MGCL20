/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGDefaultVector_HH_
#define _MGDefaultVector_HH_

#include "mg/MGCL.h"
#include "mg/Position.h"

//  MGDefault.h
//  header for class MGDefault

/** \addtogroup BASE
 *  @{
 */

// Difines default values of each class.
extern MG_DLL_DECLR const MGPosition mgNULL_Pos;
extern MG_DLL_DECLR const MGPosition mgORIGIN;
extern MG_DLL_DECLR const MGPosition mgORIGIN_2D;

extern MG_DLL_DECLR const MGVector mgX_UVEC;
extern MG_DLL_DECLR const MGVector mgY_UVEC;
extern MG_DLL_DECLR const MGVector mgZ_UVEC;

extern MG_DLL_DECLR const MGVector mgX_UVEC_2D;
extern MG_DLL_DECLR const MGVector mgY_UVEC_2D;

extern MG_DLL_DECLR const MGVector mgNULL_VEC;
extern MG_DLL_DECLR const MGVector mgZERO_VEC;
extern MG_DLL_DECLR const MGVector mgZERO_VEC_2D;

/** @}*/ //BASE
#endif
