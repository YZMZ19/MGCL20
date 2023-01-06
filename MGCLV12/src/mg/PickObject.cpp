/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
// MGPickObject.cpp : Implements MGPickObject
//
/////////////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "mg/PickObject.h"
#include "mgGL/VBO.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif


//Generate a newed clone object.
MGPickObject* MGPickObject::clone()const{
	return new MGPickObject(*this);
}

bool MGPickObject::operator<(const MGPickObject& pobj) const{
	if(static_cast<MGGelPosition>(*this)!=static_cast<MGGelPosition>(pobj))
		return static_cast<MGGelPosition>(*this)<static_cast<MGGelPosition>(pobj);
	return m_parameter<pobj.m_parameter;
}