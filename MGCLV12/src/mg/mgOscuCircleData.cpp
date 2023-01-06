/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/OscuCircleData.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
// Implement Class MGOscuCircleData.

//Constructor
MGOscuCircleData::MGOscuCircleData(
	int index,	//Index
	double radius)	//Rdius
:m_index(index),m_radius(radius){
	assert(radius>0.);
}

bool MGOscuCircleData::operator< (const MGOscuCircleData& ocd2)const{
	if(m_index==ocd2.m_index) return m_radius<ocd2.m_radius;
	return m_index<ocd2.m_index;
}

bool MGOscuCircleData::operator== (const MGOscuCircleData& ocd2)const{
	if(m_index!=ocd2.m_index) return false;
	return MGREqual(m_radius, ocd2.m_radius);
}
