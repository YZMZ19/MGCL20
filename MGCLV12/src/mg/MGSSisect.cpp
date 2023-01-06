/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Position.h"
#include "mg/SSisect.h"
#include "mg/CCisects.h"
#include "mg/Curve.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGSSisect.cc
// MGSSisect ‚ÌŽÀ‘•ƒtƒ@ƒCƒ‹

//
// Constructor
//

// Move Constructor;
MGSSisect::MGSSisect(MGSSisect&& ssi):m_rel(ssi.m_rel),
m_iline(ssi.m_iline.release()),
m_param1(ssi.m_param1.release()), m_param2(ssi.m_param2.release()){
;}

//Move Assignment
MGSSisect& MGSSisect::operator= (MGSSisect&& rhs){
	m_rel=rhs.m_rel;
	m_iline.reset(rhs.m_iline.release());
	m_param1=std::move(rhs.m_param1);
	m_param2=std::move(rhs.m_param2);
	return *this;
}

//Construct providing all the raw data.
//Copy version. Copy of the three curves will take place.
MGSSisect::MGSSisect (
	const MGCurve& iline,
	const MGCurve& param1,
	const MGCurve& param2,
	const MGSSRELATION rel)
:m_iline(iline.clone()),
m_param1(param1.clone()), m_param2(param2.clone()),m_rel(rel){;
}

bool MGSSisect::operator< (const MGSSisect& ssi2)const{
	return m_param1->start_point()<ssi2.m_param1->start_point();
}

bool MGSSisect::operator== (const MGSSisect& ssi2)const{
	if((*m_param1)!= *(ssi2.m_param1)) return false;
	if((*m_param2)!= *(ssi2.m_param2)) return false;
	return true;
}

bool MGSSisect::operator< (const MGisect& is)const{
	const MGSSisect* cis=dynamic_cast<const MGSSisect*>(&is);
	if(cis) return operator<(*cis);
	return is>(*this);
}
bool MGSSisect::operator== (const MGisect& is)const{
	const MGSSisect* cis=dynamic_cast<const MGSSisect*>(&is);
	if(!cis) return false;
	return operator==(*cis);
}

///Exchange 1st and 2nd order of the parameter line representation.
///When 1st object's manifold dimension is less than 2nd one's,
///this does nothing. Valid only when their manifold dimensions are equal.
void MGSSisect::exchange12(){
	m_param1.swap(m_param2);
}

//Test if two ssi's world curve have common parts (in line_zero()).
//Fucntion's return value is 
//		1:have commonpart.
//		0>=:no common part(except a point).
int MGSSisect::has_common(const MGSSisect& ssi2)const{
	std::vector<double> dvec;
	MGCCisects isects;
	int retcode=line().common(ssi2.line(),dvec,isects);
	if(retcode==1 || retcode==3) return 1;
	return 0;
}

//negate the direction of the intersection line.
void MGSSisect::negate(){
	m_iline->negate();	//(x,y,z)coordinate representaion of the line.
	m_param1->negate();	//(u,v) representaion of the line of the first
	m_param2->negate();	//(u,v) representaion of the line of the second
}

// Output virtual function.
std::ostream& MGSSisect::toString(std::ostream& ostrm)const{
//	ostrm.setf ( ios::scientific, ios::floatfield );
//	ostrm.precision ( 10 );
	if(m_iline){
		ostrm<<"MGSSisect::m_iline="<<(*m_iline);
		if(m_param1) ostrm<<", m_param1="<<(*m_param1);
		if(m_param2) ostrm<<", m_param2="<<(*m_param2);
		ostrm<<", m_rel="   <<(m_rel);
	}else{
		ostrm<<"MGSSisect::m_iline="<<m_iline;
	}
	return ostrm;
}

void MGSSisect::set_null(){
	m_iline=nullptr;
	m_param1=nullptr;
	m_param2=nullptr;
}
