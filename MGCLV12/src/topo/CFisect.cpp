/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Position.h"
#include "mg/FSurface.h"
#include "topo/CFisect.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGCFisect Class.
//MGCFisect is to represent an intersection of a face and a curve.
//(MGCSisect csi, MGFSurface* f) where csi consists of world point, curve parameter,
//and face(surface) parameter, and f is a face pointer.

///////Constructor////////

//Construct from all the necessary data.
MGCFisect::MGCFisect(
	const MGPosition& point,	//World coordinate point data of the isect.
	const double& t,			//curve parameter value of the isect.
	const MGPosition& uv,		//Face(Surface) parameter value of the isect.
	const MGFSurface& face)			//face.
	:MGCSisect(point,t,uv), m_face(&face){;}

///////Operator oveload///////

bool MGCFisect::operator< (const MGCFisect& fp)const{
	return param_curve()<fp.param_curve();
}

//Ordering functions.
bool MGCFisect::operator< (const MGisect& is)const{
	const MGCFisect* cis=dynamic_cast<const MGCFisect*>(&is);
	if(cis) return operator<(*cis);
	return is>(*this);
}

bool MGCFisect::operator== (const MGCFisect& fp)const{
	if(m_face!=fp.m_face) return false;
	return operator==(fp);
}

bool MGCFisect::operator== (const MGisect& is)const{
	const MGCFisect* cis=dynamic_cast<const MGCFisect*>(&is);
	if(!cis) return false;
	return operator==(*cis);
}

///////Member function///////

// Output virtual function.
std::ostream& MGCFisect::toString(std::ostream& ostrm)const{
	ostrm<<"MGCFisect::m_face="<<m_face<<";";
	MGCSisect::toString(ostrm);
	return ostrm;
}
