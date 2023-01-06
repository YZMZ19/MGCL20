/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
//////////
//A class that contains all the necessary input parameter to
//construct mgTLInputParam. This is used to serialize the tessellation parameters'

#include "StdAfx.h"
#include <math.h>
#include "mg/Ifstream.h"
#include "mg/Ofstream.h"
#include "mg/Tolerance.h"
#include "Tl2/TLInputParam.h"

mgTLInputParam::mgTLInputParam(
	double crvTol,			//バウンダリのトレランス
	double surfTol,			//平面とみなすトレランス
	double max_ratio,	//最大アスペクト比
	MGCL::fan_kind fk,
		//fk=SINGLE_TRIANGLE:   1 triangle/FAN
		//fk=MULTIPLE_TRIANGLES: as many triangles as possible/FAN
	int minimum_tri,	//Specify minimum number of triangles.
	double max_edge_len	//when max_edge_len<=0, this means no limits on an edge length.
):m_crvTol(crvTol), m_surfTol(surfTol), m_max_ratio(max_ratio),
m_fk(fk), m_minimum_tri(minimum_tri),m_max_edge_len(max_edge_len){
;}

//Construct from the object box data and the span length to draw object.
//span_length=MGOpenGLView::span_length().
mgTLInputParam::mgTLInputParam(
	const MGObject& obj,
	double span_length
):m_max_ratio(2.),m_fk(MGCL::MULTIPLE_TRIANGLES), m_minimum_tri(4),
m_max_edge_len(-1.){
	m_crvTol=compute_curve_tolerance(obj,span_length);
	m_surfTol=m_crvTol*SURFACE_TOL_BY_CURVE;
}

//////////// Operator overload.////////

std::ostream& operator<< (std::ostream& out, const mgTLInputParam& para){
	out<<"mgTLInputParam="<<(&para)<<"::m_crvTol="<<para.m_crvTol;
	out<<", m_surfTol="<<para.m_surfTol<<", m_max_ratio="<<para.m_max_ratio;
	out<<", m_fk=";
	if(para.m_fk==MGCL::SINGLE_TRIANGLE) out<<"SINGLE_TRIANGLE";
	else out<<"MULTIPLE_TRIANGLES";
	out<<", m_minimum_tri="<<para.m_minimum_tri<<", m_max_edge_len="<<para.m_max_edge_len;
	return out;
}

// Serialization fucntion.
MGOfstream& operator<< (MGOfstream& buf, const mgTLInputParam& para){
	buf<<para.m_crvTol;
	buf<<para.m_surfTol;
	buf<<para.m_max_ratio;
	buf<<int(para.m_fk);
	buf<<para.m_minimum_tri;
	buf<<para.m_max_edge_len;
	return buf;
}
MGIfstream& operator>> (MGIfstream& buf, mgTLInputParam& para){
	buf>>para.m_crvTol;
	buf>>para.m_surfTol;
	buf>>para.m_max_ratio;
	int fk;
	buf>>fk;
	para.m_fk=static_cast<MGCL::fan_kind >(fk);
	buf>>para.m_minimum_tri;
	buf>>para.m_max_edge_len;
	return buf;
}
