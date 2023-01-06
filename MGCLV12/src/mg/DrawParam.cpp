/****************************************************************/
/*   Copyright (c) 2019 by System fugen G.K.                */
/*                       All rights reserved.                   */
/****************************************************************/
//////////
//A class that contains the necessary input parameter to draw curves and surfaces

#include "StdAfx.h"
#include "mg/Box.h"
#include "mg/Object.h"
#include "mg/drawParam.h"
#include "mg/Ifstream.h"
#include "mg/Ofstream.h"
#include "mg/Tolerance.h"
#include "tl2/TlInputParam.h"
#include "mgGL/VBO.h"

void MGDrawParam::build_crv_srf_tolerance(
	double curve_tolerance,//Maximum deviation allowed to approximate curves.
	double surface_tolerance//Maximum deviation allowed to approximate surfaces.
){
	double wczero=MGTolerance::wc_zero();
	if(curve_tolerance<wczero)
		curve_tolerance=wczero;
	if(surface_tolerance<curve_tolerance)
		surface_tolerance=curve_tolerance*SURFACE_TOL_BY_CURVE;

	m_curve_tolerance_tess=curve_tolerance;
	m_surface_tolerance_tess=surface_tolerance;
}

MGDrawParam::MGDrawParam(
	double curve_tolerance,//Maximum deviation allowed to approximate curves.
	double surface_tolerance,//Maximum deviation allowed to approximate surfaces.
	double maximum_edge_length,//Maximum edge length for the tessellation.
	double span_length_wire,//length of a line segment to approximate curves.
	int line_desity_wire_face//Indicates how many wires be drawn for the inner lines of
		///the wire representation of the face.
):m_span_length_wire(span_length_wire),m_line_desity_wire_face(line_desity_wire_face),
m_maximum_edge_length_tess(maximum_edge_length){
	build_crv_srf_tolerance(curve_tolerance,surface_tolerance);
}

//Construct from the object box data and the span length to draw object.
//span_length=MGOpenGLView::span_length().
MGDrawParam::MGDrawParam(
	const MGObject& obj,
	double span_length
):m_span_length_wire(span_length),m_maximum_edge_length_tess(-1.),m_line_desity_wire_face(1){
	m_curve_tolerance_tess=compute_curve_tolerance(obj,span_length);
	m_surface_tolerance_tess=m_curve_tolerance_tess*SURFACE_TOL_BY_CURVE;
}

MGDrawParam::MGDrawParam(
	const mgTLInputParam& tlpara
):m_span_length_wire(20),m_line_desity_wire_face(1),
m_maximum_edge_length_tess(tlpara.max_edge_len()){
	build_crv_srf_tolerance(tlpara.crvTol(),tlpara.surfTol());
}

//////////// Operator overload.////////

std::ostream& operator<< (std::ostream& out, const MGDrawParam& para){
	out<<"MGDrawParam="<<(&para)<<"::m_curve_tolerance_tess="<<para.m_curve_tolerance_tess;
	out<<", m_surface_tolerance_tess="<<para.m_surface_tolerance_tess;
	out<<", m_maximum_edge_length_tess="<<para.m_maximum_edge_length_tess<<std::endl;
	return out;
}

// Serialization fucntion.
MGOfstream& operator<< (MGOfstream& buf, const MGDrawParam& para){
	buf<<para.m_curve_tolerance_tess;
	buf<<para.m_surface_tolerance_tess;
	buf<<para.m_maximum_edge_length_tess;
	return buf;
}
MGIfstream& operator>> (MGIfstream& buf, MGDrawParam& para){
	buf>>para.m_curve_tolerance_tess;
	buf>>para.m_surface_tolerance_tess;
	buf>>para.m_maximum_edge_length_tess;
	return buf;
}

double compute_curve_tolerance(
	const MGObject& obj,	///<Target object 
	double span_length		///<Span length of MGOpenglView.
){
	double object_length=obj.box().length();
	object_length*=3.;
	int div_num=int(object_length/span_length*20.+1.)+1;
	if(div_num>MAXDIVNUM)
		div_num=MAXDIVNUM;
	return object_length*(1.-cos(mgPAI/div_num));
}