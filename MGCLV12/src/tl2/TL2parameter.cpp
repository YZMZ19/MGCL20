#include "StdAfx.h"
#include "mg/Tolerance.h"
#include "mg/KnotVector.h"
#include "mg/FSurface.h"
#include "topo/Edge.h"
#include "topo/Loop.h"
#include "topo/Face.h"
#include "topo/Shell.h"
#include "Tl2/TL2parameter.h"
#include "mgGL/VBO.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

/****************************************************************/
/*   Copyright (c) 2019 by System fugen G.K.                */
/*                       All rights reserved.                   */
/****************************************************************/

//mgTL2parameter is a proprietry class for Face tessellation.
//mgTL2parameter represents one intesection of Loop and u=const(or v=const)
//line in the parameter space of a face.
//When intersection with u=const, m_t is v value, or vise versa.

///Compute maximum edge length for the tessellation from an object, twoManifold.
///twoManifold must be MGFSurface or MGShell.
double compute_max_edge_len(const MGObject& twoManifold){
	const MGBox& box=twoManifold.box();
	double s[3] = {box[0].length(), box[1].length(),box[2].length()};
	std::sort(s, s + 3);
	double max_edge_len = s[2]/EDGE_LENGTH_DENOM;
	double elenmax2 = s[1] * 1.2;
	if (max_edge_len> elenmax2)
		max_edge_len = elenmax2;
	const MGShell* sh=dynamic_cast<const MGShell*>(&twoManifold);
	if(sh){
		int n=sh->number_of_faces();
		max_edge_len/=double(1.+log(double(n)));
	}
	return max_edge_len;
}

void mgTL2parameter::build_parameter(
	const MGFSurface& srf,
	double crvTol,
	double surfTol,
	double max_edge_len
){
	if(crvTol<0.){
		double length=srf.get_box().length();
		length*=.05;
	}
	if(surfTol>0. && surfTol<crvTol)
		crvTol=surfTol;

	double wczero=MGTolerance::wc_zero();
	if(crvTol<wczero*3.)
		crvTol=wczero*3.;
	if(surfTol<=crvTol)
		surfTol=crvTol*SURFACE_TOL_BY_CURVE;

	m_tess_crvError=crvTol;
	m_tess_srfError=surfTol;

	double u0=srf.param_s_u(), u1=srf.param_e_u();
	double v0=srf.param_s_v(), v1=srf.param_e_v();
	double rczero=MGTolerance::rc_zero();
	m_puerror=rczero*(u1-u0);
	m_pverror=rczero*(v1-v0);
	m_uverror=sqrt(m_puerror*m_puerror+m_pverror*m_pverror);

	m_max_edge_len=max_edge_len;
	if(m_max_edge_len<0.){
		m_max_edge_len=compute_max_edge_len(*(srf.object_pointer()));
	}

	double crvTol10=crvTol*10.;
	if(m_max_edge_len<crvTol10)
		m_max_edge_len=crvTol10;
	m_max_edge_len_sqr=m_max_edge_len*m_max_edge_len;
}

mgTL2parameter::mgTL2parameter(
	const MGFSurface& obj,
	double crvTol,
	double surfTol,
	const std::vector<SHLL_COM_EDGES>* polylines,
		///< Input polygonized polylines for the face boundaries.
		///< polylines[i][j] is a j-th edge's polyline for face.loop(i),
		///< must be MGLBRep of order 2.
		///< polylines[i][j]=0 indicates loop i's edge j can be face's bounday and
		///< has any common edges.
		///< **polylines[i][j] must be the same direction as the faces's parameter edge.
	double max_edge_len
):m_face(obj.get_face_pointer()),m_surface(obj.get_surface_pointer())
,m_Bpolylines(polylines){
	build_parameter(obj,crvTol,surfTol,max_edge_len);
}

mgTL2parameter::mgTL2parameter(
	const MGFSurface& obj,
	const MGDrawParam& param,//parameter for the tessellation.
	const std::vector<SHLL_COM_EDGES>* polylines
		///< Input polygonized polylines for the face boundaries.
		///< polylines[i][j] is a j-th edge's polyline for face.loop(i),
		///< must be MGLBRep of order 2.
		///< polylines[i][j]=0 indicates loop i's edge j can be face's bounday and
		///< has any common edges.
		///< **polylines[i][j] must be the same direction as the faces's parameter edge.
):m_face(obj.get_face_pointer()),m_surface(obj.get_surface_pointer()),
m_Bpolylines(polylines){
	double ctol=param.curve_tolerance_tess();
	double stol=param.surface_tolerance_tess();
	double elen=param.maximum_edge_length_tess();
	build_parameter(obj,ctol,stol,elen);
}

//Copy constructor.
mgTL2parameter::mgTL2parameter(const mgTL2parameter& param2)
:m_face(param2.m_face),m_surface(param2.m_surface),
m_max_edge_len(param2.m_max_edge_len),
m_max_edge_len_sqr(param2.m_max_edge_len_sqr),
m_tess_crvError(param2.m_tess_crvError), m_tess_srfError(param2.m_tess_srfError),
m_puerror(param2.m_puerror),m_pverror(param2.m_pverror){
}

///////Operator oveload///////
mgTL2parameter& mgTL2parameter::operator=(const mgTL2parameter& param2){
	m_face=param2.m_face;
	m_surface=param2.m_surface;
	m_puerror=param2.m_puerror;
	m_pverror=param2.m_pverror;

	m_tess_crvError=param2.m_tess_crvError;
	m_tess_srfError=param2.m_tess_srfError;

	m_max_edge_len=param2.m_max_edge_len;
	m_max_edge_len_sqr=param2.m_max_edge_len_sqr;
	return *this;
}

//Update the curve error. Returned is the old error.
double mgTL2parameter::set_tess_crvError(double error) {
	double old=m_tess_crvError;
	m_tess_crvError=error;
	return old;
}

//Update the surface error. Returned is the old error.
double mgTL2parameter::set_tess_srfError(double error) {
	double old = m_tess_srfError;
	m_tess_srfError = error;
	return old;
}

ostream& operator<< (ostream& out, const mgTL2parameter& para){
	out<<"TLparam::m_face="<<para.m_face<<",m_surface="<<para.m_surface
		<<",m_tess_crvError="<<para.m_tess_crvError;
	out<<",m_tess_srfError="<<para.m_tess_srfError;
	out<<"m_Bpolylines="<<para.m_Bpolylines;
	if(para.m_Bpolylines){
		const std::vector<SHLL_COM_EDGES>& Bp=*(para.m_Bpolylines);
		size_t nBp=Bp.size();
		out<<"Number of Bpolyline="<<nBp<<endl;
		for(size_t j=0; j<nBp; j++){
			const SHLL_COM_EDGES& loopj=Bp[j];
			for(size_t k=0; k<loopj.size(); k++){
				out<<"("<<j<<","<<k<<")-th Edge::"<<loopj[k]<<endl;
			}
		}
	}
	out<<",m_puerror="<<para.m_puerror<<",m_pverror="<<para.m_pverror;
	out<<", m_max_edge_len="<<para.m_max_edge_len<<endl;
	return out;
}
