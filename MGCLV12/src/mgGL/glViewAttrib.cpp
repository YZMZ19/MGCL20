/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include <bitset>
#include "mg/Ifstream.h"
#include "mg/Ofstream.h"
#include "mgGL/glViewAttrib.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

//////////////////////////////METHOD//////////////////
MGglViewAttrib::MGglViewAttrib(bool is_perspective)
:m_viewMode(MGCL::WIREVIEW),m_perspective(is_perspective),
m_fovy(45.),m_center(0.,0.,0.),m_scale(INITIAL_SCALE),
m_cx(0.),m_cy(0.),m_diameter(1.),
m_eyeP(0.,0.,1.),m_up_vector(MGVector(0.,1.,0.)),
m_PreCenterMat(1.),m_modelViewMat(1.)
{
	if(!is_perspective)
		m_fovy=0.;
}

MGglViewAttrib::MGglViewAttrib(
	const MGBox& box, const MGColor* gridColors
):m_viewMode(MGCL::WIREVIEW),m_perspective(true),m_fovy(45.){
	m_cplane.setGridDataByBox(box,1,gridColors);//Initialize construction plane by the box.
	const MGPlane& pl=m_cplane.plane();
	const MGVector& uderi=pl.u_deriv();
	const MGVector& vderi=pl.v_deriv();
	setEyePositionUpVector(uderi*vderi, vderi);

	compute_viewing_environment(box);
}

///Compute the viewing environment the parameter box.
///compute_viewing_environment() uses (eye_position,view_up_vector) as input.
///They must be set before compute_viewing_environment.
///eye_position is used only to get the direction from the origin.
void MGglViewAttrib::compute_viewing_environment(
	const MGPosition& center,
	double diameter
){
	setHomeMatrix();

	m_center=center;
	m_diameter=diameter;
	if(m_diameter<1.)
		m_diameter=1.;
	double diam10=m_diameter*10.;

	double tan_theta2=1.;
	if(m_fovy>0.){
		tan_theta2=tan(MGCL::degree_to_radian(m_fovy*.5))*2.;
	}
	m_near=diam10/tan_theta2;
	m_far=m_near+diam10;
	double eyeLength=m_near+.5*diam10;
	MGUnit_vector eyeP(m_eyeP);
	m_eyeP=eyeP*eyeLength;
}

#define INITIAL_BOX_SCALE INITIAL_SCALE*1.1
///Compute the viewing environment the parameter box.
///compute_viewing_environment() uses (eye_position,view_up_vector) as input.
///They must be set before compute_viewing_environment.
///eye_position is used only to get the direction from the origin.
void MGglViewAttrib::compute_viewing_environment(
	const MGBox& box///<Input the box of the target scene.
){
	MGBox box2=box*INITIAL_BOX_SCALE;
	double diameter=box2.len();
	compute_viewing_environment(box.mid(),diameter);
}

void MGglViewAttrib::setHomeMatrix(){
	m_modelViewMat=m_PreCenterMat=glm::mat4(1.);
	m_scale=INITIAL_SCALE;
	m_cx=m_cy=0.;
}

///Set m_eyeP and m_up_vector(the eye position and view-up-vector of the OpenGL).
void MGglViewAttrib::setEyePositionUpVector(
	const MGPosition& eyeP,
	const MGVector& upVector
){
	m_eyeP=eyeP;
	m_up_vector=upVector;
}

///Set if this view is a perspective view(true), or orthnormal view(falsle).
void MGglViewAttrib::set_perspective(bool pers, double fovy){
	m_perspective=pers;
	if(pers)
		m_fovy=fovy;
	else
		m_fovy=0.;
}

// Serialization.
MGOfstream& operator<< (MGOfstream& buf, const MGglViewAttrib& atr){
	buf<<atr.m_cplane;

	std::bitset<32> boolData;
	boolData.set(0,atr.m_perspective);
	buf << boolData.to_ulong();

	buf<<atr.m_fovy;
	buf<<atr.m_near;
	buf<<atr.m_far;
	atr.m_eyeP.dump(buf);
	atr.m_up_vector.dump(buf);
	atr.m_center.dump(buf);
	buf<<atr.m_diameter;
	buf<<atr.m_scale;
	buf<<atr.m_cx;
	buf<<atr.m_cy;

	const float* premat=&atr.m_PreCenterMat[0][0];
	for(int i=0; i<16; i++){
		buf<<premat[i];
	}
	const float* mvmat=&atr.m_modelViewMat[0][0];
	for(int i=0; i<16; i++){
		buf<<mvmat[i];
	}
	int viewMode=atr.m_viewMode;
	buf<<viewMode;
	return buf;
}

MGIfstream& operator>> (MGIfstream& buf, MGglViewAttrib& atr){
	buf>>atr.m_cplane;

	int lbit;
	buf >> lbit;
	std::bitset<32> boolData(lbit);
	atr.m_perspective=boolData[0];

	buf>>atr.m_fovy;
	buf>>atr.m_near;
	buf>>atr.m_far;
	atr.m_eyeP.restore(buf);
	atr.m_up_vector.restore(buf);
	atr.m_center.restore(buf);
	buf>>atr.m_diameter;
	buf>>atr.m_scale;
	buf>>atr.m_cx;
	buf>>atr.m_cy;

	float* premat=&atr.m_PreCenterMat[0][0];
	for(int i=0; i<16; i++){
		buf>>premat[i];
	}
	float* mvmat=&atr.m_modelViewMat[0][0];
	for(int i=0; i<16; i++){
		buf>>mvmat[i];
	}
	int viewMode;
	buf>>viewMode;
	atr.m_viewMode=(MGCL::VIEWMODE)viewMode;
	return buf;
}

static const std::string viewMode[5]={
	"DONTCARE",
	"WIRE",		///< wire frame mode
	"SHADING",	///< surface mode
	"WIRE_AND_SHADING", ///< wire and surface mode
	"HIGHLIGHT"
};

//Debug Function.
std::ostream& operator<< (std::ostream& out, const MGglViewAttrib& atr){
	out<<"MGglViewAttrib="<<&atr<<",m_perspective="<<atr.m_perspective;
	out<<",m_fovy="<<atr.m_fovy<<",m_near="<<atr.m_near<<",m_far="<<atr.m_far;
	out<<",m_eyeP="<<atr.m_eyeP<<std::endl<<",m_up_vector="<<atr.m_up_vector<<",m_center="<<atr.m_center;
	out<<",m_diameter="<<atr.m_diameter<<",m_scale="<<atr.m_scale;
	out<<",(m_cx,m_cy)=("<<atr.m_cx<<","<<atr.m_cy<<")";
	out<<std::endl<<",m_cplane="<<atr.m_cplane<<std::endl;

	const float* premat=&atr.m_PreCenterMat[0][0];
	out<<",m_PreCenterMat=["<<premat[0];
	for(int i=1;i<16; i++) out<<","<<premat[i];
	out<<"]";

	const float* mvmat=&atr.m_modelViewMat[0][0];
	out<<",m_modelViewMat=["<<mvmat[0];
	for(int i=1;i<16; i++) out<<","<<mvmat[i];
	out<<"]";

	int vmode=atr.m_viewMode;
	out<<",m_viewMode="<<viewMode[vmode];
	out<<std::endl;
	return out;
}
