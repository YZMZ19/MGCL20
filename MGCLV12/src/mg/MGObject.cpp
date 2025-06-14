/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Box.h"
#include "mg/Ofstream.h"
#include "mg/Ifstream.h"
#include "mg/Point.h"
#include "mg/Straight.h"
#include "mg/Ellipse.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/SurfCurve.h"
#include "mg/TrimmedCurve.h"
#include "mg/CompositeCurve.h"
#include "mg/BSumCurve.h"
#include "mg/Plane.h"
#include "mg/Sphere.h"
#include "mg/Cylinder.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/BSumSurf.h"
#include "mg/MGStl.h"
#include "mg/Tolerance.h"
#include "mg/GelFactory.h"
#include "mgGL/Appearance.h"

#include "topo/Complex.h"
#include "topo/PVertex.h"
#include "topo/BVertex.h"
#include "topo/Edge.h"
#include "topo/Face.h"
#include "topo/Loop.h"
#include "topo/Shell.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGObject
// Implementation of MGObject.

//Copy constructor.
MGObject::MGObject(const MGObject& obj2)
:m_box(obj2.m_box){
	if(obj2.m_appearance)
		m_appearance=std::unique_ptr<MGAppearance>(obj2.m_appearance->clone());
}

//Move constructor.
MGObject::MGObject(MGObject && obj2)
:MGAttribedGel(std::move(obj2)),m_appearance(std::move(obj2.m_appearance)),
m_box(std::move(obj2.m_box)){
}

MGObject::~MGObject() {///Destructor.
	;
}

//Move assignment.
MGObject& MGObject::operator=(MGObject && obj2){
	MGAttribedGel::operator=(std::move(obj2));
	m_appearance=std::move(obj2.m_appearance);
	m_box=std::move(obj2.m_box);
	return *this;
}

//Return minimum box that includes whole of the geometry.
const MGBox& MGObject::box() const{
	if(m_box.is_null())
		compute_box(m_box);
	return m_box;
}
bool MGObject::box_is_null() const{
	return m_box.is_null();
}
void MGObject::invalidateBox()const{
	m_box.set_null();
	setDirty(true);
}

//Assignment.
//When the leaf object of this and gel2 are not equal, this assignment
//does nothing.
MGObject& MGObject::set_object(const MGObject& gel2){
	MGAttribedGel::operator=(gel2);
	if(gel2.m_appearance){
		m_appearance.reset(gel2.m_appearance->clone());
	}else{
		m_appearance.reset();
	}

	m_box=gel2.m_box;
	return *this;
}

//make this group has appearance and get the MGAppearance pointer.
MGAppearance* MGObject::ensure_appearance(){
	if(!m_appearance)
		m_appearance.reset(new MGAppearance());
	return m_appearance.get();
}

//set the copy of appr2 to this MGAttribedgel.
void MGObject::set_appearance(const MGAppearance& appr2){
	m_appearance.reset(appr2.clone());
}
void MGObject::set_appearance(MGAppearance* appr2){
	m_appearance.reset(appr2);
}

//Test if this and 2nd object has common area about their box(),
//taking error into account.
bool MGObject::has_common(const MGObject& obj2) const{
	double err=MGTolerance::wc_zero(); err*=.5;
	MGBox bx=box(); bx.expand(err);
	MGBox bx2=obj2.box(); bx2.expand(err);
	bx&=bx2;
	return !bx.empty();
}

//Remove the MGAppearance of this MGAttribedGel.
std::unique_ptr<MGAppearance> MGObject::remove_appearance(){
	return std::move(m_appearance);
}

//Read all member data.
void MGObject::ReadMembers(MGIfstream& buf){
	m_appearance.reset(static_cast<MGAppearance*>(buf.ReadPointer()));
	m_box.restore(buf);
}

//Write all member data
void MGObject::WriteMembers(MGOfstream& buf)const{
	buf.WritePointer(m_appearance.get());
	m_box.dump(buf);
}

MGisects MGObject::isect(const MGFSurface & f) const{
	const MGObject* obj = dynamic_cast<const MGObject*>(&f);
	return isect(*obj);
}

AUTO_GEL_REGISTER(MGPoint, MGPOINT_TID);

AUTO_GEL_REGISTER(MGStraight, MGSTRAIGHT_TID);
AUTO_GEL_REGISTER(MGEllipse, MGELLIPSE_TID);
AUTO_GEL_REGISTER(MGLBRep, MGLBREP_TID);
AUTO_GEL_REGISTER(MGRLBRep, MGRLBREP_TID);
AUTO_GEL_REGISTER(MGSurfCurve, MGSRFCRV_TID);
AUTO_GEL_REGISTER(MGTrimmedCurve, MGTRMCRV_TID);
AUTO_GEL_REGISTER(MGCompositeCurve, MGCOMPCRV_TID);
AUTO_GEL_REGISTER(MGBSumCurve, MGBSUMCRV_TID);

AUTO_GEL_REGISTER(MGPlane, MGPLANE_TID);
AUTO_GEL_REGISTER(MGSphere, MGSPHERE_TID);
AUTO_GEL_REGISTER(MGSBRep, MGSBREP_TID);
AUTO_GEL_REGISTER(MGRSBRep, MGRSBREP_TID);
AUTO_GEL_REGISTER(MGCylinder, MGCYLINDER_TID);
AUTO_GEL_REGISTER(MGBSumSurf, MGBSUMSURF_TID);

AUTO_GEL_REGISTER(MGPVertex, MGPVERTEX_TID);
AUTO_GEL_REGISTER(MGBVertex, MGBVERTEX_TID);
AUTO_GEL_REGISTER(MGEdge, MGEDGE_TID);
AUTO_GEL_REGISTER(MGFace, MGFACE_TID);
AUTO_GEL_REGISTER(MGLoop, MGLOOP_TID);
AUTO_GEL_REGISTER(MGShell, MGSHELL_TID);
AUTO_GEL_REGISTER(MGStl, MGSTL_TID);

#ifndef _CONSOLE
#include "mgGL/PlaneImage.h"
AUTO_GEL_REGISTER(MGPlaneImage, MGPLANEIMAGE_TID);
#endif
