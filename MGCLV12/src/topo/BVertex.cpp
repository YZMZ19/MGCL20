/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Point.h"
#include "topo/BVertex.h"
#include "topo/PVertex.h"
#include "topo/Edge.h"
#include "topo/Face.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGBVertex Class.
//MGBVertex is 0 manifold dimension binder cell, is an point.
//MGBVertex is a binder cell of MGPVertex.

/////// Constructor ///////

//Copy constructor.
MGBVertex::MGBVertex(const MGBVertex& rhs):m_point(rhs.m_point->clone()){
}

//Constructor that holds newPartners as partner members.
MGBVertex::MGBVertex(std::vector<const MGPCell*>& newPartners)
:MGBCell(newPartners){
}

///Construct MGBVertex of only one partner member of pv.
MGBVertex::MGBVertex(const MGPVertex& pv):MGBCell(&pv){ ; };

MGBVertex& MGBVertex::operator=(const MGBVertex& bv){
	if(this==&bv)
		return *this;
	MGBCell::operator=(bv);
	m_point.reset(bv.m_point->clone());
	return *this;
}

///Comparison.
bool MGBVertex::operator==(const MGGel& gel2)const{
	const MGBVertex* bv2 = dynamic_cast<const MGBVertex*>(&gel2);
	return bv2 ? bv2==this:false;
}
bool MGBVertex::operator<(const MGGel& gel2)const{
	const MGBVertex* bv2 = dynamic_cast<const MGBVertex*>(&gel2);
	if(bv2)
		return (*this) < (*bv2);
	return identify_type() < gel2.identify_type();
}

bool MGBVertex::operator<(const MGBVertex& bv)const{
	const MGPoint* pnt1=point();
	if(!pnt1)
		return true;

	const MGPoint* pnt2= bv.point();
	if(!pnt2)
		return false;

	return (*pnt1)<(*pnt2);
}

///////Member Function///////

///Transform the binder.
void MGBVertex::binder_tr(const MGVector& v){//Translation.
	if(m_point)
		*m_point+=v;
}
void MGBVertex::binder_tr(double s){//Scaling.
	if(m_point)
		*m_point*=s;
}
void MGBVertex::binder_tr(const MGMatrix& mat){//Matrix transform.
	if(m_point)
		*m_point*=mat;
}
void MGBVertex::binder_tr(const MGTransf& tr){//Trans transform.
	if(m_point)
		*m_point*=tr;
}

///Get extent geometry, may be null if this does not have extent.
const MGPoint* MGBVertex::extentBC() const{
	return m_point ? m_point.get():nullptr;
}
MGPoint* MGBVertex::extentBC(){
	return m_point ? m_point.get():nullptr;
}


//Make sure that this has an extent expression.
//When this did not have an extent, make the extent from the partner
//member's parameter expression and the star cell.
//This must be a binder cell that has partner members that are
//boundaries. When this is not the case or this had an extent already,
//make_extent does nothing.
void MGBVertex::make_extent() const{
	if(m_point)
		return;

	const MGPVertex* mpartner = partner_member_vertex(0);
	if(mpartner){
		const MGEdge* edge = dynamic_cast<const MGEdge*>(mpartner->star());
		if(edge){
			MGPosition P =
				mpartner->is_start_vertex() ? edge->start_point() : edge->end_point();
			m_point.reset(new MGPoint(P));
		}
	}
}

//Obtain the i-th member partner vertex.
const MGPVertex* MGBVertex::partner_member_vertex(int i)const{
	return static_cast<const MGPVertex*>(partner_member(i));
}

///Release the extent of this binder cell.
MGGeometry* MGBVertex::release_extentBC(){
	return m_point.release();
}

///Set extent of this binder cell.
void MGBVertex::set_extentBC(std::unique_ptr<MGGeometry>&& extent){
	assert(dynamic_cast<MGPoint*>(extent.get()));//must be MGPoint.
	MGPoint* pt = dynamic_cast<MGPoint*>(extent.release());
	m_point.reset(pt);
}
void MGBVertex::set_extent_as_nullBC(){
	m_point.reset();
}
void MGBVertex::set_point(const MGPosition& P){
	m_point.reset(new MGPoint(P));
}

// Output virtual function.
std::ostream& MGBVertex::toString(std::ostream& ostrm) const{
	ostrm<<"<<BV="<<(const MGGel*)this;
	MGBCell::toString(ostrm);
	if(m_point)
		ostrm<<"m_point="<<*m_point;

	ostrm<<"=BV>>";
	return ostrm;
}

//Return extent data(i.e. MGPoint*).
const MGPoint* MGBVertex::point() const{
	return dynamic_cast<const MGPoint*>(m_point.get());
}

//Get the position data of this vertex.
//Returns MGPoint data if this has extent. Otherwise obtain from partner's
//star edge by evaluating the edge's data.
MGPosition MGBVertex::position() const{
	make_extent();//Make sure that this has MGPoint extent.
	const MGPoint* pnt=point();
	return pnt->position();
}

///Clone this BCell, building the new partner membership of partners.
SharedBCell MGBVertex::cloneWithPartnerMembers(std::vector<const MGPCell*>& newPartners)const{
	SharedBVertex newBcell(new MGBVertex(newPartners));
	if(m_point)
		newBcell->set_point(m_point->position());
	return newBcell;
}
