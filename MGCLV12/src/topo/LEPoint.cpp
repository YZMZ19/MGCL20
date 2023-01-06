/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Curve.h"
#include "mg/SurfCurve.h"
#include "topo/LEPoint.h"
#include "topo/LCisect.h"
#include "topo/Edge.h"
#include "topo/Loop.h"
#include "topo/Face.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Define MGLEPoint Class.
//MGLEPoint is to represent Loop's point. This is represented as
//(l, i, t), where l is loop pointer, i is the edge's iterator in the loop, 
//and t is the parameter value of the curve of the pcell(edge).

///////Constructor////////

//Construct from the loop, edge number, and the edge's parameter value.
MGLEPoint::MGLEPoint(
	const MGLoop& lp,
	int i,			//Edge number.
	double t)			//Parameter value of the edge i.
:m_t(t), m_i(lp.pcellIterator(i)){;}

//Conversion constructor
MGLEPoint::MGLEPoint(const MGLCisect& lci)
:m_i(lci.lp().m_i), m_t(lci.lp().m_t){;}

//convert loop1's MGLEPoint lep to loop2's MGLEPoint. loop2 must be a copy of loop1.
MGLEPoint::MGLEPoint(
	const MGLoop& loop1,	//Loop's edge iterator.
	const MGLEPoint& lep,	//Parameter value of i-th pcell curve.
	const MGLoop& loop2){
	*this=MGLEPoint(loop2,lep.edge_num(),lep.param());
}

///////Operator oveload///////

bool MGLEPoint::operator< (const MGLEPoint& lp)const{
	const MGLoop* lp1=loop();
	const MGLoop* lp2=lp.loop();
	if(lp1 != lp2) return (*lp1)<(*lp2);
	if(*m_i==*(lp.m_i)) return (m_t<lp.m_t);
	return (edge_num()<lp.edge_num());
}
bool MGLEPoint::operator> (const MGLEPoint& lp)const{
	return lp<*this;
}
bool MGLEPoint::operator<= (const MGLEPoint& lp)const{
	return !(*this>lp);
}
bool MGLEPoint::operator>= (const MGLEPoint& lp)const{
	return !(*this<lp);
}
bool MGLEPoint::operator== (const MGLEPoint& lp)const{
	if(*m_i==*(lp.m_i))
		return (fabs(m_t-lp.m_t)<=(edge()->loop())->error());
	return false;
}

///////Member function///////

//Obtain the edge pointer.
const MGEdge* MGLEPoint::edge() const{
	return dynamic_cast<const MGEdge*>(m_i->get());
}
MGEdge* MGLEPoint::edge_to_update() const{
	return const_cast<MGEdge*>(dynamic_cast<const MGEdge*>(m_i->get()));
}

//return loop's edge number.
int MGLEPoint::edge_num()const{
	int i=0;
	const MGComplex* lp=(*m_i)->parent_complex();
	MGComplex::const_iterator itr=lp->pcell_begin(),	itre=lp->pcell_end();
	while(itr!=itre && itr!=m_i){
		i++; itr++;
	}
	return i;
}
bool MGLEPoint::equal_edge(const MGLEPoint& le2) const{
	return m_i->get()==le2.m_i->get(); 
}

//Test if two LEPoints are of the same position.
bool MGLEPoint::equal_position(const MGLEPoint& le2) const{
	return (eval()-le2.eval()).len()<=loop()->error();
}

//Test if this is the same position as P.
bool MGLEPoint::equal_position(const MGPosition& P) const{
	return (eval()-P).len()<=loop()->error();
}

//Evaluation of the loop at the LEPoint.
//When nderi=0, get the positional data(a parameter (u,v) of the surface)
//at the point.
MGVector MGLEPoint::eval(int nderi)const{
	const MGCurve* crv=edge()->base_curve();
	if(crv) return crv->eval(param(), nderi);
	else return MGVector();
}

//Evaluation of the star curves of the loop at the point t.
//When nderi=0, get a position of the surface at the boundary point.
//The star curve is SurfCurve(face's surface, loop's curve).
//(The star curve has the same world coordinate with the binder curve's, but
//their direction may be opposite. The star curve has always the same direction
//as the loop.)
MGVector MGLEPoint::eval_star(int nderi)const{
	assert(edge()->face()->surface());
	return edge()->eval_star(param(),nderi);
}

//Compute a vector at the MGLEPoint point that goes inside the face
//and perpendicular to the boundary loop.
MGUnit_vector MGLEPoint::inner_vector()const{
	return MGVector(0.,0.,1.)*eval(1);
}

///test if this is the end point of the edge.
bool MGLEPoint::is_Edge_end_point()const{
	return edge()->is_end_point(m_t);
}

///test if this is the start point of the edge.
bool MGLEPoint::is_Edge_start_point()const{
	return edge()->is_start_point(m_t);
}

//test if this is the end point of the loop.
bool MGLEPoint::is_end_point()const{
	return equal_position(loop()->end_point());
}

//test if this is the start point of the loop.
bool MGLEPoint::is_start_point()const{
	return equal_position(loop()->start_point());
}

//Get loop pointer.
const MGLoop* MGLEPoint::loop()const{
	return edge()->loop();
}

//Get the start point of the next edge in the loop sequence.
MGLEPoint MGLEPoint::end_of_pre_edge()const{
	const MGEdge* edg2=edge()->pre_edge();
	MGComplex::const_iterator ei2=edg2->edge_iterator();
	return MGLEPoint(ei2, edg2->param_e());
}

//Get the start point of the next edge in the loop sequence.
MGLEPoint MGLEPoint::start_of_next_edge()const{
	const MGEdge* edg2=edge()->aft_edge();
	MGComplex::const_iterator ei2=edg2->edge_iterator();
	return MGLEPoint(ei2, edg2->param_s());
}

//Set iteraror of the edge and parameter value of the edge.
void MGLEPoint::set(
	MGComplex::const_iterator i,		//Iterator of the edge.
	double t
){							//parameter of the edge.
	m_i=i; m_t=t;
}

//Debug Function
std::ostream& operator<< (std::ostream& ostrm, const MGLEPoint& lp){
	ostrm<<"MGLEPoint::m_i=";
	ostrm<<*(lp.m_i)<<", m_t="<<lp.m_t;
	return ostrm;
}
