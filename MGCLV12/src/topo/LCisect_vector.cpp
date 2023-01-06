/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Box.h"
#include "topo/LCisect_vector.h"
#include "topo/Loop.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGLCisect_vector defines linked list of MGLCisect.
// Used to represent Intersection points of Loop and curve.

// Constructor
MGLCisect_vector::MGLCisect_vector()
:m_loop(0){;}

MGLCisect_vector::MGLCisect_vector
(const MGLoop& loop):m_loop(&loop)
{
	double err=loop.error();
	m_error_square=err*err*9.;
}

//Copy Constructor.

// Destructor.

// Operator overload.

//Assignment.

// Member Function.

void MGLCisect_vector::append(
	const MGLCisect& lcis)
{
	if(!m_loop) m_loop=lcis.loop();

// Adds the MGLCisect to the end of the list.
	iterator itr=begin(), itrend=end();
	for(; itr!=itrend; itr++){
		if((*itr).distance_square(lcis)<=m_error_square) return;
	}

	push_back(lcis);
}

void MGLCisect_vector::append(
	const MGLEPoint& lp,		//loop's parameter with edge id.
	double t,				//Curve's parameter value.
	const MGPosition& uv)	//Face's parameter value(u,v) data.
{
	append(MGLCisect(lp,t,uv));
}

void MGLCisect_vector::append(const MGLCisect_vector& list){
// Adds the MGLCisect_vector to the end of the list.
	const_iterator i;
	for(i=list.begin(); i!=list.end(); i++) append(*i);
}

//Debug Function
std::ostream& operator<< (std::ostream& out, const MGLCisect_vector& list){
	out<<"MGLCisect_vector::m_loop="<<(list.m_loop)
		<<" ,m_error_square="<<list.m_error_square;
	int n=list.entries();
	out<<", number of isect="<<n<<std::endl;
	MGLCisect_vector::const_iterator itr; int i=0;
	for(itr=list.begin(); itr!=list.end(); itr++)
		out<<i++<<":"<<(*itr);
	return out;
}

//Update MGLEPoint in this LCisect_vector.
//This is to update MGLEPoints in this obtained before MGLoop::make_vertex
//and the loop is updated by MGLoop::make_vertex.
//Generally speaking, when make_vertex is invoked after MGLCisect_vector is obtaind,
//the MGLEPoints in the MGLCisect_vector do not contain correct values, since
//a new edge is inserted int the MGComplex's cell vector.
void MGLCisect_vector::update_lepoint(
	const MGLEPoint& lep	//MGLEPoint used for MGLoop::make_vertex.
){
	//
	MGComplex::const_iterator lep_itr=lep.iterator();
	MGComplex::const_iterator lep_itr_new=lep_itr; lep_itr_new++;
	double t_vertex=lep.param();
	for(iterator i=begin(); i!=end(); i++){
		MGLCisect& lci=*i;
		MGLEPoint lepi=lci.lp();
		if(lepi.iterator()!=lep_itr)
			continue;
		double t=lepi.param();
		if(t_vertex<t){
			lci.set_lepoint(MGLEPoint(lep_itr_new,t));
		}
	}
}
