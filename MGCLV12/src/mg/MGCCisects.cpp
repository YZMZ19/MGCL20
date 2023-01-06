/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Position.h"
#include "mg/Curve.h"
#include "mg/CCisects.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGCCisects defines linked list of MGCCisect.
// Used to represent intersection points of two curves.

// Constructor
MGCCisects::MGCCisects(const MGCurve* c1, const MGCurve*c2)
: MGisects(c1,c2){
	double er1=0., er2=0.;
	if(c1) er1=c1->param_error();
	if(c2) er1=c2->param_error();
	m_error=(er1*er1+er2*er2)*9.;	//set error.
}

/*
void MGCCisects::append(const MGCCisect& isect){
// Adds the MGCCisect to the end of the list.
	iterator itr; double dif1,dif2;
	for(itr=begin(); itr!=end(); itr++){
		MGCCisect& i=*static_cast<MGCCisect*>(itr->get());
		dif1=i.param1()-isect.param1();
		dif2=i.param2()-isect.param2();
		if((dif1*dif1+dif2*dif2)<=m_error) return;
	}
	append(isect);
}
*/

// 交点の全てのコンポーネントを指定して，交点リストに追加
//Add one intersection point to the list.
void MGCCisects::append(
		const MGPosition& point,	//Intesection point(x,y,)
		double t1,					//parameter value of curve 1.
		double t2,					//parameter value of curve 2.
		const MGCCRELATION r1){
	emplace_back<MGCCisect>(point,t1,t2,r1);
}

///Return the pointer to curve1.
const MGCurve* MGCCisects::curve1() const{
	return dynamic_cast<const MGCurve*>(object1()); 
}

///Return the pointer to curve2.
const MGCurve* MGCCisects::curve2()const{
	return dynamic_cast<const MGCurve*>(object2());
}

void MGCCisects::append(MGCCisects&& list){
// Adds the MGCCisects to the end of the list.
	push_back(std::move(list));
}

/// Adds the MGCCisect to the beginning of the list.
void MGCCisects::prepend(MGCCisect&& isect){
// Adds the MGCCisects to the end of the list.
	m_list.emplace_front(new MGCCisect(std::move(isect)));
}

std::unique_ptr<MGCCisect> MGCCisects::removeFirst(){
//Remove the first MGCCisect int the list and return the MGCCisect.
//If i is not valid, behavior is undefined.
	std::unique_ptr<MGCCisect> isect(static_cast<MGCCisect*>(front().release()));
	m_list.pop_front();
	return isect;
}

std::unique_ptr<MGCCisect> MGCCisects::removeLast(){
//Remove the first MGCCisect int the list and return the MGCCisect.
//If i is not valid, behavior is undefined.
	std::unique_ptr<MGCCisect> isect(static_cast<MGCCisect*>(back().release()));
	m_list.pop_back();
	return isect;
}
