/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Position.h"
#include "mg/Curve.h"
#include "mg/FSurface.h"
#include "mg/CSisects.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGCSisects defines linked list of MGCSisect.
// Used to represent Intersection points of two curves.

#define ERROR_EXPAND 10.
// Constructor
MGCSisects::MGCSisects(const MGCurve* crv, const MGFSurface* srf)
:MGisects(crv,dynamic_cast<const MGObject*>(srf)){
	m_errort=m_erroru=m_errorv=0.;
	if(crv) m_errort=crv->param_error()*ERROR_EXPAND;
	if(srf){
		m_erroru=srf->param_error_u()*ERROR_EXPAND;
		m_errorv=srf->param_error_v()*ERROR_EXPAND;
	}
}

///Return the pointer to the curve.
const MGCurve* MGCSisects::curve() const{
	return dynamic_cast<const MGCurve*>(object1()); 
}

///Return the pointer to the surface.
const MGFSurface* MGCSisects::surface()const{
	return dynamic_cast<const MGFSurface*>(object2());
}

// Adds the MGCSisect to the end of the list.
//End points will be prefered.
void MGCSisects::append(MGCSisect&& isect){
	const MGCurve& curve=*static_cast<const MGCurve*>(object1());
	double t0=curve.param_s(), t1=curve.param_e(), ter;
	double uer,ver;

	iterator itr;
	for(itr=begin(); itr!=end(); itr++){
		MGCSisect& i=isectCast<MGCSisect>(itr);
		isect.distance(i,ter,uer,ver);
		if(ter<=m_errort && uer<=m_erroru && ver<=m_errorv){
			double s=i.param_curve(), t=isect.param_curve();
			double slen, tlen;
			if(s<=t){
				slen=s-t0; tlen=t1-t;
			}else{
				slen=t1-s; tlen=t-t0;
			}
			if(slen>tlen)
				itr->reset(new MGCSisect(isect));
			return;
		}
	}
	emplace_back(new MGCSisect(std::move(isect)));
}

// 全てのコンポーネントを指定して交点を生成する
void MGCSisects::append(
		const MGPosition& point,		//intersection point.
		double t,				//Curve's parameter value.
        const MGPosition& uv,	//Surface's parameter values.
		const MGCSRELATION rl	//Curve and Surface relation
		)
{
	append(MGCSisect(point,t,uv,rl));
}

// Adds the MGCSisects to the end of the list.
void MGCSisects::append(MGCSisects&& list){
	assert(object1()==list.object1() && object2()==list.object2());
	m_list.splice(m_list.end(),std::move(list.m_list));
}

std::unique_ptr<MGCSisect> MGCSisects::removeAt(iterator i){
//Remove the MGCSisect and return the MGCSisect. If i is no valid, 
// behavior is undefined.
	return release<MGCSisect>(i);
}

std::unique_ptr<MGCSisect> MGCSisects::removeFirst(){
//Remove the first MGCSisect int the list and return the MGCSisect.
//If i is not valid, behavior is undefined.
	return releaseFront<MGCSisect>();
}

std::unique_ptr<MGCSisect> MGCSisects::removeLast(){
//Remove the first MGCSisect int the list and return the MGCSisect.
//If i is not valid, behavior is undefined.
	return releaseBack<MGCSisect>();
}
