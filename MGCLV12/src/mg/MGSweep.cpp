/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Box.h"
#include "mg/BPointSeq.h"
#include "mg/Position.h"
#include "mg/Ellipse.h"
#include "mg/LBRep.h"
#include "mg/Straight.h"
#include "mg/SBRep.h"
#include "mg/RLBRep.h"
#include "mg/RSBRep.h"
#include "mg/TrimmedCurve.h"
#include "mg/SurfCurve.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//The sweep surface is defined as:
//rail(say c(u)) is the rail and the straight line segments
//from C(u)+start_dist*uvec to C(u)+end_dist*uvec are the generatrix.
//The surface is expressed as: S(u,v)=c(u)+uvec*v,
//for rail.param_s()<=u<=rail.param_e(), start_dist<=v<=end_dist.
void MGSBRep::buildSweep(
	const MGLBRep& rail,		//Sweep(rail) crv.
	const MGUnit_vector& uvec,	//Sweep Direction.
	double start_dist,			//distance to start edge.
	double end_dist			//distance to end edge.
){
	MGVector N = uvec;
	if (start_dist > end_dist) {
		double save = end_dist;
		end_dist = start_dist;
		start_dist = save;
		N *= -1.;
	}

	int n = rail.bdim();
	m_surface_bcoef=MGSPointSeq(n, 2, 3);	//面のB表現
	MGVector SPoint = N * start_dist, EPoint = N * end_dist;

	//面B表現を作成する
	const MGBPointSeq& bpnt1 = rail.line_bcoef();//曲線のB表現
	for(int i=0; i < n; i++){
		MGVector lbi = bpnt1(i);		
		m_surface_bcoef.store_at(i, 0, lbi+ SPoint);//スイープ始点
		m_surface_bcoef.store_at(i, 1, lbi+ EPoint);//スイープ終点
	}

	m_uknot = rail.knot_vector();
	m_vknot = MGKnotVector(2, 2, start_dist, end_dist);

	invalidateBox();
	copy_appearance(rail);
}

//The sweep surface is defined as:
//rail(say c(u)) is the rail and the straight line segments
//from C(u)+start_dist*uvec to C(u)+end_dist*uvec are the generatrix.
//The surface is expressed as: S(u,v)=c(u)+uvec*v,
//for rail.param_s()<=u<=rail.param_e(), start_dist<=v<=end_dist.
void MGSBRep::buildSweep(
	const MGStraight& rail,		//rail crv.
	const MGUnit_vector& uvec,	//Sweep Direction.
	double start_dist,			//distance to start edge.
	double end_dist			//distance to end edge.
){
	MGVector N = uvec;
	if (start_dist > end_dist) {
		double save = end_dist;
		end_dist = start_dist;
		start_dist = save;
		N *= -1.;
	}

	int n = 2;
	m_surface_bcoef=MGSPointSeq(n, 2, 3);	//面のB表現
	MGVector SPoint = N * start_dist, EPoint = N * end_dist;

	//直線の始終点を与えられたベクトル、長さで移動し、面B係数に入力する
	MGPosition railS = rail.start_point();
	m_surface_bcoef.store_at(0, 0, railS+SPoint);
	m_surface_bcoef.store_at(0, 1, railS+EPoint);

	MGPosition railE = rail.end_point();
	m_surface_bcoef.store_at(1, 0, railE+SPoint);
	m_surface_bcoef.store_at(1, 1, railE+EPoint);

	m_uknot = MGKnotVector(2, 2, rail.param_s(), rail.param_e());
	m_vknot = MGKnotVector(2, 2, start_dist, end_dist);
	invalidateBox();
	copy_appearance(rail);
}

//The sweep surface is defined as:
//rail(say c(u)) is the rail and the straight line segments
//from C(u)+start_dist*uvec to C(u)+end_dist*uvec are the generatrix.
//The surface is expressed as: S(u,v)=c(u)+uvec*v,
//for rail.param_s()<=u<=rail.param_e(), start_dist<=v<=end_dist.
void MGRSBRep::buildSweep(
	const MGRLBRep& rail,		//Sweep crv.
	const MGUnit_vector& uvec,	//Sweep Direction.
	double start_dist,			//distance to start edge.
	double end_dist			//distance to end edge.
){
	MGVector N = uvec;
	if (start_dist > end_dist) {
		double save = end_dist;
		end_dist = start_dist;
		start_dist = save;
		N *= -1.;
	}

	int n = rail.bdim();
	const MGBPointSeq& rbpnt1 = rail.homogeneous().line_bcoef();
	MGSPointSeq& spnt1 = surface_bcoef(); spnt1.resize(n, 2, 4);//面のB表現(weight含むからsdim+1)

	//面B表現を作成する
	MGBPointSeq staBpnt = rbpnt1, endBpnt = rbpnt1;
	staBpnt.homogeneous_transform(N * start_dist);
	endBpnt.homogeneous_transform(N * end_dist);

	for(int i=0; i < n; i++){
		spnt1.store_at(i, 0, staBpnt(i));
		spnt1.store_at(i, 1, endBpnt(i));
	}
	knot_vector_u() = rail.knot_vector();
	knot_vector_v() = MGKnotVector(2, 2, start_dist, end_dist);

	invalidateBox();
	copy_appearance(rail);
}

//Return sweep surface from crv
//Returned is a newed MGSurface, must be deleted.
//The sweep surface is defined as:
//This curve(say c(t)) is the rail and the straight line segments from
//C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGSurface* MGLBRep::sweep(
	const MGUnit_vector& uvec,		//Sweep Direction.
	double start_dist,				//distance to start edge.
	double end_dist) const			//distance to end edge.
{
	auto srf = new MGSBRep;
	srf->buildSweep(*this, uvec, start_dist, end_dist);
	return srf;
}

//Return sweep surface from crv
//Returned is a newed MGSurface, must be deleted.
//The sweep surface is defined as:
//This straight(say c(t)) is the rail and the straight line segments from
//C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGSurface* MGStraight::sweep(
	const MGUnit_vector& uvec,		//Sweep Direction.
	double start_dist,				//distance to start edge.
	double end_dist) const			//distance to end edge.
{
	auto srf = new MGSBRep;
	srf->buildSweep(*this, uvec, start_dist, end_dist);
	return srf;
}

//Return sweep surface from crv
//Returned is a newed MGSurface, must be deleted.
//The sweep surface is defined as:
//This curve(say c(t)) is the rail and the straight line segments from
//C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGSurface* MGRLBRep::sweep(
	const MGUnit_vector& uvec,		//Sweep Direction.
	double start_dist,				//distance to start edge.
	double end_dist) const			//distance to end edge.
{
	auto srf = new MGRSBRep;
	srf->buildSweep(*this, uvec, start_dist, end_dist);
	return srf;
}

//Return sweep surface from crv
//Returned is a newed MGSurface, must be deleted.
//The sweep surface is defined as:
//This curve(say c(t)) is the rail and the straight line segments from
//C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGSurface* MGEllipse::sweep(
	const MGUnit_vector& uvec,		//Sweep Direction.
	double start_dist,				//distance to start edge.
	double end_dist) const			//distance to end edge.
{
	auto srf = new MGRSBRep;
	srf->buildSweep((MGRLBRep)*this, uvec, start_dist, end_dist);
	return srf;
}

//Return sweep surface from crv
//Returned is a newed MGSurface, must be deleted.
//The sweep surface is defined as:
//This curve(say c(t)) is the rail and the straight line segments from
//C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGSurface* MGTrimmedCurve::sweep(
	const MGUnit_vector& uvec,		//Sweep Direction.
	double start_dist,				//distance to start edge.
	double end_dist) const			//distance to end edge.
{
	MGCurve *tempCrv = clone();
	MGSurface *rtnSrf = tempCrv->sweep(uvec, start_dist, end_dist);
	delete tempCrv;
	return rtnSrf;
}

//Return sweep surface from crv
//Returned is a newed MGSurface, must be deleted.
//The sweep surface is defined as:
//This curve(say c(t)) is the rail and the straight line segments from
//C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGSurface* MGSurfCurve::sweep(
	const MGUnit_vector& uvec,		//Sweep Direction.
	double start_dist,				//distance to start edge.
	double end_dist) const			//distance to end edge.
{
	MGLBRep tempCrv(*this);
	MGSurface *rtnSrf = tempCrv.sweep(uvec, start_dist, end_dist);
	return rtnSrf;
}
