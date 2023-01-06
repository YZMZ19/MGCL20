/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "cskernel/Blgi2d.h"
#include "mg/Tolerance.h"
#include "mg/Curve.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/SPointSeq.h"
#include "mg/Surface.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//同じカーブの種類で同じノットベクトルを持っているスプライン曲線列かどうか調べる
bool isSameKnotVector(const std::vector<const MGCurve*>& pvCurves){
	std::vector<const MGCurve*>::const_iterator citer = pvCurves.begin();
	const MGCurve& curve0=**citer;
	MGCURVE_TYPE CrvId =curve0.type();
	//スプラインカーブ以外が含まれていたらfalseを返す
	if(CrvId != MGCURVE_TYPE::MGCURVE_SPLINE && CrvId != MGCURVE_TYPE::MGCURVE_RSPLINE)
		return false;

	std::vector<const MGCurve*>::const_iterator cend = pvCurves.end();
	const MGKnotVector& t0 = curve0.knot_vector();//1st knot vector.
	for(citer++; citer!=cend; citer++){
		const MGCurve& curvei=**citer;
		//カーブの種類が違っていたり、ノットベクトルが違っていたらfalseを返す
		if(curvei.type()!=CrvId || curvei.knot_vector()!=t0)
			return false;
	}
	return true;
}

//リブ曲線列から面を作成する
//リブ曲線の全てのノットが同じスプラインの時MGSBRepかMGRSBRepが返却される
//それ以外の場合は、曲線をLBRepで再構成して面を作成するのでMGSBRepが返却される
//作成する面のノットベクトルはリブ曲線の向きをu,リブ列方向をvとする
//Let v0=start parameter value, v1=terminate parameter value along v, then
//v=v0 const parameter line is curves[0], and v=v1 const parameter line is
//curves[n-1], where n=curves.size(). n must be greater or equal to 2.
//When n==2, the surface is a ruled surface(that is, order_u() is 2).
std::unique_ptr<MGSurface> MGCL::createSurfaceFromRibs(
	const std::vector<const MGCurve*>& curves,//リブ曲線列
	bool direction_adjustment	//=true, curves[.] direction are adjusted to line
								//to the same direction.
){
	//リブ列のノットベクトルがそろえられているか調べる(Bスプライン以外はfalseが返る)
	if(isSameKnotVector(curves)){
		const MGCurve& crv0 = *(curves.front());
		if(crv0.type() == MGCURVE_TYPE::MGCURVE_RSPLINE){//Rational or non Rationalを調べる
			size_t nRibNum = curves.size();
			std::vector<const MGRLBRep*> vecPtrRLBReps(nRibNum);
			for(size_t i=0; i<nRibNum; i++)
				vecPtrRLBReps[i]=static_cast<const MGRLBRep*>(curves[i]);
			std::unique_ptr<MGRSBRep> rsb(new MGRSBRep);
			rsb->buildByRibRLBRep(vecPtrRLBReps, direction_adjustment);
			return rsb;
		}
	}

	//ノットベクトルが異なったリブの場合曲線をリビルドする
	std::unique_ptr<MGSBRep> srf(new MGSBRep);
	srf->buildByRibCurves(curves, direction_adjustment);
	return srf;
}

//リブ曲線列の両端点、中点の距離を平均してdata points tauを作成する
void createRibDataPoints(const std::vector<UniqueLBRep>& curves, MGNDDArray& tau){
	int nBdim = (int)curves.size();
	tau.resize(nBdim);
	double taui=0.0;
	int nBdimM1=nBdim-1; 
	double error=MGTolerance::wc_zero()*10.;
	for(int i=0; i<nBdimM1; i++){
		tau[i]=taui;
		const MGCurve& curvei=*(curves[i]); const MGCurve& curveiP1=*(curves[i+1]);
		double dStartDist = curvei.start_point().distance(curveiP1.start_point());
		double dEndDist = curvei.end_point().distance(curveiP1.end_point());
		double dCenterDist = curvei.center().distance(curveiP1.center());
		double dist=(dStartDist+dEndDist+dCenterDist)/3.0;
		if(dist<error)
			dist=error;
		taui += dist;
	}
	tau[nBdimM1]=taui;
}

// Creates a ruled surface.
std::unique_ptr<MGSurface> MGCL::create_ruled_surface(
	const MGCurve& cross1,  // a curve as Edge No.0.
	const MGCurve& cross2,   // another curve as Edge No.2.
	bool direction_adjustment	//=true, curves[.] direction are adjusted to line
								//to the same direction.
){
	std::vector<const MGCurve*> curves(2);
	curves[0]=&cross1; curves[1]=&cross2;
	return MGCL::createSurfaceFromRibs(curves,direction_adjustment);
}

///リブ曲線列から面を作成する
///作成する面のノットベクトルはリブ曲線の向きをu,リブ列方向をvとする
///This constructor only generates MGSBRep even if curves are MGRLBRep of the same
///knot configuration. To avoid this, use createSurfaceFromRibs() that generates
///MGRSBRep when curves are MGRLBRep of the same knot configuration.
///
///Let v0=start parameter value, v1=terminate parameter value along v, then
///v=v0 const parameter line is curves[0], and v=v1 const parameter line is
///curves[n-1], where n=curves.size(). n must be greater or equal to 2.
///When n==2, the surface is a ruled surface(that is, this->order_u() is 2).
///When direction_adjustment=true, curves[i] directions are adjusted to line up
///to the same direciton as the curves[i-1]'s.
///When direction_adjustment=false, no curves of curves[i] are negated and
///the directions of curves[i] are the direcitons of ribs.
void MGSBRep::buildByRibCurves(
	const std::vector<const MGCurve*>& ribCurves,
	bool direction_adjustment//=true, curves[.] direction are adjusted to line
								//to the same direction.
){
	invalidateBox();

	int i=0, nRibNum = (int)ribCurves.size();
	assert(nRibNum>=2);

	std::vector<UniqueCurve> crvs(nRibNum);
	for(size_t i=0; i<nRibNum; i++)
		crvs[i].reset(ribCurves[i]->clone());
	updateRibsDirection(crvs.begin(), crvs.end());//crvs have the same directions.

	std::vector<const MGCurve*> crvs2(nRibNum);
	extractConstPointerVec(crvs.begin(), crvs.end(), crvs2.begin());
	std::vector<UniqueLBRep> lbs=rebuildAsSameKnotVector(crvs2);	//リビルド実行
	m_uknot = lbs[0]->knot_vector();
	
	MGNDDArray vtau;
	createRibDataPoints(lbs,vtau);
	int lenu = m_uknot.bdim(), lenv = vtau.length();
	int nvOrder = 4;//Default order is 4.
	if(lenv < nvOrder)
		nvOrder = lenv;
	m_vknot=MGKnotVector(vtau,nvOrder);

	//面作成準備
	int len=lenu;
	if(lenv > len)
		len = lenv;	//lenはlenu,lenvの大きい方の値が入る
	int lenuv = lenu * lenv;
	int nOrder = m_uknot.order(), nSdim = lbs[0]->sdim();
	if(nOrder<nvOrder)
		nOrder=nvOrder;
	m_surface_bcoef = MGSPointSeq(lenu, lenv, nSdim);
	double *wk = new double[len*2+lenuv*nSdim+len*(2*nOrder-1)];
	double *wk2 =wk+len*2;
	double *q =wk2+lenuv*nSdim;

	//曲線から元となるコントロールポイントを作成する
	for(i=0; i<lenu; i++){
		for(int j=0; j<lenv; j++){
			for(int k=0; k<nSdim; k++){
				wk2[j+lenv*i+lenuv*k] = lbs[j]->coef(i,k);
			}
		}
	}
	int blg4spError = 2;
	for(i=0; i<nSdim; i++){
		blgi2d_(&blg4spError,vtau.data(),wk2+lenuv*i,m_vknot.data(),
			nvOrder,lenv,lenu,lenv,lenu,wk,q,&m_surface_bcoef(0,0,i));
	}
	delete[] wk;
}

//リブ曲線列から面を作成する
//作成する面のノットベクトルはリブ曲線の向きをu,リブ列方向をvとする
//Let v0=start parameter value, v1=terminate parameter value along v, then
//v=v0 const parameter line is curves[0], and v=v1 const parameter line is
//curves[n-1], where n=curves.size(). n must be greater or equal to 2.
//When n==2, the surface is a ruled surface(that is, this->order_u() is 2).
void MGRSBRep::buildByRibRLBRep(
	const std::vector<const MGRLBRep*>& ribCurves,
	bool direction_adjustment//=true, curves[.] direction are adjusted to line
								//to the same direction.
){
	size_t i=0, nRibNum= ribCurves.size();
	std::vector<const MGCurve*> vecPtrRibLBReps(nRibNum);
	for(i=0; i<nRibNum; i++){
		vecPtrRibLBReps[i] = &ribCurves[i]->homogeneous();
	}
	m_surface.buildByRibCurves(vecPtrRibLBReps, direction_adjustment);
}
