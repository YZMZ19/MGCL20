/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/CCisects.h"
#include "mg/CompositeCurve.h"
#include "mg/LBRep.h"
#include "mg/Plane.h"
#include "mg/SurfCurve.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGCompositeCurve Class.
//MGCompositeCurve is a composite of other leaf curves.
//Assumedly they are connected as C0 continuity. However, MGCompositeCurve
//does not check their continuity, but only put first or last as the user says
// (in connect_to_end or in connect_to_start), except two curves connecting
//are both MGLBRep or MGRLBRep. When the two connecting curves are both
//MGLBRep or MGRLBRep and they have more than C0 continuity,
//they are changed to one MGLBRep orMGRLBRep representation.
//Parameter ranges of the member curves are always continuous, from param_s() of the
//1st curve to param_e() of the last form MGCompositeCurve's paramter range.

//関数名：
//		common
//目的：
//		与えられた曲線と自身の交点もしくは共通部分があるかどうか調べる。
//引数：
//		const MGCurve&			crv2,		(I/ )	与えられる曲線
//		std::vector<double>&	vec_param	( /O)	共通部分のパラメータ範囲
//		MGCCisects&			isect		( /O)	交点
//				 4nの配列で、t(4*i+0),t(4*i+1)が自身のパラメータ範囲(t(4*i+0) < t(4*i+1))、
//							 t(4*i+2),t(4*i+3)が与曲線のパラメータ範囲(f(t(4*i+0))=f(t(4*i+2))
//戻り値：
//		3:交点も共通部分も求まった
//		2:交点のみが求まった
//		1:共通部分のみが求まった
//		0:交点も共通部分もなかった
//		-1:共通エッジの収束計算エラー
//		-2:共通エッジが４個以上求まった(のっていないと見なす)
//追記：
//	曲線が共通かどうかの誤差にはline_zero()、をパラメータ範囲の収束計算の
//	誤差には、パラメータ範囲*rc_zero()を使用した
int MGCompositeCurve::common(
	const MGCurve& crv2,
	std::vector<double>& vecComSpan,
	MGCCisects& isect
)const{
	int obtainedIsect=0, obtainedCommon=0;
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		std::vector<double> vecComSpan2;
		MGCCisects isect2;
		int obtained2=(**i).common(crv2,vecComSpan2,isect2);
		if(obtained2<0)
			return obtained2;
		size_t ncom=vecComSpan2.size();
		if(ncom){
			obtainedCommon=1;
			for(size_t j=0; j<ncom; j++)
				vecComSpan.push_back(vecComSpan2[j]);
		}
		size_t nisect=isect2.size();
		if(nisect){
			obtainedIsect=2;
			isect.append(std::move(isect2));
		}
	}
	return obtainedCommon+obtainedIsect;
}

//関数名：
//		common
//目的：
//		与えられた曲線と自身の共通部分があるかどうか調べる。
//引数：
//		const MGCurve&			crv2,		(I/ )	与えられる曲線
//		std::vector<double>&	vec_param	( /O)	共通部分のパラメータ範囲
//				 4nの配列で、t(4*i+0),t(4*i+1)が自身のパラメータ範囲(t(4*i+0) < t(4*i+1))、
//							 t(4*i+2),t(4*i+3)が与曲線のパラメータ範囲(f(t(4*i+0))=f(t(4*i+2))
//戻り値：
//		共通部分の数:	共通部分が求まった
//		0:				共通部分がなかった
//		-1:				共通エッジの収束計算エラー
//		-2:				共通エッジが４個以上求まった(のっていないと見なす)
//追記：
//	曲線が共通かどうかの誤差にはline_zero()を、パラメータ範囲の収束計算の誤差には、
//  パラメータ範囲*rc_zero()を使用した
int MGCompositeCurve::common(
	const MGCurve& crv2,
	std::vector<double>& vecComSpan
)const{
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		std::vector<double> vecComSpan2;
		int obtained2=(**i).common(crv2,vecComSpan2);
		if(obtained2<0)
			return obtained2;
		size_t ncom=vecComSpan2.size();
		for(size_t j=0; j<ncom; j++)
			vecComSpan.push_back(vecComSpan2[j]);

	}
	return (int)(vecComSpan.size()/4);
}

/// Offset of constant deviation from this curve.
/// The offset value must be less than radius of curvature.
/// When this curve is not C1 continuous, this is divided into C1 curves,
/// and more than one offset curves are obtained.
/// line_zero() is used to approximate curves of the offset.
std::vector<UniqueCurve> MGCompositeCurve::offset(
	double ofs_value, //オフセット量
	bool principalNormal/// true: Offset direction is to principal normal
						/// false: to binormal

)const{
	std::vector<UniqueCurve> offsetCrvs;
	for (auto crvi : m_composite) {
		for(auto& crvj: crvi->offset(ofs_value, principalNormal))
			offsetCrvs.emplace_back(crvj.release());
	}
	return offsetCrvs;
}

/// Offset of variable deviation from this curve.
/// When this curve is not C1 continuous, divided into C1 curves,
/// and more than one offset curves are obtained.
/// The direction of offset is toward the principal normal,
/// or to the direction to center of curvature.
/// line_zero() is used approximate the offset curve.
std::vector<UniqueCurve> MGCompositeCurve::offset(
	const MGLBRep& ofs_value_lb,	//空間次元１の線B表現で示したオフセット量
	bool principalNormal/// true: Offset direction is to principal normal
						/// false: to binormal
)const{
	std::vector<UniqueCurve> offsetCrvs;
	for (auto crvi : m_composite) {
		for (auto& crvij : crvi->offset(ofs_value_lb, principalNormal))
			offsetCrvs.emplace_back(crvij.release());
	}
	return offsetCrvs;
}

//Return sweep surface from crv
//Returned is a newed MGSurface, must be deleted.
MGSurface* MGCompositeCurve::sweep(
	const MGUnit_vector& uvec,			//Sweep Direction.
	double start_dist,					//distance to start edge.
	double end_dist			//distance to end edge.
)const{
	return new MGPlane();
}

//Obtain the projected curve of this onto the surface srf.
//The direction of the projection is along the vector vec if the vec is not
//NULL, and normal to the surface if the vec is NULL.
//Output of 'project' is two kind of curves:
//one is general world coordinate curves('vec_crv'), and the other is (u,v) curves of
//the parameter space of the surfaces(vec_crv_uv).
//vec_crv_uv.size() is equal to vec_crv.size(). Let the size be n, then
// (vec_crv_uv[i], vec_crv[i]) is one pair for 0<=i<n.
size_t MGCompositeCurve::project_onto_surface(
	const MGFSurface& srf,
	std::vector<UniqueCurve>& vec_crv_uv,
		//Projected curve(surface parameter (u,v) representation) will be appended.
	std::vector<UniqueCurve>& vec_crv,
		//Projected curve(world coordinate(x,y,z) representation) will be appended.
	const MGVector& vec
)const{
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++)
		srf.project(**i,vec_crv_uv,vec_crv,vec);
	return vec_crv.size();
}

///Approximate this curve by a polyline and output to lb2.
///The tolerance of the approximation is error.
void MGCompositeCurve::polygonize(
	double error,	///<tolerance allowed for the approximation
	MGLBRep& lb2	///<Obtained polyline will be output as an MGLBRep of order2.
)const{
	const_iterator i=begin(), ie=end();
	if(i!=ie)
		(**i++).polygonize(error,lb2);
	for(; i!=ie; i++){
		MGLBRep lbi;
		(**i).polygonize(error,lbi);
		lb2.connect(0,2,lbi);
	}
}

///Approximate this curve as a MGLBRep curve
///within the tolerance MGTolerance::line_zero().
///When parameter_normalization=0, reparameterization will not done, and
///the evaluation at the same parameter has the same values before and after
///of approximate_as_LBRep.
void MGCompositeCurve::approximate_as_LBRep(
	MGLBRep& lb,		///<Approximated obrep will be set.
	int ordr,		///<new order
	int parameter_normalization,
		//Indicates how the parameter normalization be done:
		//=0: no parameter normalization.
		//=1: normalize to range=(0., 1.);
		//=2: normalize to make the average length of the 1st derivative 
		//    is as equal to 1. as possible.
	bool neglectMulti///<Indicates if multiple knots be kept.
		///< true: multiplicity is removed.
		///< false: multiplicity is kept.
)const{
	const_iterator i=begin(), ie=end();
	if(i!=ie)
		(**i++).approximate_as_LBRep(lb,ordr,parameter_normalization,neglectMulti);
	for(; i!=ie; i++){
		MGLBRep lbi;
		(**i).approximate_as_LBRep(lbi,ordr,parameter_normalization,neglectMulti);
		int which;
		double ratio;
		int cntnty=lb.continuity(lbi,which,ratio);
		lb.connect(cntnty,2,lbi);
	}
}
