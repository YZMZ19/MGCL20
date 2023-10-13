/********************************************************************/
/* Copyright (c) 2023 System fugen G.K. and Yuzi Mizuno             */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Tolerance.h"
#include "mg/Box.h"
#include "mg/BPointSeq.h"
#include "mg/Curve.h"
#include "mg/Ellipse.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/Straight.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// Implementation of curve offset.


///曲線をオフセットするのに十分分割したノットベクトルを返却する
///オフセット量曲線も考慮に入れ、分割数の多い方にあわせている。
MGKnotVector offset_make_knotvector(const MGCurve* lb1, const MGLBRep* lb2){
	//元となるノットベクトルを生成する
	const MGKnotVector& t1 = lb1->knot_vector();
	MGKnotVector tNew(t1, 4);	//オーダー4のノットベクトルに作り替える
	double error = lb1->param_error();
	for (int i = t1.order() - 1; i < t1.bdim(); i++) {
		double spara = t1(i), epara = t1(i + 1);	//スパンの始終終点
		if (epara - spara < error)
			continue;	//マルチノットのときの処理(RLBRepのみ)

		//1スパンの分割数を決定する
		MGInterval interval(spara, epara);
		int ndiv = lb1->divideNum(interval);//オフセット曲線の分割数
		int tmp_ndiv = lb2 ?  lb2->divideNum(interval) : -1;//オフセット量曲線の分割数
		if (tmp_ndiv > ndiv)
			ndiv = tmp_ndiv;//分割数の多い方を用いる

		int id=tNew.locate(spara)+1;
		int ndivm1 = ndiv - 1;
		tNew.reserveInsertArea(ndivm1, id);

		double delta = (epara - spara) / ndiv;
		double tpara = spara + delta;
		for (int j = 0; j < ndivm1; j++, tpara += delta)
			tNew(id + j) = tpara;
	}
	return tNew;
}

///曲線 crvを折れ(C0 continuous)で分割する
size_t divide_multi_ofs(
	const MGCurve& crv,
	std::vector<UniqueCurve>& crv_list ///divided curves are output.
) {
	const MGKnotVector& t = crv.knot_vector();
	int k = t.order();
	int km1 = k - 1;
	int	start_index = km1, index = 0, multi = 1, vbdim = t.bdim();
	do {
		if (km1 <= 1)
			index = start_index + 1;
		else
			multi = t.locate_multi(start_index, km1, index);
		crv_list.emplace_back(crv.part(t(start_index), t(index)));
		start_index = index + multi - 1;
	} while (index < vbdim);	//多重度が見つからなかったら終わり
	return crv_list.size();
}

//オフセット量固定のC1連続曲線のオフセット
//Curve offset funciont when offset value is constant.
//original curve is supposed to be C1 continuous everywhere.
std::unique_ptr<MGLBRep> C1CurveConstantOffset(
	const MGCurve& original,
	double  ofs_value,//固定オフセット量
	bool pricipalNormal//true if offset direction is toward principal normal,
						//false if to binormal.
){
	MGKnotVector t = offset_make_knotvector(&original, nullptr);//十分分割したノットベクトルを求める
	MGNDDArray tau; tau.buildByKnotVector(t);

	//制御点を生成する
	MGVector T, N, B, Nold, Bold;
	double curvature, torsion;
	original.Frenet_frame(tau(0), T, Nold, Bold, curvature, torsion);
		//save the 1st N, or B data.
	MGVector& dirOld = pricipalNormal ? Nold : Bold;
	MGVector& dir = pricipalNormal ? N : B;

	int n = tau.length();
	MGBPointSeq bp1(n, original.sdim());
	for (int i = 0; i < n; i++) {
		double taui = tau(i);
		original.Frenet_frame(taui, T, N, B, curvature, torsion);
		if (dir % dirOld < 0.)
			dir.negate();
		MGPosition pos = original.eval(taui) + dir * ofs_value;
		bp1.store_at(i, pos);
		dirOld = dir;
	}
	//ノットベクトルを求めてオフセット曲線を生成する(精度十分の曲線を生成する)
	std::unique_ptr<MGLBRep> newLB(new MGLBRep);
	if(newLB->buildByInterpolationDataPoints(tau, bp1))
		return nullptr;//if error detected.
	newLB->remove_knot();
	return newLB;
}

//オフセット量可変のC1連続曲線のオフセット
//Curve offset funciont when offset value is variable.
//original curve is supposed to be C1 continuous everywhere.
UniqueLBRep C1CurveVariableOffset(
	const MGCurve& crv,
	const MGLBRep& ofs_value_lb,	//Variable offset value input(LBRep of space dim =1)
	bool pricipalNormal//true if offset direction is principal normal,
					//false if binormal.
){
	MGKnotVector knotVector = offset_make_knotvector(&crv, &ofs_value_lb);	//十分分割したノットベクトルを求める
	MGNDDArray tau;
	tau.buildByKnotVector(knotVector);

	MGVector T, N, B, Nold, Bold;
	double curvature, torsion;
	crv.Frenet_frame(tau(0), T, Nold, Bold, curvature, torsion);
	MGVector& dirOld = pricipalNormal ? Nold : Bold;
	MGVector& dir = pricipalNormal ? N : B;

	//Generate data of interpolation at tau.
	int n = tau.length();
	MGBPointSeq bp1(n, crv.sdim());
	for (int i = 0; i < n; i++) {
		double taui = tau(i);
		double delta = ofs_value_lb.eval_position(taui).ref(0);//offset value at tau i.
		crv.Frenet_frame(taui, T, N, B, curvature, torsion);
		if (dir % dirOld < 0.)
			dir.negate();

		MGPosition pos = crv.eval(taui) + dir * delta;
		bp1.store_at(i, pos);
		dirOld = dir;
	}

	//ノットベクトルを求めてオフセット曲線を生成する
	UniqueLBRep offsetCrv(new MGLBRep);
	offsetCrv->setKnotVector(std::move(knotVector));
	if(offsetCrv->buildByInterpolationWithKTV(tau, bp1))	//精度十分の曲線を生成する
		return nullptr;//When error detected in buildByInterpolationWithKTV().
	offsetCrv->remove_knot();
	return offsetCrv;
}

///Offset of costant deviation from this curve that is C1 cotinuous everywhere.
///The offset value ofs_value must be less than radius of curvature.
///The offset direction is Normal(obtained by Frenet_frame()) at each point.
UniqueLBRep MGCurve::offsetC1(
	double ofs_value, ///offset valuem may be negative.
	bool principalNormal/// true: Offset direction is to principal normal
						/// false: to binormal
) const {
	return C1CurveConstantOffset(*this, ofs_value, principalNormal);
}

//一定オフセット関数
//オフセット方向は、法線方向から見て入力曲線の進行方向左側を正とする。
//法線ベクトルがヌルの場合、始点において曲率中心方向を正とする。
//ただし、曲率中心へ曲率半径以上のオフセットは行わない。トレランスはline_zero()を使用している。
//戻り値は、オフセット曲線リストが返却される。
std::vector<UniqueCurve> MGCurve::offset(
	double ofs_value,
	bool principalNormal/// true: Offset direction is to principal normal
						/// false: to binormal
)const{
	//曲線を折れ(C0 continuity)で分割する
	std::vector<UniqueCurve> crv_list, offseted;
	divide_multi_ofs(*this,crv_list);
	for (auto& crv : crv_list) {
		UniqueCurve crvOffset = C1CurveConstantOffset(*crv, ofs_value, principalNormal);
		if(crvOffset)
			offseted.push_back(std::move(crvOffset));
	}	
	return offseted;
}

//可変オフセット関数
//オフセット量は空間次元1の線B表現で与えられる。
//オフセット方向は、法線方向から見て入力曲線の進行方向左側を正とする。
//法線ベクトルがヌルの場合、始点において曲率中心方向を正とする。
//ただし、曲率中心へ曲率半径以上のオフセットは行わない。トレランスはline_zero()を使用している。
//戻り値は、オフセット曲線リストが返却される。
std::vector<UniqueCurve> MGCurve::offset(
	const MGLBRep& ofs_value_lb,	//空間次元１の線B表現で示したオフセット量
	bool principalNormal/// true: Offset direction is to principal normal
						/// false: to binormal
)const{
	//曲線を折れ（multiple knot）で分割する
	std::vector<UniqueCurve> crv_list;
	divide_multi_ofs(*this, crv_list);

	std::vector<UniqueCurve> tmp_list;
	for (auto& crv : crv_list) {
		UniqueCurve crvOffset = C1CurveVariableOffset(*crv, ofs_value_lb, principalNormal);
		if (crvOffset)
			tmp_list.push_back(std::move(crvOffset));
	}

	return tmp_list;
}

/// Offset of constant deviation from this curve.
/// The offset value must be less than radius of curvature.
/// When this curve is not C1 continuous, this is divided into C1 curves,
/// and more than one offset curves are obtained.
/// line_zero() is used to approximate curves of the offset.
std::vector<UniqueCurve> MGRLBRep::offset(
	double ofs_value,
	bool principalNormal/// true: Offset direction is to principal normal
						/// false: to binormal
) const {
	std::vector<UniqueCurve> curves;
	curves.emplace_back(C1CurveConstantOffset(*this, ofs_value, principalNormal).release());
	return curves;
}

/// Offset of variable deviation from this curve.
/// When this curve is not C1 continuous, divided into C1 curves,
/// and more than one offset curves are obtained.
/// The direction of offset is toward the principal normal,
/// or to the direction to center of curvature.
/// line_zero() is used approximate the offset curve.
std::vector<UniqueCurve> MGRLBRep::offset(
	const MGLBRep& ofs_value_lb,			///<空間次元１の線B表現で示したオフセット量
	bool principalNormal/// true: Offset direction is to principal normal
						/// false: to binormal
) const {
	std::vector<UniqueCurve> curves;
	curves.emplace_back(C1CurveVariableOffset(*this, ofs_value_lb, principalNormal).release());
	return curves;

}
