/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Box.h"
#include "mg/Straight.h"
#include "mg/Surface.h"
#include "mg/Plane.h"
#include "mg/SPointSeq.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// Implementation of surface offset.

#define POW_LOW 0.34   //最初にオフセットさせるポイント数を決める係数(速度重視)
#define POW_HIGH 0.50  //最初にオフセットさせるポイント数を決める係数(精度重視)
#define NUM_DIV 25     //1パッチの分割数がこれを越えたときPOW_LOWを使用するようにする

//一定オフセット関数
//オフセット方向は、ノーマル方向を正とする。曲率半径より大きいオフセットは行わない
//戻り値は、オフセット曲面リストが返却される。エラーのとき長さ0の曲面リストが返る。
//トレランスはline_zero()を使用している。
std::vector<UniqueSurface> MGSurface::offset(
	double ofs_value,			//オフセット量
	int& error) const			//エラーコード 0:成功 -2:曲率半径以上のオフセット不可 -3:面生成コンストラクタエラー
{
	//先にu方向の折れの部分で分割する
	std::vector<UniqueSurface> vecSrf, ofs_srfl, null_srf_vec;
	int div = divide_multi_knot(vecSrf);
	for(int i = 0; i < div; i++){
		UniqueSurface pofs_srf = vecSrf[i]->offset_c1(ofs_value, error);	//オフセットを行う
		if(error < 0)return null_srf_vec;
		ofs_srfl.emplace_back(pofs_srf.release());
	}
	return ofs_srfl;
}

//Offset.
//distance is plus value if the direction is toward normal vector of the
//FSurface. Minus if against the normal vector.
//エラーコード 0:成功 -1:曲率半径以上のオフセット不可 -3:面生成コンストラクタエラー
int MGSurface::offset_fs(
	double distance,
	std::vector<UniqueFSurface>& vecOfsFSurface	//Offset MGFSurfaces are appended.
)const{
	int error;
	std::vector<UniqueSurface> srfs=offset(distance,error);
	if(error)
		return error;

	std::move(srfs.begin(), srfs.end(), std::back_inserter(vecOfsFSurface));
	return 0;
}

//C1連続曲面の一定オフセット関数
//オフセット方向は、ノーマル方向を正とする。曲率半径より大きいオフセットは行わない。
//戻り値は、オフセットした曲面のオートポインタである。エラーのときヌルが返る。
//トレランスはline_zero()を使用している。
std::unique_ptr<MGSurface> MGSurface::offset_c1(
	double ofs_value,		//オフセット量
	int& error) const		//エラーコード 0:成功 -1:面におれがある
 					// -2:曲率半径以上のオフセット不可 -3:面生成コンストラクタエラー
{
	//C0連続のときエラー
	std::unique_ptr<MGSurface> pofs_srf;
	const MGKnotVector &knotu = knot_vector_u(), &knotv = knot_vector_v();	//SBRep, RSBRepのみだから
	int startu = knotu.order() - 1, startv = knotv.order() - 1, index = 0;
	if(type() == MGSURFACE_TYPE::MGSURFACE_SPLINE){	//ノンラショナルのときだけおれのチェックを行う
		if(knotu.order() != 2){
			if(knotu.locate_multi(startu, knotu.order() - 1, index)){error = -1; return pofs_srf;}
		}else if(knotu.bdim() > 2){error = -1; return pofs_srf;}
		if(knotv.order() != 2){
			if(knotv.locate_multi(startv, knotv.order() - 1, index)){error = -1; return pofs_srf;}
		}else if(knotv.bdim() > 2){error = -1; return pofs_srf;}
	}

	//オフセットが曲率半径以内かどうか調べる
	if(!offset_check_curva(ofs_value)){error = -2; return pofs_srf;}

	//面を分割するノットベクトルとデータポイントを計算する
	MGKnotVector knotVec_u, knotVec_v;
	MGNDDArray data_point_u, data_point_v;
	offset_calc_knot_vec(knotVec_u, knotVec_v, data_point_u, data_point_v);

	//オフセット面を作成する
	int	ulen = knotVec_u.bdim(), vlen = knotVec_v.bdim();
	MGSPointSeq sb1(ulen, vlen, sdim());
	for(int i = 0; i < ulen; i++){
		for(int j = 0; j < vlen; j++){
			double data_u = data_point_u(i), data_v = data_point_v(j);
			MGPosition pos, ofs_pos;
			pos = eval(data_u, data_v, 0, 0);
			MGVector N(unit_normal(data_u, data_v));
			ofs_pos = pos + (N* ofs_value);
			sb1.store_at(i, j, ofs_pos);
		}
	}
	MGSBRep* ofs_sbrep=new MGSBRep;
	ofs_sbrep->setKnotVector(std::move(knotVec_u), std::move(knotVec_v));
	ofs_sbrep->buildByInterpolationWithKTV(data_point_u, data_point_v, sb1);
	ofs_sbrep->remove_knot();
	ofs_sbrep->copy_appearance(*this);
	pofs_srf.reset(ofs_sbrep);
	error = 0;
	return pofs_srf;
}

//C1連続曲面の一定オフセット関数
//オフセット方向は、ノーマル方向を正とする。曲率半径より大きいオフセットは行わない
//戻り値は、オフセットした曲面のオートポインタである。エラーのときはヌルが返る。
//トレランスはline_zero()を使用している。
std::unique_ptr<MGSurface> MGPlane::offset_c1(
	double ofs_value,				//オフセット量
	int& error) const				//エラーコード 0:成功 -1:面におれがある
									// -2:曲率半径以上のオフセット不可 -3:面生成コンストラクタエラー
{
	MGPlane ofs_plane = *this + (this->normal() * ofs_value);
	std::unique_ptr<MGSurface> pofs_srf(new MGPlane(ofs_plane));
	error = 0;
	return pofs_srf;
}

//C1連続曲面の一定オフセット関数で使用するノットベクトルとデータポイントを求める関数
void MGSurface::offset_calc_knot_vec(
	MGKnotVector& knot_vec_u,			//求まったuノットベクトル
	MGKnotVector& knot_vec_v,			//求まったvノットベクトル
	MGNDDArray& data_point_u,			//求まったuデータポイント
	MGNDDArray& data_point_v) const		//求まったvデータポイント
{
	int num_div = offset_div_num();		//2次微分値に応じて分割数を求める
	int ord_u = order_u(), ord_v = order_v();
	int knotNumU = (bdim_u() - ord_u + 1) * num_div, knotNumV = (bdim_v() - ord_v + 1) * num_div;
	if(ord_u < 4)ord_u = 4; if(ord_v < 4)ord_v = 4;
	knot_vec_u = MGKnotVector(knot_vector_u(), ord_u);
	knot_vec_v = MGKnotVector(knot_vector_v(), ord_v);
	knot_vec_u.change_knot_number(knotNumU);
	knot_vec_v.change_knot_number(knotNumV);
	data_point_u.buildByKnotVector(knot_vec_u);
	data_point_v.buildByKnotVector(knot_vec_v);
}

//曲率半径とオフセット値の判定を行う
//errorのときfalseが返却される
int MGSurface::offset_check_curva(
	double ofs_value) const		//オフセット量
{
	offset_check_curva_one(ofs_value);
	MGSurface *temp_srf = this->copy_surface();
	temp_srf->exchange_uv();

	int rtn = temp_srf->offset_check_curva_one(-ofs_value);	//uv入れ替えたのでnormalが逆になる
	delete temp_srf;
	return rtn;
}

//曲率半径とオフセット値の判定を行う(u=constのパラメータ曲線を使用する)
//errorのときfalseが返却される
int MGSurface::offset_check_curva_one(
	double ofs_value) const		//オフセット量
{
	//u方向スパンを2*orderで分割するパラメータ曲線の曲率半径がオフセット値よりも大きいかどうか調べる
	int i, j, k, l;
	for(i = order_u() - 1; i < bdim_u(); i++){
		double uspan = (knot_u(i + 1) - knot_u(i)) / (2. * order_u());
		int ndiv_u = 2 * order_u();
		if(i == bdim_u() - 1)ndiv_u++;
		for(k = 0; k < ndiv_u; k++){
			double uparam = knot_u(i) + (uspan * k);
			MGCurve* pcrv = parameter_curve(1, uparam);	//u = constのパラメータ曲線
			for(j = order_v() - 1; j < bdim_v(); j++){
				double vspan = (pcrv->knot(j + 1) - pcrv->knot(j)) / (2. * pcrv->order());
				int ndiv_v = 2 * pcrv->order();
				if(j == bdim_v() - 1)ndiv_v++;
				for(l = 0; l < ndiv_v; l++){
					double vparam = pcrv->knot(j) + (vspan * l), curva, torsion, radius;
					MGUnit_vector T, N, B, norm;
					pcrv->Frenet_frame(vparam, T, N, B, curva, torsion);
					norm = unit_normal(uparam, vparam);
					if((N % norm) * ofs_value <= 0 || MGMZero(curva))continue;
					radius = 1. / curva;
					if((radius * radius) < (ofs_value * ofs_value)){delete pcrv; return false;}
				}
			}
			delete pcrv;
		}
	}
	return true;
}

//オフセットするサンプルポイントの1パッチごとの分割数を求める
//全てのパッチ中の分割数で最大の値を返す
int MGSurface::offset_div_num() const{
	int max_div = 0;
	int bdu=bdim_u(), bdv=bdim_v();
	for(int i = order_u() - 1; i < bdu; i++){
		MGInterval itvl_u(knot_u(i), knot_u(i + 1));
		for(int j = order_v() - 1; j < bdv; j++){
			MGInterval itvl_v(knot_v(j), knot_v(j + 1));
			int div = offset_div_num_one(MGBox(itvl_u, itvl_v));
			if(div > max_div)max_div = div;
		}
	}
	return max_div;
}

//オフセットするサンプルポイントの1パッチごとの分割数を求める
//全てのパッチ中の分割数で最大の値を返す
int MGPlane::offset_div_num() const
{
	return 1;
}

//オフセットするサンプルポイントの1パッチごとの分割数を求める
//分割数n = sqrt(1 / tol) * sqrt((M1 + 2 *M2 + M3) / 8)は以上の式で求まる
int MGSurface::offset_div_num_one(
	const MGBox& param_range) const{
	int orderu=order_u(), orderv=order_v();

	//calculate maximum second derivative
	int div_u = 2 * orderu, div_v = 2 * orderv;
	double span_u = fabs(param_range.high().ref(0) - param_range.low().ref(0)) / double(div_u);
	double span_v = fabs(param_range.high().ref(1) - param_range.low().ref(1)) / double(div_v);
	double start_u = param_range.low().ref(0), start_v = param_range.low().ref(1);
	double M1 = 0.0, M2 = 0.0, M3 = 0.0;	//M1: Maximum Second derivative of Suu, M2: Suv, M3: Svv
	for(int i = 0; i <= div_u; i++){
		for(int j = 0; j <= div_v; j++){
			double param_u = start_u + (span_u * i), param_v = start_v + (span_v * j);
			double d1 = eval(param_u, param_v, 2, 0).len() * span_u * span_u;	//Second derivative of Suu
			double d2 = eval(param_u, param_v, 2, 2).len() * span_v * span_v;	//Second derivative of Suv
			double d3 = eval(param_u, param_v, 0, 2).len() * span_u * span_v;	//Second derivative of Svv
			if(d1 > M1)M1 = d1;
			if(d2 > M2)M2 = d2;
			if(d3 > M3)M3 = d3;
		}
	}
	int n = int(pow((1. / MGTolerance::line_zero()), POW_HIGH) * sqrt((M1 + (2. * M2) + M3) / 8.));
	if(n > NUM_DIV)n = int(pow((1. / MGTolerance::line_zero()), POW_LOW) * sqrt((M1 + (2. * M2) + M3) / 8.));
	if(n < 2*orderu) n=2*orderu;
	if(n < 2*orderv) n=2*orderv;
	return n;
}

//uまたはv方向に折れ(マルチノット)があるとき面を分割する
//戻り値は、分割数を返却する
int MGSBRep::divide_multi_knot(
	std::vector<UniqueSurface>& srfl	///Divided objects are appended.
) const {
	//先にu方向の折れの部分で分割する
	std::vector<UniqueSBRep> srf_u;
	int udiv = divide_multi_knot_u(srf_u), num = 0;
	for(int i = 0; i < udiv; i++){
		//次にv方向の折れで分割する
		std::vector<UniqueSBRep> srf_v;
		num += srf_u[i]->divide_multi_knot_v(srf_v);
		std::move(srf_v.begin(), srf_v.end(), std::back_inserter(srfl));
	}
	return num;
}

//u方向に折れ(マルチノット)があるとき面を分割する
//戻り値は、分割数を返却する
int MGSBRep::divide_multi_knot_u(
	std::vector<UniqueSBRep>& srfl) const		//分割した曲面リスト
{
	int start_index = 0, index = 0, count = 0, multi = 0, bdim = 0, order = 0;
	MGInterval intv = param_range().ref(1);
	const MGKnotVector &knot_vector = knot_vector_u();
	start_index = order_u() - 1;
	bdim = bdim_u();
	order = order_u();
	do{		//u,v方向に折れ(C0連続)があるとき面を分割する
		if(order == 2){		//オーダー2のときの処理
			index = start_index + 1; multi = 1;
		}else{
			multi = knot_vector.locate_multi(start_index, order - 1, index);
		}
		MGBox uv_range = MGBox(MGInterval(knot_u(start_index), knot_u(index)), intv);
		MGSBRep* sb=new MGSBRep;
		shrinkToParameters(uv_range, *sb);
		srfl.emplace_back(sb);
		start_index = index + multi - 1;
		count++;
	}while(index != bdim);	//多重度が見つからなかったら終わり
	return count;
}

//v方向に折れ(マルチノット)があるとき面を分割する
//戻り値は、分割数を返却する
int MGSBRep::divide_multi_knot_v(
	std::vector<UniqueSBRep>& srfl///Divided MGRSBRep are appended.
)const{
	int start_index = 0, index = 0, count = 0, multi = 0, bdim = 0, order = 0;
	MGInterval intv = param_range().ref(0);
	const MGKnotVector &knot_vector = knot_vector_v();
	start_index = order_v() - 1;
	bdim = bdim_v();
	order = order_v();
	do{		//u,v方向に折れ(C0連続)があるとき面を分割する
		if(order == 2){		//オーダー2のときの処理
			index = start_index + 1; multi = 1;
		}else{
			multi = knot_vector.locate_multi(start_index, order - 1, index);
		}
		MGBox uv_range = MGBox(intv, MGInterval(knot_v(start_index), knot_v(index)));
		MGSBRep* sb=new MGSBRep;
		shrinkToParameters(uv_range, *sb);
		srfl.emplace_back(sb);
		start_index = index + multi - 1;
		count++;
	}while(index != bdim);	//多重度が見つからなかったら終わり
	return count;
}

//uまたはv方向に折れ(マルチノット)があるとき面を分割する
//戻り値は、分割数を返却する
int MGRSBRep::divide_multi_knot(
	std::vector<UniqueSurface>& srfl	///Divided objects are appended.
) const {
	//先にu方向の折れの部分で分割する
	std::vector<UniqueRSBRep> srf_u;
	int udiv = divide_multi_knot_u(srf_u), num = 0;
	for(int i = 0; i < udiv; i++){
		//次にv方向の折れで分割する
		std::vector<UniqueRSBRep> srf_v;
		num += srf_u[i]->divide_multi_knot_v(srf_v);
		std::move(srf_v.begin(), srf_v.end(), std::back_inserter(srfl));
	}
	return num;
}

//u方向に折れ(マルチノット)があるとき面を分割する
//戻り値は、分割数を返却する
int MGRSBRep::divide_multi_knot_u(
	std::vector<UniqueRSBRep>& srfl///Divided MGRSBRep are appended.
)const{
	MGInterval intv = param_range().ref(1);
	const MGKnotVector &knot_vector = knot_vector_u();
	int start_index = order_u() - 1, index = 0, count = 0, multi = 0, bdim = bdim_u(), order = order_u();
	do{		//u方向に折れ(マルチノット)があるとき面を分割する
		if(order == 2){		//オーダー2のときの処理
			index = start_index + 1; multi = 1;
		}else{
			multi = knot_vector.locate_multi(start_index, order - 1, index);
		}
		MGBox uv_range = MGBox(MGInterval(knot_u(start_index), knot_u(index)), intv);
		MGRSBRep* rsb=new MGRSBRep;
		shrinkToParameters(uv_range, *rsb);
		srfl.emplace_back(rsb);
		start_index = index + multi - 1;
		count++;
	}while(index != bdim);	//多重度が見つからなかったら終わり
	return count;
}

//v方向に折れ(マルチノット)があるとき面を分割する
//戻り値は、分割数を返却する
int MGRSBRep::divide_multi_knot_v(
	std::vector<UniqueRSBRep>& srfl	///Divided MGRSBRep are appended.
) const{
	MGInterval intv = param_range().ref(0);
	const MGKnotVector &knot_vector = knot_vector_v();
	int start_index = order_v() - 1, index = 0, count = 0, multi = 0, bdim = bdim_v(), order = order_v();
	do{		//u方向に折れ(マルチノット)があるとき面を分割する
		if(order == 2){		//オーダー2のときの処理
			index = start_index + 1; multi = 1;
		}else{
			multi = knot_vector.locate_multi(start_index, order - 1, index);
		}
		MGBox uv_range = MGBox(intv, MGInterval(knot_v(start_index), knot_v(index)));
		MGRSBRep* rsb=new MGRSBRep;
		shrinkToParameters(uv_range, *rsb);
		srfl.emplace_back(rsb);
		//first_index = start_index = index + multi - 1;
		start_index = index + multi - 1;
		count++;
	}while(index != bdim);	//多重度が見つからなかったら終わり
	return count;
}

//uまたはv方向に折れ(C0連続)があるとき面を分割する
//戻り値は、分割数を返却する
int MGSurface::divide_multi_knot(
	std::vector<UniqueSurface>& srfl	///Divided objects are appended.
) const {
	srfl.emplace_back(copy_surface());
	return 1;
}

