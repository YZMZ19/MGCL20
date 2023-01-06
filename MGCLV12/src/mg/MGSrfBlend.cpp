/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#include "StdAfx.h"
#include "mg/Tolerance.h"
#include "mg/Box.h"
#include "mg/KnotArray.h"
#include "mg/SPointSeq.h"
#include "mg/Straight.h"
#include "mg/CCisects.h"
#include "mg/LBRepEndC.h"
#include "mg/SBRepEndC.h"
#include "mg/SBRepTP.h"
#include "mg/Coons.h"
#include "mg/SBRep.h"

#include "cskernel/Bvstan.h"
using namespace std;

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void buildGRuledSurfaceWithTP(
	bool tpIsSpecified[4],//Indicates if a tangent plane is specified for edge[].
	const MGLBRep*	edge[4],//境界線リスト(vmin,umax,vmax,uminの順,辺番号0,1,2,3の順)
	MGSBRep& sb	//generalized ruled surface is output.
);

//Compute 2nd derivative.
MGVector get_2nd(const MGVector& first,const MGVector& second_guess,double curvature){
	if(MGRZero(second_guess.len()))
		return second_guess;

	MGUnit_vector e2(second_guess);
	double len=first.len()*curvature;
	return e2*len;
}

//ノットをあわせた後の境界線、ブレンド関数、データポイントから点列を求める。
void bilinear_spoint_proc(
	const MGLBRep* brepl[4],//ノットベクトルをあわせた後の境界線
	const MGCurve& blendCrvU,	//空間次元1のu方向のブレンド曲線(パラメータ、値域ともに0,1)
	const MGCurve& blendCrvV,	//空間次元1のv方向のブレンド曲線(パラメータ、値域ともに0,1)
	const MGNDDArray& utau,		//u方向のデータポイント
	const MGNDDArray& vtau,		//v方向のデータポイント
	MGSPointSeq& spoint		//点列
){
	spoint = MGSPointSeq(utau.length(), vtau.length(), brepl[0]->sdim());

	//4隅の点は境界線の端点の中間とする
	MGPosition	p00 = (brepl[0]->start_point() + brepl[3]->start_point()) / 2.0,
				p01 = (brepl[2]->start_point() + brepl[3]->end_point()) / 2.0,
				p10 = (brepl[0]->end_point() + brepl[1]->start_point()) / 2.0,
				p11 = (brepl[1]->end_point() + brepl[2]->end_point()) / 2.0;

	//境界線上の点を求める
	int i, nu=utau.length(), nv=vtau.length();
	int num1=nu-1, nvm1=nv-1;
	for(i=0; i<nu; i++){
		spoint.store_at(i, 0, brepl[0]->eval(utau(i)));
		spoint.store_at(i, nvm1, brepl[2]->eval(utau(i)));
	}
	for(i=0; i<nv; i++){
		spoint.store_at(0, i, brepl[3]->eval(vtau(i)));
		spoint.store_at(num1, i, brepl[1]->eval(vtau(i)));
	}

	//内部点を求める
	double spanu = brepl[0]->param_span(), spanv = brepl[1]->param_span();
	for(i=1; i<num1; i++){
		double blendu = (blendCrvU.eval(utau(i) / spanu))(0);
		double onembu=1.-blendu;
		MGPosition	pu0(spoint(i, 0)), pu1(spoint(i, nvm1));
		for(int j=1; j<nvm1; j++){
			double blendv = (blendCrvV.eval(vtau(j) / spanv))(0);
			double onembv=1.-blendv;
			MGPosition	p0v(spoint(0, j)), p1v(spoint(num1, j));
			spoint.store_at(i,j,onembu*p0v+blendu*p1v+onembv*pu0+blendv*pu1
				- (onembv*(onembu*p00+blendu*p10)+blendv*(onembu*p01+blendu*p11)));
		}
	}
}

//境界線、ブレンド関数、接続面を与え、対辺が同じノットベクトルのスプライン曲線
//になるように再作成して、点列、ノットベクトル、データポイントを求め、境界線の
//パラメータにあわせてあった接続面をリビルドした後のパラメータに変更する。
//境界線はC1連続であり、vmin,umax,vmax,uminの順で、vmin,vmaxの向きをuminから
//umaxの方向にumin,umaxの向きをvminからvmaxの方向になっているものとする。
//境界線のノットベクトルをあわせるときの誤差はline_zero()を使用している。
void bilinear_spoint(
	const MGCurve& blendCrvU,//空間次元1のu方向のブレンド曲線(パラメータ、値域ともに0,1)
	const MGCurve& blendCrvV,//空間次元1のv方向のブレンド曲線(パラメータ、値域ともに0,1)
	const MGLBRep* perimeters[4],
	MGSPointSeq& spoint,//点列
	MGNDDArray& utau,	//u方向のデータポイント
	MGNDDArray& vtau,	//v方向のデータポイント
	const MGSBRepTP* tp//接続面(パラメータ範囲は境界線と同じ)
){
	//端末ベクトル、データポイントを求める。
	MGENDCOND condk[4]={ MGENDCOND::MGENDC_NO,MGENDCOND::MGENDC_NO,MGENDCOND::MGENDC_NO,MGENDCOND::MGENDC_NO};
	for(int i=0; i<4; i++)
		if(tp && tp->specified(i))
			condk[i] = MGENDCOND::MGENDC_1D;

	//点列を求める
	const MGKnotVector& knotu = perimeters[0]->knot_vector();
	const MGKnotVector& knotv = perimeters[1]->knot_vector();
	utau = MGNDDArray(condk[3], condk[1], knotu);
	vtau = MGNDDArray(condk[0], condk[2], knotv);
	bilinear_spoint_proc(perimeters, blendCrvU, blendCrvV, utau, vtau, spoint);
}

//境界線、接続面、点列データポイント、接ベクトルから接続条件(MGSBRepEndC)を求める
bool bilinear_endc(
	MGSPointSeq&		spoint,		//点列
	const MGSBRepTP&		tp,			//接続面(パラメータ範囲は境界線と同じ)
	const MGLBRep* perimeters[4],
	const MGCurve&			blendCrvU,		//空間次元1のu方向のブレンド曲線(パラメータ、値域ともに0,1)
	const MGCurve&			blendCrvV,		//空間次元1のv方向のブレンド曲線(パラメータ、値域ともに0,1)
	const MGNDDArray&		utau,		//u方向のデータポイント
	const MGNDDArray&		vtau,		//v方向のデータポイント
	MGSBRepEndC&			endc)	//端末条件
{
	if(utau.length() != spoint.length_u() || vtau.length() != spoint.length_v() || spoint.sdim() != 3)return false;
	const MGLBRep& peri0=*(perimeters[0]);
	const MGLBRep& peri1=*(perimeters[1]);
	const MGLBRep& peri2=*(perimeters[2]);
	const MGLBRep& peri3=*(perimeters[3]);

	MGENDCOND condk[4]={ MGENDCOND::MGENDC_NO,MGENDCOND::MGENDC_NO,MGENDCOND::MGENDC_NO,MGENDCOND::MGENDC_NO};
	const int dim=3; int m;
	int i,p;
	int lenu=utau.length(), lenv=vtau.length();
	int lenum1=lenu-1, lenvm1=lenv-1;
	double ur[2]={utau(0),utau(lenum1)}, vr[2]={vtau(0),vtau(lenvm1)};
	const int ipse[2]={1,2};

	int ktp, ntp, n;
	const double* ttp; const double* rtp;
	int itanse[2]={1,1}; double tanse[6];
	int irtp, ip1, ip2;
	double work[105], tangen[3];

	int tps[4]={0,0,0,0};
	for(i=0; i<4; i++) if(tp.specified(i)) tps[i]=1;

	if(tps[0] || tps[2]){

	const MGKnotVector& knotv = peri1.knot_vector();
	MGNDDArray vtau2(condk[0], condk[2], knotv);
	MGSPointSeq spoint2v;
	bilinear_spoint_proc(perimeters, blendCrvU, blendCrvV, utau, vtau2, spoint2v);
	spoint2v.capacity(ip1, ip2);
	n=vtau2.length();
	double spanu = peri0.param_span();
	for(m=0; m<2; m++){
		//Process of perimeter num 0 and 2(v=min and max line)
		int perimeter=m*2;
		if(tps[perimeter]){
			const MGLBRep& tpm = tp.TP(perimeter);
			ktp=tpm.order(); ntp=tpm.bdim();
			ttp=tpm.knot_data(); rtp=tpm.coef_data();
			irtp=tpm.line_bcoef().capacity();
		//We have to generate first and 2nd derivative data for this perimeter.
			double v=vr[m];
			MGVector deri1[2]={peri3.eval(v,1),peri1.eval(v,1)};

			MGBPointSeq first(lenu,dim);
			for(p=0; p<dim; p++) tanse[p]=first(0,p)=deri1[0][p];
			for(p=0; p<dim; p++) tanse[p+3]=first(lenum1,p)=deri1[1][p];
			for(i=1; i<lenum1; i++){
				double tptau=utau(i);
				bvstan_(ur,ktp,ntp,ttp,rtp,tptau,n,vtau2.data(),
					spoint2v.data(i,0,0),ipse[m],itanse,tanse,irtp,
					ip1,ip2,work,tangen);
				MGVector deri1ati(dim,tangen);
				first.store_at(i,deri1ati);
			}
			endc.set_1st(perimeter,std::move(first));
		}
	}

	}

	int ntemp=lenv-tps[0]-tps[2];
	if((tps[0] || tps[2]) && ntemp>=2){

	MGVector deri2sv[2]={peri3.eval(vr[0],2),peri1.eval(vr[0],2)};
	MGVector deri2ev[2]={peri3.eval(vr[1],2),peri1.eval(vr[1],2)};
	double crvtr0s=peri3.curvature(vr[0]), crvtr1s=peri1.curvature(vr[0]);
	double crvtr0e=peri3.curvature(vr[1]), crvtr1e=peri1.curvature(vr[1]);
	double spanu = peri0.param_span();
	MGNDDArray taut(ntemp);
	taut(0)=vtau(0);
	for(int j=1 ; j<ntemp-1; j++) taut(j)=vtau(j+tps[0]);
	taut(ntemp-1)=vtau(lenv-1);
	for(i=1; i<lenum1; i++){
		double tptau=utau(i);
		MGLBRepEndC endc0, endc2;
		MGBPointSeq bpt(ntemp,3);
		bpt.store_at(0,spoint(i,0));
		double blendu = (blendCrvU.eval(tptau / spanu))(0);
		if(tps[0]){
			MGVector deri1=MGVector(3,endc.first(0)(i));
			endc0.set_1st(deri1);
			MGVector deri2ati0=deri2sv[0].interpolate_by_rotate(blendu,deri2sv[1]);
			deri2ati0=get_2nd(deri1,deri2ati0,crvtr0s+(crvtr1s-crvtr0s)*blendu);
			endc0.set_2nd(deri2ati0);
		}
		for(int j=1 ; j<ntemp-1; j++)
			bpt.store_at(j,spoint(i,j+tps[0]));
		if(tps[2]){
			MGVector deri1=MGVector(3,endc.first(2)(i));
			endc2.set_1st(deri1);
			MGVector deri2ati0=deri2ev[0].interpolate_by_rotate(blendu,deri2ev[1]);
			deri2ati0=get_2nd(deri1,deri2ati0,crvtr0e+(crvtr1e-crvtr0e)*blendu);
			endc2.set_2nd(deri2ati0);
		}
		bpt.store_at(ntemp-1,spoint(i,lenv-1));
		int k=6;
		int new_bdim=ntemp+tps[0]*2+tps[2]*2;
		if(new_bdim<k) k=new_bdim;
		MGLBRep lbt;lbt.buildByInterpolationEC(endc0, endc2, taut, bpt, k);
		if(tps[0]) spoint.store_at(i,1,lbt.eval(vtau(1)));
		if(tps[2]) spoint.store_at(i,lenv-2,lbt.eval(vtau(lenv-2)));
	}

	}

	if(tps[3] || tps[1]){

	const MGKnotVector& knotu = peri0.knot_vector();
	MGNDDArray utau2(condk[3], condk[1], knotu);
	MGSPointSeq spoint2u;
	bilinear_spoint_proc(perimeters, blendCrvU, blendCrvV, utau2, vtau, spoint2u);
	int psizeu, psizev;
	spoint2u.capacity(psizeu, psizev);
	n=utau2.length();
	ip1=1; ip2=psizeu*psizev;
	const int iv[2]={3,1};
	double spanv = peri1.param_span();
	for(m=0; m<2; m++){
		//Process of perimeter num 3 and 1(u=min and max line)
		int perimeter=iv[m];
		if(tps[perimeter]){
			const MGLBRep& tpm = tp.TP(perimeter);
			ktp=tpm.order(); ntp=tpm.bdim();
			ttp=tpm.knot_data(); rtp=tpm.coef_data();
			irtp=tpm.line_bcoef().capacity();
		//We have to generate first and 2nd derivative data for this perimeter.
			double u=ur[m];
			MGVector deri1[2]={peri0.eval(u,1),peri2.eval(u,1)};

			MGBPointSeq first(lenv,dim);
			for(p=0; p<dim; p++) tanse[p]=first(0,p)=deri1[0][p];
			for(p=0; p<dim; p++) tanse[p+3]=first(lenvm1,p)=deri1[1][p];
			for(i=1; i<lenvm1; i++){
				double tptau=vtau(i);
				bvstan_(vr,ktp,ntp,ttp,rtp,tptau,n,utau2.data(),
					spoint2u.data(0,i,0),ipse[m],itanse,tanse,irtp,
					ip1,ip2,work,tangen);
				MGVector deri1ati(dim,tangen);
				first.store_at(i,deri1ati);
			}
			endc.set_1st(perimeter,std::move(first));
		}
	}

	}

	ntemp=lenu-tps[3]-tps[1];
	if((tps[3] || tps[1]) && ntemp>=2){

	MGVector deri2su[2]={peri0.eval(ur[0],2),peri2.eval(ur[0],2)};
	MGVector deri2eu[2]={peri0.eval(ur[1],2),peri2.eval(ur[1],2)};
	double crvtr0s=peri0.curvature(ur[0]), crvtr1s=peri2.curvature(ur[0]);
	double crvtr0e=peri0.curvature(ur[1]), crvtr1e=peri2.curvature(ur[1]);

	double spanv = peri1.param_span();
	MGNDDArray taut(ntemp);
	taut(0)=utau(0);
	for(int j=1 ; j<ntemp-1; j++) taut(j)=utau(j+tps[3]);
	taut(ntemp-1)=utau(lenu-1);
	for(i=1; i<lenvm1; i++){
		double tptau=vtau(i);
		MGLBRepEndC endc0, endc2;
		MGBPointSeq bpt(ntemp,3);
		bpt.store_at(0,spoint(0,i));
		double blendv = (blendCrvV.eval(tptau / spanv))(0);
		if(tps[3]){
			MGVector deri2ati0=deri2su[0].interpolate_by_rotate(blendv,deri2su[1]);
			endc0.set_1st(MGVector(3,endc.first(3)(i)));
			endc0.set_2nd(deri2ati0);
		}
		for(int j=1 ; j<ntemp-1; j++) bpt.store_at(j,spoint(j+tps[3],i));
		if(tps[1]){
			MGVector deri2ati0=deri2eu[0].interpolate_by_rotate(blendv,deri2eu[1]);
			endc2.set_1st(MGVector(3,endc.first(1)(i)));
			endc2.set_2nd(deri2ati0);
		}
		bpt.store_at(ntemp-1,spoint(lenu-1,i));
		int k=6;
		int new_bdim=ntemp+tps[3]*2+tps[1]*2;
		if(new_bdim<k) k=new_bdim;
		MGLBRep lbt;lbt.buildByInterpolationEC(endc0, endc2, taut, bpt, k);
		if(tps[3]) spoint.store_at(1,i,lbt.eval(utau(1)));
		if(tps[1]) spoint.store_at(lenu-2,i,lbt.eval(utau(lenu-2)));
	}

	}
	return true;
}

double get_deri_coef(double t0, double t1,double alpha, double t){
	if(t<=.5) return (2.*t*(1.-alpha)+2.*alpha*t0-(t0+t1))/(t0-t1);
	return (2.*t*(1.-alpha)+2.*alpha*t1-(t0+t1))/(t1-t0);
}

//4本の境界線、ブレンド関数、接続面を与えて面を生成する。
//境界線はvmin,umax,vmax,uminの順で、vmin,vmaxの向きをuminからumaxの方向に
//umin,umaxの向きをvminからvmaxの方向になっているものとする。境界線のノットベクトル
//をあわせるときの誤差はline_zero()を使用している。ブレンド曲線はパラメータ範囲0,1
//で値域も0,1である。接続面(MGSBRepTP)のパラメータ範囲は各境界線と同じとする。
//		エラーコード：
//		0:			正常終了
//		-1:			境界線がC1連続でなかった
//		-2:			接続面のパラメータ範囲が境界線と違った
//		-3:			境界線のノットベクトルをあわせられなかった
//		-4:			接続面のパラメータ範囲を変更できなかった
//		-5:			点列が求まらなかった
//		-6:			端末条件が求まらなかった
//		-7:			面が生成できなかった
int MGSBRep::buildByBlendWithTP(
	const MGLBRep * perimeters[4],
	MGSBRepTP&& tp,
	const MGCurve * blendCrvU,
	const MGCurve * blendCrvV
){
	//1. Define blending function.
	MGPosition zero(1); zero(0) = 0.;
	MGPosition one(1); one(0) = 1.;
	MGStraight zero_one(one, zero);
	const MGCurve& blendu = blendCrvU ? *blendCrvU : zero_one;
	const MGCurve& blendv = blendCrvV ? *blendCrvV : zero_one;;

	//2. Build spoint data.
	MGSPointSeq spoint;
	MGNDDArray utau, vtau;
	bilinear_spoint(blendu, blendv, perimeters, spoint, utau, vtau, &tp);

	//3. 接続条件(MGSBRepEndC)を求める
	MGSBRepEndC endc;
	if(!bilinear_endc(spoint,tp,perimeters,blendu,blendv,utau,vtau,endc))
		return -6;

	//4. 点列、データポイント、ノットベクトル、接続条件より面を生成する
	m_uknot = std::move(perimeters[0]->knot_vector());
	m_vknot = std::move(perimeters[1]->knot_vector());
	return  buildByInterpolationECWithKTV(endc, utau, vtau, spoint);
}
int MGSBRep::buildByBlendWithTP(
	const MGCurve*	edge[4],//境界線リスト(vmin,umax,vmax,uminの順,辺番号0,1,2,3の順)
	const MGSBRepTP& tp,	//接続面(パラメータ範囲は境界線と同じ)
	const MGCurve*	blendCrvU,//空間次元1のu方向のブレンド曲線(パラメータ、値域ともに0,1)
	const MGCurve*	blendCrvV//空間次元1のv方向のブレンド曲線(パラメータ、値域ともに0,1)
){

	//1. Adjust parameter range of the tp and edge.
	MGSBRepTP tempTP(tp);
	std::unique_ptr<MGLBRep> perimeters[4];
	rebuildAsSurfacePerimeters(edge,perimeters, &tempTP);
	const MGLBRep* peris[4];
	extractConstPointerVec(perimeters, perimeters+4, peris);
	return buildByBlendWithTP(peris, std::move(tempTP), blendCrvU, blendCrvV);
}

//Easy to use version of buildByBlendWithTP.
//When blendCrvU,V were null, straight line from 0. to 1. will be used.
int MGSBRep::buildByBlend(
	const MGCurve*	edges[4],//境界線リスト(vmin,umax,vmax,uminの順,辺番号0,1,2,3の順)
	const MGCurve*	blendCrvU,//空間次元1のu方向のブレンド曲線(パラメータ、値域ともに0,1)
	const MGCurve*	blendCrvV //空間次元1のv方向のブレンド曲線(パラメータ、値域ともに0,1)
){
	UniqueLBRep peris[4];
	rebuildAsSurfacePerimeters(edges, peris);
	const MGLBRep*	perimeters[4];
	extractConstPointerVec(peris, peris+4, perimeters);
	return buildByBlend(perimeters, blendCrvU, blendCrvV);
}

int MGSBRep::buildByBlend(
	const MGLBRep*	perimeters[4],///<境界線リスト(vmin,umax,vmax,uminの順,辺番号0,1,2,3の順)
	const MGCurve*	blendCrvU,///<空間次元1のu方向のブレンド曲線(パラメータ、値域ともに0,1)
	const MGCurve*	blendCrvV ///<空間次元1のv方向のブレンド曲線(パラメータ、値域ともに0,1)
){
	assert(is_valid_perim(perimeters));

	MGPosition zero(1); zero(0)=0.;
	MGPosition one(1); one(0)=1.;
	MGStraight zero_one(one,zero);
	const MGCurve& blendu = blendCrvU ? *blendCrvU : zero_one;
	const MGCurve& blendv = blendCrvV ? *blendCrvV : zero_one;;

	//2. Build spoint data.
	MGSPointSeq spoint;
	MGNDDArray utau, vtau;
	bilinear_spoint(blendu,blendv, perimeters,spoint,utau,vtau,0);

	//3. 点列、データポイント、ノットベクトル、接続条件より面を生成する
	m_uknot = std::move(perimeters[0]->knot_vector());
	m_vknot = std::move(perimeters[1]->knot_vector());
	
	MGSBRepEndC endc;
	return buildByInterpolationECWithKTV(endc, utau, vtau, spoint);
}

void get_derivatives(
	bool along_u,			//true if for perimeter 0 and 2.
							//false if for perimeter 3 and 1.
	const MGSBRep& surf,	//Original surface.
	const MGSBRepTP& tp,	//接続面(パラメータ範囲は境界線と同じ)
	std::unique_ptr<MGLBRep> derivatives[4]//array of derivatives.
){
	const MGKnotVector& t=along_u ? surf.knot_vector_u():surf.knot_vector_v();
	MGNDDArray tau; tau.buildByKnotVector(t);

	int n=t.bdim();
	int sd=surf.sdim();
	MGBPointSeq deris(n, sd);

	const MGKnotVector& tOther=along_u ? surf.knot_vector_v():surf.knot_vector_u();
	double tS=t.param_s(), tE=t.param_e();

	double paraPeri[2]={tOther.param_s(), tOther.param_e()};
	int peri[2]={3, 1};
	if(along_u){
		peri[0]=0; peri[1]=2;
	}

	for(int m=0; m<2; m++){//for perimeter 3 and 1.
		int perim=peri[m];
		double paraPerim=paraPeri[m];
		for(int j=0; j<n; j++){
			double tauj=tau[j];
			MGVector deri=along_u ? 
				surf.eval(tauj,paraPerim,0,1):surf.eval(paraPerim, tauj, 1, 0);
			if(tp.specified(perim)){
				MGVector N=tp.TP(perim).eval(tauj);
					//Normal of tangent plane at surf(um,vtau[j]).
				double vlen=deri.len();
				deri=MGUnit_vector(N*deri*N)*vlen;
			}
			deris.store_at(j, deri);
		}
		MGLBRep* deriLB=new MGLBRep;
		deriLB->setKnotVector(t);
		deriLB->buildByInterpolationWithKTV(tau, deris);
		derivatives[perim].reset(deriLB);
	}
}



///4本の境界線、ブレンド関数、接続面を与えて面を生成する。
///境界線はvmin,umax,vmax,uminの順で、vmin,vmaxの向きをuminからumaxの方向に
///umin,umaxの向きをvminからvmaxの方向になっているものとする。境界線のノットベクトル
///をあわせるときの誤差はline_zero()を使用している。
///接続面(MGSBRepTP)のパラメータ範囲は各境界線と同じとする。
void MGSBRep::buildFromSidesCoonsWithTP(
	const MGCurve*	edge[4],//境界線リスト(vmin,umax,vmax,uminの順,辺番号0,1,2,3の順)
	const MGSBRepTP& tp	//接続面(パラメータ範囲は境界線と同じ)
){
	invalidateBox();

	MGSBRepTP tp2(tp);
	std::unique_ptr<MGLBRep> perimeters[4];
	rebuildAsSurfacePerimeters(edge,perimeters, &tp2);

	//Save the original parameter range.
	MGKnotVector& tu=perimeters[0]->knot_vector();
	double u0=tu.param_s(), u1=tu.param_e();
	MGKnotVector& tv=perimeters[1]->knot_vector();
	double v0=tv.param_s(), v1=tv.param_e();

	//Change parameter range to (0., 1.).
	int i;
	for(i=0; i<4; i++){
		perimeters[i]->change_range(0., 1.);
		if(tp2.specified(i))
			tp2.TP(i).change_range(0., 1.);
	}

//Method 1. Get uled surfaces according to the specification of TP[i].
	bool tpSpecified[4] = { tp2.specified(0), tp2.specified(1),
						tp2.specified(2),tp2.specified(3) };
	const MGLBRep* peris[4];
	extractConstPointerVec(perimeters, perimeters+4, peris);
	MGSBRep ruled01;
	buildGRuledSurfaceWithTP(tpSpecified, peris, ruled01);

	std::unique_ptr<MGLBRep> derivatives[4];
	bool along_u=true;
	get_derivatives(along_u,ruled01,tp2,derivatives);
	get_derivatives(!along_u,ruled01,tp2,derivatives);

//Method 2.
	//Save the old knot vector.
	MGKnotVector& tu2 = m_uknot = tu;
	MGKnotVector& tv2 = m_vknot = tv;
	//Construct a coons' patch.
	MGCoons coons(perimeters,derivatives);

	MGNDDArray utau(3,tu2.bdim()-2,tu2);
	//tv2.change_knot_number(tv2.bdim()*3);
	MGNDDArray vtau(3,tv2.bdim()-2,tv2);
	MGSPointSeq spoint(utau.length(),vtau.length(),3);
	coons.eval(utau,vtau,spoint);
	MGSBRepEndC endc(utau,vtau,coons);
	buildByInterpolationECWithKTV(endc, utau, vtau, spoint);
	change_range(1,u0,u1);
	change_range(0,v0,v1);
}

///Build MGSBRep, given the opposing 2 side curves
///and optionally their tangent plane MGLBRep.
///edge[0] makes perimeter 0.
///When two edge edge[] have opposite direction(at the middle point), edge[1] is negated.
///When tangent plane's direction is opposite, they are negated.
///Tangetn agnitude along v parameter can be input through tangentMagnitude[2].
///When tangentMagnitude is specified, the derivative vector 
///at start(tangentMagnitude[0]) and at end([1]) is multiplied by tangentMagnitude[].
void MGSBRep::buildFrom2Sides(
	const MGCurve*	edge[2],///<Two of the 4 perimeters of the building MGSBRep.
		///edge[0] makes perimeter 0 and edge[1] makes perimeter 2.
	const MGLBRep* tpIn[2],	///<Tangent plane LBRep of edge[0] and [1], may be null.
	const double* tangentMagnitude
){
	int k=4;	//オーダーは４とする 

	//Make edge[] parameter ranges are the same.
	const MGCurve& edge0=*edge[0];
	double t0=edge0.param_s(), t1=edge0.param_e();
	UniqueCurve edge1Temp(edge[1]->clone());//Make temporal curve to update the parameter range.
	edge1Temp->change_range(t0, t1);

	//Make tp[] parameter ranges are the same.
	UniqueLBRep tpTemp[2];
	MGLBRep* tp[2]{nullptr, nullptr};//To input rebuildAsSameKnotVector.
	for(int i=0; i<2; i++){
		if(tpIn[i]){
			MGLBRep* tpi=tpIn[i]->clone();
			tpTemp[i].reset(tpi);
			tp[i]=tpi;
			tpi->change_range(t0, t1);
		}
	}

	//Adjust parameter direction of edge[1] and tp[].
	double tmid=(t0+t1)*.5;
	MGVector TPeri2=edge1Temp->eval(tmid, 1);
	if(edge0.eval(tmid, 1)%TPeri2<0.){
		edge1Temp->negate();
		if(tp[1]) tp[1]->negate();
	}

	//Adjust normal direction of tp[].
	std::vector<const MGCurve*> temp_crvl{&edge0, edge1Temp.get()};
	MGVector T=edge1Temp->eval(tmid)-edge0.eval(tmid);//v=const parameter direction of this.
	for(int i=0; i<2; i++){
		MGVector B=temp_crvl[i]->eval(tmid, 1);
		if(tp[i] && (B*T)%tp[i]->eval(tmid)<0.)
			*tp[i] *=-1.;
	}

	//Rebuild temp_crvl[](edge[]) to have the same knot configuration.
	std::vector<UniqueLBRep> lbreps = rebuildAsSameKnotVector(temp_crvl, k, tp);
	MGLBRep& peri0=*lbreps[0]; MGLBRep& tp0=*tp[0];
	MGLBRep& peri2=*lbreps[1]; MGLBRep& tp2=*tp[1];

	MGKnotVector& tu=peri0.knot_vector();
	MGNDDArray utau; utau.buildByKnotVector(tu);
	int nu=utau.length();
	int nv=2;
	int sd=3;
	MGSPointSeq sp(nu, nv,sd);

	double mag0=1.,mag2=1.;
	if(tangentMagnitude){
		mag0=tangentMagnitude[0];
		mag2=tangentMagnitude[1];
	}

	//Define v data point.
	double utau3[3]={utau[0], utau[nu/2], utau[nu-1]};
	double vlen=0.;
	for(int j=0; j<3; j++){
		double utauj=utau3[j];
		MGPosition P0=peri0.eval(utauj), P2=peri2.eval(utauj);
		vlen+=(P2-P0).len();
	}
	vlen/=3.;
	MGNDDArray vtau(2); vtau(0)=0.; vtau(1)=vlen;
	mag0/=vlen; mag2/=vlen;

	MGBPointSeq firstD0(nu,sd);//1st derivatives along v direction at peri0  are stored.
	MGBPointSeq firstD2(nu,sd);//at peri2  are stored.
	bool peri0TP=tp[0], peri2TP=tp[1];//true if tp is specified.
	for(int i=0; i<nu; i++){
		double utaui=utau[i];
		MGPosition P0=peri0.eval(utaui), P2=peri2.eval(utaui);
		sp.store_at(i, 0, P0);
		sp.store_at(i, 1, P2);
		if(!peri0TP && !peri2TP)
			continue;

		MGVector V02={P2-P0};
		double vleni=V02.len();
		if(peri0TP){
			MGVector N=tp0.eval(utaui);
			MGVector T=(N*V02*N).normalize();
			double mag=mag0*vleni;
			T*=mag;
			firstD0.store_at(i, T);
		}
		if(peri2TP){
			MGVector N=tp2.eval(utaui);
			MGVector T=(N*V02*N).normalize();
			double mag=mag2*vleni;
			T*=mag;
			firstD2.store_at(i, T);
		}
	}

	MGSBRepEndC sEndc;
	if(peri0TP)
		sEndc.set_1st(0, std::move(firstD0));
	if(peri2TP)
		sEndc.set_1st(2, std::move(firstD2));

	m_uknot=std::move(tu);
	makeKnotVectorFromEC1Dire(false,sEndc, vtau);
	buildByInterpolationECWithKTV(sEndc, utau, vtau, sp);
}