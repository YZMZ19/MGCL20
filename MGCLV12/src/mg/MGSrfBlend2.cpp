/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#include "StdAfx.h"
#include "mg/SPointSeq.h"
#include "mg/LBRepEndC.h"
#include "mg/SBRep.h"
#include "mg/SBRepEndC.h"
#include "mg/SBRepTP.h"
#include "mg/SBRepVecTP.h"
#include "mg/BSumSurf.h"
#include "mg/Tolerance.h"

#include "cskernel/blg4sp2.h"

using namespace std;

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void get_all_derivatives(
	const MGLBRep* perimeters[4],//perimeters.
	const MGSBRep& surf,	//Original surface.
	const MGSBRepVecTP& vectp,	//接続面(パラメータ範囲は境界線と同じ)
	std::unique_ptr<MGLBRep> derivatives[4]//array of derivatives.
) {
	int sd = surf.sdim();
	const MGKnotVector& tu = surf.knot_vector_u();
	const MGKnotVector& tv = surf.knot_vector_v();

	MGVector N;
	///////////perimeter 0 and 2.
	int nu = surf.bdim_u(); int num1 = nu - 1;
	MGNDDArray utau2; utau2.buildByKnotVector(tu);
	MGBPointSeq deris_u(nu, sd);
	double v[2] = { tv.param_s(),tv.param_e() };

	MGNDDArray utau(nu); utau(0) = utau2[0];
	for (int m = 0; m < 2; m++) {//for perimeter 0 and 2.
		int perim = m * 2;
		if (!vectp.specified(perim)) continue;

		double vm = v[m];
		int i2 = 0;
		MGVector deri = perimeters[3]->eval(vm, 1);
		deris_u.store_at(i2++, deri);
		for (int i = 1; i < num1; i++) {
			double utaui = utau2[i];
			if (vectp.eval(perim, utaui, N)) {//get Normal of tangent plane.
				deri = surf.eval(utaui, vm, 0, 1);
				double vlen = deri.len();
				deri = MGUnit_vector(N*deri*N)*vlen;
				deris_u.store_at(i2, deri);
				utau(i2++) = utaui;
			}
		}
		deri = perimeters[1]->eval(vm, 1);
		deris_u.store_at(i2, deri);
		utau(i2++) = utau2[num1];
		utau.set_length(i2);
		deris_u.set_length(i2);
		MGLBRep* lb = new MGLBRep;
		lb->buildByInterpolationDataPoints(utau, deris_u);
		derivatives[perim].reset(lb);
	}

	////////// perimeter 1 and 3.
	int nv = surf.bdim_v(); int nvm1 = nv - 1;
	MGNDDArray vtau2; vtau2.buildByKnotVector(tv);
	MGBPointSeq deris_v(nv, sd);
	double u[2] = { tu.param_s(),tu.param_e() };
	int peri[2] = { 3,1 };

	MGNDDArray vtau(nv); vtau(0) = vtau2[0];
	for (int m = 0; m < 2; m++) {//for perimeter 3 and 1.
		int perim = peri[m];
		if (!vectp.specified(perim))
			continue;

		double um = u[m];
		int j2 = 0;
		MGVector deri = perimeters[0]->eval(um, 1);
		deris_v.store_at(j2++, deri);
		for (int j = 1; j < nvm1; j++) {
			double vtauj = vtau2[j];
			if (vectp.eval(perim, vtauj, N)) {//Normal of tangent plane at surf(um,vtau[j]).
				deri = surf.eval(um, vtauj, 1, 0);
				double vlen = deri.len();
				deri = MGUnit_vector(N*deri*N)*vlen;
				deris_v.store_at(j2, deri);
				vtau(j2++) = vtauj;
			}
		}
		deri = perimeters[2]->eval(um, 1);
		deris_v.store_at(j2, deri);
		vtau(j2++) = vtau2[nvm1];
		vtau.set_length(j2);
		deris_v.set_length(j2);
		MGLBRep* lb = new MGLBRep;
		lb->buildByInterpolationDataPoints(vtau, deris_v);
		derivatives[perim].reset(lb);
	}

	int i, l;
	//*****Twist adjustment. Here dudv and dvdu are only averaged, may be able to improve.
	//Adjust the twist vector at the corners.
	MGVector twist[4];
	for (i = 0; i < 4; i++) {
		int im1 = (i + 3) % 4;
		if (vectp.specified(im1) && vectp.specified(i)) {
			double tim1 = derivatives[im1]->param_e(), ti = derivatives[i]->param_s();
			if (im1 >= 2) tim1 = derivatives[im1]->param_s();
			if (i >= 2) ti = derivatives[i]->param_e();
			MGVector deriim1 = derivatives[im1]->eval(tim1, 1), derii = derivatives[i]->eval(ti, 1);
			twist[i] = (deriim1 + derii)*.5;
		}
	}
	MGLBRepEndC ECS, ECE;
	MGNDDArray utau3(3, tu.bdim() - 2, tu);
	MGNDDArray vtau3(3, tv.bdim() - 2, tv);
	int ntau2, nutau2 = utau3.length(), nvtau2 = vtau3.length();
	MGBPointSeq dyu(nutau2, 1), dyv(nvtau2, 1);
	for (l = 0; l < nutau2; l++) dyu(l, 0) = 1.;
	dyu(0, 0) = .01; dyu(nutau2 - 1, 0) = .01;
	for (l = 0; l < nvtau2; l++) dyv(l, 0) = 1.;
	dyv(0, 0) = .01; dyv(nvtau2 - 1, 0) = .01;
	MGNDDArray* tau2;
	const MGKnotVector* t;
	double* dy;
	for (i = 0; i < 4; i++) {
		if (!vectp.specified(i))continue;
		double ts, te;
		if (i % 2) {
			tau2 = &vtau3; t = &tv; dy = dyv.data(); ntau2 = nvtau2;
		}
		else {
			tau2 = &utau3; t = &tu; dy = dyu.data(); ntau2 = nutau2;
		}
		ts = (*tau2)[0]; te = (*tau2)[tau2->length() - 1];
		int ids = i, ide = i + 1; ide %= 4;
		if (i >= 2) { ids = ide; ide = i; }
		if (twist[ids].sdim() == 0)
			ECS.set_1st(derivatives[i]->eval(ts, 1));
		else
			ECS.set_1st(twist[ids]);
		if (twist[ide].sdim() == 0)
			ECE.set_1st(derivatives[i]->eval(te, 1));
		else
			ECE.set_1st(twist[ide]);
		MGBPointSeq bp;
		derivatives[i]->eval_line(*tau2, bp);
		double dev = bp(0).len() + bp(ntau2 / 2).len() + bp(ntau2 - 1).len();
		dev /= 3.;
		dev *= MGTolerance::angle_zero()*.7;
		MGLBRep* lb = new MGLBRep;
		lb->buildSRSmoothedLB_of_1stDeriv(ECS, ECE, *tau2, bp, dy, dev, false);
		derivatives[i].reset(lb);
	}
}
MGLBRep* get_perriderisub(
	const MGLBRep& lb,
	const MGLBRep* deriS,
	const MGLBRep* deriE
){
	double t0=lb.param_s(), t1=lb.param_e();
	MGNDDArray tau(2); tau(0)=t0; tau(1)=t1;
	MGBPointSeq bp(2,3);
	MGBPointSeq bp1(2,3);
	MGLBRepEndC endcS, endcE;
	bp.store_at(0,lb.start_point());
	if(deriS)
		endcS.set_1st(lb.eval(t0,1));
	bp.store_at(1,lb.end_point());
	if(deriE)
		endcE.set_1st(lb.eval(t1,1));
	MGLBRep* lbNew=new MGLBRep;
	lbNew->buildByInterpolationEC(endcS,endcE,tau,bp);
	return lbNew;
}

//compute perimeters and derivatives from only 4 corner point data
//of input perimeters and derivatives.
void get_peri2_deri2(
	const MGLBRep* perimeters[4],
	const std::unique_ptr<MGLBRep> derivatives[4],
	std::unique_ptr<MGLBRep> perimeters2[4],
	std::unique_ptr<MGLBRep> derivatives2[4]
				//array of derivatives[4].
){
	int m;
///////////perimeter 0 and 2.
	const MGLBRep* deriS=derivatives[3].get();
	const MGLBRep* deriE=derivatives[1].get();
	for(m=0; m<=2; m+=2){//for perimeter 0 and 2.
		perimeters2[m].reset(get_perriderisub(*(perimeters[m]),deriS,deriE));
		if(derivatives[m]){
			derivatives2[m].reset(get_perriderisub(*(derivatives[m]),deriS,deriE));
		}else derivatives2[m].reset();
	}

////////// perimeter 1 and 3.
	deriS=derivatives[0].get();
	deriE=derivatives[2].get();
	for(m=1; m<=3; m+=2){//for perimeter 1 and 3.
		perimeters2[m].reset(get_perriderisub(*(perimeters[m]),deriS,deriE));
		if(derivatives[m]){
			derivatives2[m].reset(get_perriderisub(*(derivatives[m]),deriS,deriE));
		}else derivatives2[m].reset();
	}
}

MGLBRep get1deri_of_peri(
	int peri,	//perimeter number
	const MGSBRep& surf
){
	bool alongu=true; if(peri%2) alongu=false;
	const MGKnotVector* t;
	const MGKnotVector* t_other;
	if(alongu){
		t=&(surf.knot_vector_u());
		t_other=&(surf.knot_vector_v());
	}else{
		t=&(surf.knot_vector_v());
		t_other=&(surf.knot_vector_u());
	}
	double t0=t_other->param_s();
	if(peri==1 || peri==2)
		t0=t_other->param_e();
	int len=t->bdim();
	MGBPointSeq dbp(len,3);
	const MGSPointSeq& spoint=surf.surface_bcoef();
	for(int j=0; j<len; j++){
		MGLBRep l0;
		l0.buildLBRepFromMemberData(*t_other, MGBPointSeq(!alongu, j, spoint));
		dbp.store_at(j,l0.eval(t0,1));
	}
	MGLBRep lb;
	lb.buildLBRepFromMemberData(*t, std::move(dbp));
	return lb;
}

MGSBRep* get_1DireSurf(
	bool udire,	//indicates if perimetes[0],[2] or [3],[1] should be used to construct
				//the surface. if udire=true, [0] and [2] are used.
	const MGLBRep* perimeters[4],
	const std::unique_ptr<MGLBRep> derivatives[4]
){
	int ncd=3;

	int ids,ide,otherS, otherE;
	if(udire){
		ids=0; ide=2;
		otherS=3; otherE=1;
	}else{
		ids=3; ide=1;
		otherS=0; otherE=2;
	}

//Herafter, u and v directions are temporal. When udire==true, u is v, and v is u.
//Construct u-direction data point and knot vector.
	double u0=perimeters[otherS]->param_s(),u1=perimeters[otherS]->param_e();
	MGENDCOND ecu0= MGENDCOND::MGENDC_1D, ecu1= MGENDCOND::MGENDC_1D;
	int lenu=2;
	if(!derivatives[ids]) ecu0= MGENDCOND::MGENDC_NO; else lenu++;
	if(!derivatives[ide]) ecu1= MGENDCOND::MGENDC_NO; else lenu++;
	int orderu=4; if(orderu>lenu) orderu=lenu;
	MGNDDArray utau(lenu);
	size_t i=0;
	utau((int)i++)=u0;
	if(ecu0== MGENDCOND::MGENDC_1D) utau((int)i++)=u0;
	if(ecu1== MGENDCOND::MGENDC_1D) utau((int)i++)=u1;
	utau((int)i++)=u1;
	MGKnotVector tu(utau,orderu);

//Construct v-direction data point and knot vector.
	const MGKnotVector& tv=perimeters[ids]->knot_vector();
	MGNDDArray vtau(3,tv.bdim()-2,tv);
	if(!derivatives[otherS]) vtau.add_data((vtau(0)+vtau(1))*.5);
	if(!derivatives[otherE]){
		int nvtau2=vtau.length();
		vtau.add_data((vtau(nvtau2-1)+vtau(nvtau2-2))*.5);		
	}

	double* wk=new double[lenu*2+lenu*(2*orderu-1)];
	double* q=wk+lenu*2;

	int ecu0i = static_cast<int>(ecu0);
	int ecu1i = static_cast<int>(ecu1);
	int lenum1=lenu-1, lenum2=lenu-2;
	MGBPointSeq bp(lenu,ncd);
	int nvtau=vtau.length();
	std::vector<MGLBRep*> lines;
	int error=2;
	for(int j=0; j<nvtau; j++){
		double vtauj=vtau(j);
		bp.store_at(0,perimeters[ids]->eval(vtauj));
		if(ecu0== MGENDCOND::MGENDC_1D){
			bp.store_at(1,derivatives[ids]->eval(vtauj));
			MGVector der=derivatives[ids]->eval(vtauj);
		}
		if(ecu1== MGENDCOND::MGENDC_1D){
			bp.store_at(lenum2,derivatives[ide]->eval(vtauj));
			MGVector der=derivatives[ide]->eval(vtauj);
		}
		bp.store_at(lenum1,perimeters[ide]->eval(vtauj));
		MGLBRep* lb=new MGLBRep(lenu,orderu,ncd);
		MGBPointSeq& obp=lb->line_bcoef();
		lb->knot_vector()=tu;
		for(int k=0; k<ncd; k++){
			blg4sp2_(orderu,&error,ecu0i,ecu1i,utau.data(),bp.data(0,k)
				,lenu,lenu,1,tu.data(),1,wk,wk+lenu,q,obp.data(0,k));
		}
		lines.push_back(lb);
	}

	std::vector<UniqueLBRep> derivatives2(2);
	int id[2]={otherS, otherE};
	for(int m=0; m<=1; m++){
	const MGLBRep* deri=derivatives[id[m]].get();
	if(deri){
		bp.store_at(0,deri->eval(u0));
		if(ecu0== MGENDCOND::MGENDC_1D) bp.store_at(1,deri->eval(u0,1));
		if(ecu1== MGENDCOND::MGENDC_1D) bp.store_at(lenum2,deri->eval(u1,1));
		bp.store_at(lenum1,deri->eval(u1));
		MGLBRep* deri2=new MGLBRep(lenu,orderu,ncd);
		MGBPointSeq& obp=deri2->line_bcoef();
		deri2->knot_vector()=tu;
		for(int k=0; k<ncd; k++){
			blg4sp2_(orderu,&error,ecu0i,ecu1i,utau.data(),bp.data(0,k)
				,lenu,lenu,1,tu.data(),1,wk,wk+lenu,q,obp.data(0,k));
		}
		derivatives2[m].reset(deri2);
	}
	}

	MGSBRep* sb=new MGSBRep;
	sb->buildByRibCurvesTangent(vtau, lines, derivatives2[0].get(), derivatives2[1].get());
	if(udire) sb->exchange_uv();

	delete[] wk;
	size_t nlines=lines.size();
	for(i=0; i<nlines; i++)
		delete lines[i];
	return sb;
}

//Construct Surface B-rep from lines and derivatives.
//Interpolation will be done only along one parameter direction,
//along v.
//tau provides data point sequence along v direction of surface (u,v) parameter
//configuration. deriS and deriE are used to provide the 1st derivative
//B-representation along the perimeter 0 and 2, may be null
//if 1st derivative B-rep is not provided. If derivative
//B-rep is provided, deriS and deriE must have the same knot configuration
//as the one of lines which makes u kont configuration of this surface (u,v)
//parameter. tau[i] is the parameter for lines[i].
int MGSBRep::buildByRibCurvesTangent(
	const MGNDDArray& tau,
	const std::vector<MGLBRep*>& lines,
	const MGLBRep* deriS,
	const MGLBRep* deriE
){
	invalidateBox();

	int i,j,k;
	const int ncd=3;

	m_uknot=lines[0]->knot_vector();
	MGENDCOND ec0= MGENDCOND::MGENDC_1D, ec1= MGENDCOND::MGENDC_1D;
	if(!deriS) ec0= MGENDCOND::MGENDC_NO;
	if(!deriE) ec1= MGENDCOND::MGENDC_NO;

// Compute b-rep dimension along v in lenv.
	//v=min and max condition.
	int lenu=m_uknot.bdim(), lenv1,lenv;
	lenv1=lenv=(int)lines.size();
	int ivs=0;
	if(ec0== MGENDCOND::MGENDC_1D){lenv+=1; ivs=1;}
	if(ec1== MGENDCOND::MGENDC_1D) lenv+=1;
	int orderv=4; if(orderv>lenv) orderv=lenv;
	MGSPointSeq surface_bcoef(lenv,lenu,ncd);//Temporal spoint seq.

// Prepare data point ordinate.
	// 1. Copy original data.
	int js;
	for(j=0; j<lenv1; j++){
		js=ivs+j;
		const MGBPointSeq& bcj=lines[j]->line_bcoef();
		for(k=0; k<ncd; k++){
			for(i=0; i<lenu; i++){
				surface_bcoef(js,i,k)=bcj(i,k);
			}
		}
	}

	// 2. Copy derivative data along perimeter from endc.
	//v=min condition.
	if(ec0== MGENDCOND::MGENDC_1D){
		const MGBPointSeq& bp1S=deriS->line_bcoef();
		for(k=0; k<ncd; k++){
			for(i=0; i<lenu; i++){
				surface_bcoef(0,i,k)=bp1S(i,k);
			}
		}
	}

	//v=max condition.
	int lenvm1=lenv-1, lenvm2=lenv-2;
	if(ec1== MGENDCOND::MGENDC_1D){
		const MGBPointSeq& bp1E=deriE->line_bcoef();
		for(k=0; k<ncd; k++){
			for(i=0; i<lenu; i++){
				surface_bcoef(lenvm1,i,k)=bp1E(i,k);
			}
		}
	}

	//Exchange positional data for blg4sp2_.
	//Perimeter 0.
	double save;
	if(ec0== MGENDCOND::MGENDC_1D){
		for(i=0; i<lenu; i++)
			for(k=0; k<ncd; k++){
				save=surface_bcoef(1,i,k);
				surface_bcoef(1,i,k)=surface_bcoef(0,i,k);
				surface_bcoef(0,i,k)=save;
			}
	}
	//Perimeter 2.
	if(ec1== MGENDCOND::MGENDC_1D){
		for(i=0; i<lenu; i++)
			for(k=0; k<ncd; k++){
				save=surface_bcoef(lenvm2,i,k);
				surface_bcoef(lenvm2,i,k)=surface_bcoef(lenvm1,i,k);
				surface_bcoef(lenvm1,i,k)=save;
			}
	}

	int lenuv=lenu*lenv;
	double* wk=new double[lenv*2+lenv*(2*orderv-1)];
	double* q=wk+lenv*2;
	MGNDDArray vtau(lenv);
	for(i=0; i<ivs; i++) vtau(i)=tau(0);
	for(j=0; j<lenv1; j++) vtau(i++)=tau(j);
	for(;i<lenv; i++) vtau(i)=tau(lenv1-1);
	vtau.set_length(lenv);
	m_vknot=MGKnotVector(vtau,orderv);

	int ec0i = static_cast<int>(ec0), ec1i = static_cast<int>(ec1);
	int error=2;
	m_surface_bcoef.resize(lenu,lenv,ncd);
	for(k=0; k<ncd; k++){
		blg4sp2_(orderv,&error,ec0i,ec1i,vtau.data(),surface_bcoef.data(0,0,k)
			,lenv,lenv,lenu,m_vknot.data(),lenu,wk,wk+lenv,q,m_surface_bcoef.data()+lenuv*k);
		if(error!=1) break;
	}
	if(error==1) error=0;
	delete[] wk;
	return error;
}

//Auxiliary fucntion of buildFromSidesCoonsWithTP and buildByBlendCoonsWithTP,
//gets a tenmporal surface to compute tangent at perimeters. The surface is obtained by
//buildGeneralizedRuledSurface according to the specification of tangent plane
//at the perimeter.
//buildGeneralizedRuledSurface() mixes two opposing edge's.
void buildGRuledSurfaceWithTP(
	bool tpIsSpecified[4],//Indicates if a tangent plane is specified for edge[].
	const MGLBRep*	edge[4],//境界線リスト(vmin,umax,vmax,uminの順,辺番号0,1,2,3の順)
	MGSBRep& sb
){
	sb=MGSBRep();
	bool rail13 = (!tpIsSpecified[0] && !tpIsSpecified[2]);
	bool rail02 = (!tpIsSpecified[1] && !tpIsSpecified[3]);
	bool hasTP = !rail13 || !rail02;//at least one tp is specified.

	if (hasTP){
		if(rail02)
			sb.buildGeneralizedRuledSurface(edge);
		else if(rail13)
			sb.buildGeneralizedRuledSurface(edge, false);
	}

	if (sb.is_null()) {
		MGSBRep ruled0; ruled0.buildGeneralizedRuledSurface(edge);
		MGSBRep ruled1; ruled1.buildGeneralizedRuledSurface(edge, false);
		MGSPointSeq sp(ruled0.surface_bcoef() + ruled1.surface_bcoef());
		sp *= .5;
		const MGKnotVector& tu = edge[0]->knot_vector();
		const MGKnotVector& tv = edge[1]->knot_vector();
		sb.buildSBRepFromMemberData(std::move(sp), tu, tv);
	}
}

///4本の境界線、ブレンド関数、接続面を与えて面を生成する。
///境界線はvmin,umax,vmax,uminの順で、vmin,vmaxの向きをuminからumaxの方向に
///umin,umaxの向きをvminからvmaxの方向になっているものとする。境界線のノットベクトル
///をあわせるときの誤差はline_zero()を使用している。
///接続面(MGSBRepTP)のパラメータ範囲は各境界線と同じとする。
void MGSBRep::buildFromSidesBoolSumWithTP(
	const MGCurve*	edge[4],///<境界線リスト(vmin,umax,vmax,uminの順,辺番号0,1,2,3の順)
	const MGSBRepTP&tp		///<接続面(パラメータ範囲は境界線と同じ)
){
	std::unique_ptr<MGLBRep> perimeters[4];
	rebuildAsSurfacePerimeters(edge, perimeters);
	buildFromSidesBoolSumWithTP(perimeters, MGSBRepVecTP(tp));
}

void MGSBRep::buildFromSidesBoolSumWithTP(const MGCurve * edge[4], const MGSBRepVecTP & tp){
	std::unique_ptr<MGLBRep> perimeters[4];
	rebuildAsSurfacePerimeters(edge, perimeters);
	buildFromSidesBoolSumWithTP(perimeters,MGSBRepVecTP(tp));
}

///4本の境界線、ブレンド関数、接続面を与えて面を生成する。
///境界線はvmin,umax,vmax,uminの順で、vmin,vmaxの向きをuminからumaxの方向に
///umin,umaxの向きをvminからvmaxの方向になっているものとする。境界線のノットベクトル
///をあわせるときの誤差はline_zero()を使用している。
///接続面(MGSBRepTP)のパラメータ範囲は各境界線と同じとする。
void MGSBRep::buildFromSidesBoolSumWithTP(
	const MGCurve * edge[4], MGSBRepVecTP && vtp
){
	std::unique_ptr<MGLBRep> perimeters[4];
	rebuildAsSurfacePerimeters(edge, perimeters);
	buildFromSidesBoolSumWithTP(perimeters, std::move(vtp));
}
void MGSBRep::buildFromSidesBoolSumWithTP(
	const UniqueLBRep	perimeters[4],///<境界線リスト(vmin,umax,vmax,uminの順,辺番号0,1,2,3の順)
	MGSBRepVecTP && vtp
){
	const MGLBRep*	peris[4];
	extractConstPointerVec(perimeters, perimeters+4, peris);
	buildFromSidesBoolSumWithTP(peris, std::move(vtp));
}
void MGSBRep::buildFromSidesBoolSumWithTP(
	const MGLBRep*	perimeters[4],///<境界線リスト(vmin,umax,vmax,uminの順,辺番号0,1,2,3の順)
	MGSBRepVecTP && vtp
){
	invalidateBox();

	MGKnotVector&  tu = knot_vector_u() = perimeters[0]->knot_vector();
	MGKnotVector&  tv = knot_vector_v() = perimeters[1]->knot_vector();
	vtp.change_range(true, tu.param_s(), tu.param_e());
	vtp.change_range(false, tv.param_s(), tv.param_e());

	bool tpSpecified[4] = { vtp.specified(0), vtp.specified(1),
						vtp.specified(2),vtp.specified(3) };
	MGSBRep ruled01;
	buildGRuledSurfaceWithTP(tpSpecified, perimeters, ruled01);
	std::unique_ptr<MGLBRep> derivatives[4];
	get_all_derivatives(perimeters, ruled01, vtp, derivatives);

	//Construct the boolian sum surface.
	MGSBRep* g1 = get_1DireSurf(true, perimeters, derivatives);
	MGSBRep* g2 = get_1DireSurf(false, perimeters, derivatives);

	std::unique_ptr<MGLBRep> derivatives2[4];
	std::unique_ptr<MGLBRep> perimeters2[4];
	get_peri2_deri2(perimeters, derivatives, perimeters2, derivatives2);

	std::vector<MGLBRep*> lines(2);
	MGNDDArray vtau0(2); vtau0(0) = tv.param_s(); vtau0(1) = tv.param_e();
	lines[0] = perimeters2[0].get(); lines[1] = perimeters2[2].get();
	MGSBRep* g12 = new MGSBRep;
	g12->buildByRibCurvesTangent(vtau0, lines, derivatives2[0].get(), derivatives2[2].get());

	MGBSumSurf g(g1, g2, g12);

	MGNDDArray utau(3, tu.bdim() - 2, tu);
	MGNDDArray vtau(3, tv.bdim() - 2, tv);
	MGSPointSeq spoint(utau.length(), vtau.length(), 3);
	g.eval_spoint(utau, vtau, spoint);
	MGSBRepEndC endc(utau, vtau, g);
	buildByInterpolationECWithKTV(endc, utau, vtau, spoint);
}

template<class InputIterator>
bool is_conerCoincide(
	InputIterator first
){
	MGPosition S[4], E[4];
	for(int i = 0; i<4; i++, first++){
		S[i] = (*first)->start_point();
		E[i] = (*first)->end_point();
	}
	return (S[0]==S[3] && E[0]==S[1] && E[1]==E[2] && S[2]==E[3]);
};

bool is_valid_perim(const MGCurve* perims[4]){
	if(is_conerCoincide(perims))
		if(isMGLBRep(perims, perims+4)){
			if(perims[0]->knot_vector()==perims[2]->knot_vector())
				return (perims[3]->knot_vector()==perims[1]->knot_vector());
		}
	return false;
}

bool is_valid_perim(std::unique_ptr<MGLBRep> perims[4]){
	if(is_conerCoincide(perims)){
		if(perims[0]->knot_vector()==perims[2]->knot_vector())
		return (perims[3]->knot_vector()==perims[1]->knot_vector());
	}
	return false;
}
bool is_valid_perim(const MGLBRep * perims[4]){
	if(is_conerCoincide(perims)){
		if(perims[0]->knot_vector()==perims[2]->knot_vector())
			return perims[3]->knot_vector()==perims[1]->knot_vector();
	}
	return false;
}
bool is_valid_perim(const std::unique_ptr<MGCurve> perims[4]) {
	const MGCurve* peri2[4];
	extractConstPointerVec(perims, perims + 4, peri2);
	return is_valid_perim(peri2);
}

bool is_valid_perim(MGCurve* perims[4]) {
	const MGCurve* peri2[4];
	extractConstPointerVec(perims, perims + 4, peri2);
	return is_valid_perim(peri2);
}

//Trim perimes[] at their corner points,
//which are obtained by neighbor's two perim's closest().
//Function's return value is:
//=0 trimmed successfully.
//!=0  the corner (ret-1) is degenerated(perimeter (ret-1) beccame a point).
int trimPerimeters(
	std::unique_ptr<MGCurve> perims[4]
) {
	MGCurve& vmin = *perims[0];
	MGCurve& umax = *perims[1];
	MGCurve& vmax = *perims[2];
	MGCurve& umin = *perims[3];

	MGPosition vmin_umax = vmin.closest(umax);
	MGPosition vmin_umin = vmin.closest(umin);
	MGPosition vmax_umax = vmax.closest(umax);
	MGPosition vmax_umin = vmax.closest(umin);

	if (MGRZero(vmin_umax(0) - vmin_umin(0)))
		return 1;
	vmin.limit(vmin_umax(0), vmin_umin(0));

	if (MGRZero(vmin_umax(1) - vmax_umax(1)))
		return 2;
	umax.limit(vmin_umax(1), vmax_umax(1));

	if (MGRZero(vmax_umax(0) - vmax_umin(0)))
		return 3;
	vmax.limit(vmax_umax(0), vmax_umin(0));

	if (MGRZero(vmin_umin(1) - vmax_umin(1)))
		return 4;
	umin.limit(vmin_umin(1), vmax_umin(1));
	return 0;
}

//Given 4 perimeters, update each curve's direction so that:
//(1) perim[i] makes perimeter i(0-3).
//(2) perim[0]'s end point coincides to perim[1]'s start.
//(3) The direction of [0] and [2], and of [3] and [1] are the same.
void updateDirection(std::unique_ptr<MGCurve> perims[4]) {
	MGCurve& peri0 = *perims[0]; MGCurve& peri1 = *perims[1];
	MGCurve& peri2 = *perims[2]; MGCurve& peri3 = *perims[3];
	MGPosition peri0S = peri0.start_point(), peri0E = peri0.end_point();
	MGPosition peri1S = peri1.start_point(), peri1E = peri1.end_point();

	double len[4];
	MGVector peri0Eperi1S = peri0E - peri1S;
	len[0] = peri0Eperi1S % peri0Eperi1S;

	MGVector peri0Eperi1E = peri0E - peri1E;
	len[1] = peri0Eperi1E % peri0Eperi1E;

	MGVector peri0Speri1S = peri0S - peri1S;
	len[2] = peri0Speri1S % peri0Speri1S;

	MGVector peri0Speri1E = peri0S - peri1E;
	len[3] = peri0Speri1E % peri0Speri1E;

	int minID=0;
	for(int i=1; i<4; i++)
		if(len[minID]>len[i])
			minID=i;

	if (minID == 2 || minID == 3)
		peri0.negate();
	if (minID == 1 || minID == 3)
		peri1.negate();

	double  t0m = (peri0.param_s() + peri0.param_e())*.5;
	double  t2m = (peri2.param_s() + peri2.param_e())*.5;
	if (peri0.eval(t0m, 1) % peri2.eval(t2m, 1) < 0.)
		peri2.negate();
	double  t1m = (peri1.param_s() + peri1.param_e())*.5;
	double  t3m = (peri3.param_s() + peri3.param_e())*.5;
	if (peri1.eval(t1m, 1) % peri3.eval(t3m, 1) < 0.)
		peri3.negate();
}

//Rebuild perimeters periIn to input to MGSBRep::buildByBlendXXX().
//rebuildCurveTrimDirectionUpdate() does:
//(1) Trim periIn at their corner points.
//(2) adjust for the opposite perimeters to have the same directions.
//(3) rebuild periIn for the opposite perimeters to have the same knot configuration.
//    At this rebuild, the corner points are updated to be the same.
int rebuildCurveTrimDirectionUpdate(
	const MGCurve*	periIn[4],//境界線リスト(vmin,umax,vmax,uminの順)
	std::unique_ptr<MGLBRep> perimeters[4]
) {
	UniqueCurve peri[4];
	for (int i = 0; i < 4; i++)
		peri[i].reset(periIn[i]->clone());

	// Do trim perims at first.
	int err = trimPerimeters(peri);
	if (!err) {
		// Arrange directions each of perims.
		updateDirection(peri);

		const MGCurve* edges[4];
		extractConstPointerVec(peri, peri + 4, edges);
		rebuildAsSurfacePerimeters(edges, perimeters);
	}
	return err;
}
