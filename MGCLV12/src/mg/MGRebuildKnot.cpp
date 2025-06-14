/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Box.h"
#include "mg/BPointSeq.h"
#include "mg/Position.h"
#include "mg/Curve.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/Tolerance.h"
#include "mg/SBRepTP.h"
using namespace std;

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// Implementation of knot rebuild.

//複数カーブの共通で削除できるノットを削除する。
//ただし、入力カーブは同じノットベクトルを持つものとする。
void remove_knot_curves(
	std::vector<UniqueLBRep>& brepList,		//曲線列
	MGLBRep**		tp,	//接続面	input and output.
		//if tp[i] for crvl[i] was not null, converted new tp will be output.
	double tp_length///<ratio of angle zero. When =zero, removal is done unconditionally.
		///<The actual error of the angle is set to tp_length*MGTolerance::angle_zero().
){
	MGLBRep& lb0 = *(brepList[0]);
	int nCrv = (int)brepList.size();	//曲線数
	int k = lb0.order(), n = lb0.bdim();//次元数&B表現次元
	int km1 = k - 1;
	int mid = (km1 + n) / 2;
	if (mid <= km1) mid = k;
	const double line0 = MGTolerance::line_zero();
	const double rc0 = tp_length * MGTolerance::angle_zero();
	int num_remained;
	int preDel, nDel, nDel2;

	std::vector<double> totalTol(nCrv, 0.);
	std::vector<double> totalTolTP(nCrv, 0.);
	std::vector<MGLBRep> lb(nCrv), tpSave(nCrv);
	std::vector<bool> tp_specified(nCrv);

	int i;
	for (i = n - 1; i >= mid;) {
		int j = 0;
		for (; j < nCrv; j++) {
			tp_specified[j] = tp && tp[j];
			MGLBRep& lbj = *(brepList[j]); lb[j] = lbj; //曲線は保存しておく
			nDel2 = nDel = lbj.remove_knot_one(line0, i, totalTol[j], num_remained);
			if (tp_specified[j]) {
				tpSave[j] = *tp[j];
				nDel2 = tp[j]->remove_knot_one(rc0, i, totalTolTP[j], num_remained);
			}
			if (j == 0) preDel = nDel;
			if (!nDel || nDel != nDel2 || preDel != nDel) {//ノット削除に失敗したときの処理
				std::fill(totalTol.begin(), totalTol.end(), 0.0);	//誤差合計をクリアする
				std::fill(totalTolTP.begin(), totalTolTP.end(), 0.0);	//誤差合計をクリアする
				int j2 = j;
				while (j2 >= 0) {
					MGLBRep& lbj2 = *(brepList[j2]); lbj2 = lb[j2];	//曲線を元に戻す
					if (tp_specified[j2]) *(tp[j2]) = tpSave[j2];
					j2--;
				}
				break;
			}
		}
		if (j >= nCrv)
			i -= (nDel + num_remained);
		else i--;
	}

	std::fill(totalTol.begin(), totalTol.end(), 0.0);	//誤差合計をクリアする
	std::fill(totalTolTP.begin(), totalTolTP.end(), 0.0);//誤差合計をクリアする
	i = k;
	int nt = lb0.bdim() - mid;
	while (lb0.bdim() - i > nt) {
		int j = 0;
		for (; j < nCrv; j++) {
			tp_specified[j] = false; if (tp) { if (tp[j]) tp_specified[j] = true; };
			MGLBRep& lbj = *(brepList[j]); lb[j] = lbj;//曲線は保存しておく
			nDel2 = nDel = lbj.remove_knot_one(line0, i, totalTol[j], num_remained);
			if (tp_specified[j]) {
				tpSave[j] = *tp[j];
				nDel2 = tp[j]->remove_knot_one(rc0, i, totalTolTP[j], num_remained);
			}
			if (j == 0) preDel = nDel;
			if (!nDel || nDel != nDel2 || preDel != nDel) {//ノット削除に失敗したときの処理
				std::fill(totalTol.begin(), totalTol.end(), 0.0);	//誤差合計をクリアする
				std::fill(totalTolTP.begin(), totalTolTP.end(), 0.0);	//誤差合計をクリアする
				int j2 = j;
				while (j2 >= 0) {
					MGLBRep& lbj2 = *(brepList[j2]); lbj2 = lb[j2];	//曲線を元に戻す
					if (tp_specified[j2]) *(tp[j2]) = tpSave[j2];
					j2--;
				}
				break;
			}
		}
		if (j >= nCrv)
			i += num_remained;
		else i++;
	}
}

//Rebuild input curves as:
//(1) Change all curves to MGLBRep's of input order that have the same parameter range.
//(2) When non MGLBRep curves are included, or different knot cinfiguration MGLBRep's are
//    mixed, they are rebuilt to MGLBRep within the line_zero() tolerance.
//入力された複数曲線列を指定オーダーで再構築する。トレランスはline_zero()を使用している。
//オーダーが指定されていないとき曲線列のうちで最も大きいオーダーを使用する。このとき、
//Ellipse, Straightのオーダーは4として考える。
//パラメータ範囲は1次微分値の大きさが１になるようにしたときの長さの平均を使用している。
//戻り値は再構築後の曲線列が返却される。
std::vector<UniqueLBRep> rebuildAsSameKnotVector(
	const std::vector<const MGCurve*>& crvl,//入力曲線列
	int ordr ,			//指定オーダー
	MGLBRep**		tp	//接続面	input and output.
		//if tp[i] for crvl[i] was not null, tp is converted.
){
	//使用するオーダーと空間次元を求める
	int nCrv = (int)crvl.size(), maxOrder = 0, maxSdim = 0, i = 0;
	int idMaxOrder=0;
	double allParam = 0.0;
	for(i = 0; i < nCrv; i++){
		int ord = crvl[i]->order(), sdim = crvl[i]->sdim();
		if(ord > maxOrder){
			maxOrder = ord; idMaxOrder=i;
		}
		if(sdim > maxSdim)
			maxSdim = sdim;
	}
	if(!ordr)
		ordr = maxOrder;		//オーダーが指定されていないときの処理
	if(ordr <= 3)
		ordr = 4;	//Ellipse, Straightのときオーダー4にする
					//Also we avoid order 3.

	//LBRepにrebuidしてリストに入れる
	std::vector<UniqueLBRep> rtnBrepList(nCrv);
	//bool fsame = true;
	double errorLine = MGTolerance::line_zero()*.3;
	for(i = 0; i < nCrv; i++){
		UniqueLBRep& lbi = rtnBrepList[i];
		const MGCurve& crvi = *crvl[i];
		std::unique_ptr<MGCurve> crviRebuilt = crvi.rebuild(2, 2, errorLine, 4);//rebuild as MGLBRep.
		lbi.reset(static_cast<MGLBRep*>(crviRebuilt.release()));

		//1次微分値の大きさが１になるようなパラメータ範囲を求める
		double t0= lbi->param_s(), t1= lbi->param_e();
		double t2=(t0+t1)*.5;
		double mag= lbi->eval(t0,1).len();
		mag+= lbi->eval(t1,1).len();
		mag+= lbi->eval(t2,1).len();
		mag /= 3.;
		allParam += lbi->param_span() * mag;
	}

	//パラメータ範囲を合わせたノットベクトルを作成する。
	double avgSpan = allParam / nCrv;
	for(i = 0; i< nCrv; i++){
		rtnBrepList[i]->change_range(0.0, avgSpan);
		if(tp && tp[i])
			tp[i]->change_range(0.0, avgSpan);
	}

	//全てのノットベクトルを足し合わせる。
	MGNDDArray tau_temp;
	tau_temp.buildByKnotVector(rtnBrepList[idMaxOrder]->knot_vector());
	if(tau_temp.length()<ordr)
		tau_temp.change_number(ordr);
	MGKnotVector mixedKnotVector(tau_temp,ordr);

	for(i = 1; i < nCrv; i++){
		MGKnotVector& ti=rtnBrepList[i]->knot_vector();
		int orderti=ti.order();
		int orderi= ordr<=orderti ? ordr: orderti;
		mixedKnotVector.mix_knot_vector(MGKnotVector(tau_temp.buildByKnotVector(ti),orderi));
	}

	//ノットベクトルを十分に増やす。divideNum()を使用する。
	MGKnotVector tempKnotVector = mixedKnotVector;
	int bd=tempKnotVector.bdim();
	for(i=ordr-1; i<bd; i++){
		int maxNumDiv = 0, j = 0;
		double	spara = tempKnotVector(i),
				epara = tempKnotVector(i + 1);
		if(MGRZero((epara - spara) / avgSpan))
			continue;	//マルチノットのときの処理

		for(; j < nCrv; j++){	//各スパンの最大分割数を求める
			int ndiv = rtnBrepList[j]->divideNum(MGInterval(spara, epara));
			if(ndiv>maxNumDiv)
				maxNumDiv = ndiv;
		}
		double span = (epara - spara) / maxNumDiv, tParam = tempKnotVector(i);
		for(j = 0; j < maxNumDiv - 1; j++){
			tParam += span;
			mixedKnotVector.add_data(tParam, ordr - 1);	//ノットベクトルを増やす
		}
	}

	//増やしたノットで曲線を再作成する
	MGNDDArray tau;
	tau.buildByKnotVector(mixedKnotVector);
	int n = tau.length();
	double tplen=0.;
	for(i = 0; i < nCrv; i++){
		UniqueLBRep& lbi = rtnBrepList[i];
		MGBPointSeq bp(n, maxSdim);
		for(int j = 0; j < n; j++)
			bp.store_at(j, lbi->eval(tau(j)));
		MGLBRep* lb=new MGLBRep;
		lb->setKnotVector(mixedKnotVector);
		lb->buildByInterpolationWithKTV(tau, bp);
		rtnBrepList[i].reset(lb);
		if(tp && tp[i]){
			double tpleni=0.;
			MGLBRep& tpi=*tp[i];
			MGBPointSeq bptp(n, maxSdim);
			for(int j = 0; j < n; j++){
				MGVector tpatj=tpi.eval(tau(j));
				bptp.store_at(j, tpatj);
				tpleni+=tpatj.len();
			}
			tpi.setKnotVector(mixedKnotVector);
			tpi.buildByInterpolationWithKTV(tau, bptp);
			tpleni/=double(n);
			tplen+=tpleni;
		}
	}
	tplen/=double(nCrv);

	//曲線列のノットベクトルを精度を保持して削除を行う
	mgTolSetLineZero setLzero(errorLine);
	remove_knot_curves(rtnBrepList,tp,tplen);
	return rtnBrepList;
}

///Rebuild this curve.
std::unique_ptr<MGCurve> MGCurve::rebuild(
	int how_rebuild,
		//intdicates how rebuild be done.
		// =0: no approximation(only parameter change)
		// =1: if this is rational spline(MGRLBRep), reconstructed with new knot configuration
		//     as rational spline(MGRLBRep).
		//     Otherwise approximated by non-rational spline(MGLBRep) with new knot configuration.
		// =2: approximated by non-rational spline(MGLBRep) with new knot configuration,
		//     if this is rational spline. If this is not rational spline, same as =1.
	int parameter_normalization,
		//Indicates how the parameter normalization be done:
		//=0: no parameter normalization.
		//=1: normalize to range=(0., 1.);
		//=2: normalize to make the average length of the 1st derivative 
		//    is as equal to 1. as possible.
		//=3: specify parameter range in param_range.
	double tol,	///<tolerance allowed for the approximation
		///When tol<=0., MGTolerance::line_zero() will be employed.
	int ordr,	///<order of the new MGLBRep, >=4 is recommended.
		///When order=0 is input, the original order is unchanged if this curve is
		///MGLBRep or MGRLBRep. Otherwise order is set to 4.
	const double* param_range
)const{
	if(parameter_normalization>=3 && param_range==0)
		parameter_normalization=2;

	double tol_old, err=tol;
	if(tol>0.)
		tol_old=MGTolerance::set_line_zero(tol);
	else
		err=MGTolerance::line_zero();

	int neworder=ordr;
	if(neworder==0){
		neworder=order();
	}
	int reparaType=parameter_normalization;
	int reparaType2=reparaType;
		
	double tspan[2]={0., 1.};
	double& t1=tspan[0];
	double& t2=tspan[1];
	if(reparaType>=2)
		reparaType2=1;
	if(reparaType==3){
		t1=param_range[0];
		t2=param_range[1];
		if(t1>=t2){
			t2=t1+1.;
		}
	}

	std::unique_ptr<MGCurve> curvenew;
	if(how_rebuild){
		const MGRLBRep* rlb=dynamic_cast<const MGRLBRep*>(this);
		if(rlb && how_rebuild==1){
			std::unique_ptr<MGRLBRep> rlb2=rlb->rebuild_with_new_knot_configuration(err,reparaType2);
			curvenew.reset(rlb2.release());
		}else{
			MGLBRep* lbrep=new MGLBRep;
			approximate_as_LBRep(*lbrep,neworder,reparaType2,true);
			curvenew.reset(lbrep);
		}
	}else{
		curvenew=std::unique_ptr<MGCurve>(clone());
		if(reparaType)
			curvenew->change_range(0.,1.);
	}
	if(reparaType==2)
		t2=curvenew->get_average_tangent_length();

	if(reparaType>=2)
		curvenew->change_range(t1,t2);

	if(tol>0.)
		MGTolerance::set_line_zero(tol_old);
	curvenew->copy_appearance(*this);
	return curvenew;
}

///Approximate this curve as a MGLBRep curve
///within the tolerance MGTolerance::line_zero().
///When parameter_normalization=0, reparameterization will not be done, and
///the evaluation at the same parameter has the same values before and after
///of approximate_as_LBRep.
void MGCurve::approximate_as_LBRep(
	MGLBRep& lb,	///<Approximated LBRep will be set.
	int ordr,		///<new order
	int parameter_normalization,
		//Indicates how the parameter normalization be done:
		//=0: no parameter normalization.
		//=1: normalize to range=(0., 1.);
		//=2: normalize to make the average length of the 1st derivative 
		//    is as equal to 1. as possible.
	bool neglectMutli
)const{
	lb.copy_appearance(*this);
	lb.invalidateBox();

	const MGKnotVector& t=knot_vector();
	double ts=param_s(), te=param_e();
	int k=t.locate(ts)+1, n=t.locate(te)+1;
		//k!=t.order() or n!=t.bdim() of t occurs when this is a trimmedcurve.
	int km1=t.order()-1, start=k;
	int index, multi_found;
	MGInterval pspan(ts,te);
	MGKnotVector& tnew=lb.knot_vector();
	int norder=ordr ? ordr:4;

	do{	//Approximation by dividing to parts of continuity>=C0.
		multi_found=t.locate_multi(start,km1,index);//Locate C0 continuity point.
		if(start==k){//For the 1st span.
			approximate_as_LBRep2(lb,norder,start-1,index,neglectMutli);//1st approximation.
			if(ts>lb.param_s() || te<lb.param_e())
				lb.limit(pspan);

			double te2=lb.param_e();
			if(parameter_normalization){
				double vlen=lb.eval(te2,1).len();
				double tsnew=tnew.param_s();
				tnew-=tsnew;
				tnew*=vlen;
			}
		}else{//For the span from the 2nd.
			MGLBRep lbt;
			approximate_as_LBRep2(lbt,norder,start-1,index,neglectMutli);//from the 2nd approximation.
			if(te<lbt.param_e())
				lbt.limit(pspan);
			int which=2;
			int cn=0;
			if(parameter_normalization){
				double ratio;
				cn=lb.continuity(lbt,which,ratio);
			}
			lb.connect(cn,which,lbt);
		}
		start=index+multi_found;
	}while(index<n);
	if(parameter_normalization==1)
		tnew.change_range(0.,1.);
}

#define INCREASE_NUM 10
///Get data points for approximate_as_LBRep2.
void MGCurve::data_points_for_approximate_as_LBRep2(
	int is, int ie,//approximation parameter range, from knot_vector()[is] to [ie].
	MGKnotVector& t,//New knot configuration will be output.
				//t's order is input. other information of t will be updated.
	MGNDDArray& tau,//Data point for t will be output.
	bool neglectMulti///<Indicates if multiple knots be kept.
		///< true: multiplicity is removed.
		///< false: multiplicity is kept.
)const{
	const MGKnotVector& told=knot_vector();
	int kold=told.order();
	int knew=t.order();//new order.
	int kdiff=knew-kold;
	int nnew=knew;//nnew will be new B-Rep dimension.
	if(neglectMulti){
		tau.buildByKnotVector(told);
		tau.change_number(told.bdim()*INCREASE_NUM);
		t=MGKnotVector(tau,knew);
		return;
	}

	int index, multi_found,multi_new;
	int start=is+1;
	do{	//Count the new B-Rep dimension.
		multi_found=told.locate_multi(start,1,index);//Get the next multiplicity.
		multi_new=multi_found+kdiff;
		if(multi_new<1 || multi_found==1)
			multi_new=1;

		nnew+=multi_new+INCREASE_NUM;
		start+=multi_found;
	}while(index<ie);
	nnew-=multi_new;
	t.size_change(knew,nnew);
	
	const double& ts=told(is);
	int i=0;//is the index of t to store the data t(i).
	for(; i<knew; i++)
		t(i)=ts;
	double incNp1=double(INCREASE_NUM+1);

	start=is+1;
	do{	//Count the new bspline rep dimension.
		multi_found=told.locate_multi(start,1,index);//Get the next multiplicity.
		multi_new=multi_found+kdiff;
		if(multi_new<1 || multi_found==1)
			multi_new=1;

		double dif=(told(start)-told(start-1))/incNp1;
		for(int j=0; j<INCREASE_NUM; j++,i++)
			t(i)=t(i-1)+dif;
		for(int j=0; j<multi_new; j++,i++)
			t(i)=told(index);
		start+=multi_found;
	}while(index<ie);
	const double& te=told(ie);
	for(int j=0; j<knew-multi_new; j++)
		t(i++)=te;
	tau.buildByKnotVector(t);assert(tau.length()==nnew);
}

//Approximate this curve as a MGLBRep curve from knot_vector[is] to [ie].
//This is an internal program of MGLBRep constructor.
void MGCurve::approximate_as_LBRep2(
	MGLBRep& lb,		//Approximated LBRep will be set.
	int ordr,		//new order
	int is, int ie,//approximation parameter range, from knot_vector()[is] to [ie].
	bool neglectMulti///<Indicates if multiple knots be kept.
		///< true: multiplicity is removed.
		///< false: multiplicity is kept.
)const{
	//精度十分の曲線を生成する
	MGNDDArray tau;
	MGKnotVector& t=lb.knot_vector();
	t.change_order(ordr);
	data_points_for_approximate_as_LBRep2(is,ie,t,tau,neglectMulti);
	MGBPointSeq bp1;
	eval_line(tau,bp1);

	//制御点を生成する
	lb.buildByInterpolationWithKTV(tau, bp1);
	lb.remove_knot();
	lb.copy_appearance(*this);
}

//境界線 edge[]、接続面tpを与え(optional)、rebuildされた境界をperimetersに出力、また
//tpのパラメータ範囲も整える。
//(1) 対辺が同じノットベクトルのスプライン曲線になるように再作成
//(edge[0]と[2]は同じ方向でu方向, [1]と[3]も同じ方向でv方向を入力.方向の調整は行わない)
//(2) 相対する境界線、接続面が同じノットベクターを持つようrebuildする。
//(3) コーナーの点が同じとなるよう無条件に調整する。
//
//境界線はC1連続であり、vmin,umax,vmax,uminの順で、vmin,vmaxの向きをuminから
//umaxの方向にumin,umaxの向きをvminからvmaxの方向になっているものとする。
//境界線のノットベクトルをあわせるときの誤差はline_zero()を使用している。
void rebuildAsSurfacePerimeters(
	const MGCurve*	edge[4],//境界線リスト(vmin,umax,vmax,uminの順)
	std::unique_ptr<MGLBRep> perimeters[4],
	MGSBRepTP*		tp		//接続面(パラメータ範囲は境界線と同じ)
){
	//接続面のパラメータ範囲が境界線と同じかどうか調べる
	if(tp)
		for(int i=0; i<4; i++){
			if(!tp->specified(i))
				continue;
			MGLBRep& tpi= tp->TP(i);
			MGInterval cirange=edge[i]->param_range();
			if(cirange != tpi.param_range())
				tpi.change_range(cirange.low_point(), cirange.high_point());
		}

	//(1) 対辺が同じノットベクトルのスプライン曲線(MGLBRep)になるように再作成
	int k=4;	//オーダーは４とする 
	std::vector<const MGCurve*> temp_crvl(2);
	for(int j=0; j<2; j++){
		int jp2 = j + 2;
		temp_crvl[0] = edge[j];	temp_crvl[1] = edge[jp2];
		MGLBRep* tplb[2]={nullptr, nullptr };
		if(tp){
			if (tp->specified(j))
				tplb[0] = &(tp->TP(j));
			if (tp->specified(jp2))
				tplb[1] = &(tp->TP(jp2));
		}
		std::vector<UniqueLBRep> temp_brepl = rebuildAsSameKnotVector(temp_crvl, k, tplb);
		perimeters[j].reset(temp_brepl[0].release());
		perimeters[j+2].reset(temp_brepl[1].release());
	}

	//(2) コーナーの点が同じとなるよう無条件に調整する。
	MGPosition P[4];
	double tse[4];
	for(int i=0; i<4; i++){
		MGLBRep& peri=*perimeters[i];

		int inxt=(i+1)%4;
		int id0=inxt, id1=i;
		if(i<2){
			id0=i, id1=inxt;
			tse[i*2]=peri.param_s();
			tse[1+i*2]=peri.param_e();
		}
		P[id0]+=peri.start_point();
		P[id1]+=peri.end_point();
	}

	for(int i=0; i<4; i++)
		P[i]*=.5;

	//P0
	perimeters[3]->move(2,tse[2],P[0],tse+3);
	perimeters[0]->move(2,tse[0],P[0], tse+1);
	//P1
	perimeters[0]->move(2,tse[1],P[1],tse);
	perimeters[1]->move(2,tse[2],P[1], tse+3);
	//P2
	perimeters[1]->move(2,tse[3],P[2],tse+2);
	perimeters[2]->move(2,tse[1],P[2], tse);
	//P3
	perimeters[2]->move(2,tse[0],P[3],tse+1);
	perimeters[3]->move(2,tse[3],P[3], tse+2);
}
