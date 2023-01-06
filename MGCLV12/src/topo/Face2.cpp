/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Tolerance.h"
#include "mg/CParam_list.h"
#include "mg/Position_list.h"
#include "mg/CSisects.h"
#include "mg/SSisects.h"
#include "mg/Straight.h"
#include "mg/SurfCurve.h"
#include "mg/LBRep.h"
#include "mg/TrimmedCurve.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "topo/Edge.h"
#include "topo/Loop.h"
#include "topo/LPoint.h"
#include "topo/Face.h"
#include "topo/Shell.h"
#include "topo/FOuterCurve.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
using namespace std;

//
//Implements MGFace Class.
//MGFace is an instance of MGCell.

///////Member Function///////

//Extract all the loops of this face.
void MGFace::extract_loops(std::vector<const MGLoop*>& loops)const{
	loops.clear();
	int n=(int)m_boundaries.size();
	loops.resize(n);
	for(int i=0; i<n; i++){
		loops[i]=m_boundaries[i].get();
	}
}

//Proprietry routine for make_outer_boundary(), will append a newed straight line
//edge to the end of lp if param (uv1 and uv2) are far away enough compared with error
//err.
void MGFace::sl_edge_append(
	MGLoop*& lp,	//lp that the edge should be append.
					//If lp is null, new lp will be generated.
	int id,			//perimeter num of the staight line edge.
	const MGPosition& uv2,	//face parameter of the end point of sl.
	const MGPosition& uv1,	//face parameter of the start point of sl.
	double err_sqr			//Error square allowed to regard as same points.
){
	double dist1=uv1[0]-uv2[0]; dist1*=dist1;
	double dist2=uv1[1]-uv2[1]; dist2*=dist2;
	if((dist1+dist2)>err_sqr){
		MGEdge* e=new MGEdge(new MGStraight(uv2,uv1));
		MGCurve* crv=surface()->perimeter_curve(id);
		int is_v=id%2;
		double ts=uv1[is_v], te=uv2[is_v];
		if(id>=2){
			ts=crv->negate_param(ts); te=crv->negate_param(te);
			crv->negate();
		}
		e->resetBinder(new MGEdge(crv,MGInterval(ts,te)));
		if(lp) lp->append(e);
		else{
			lp=new MGLoop(e);
			prepend_boundary(lp);
		}
	}
}

//This is a newed MGFace or MGSurface object.
//If this is a MGFace, returns this pointer.
//If this is a MGSurface, construct a newed MGFace using this newed MGSurface,
//and returns the MGFace*.
MGFace* MGSurface::make_face(){
	MGFace* f=new MGFace(this);
	return f;
}

//Make outer boundary if not existed.
void MGFace::make_outer_boundary(){
	const MGSurface& srf=*(surface());
	if(srf.perimeter_num()<=0) return;

	MGBox prange=srf.param_range();
	double u0=prange[0].low_point(), u1=prange[0].high_point();
	double v0=prange[1].low_point(), v1=prange[1].high_point();
	MGPosition P[4]={MGPosition(u0,v0),MGPosition(u1,v0),
						MGPosition(u1,v1),MGPosition(u0,v1)};
	if(no_outer_boundaries()){
	//If no outer boundaries, make from surface perimeter.
		MGEdge *e0, *ei; MGCurve* binder;
		//Perimeter 0.
		e0=new MGEdge(new MGStraight(P[1],P[0]));
		e0->resetBinder(new MGEdge(srf.perimeter_curve(0)));
		MGLoop* lp=new MGLoop(e0);
		//Perimeter 1.
		ei=new MGEdge(new MGStraight(P[2],P[1]));
		ei->resetBinder(new MGEdge(srf.perimeter_curve(1)));
		lp->append(ei);
		//Perimeter 2.
		ei=new MGEdge(new MGStraight(P[3],P[2]));
		binder=srf.perimeter_curve(2); binder->negate();
		ei->resetBinder(new MGEdge(binder));
		lp->append(ei);
		//Perimeter 3.
		ei=new MGEdge(new MGStraight(P[0],P[3]));
		binder=srf.perimeter_curve(3); binder->negate();
		ei->resetBinder(new MGEdge(binder));
		lp->append(ei);
		//Make the loop closed.
		ei->join(false,e0);
		prepend_boundary(lp);
		return;
	}
	MGLoop* lp=loop(int(0)).get();
	if((*lp).is_outer_boundary()) return;

//Case that perimeter boundary.
	int id, idp1, ide,i, id1,id2;
	//id: surface perimeter id of uv.
	//ids: perimeter id of uvstart;
	//i: loop id.
	MGPosition uv,uve;
	//uv: last point processed.
	//uve: start point of the lp(and will be the end point of the outer_boundary).

	lp->both_end_on_perimeter(id1,id2);
	MGPosition uv1=lp->start_point(), uv2=lp->end_point();
	double t1=uv1[id1%2], t2=uv2[id2%2];
	if(id1>=2) t1*=-1.; if(id2>=2) t2*=-1.;
	bool reverse_on_1peri=(id1==id2 && t2<t1);
	if(id1>id2 || reverse_on_1peri){
		id=id2; uv=uv2;
		ide=id1; uve=uv1;
		i=1;	//To process next loop.
	}else{
		uv=P[0];id=0;
		uve=P[0]; ide=3;
		lp=0; i=0;	//To process first loop.
	}

	double err_sqr=parameter_error(); err_sqr*=err_sqr;
	int n=number_of_loops();
	MGLoop* tlp;
	while(i<n){//Process all the perimeter boudary loops.
		tlp=loop(i++).get();
		if(tlp->both_end_on_perimeter(id1,id2)){
			while(id1>id){
				idp1=id+1;
				sl_edge_append(lp,id,P[idp1],uv,err_sqr);
				uv=P[idp1];id=idp1;
			}
			uv2=tlp->start_point();
			sl_edge_append(lp,id,uv2,uv,err_sqr);
			if(!lp)
				lp=tlp;//tlp is the 1st loop in this and will be the outer loop.
			else{
				free_boundary(tlp);//Release tlp from this Face.
				lp->join(false, tlp);
			}
			uv=lp->end_point(); id=id2;
		}else{break;}
	}
	while(ide>id){
		idp1=id+1;
		sl_edge_append(lp,id,P[idp1],uv,err_sqr);
		uv=P[idp1]; id=idp1;
	}
	sl_edge_append(lp,id,uve,uv,err_sqr);
	//Make the loop closed.
	lp->last_edge()->join(false,lp->first_edge());
}

///Obtain outer_boundary curves(world coordinates representation) of the FSurface.
///Let the output of outer_boundary() be wcurves and of outer_boundary_param()
///be pcurves, then wcurves[i] corresponds to pcurves[i] one to one.
///The output curves can be considered as a continuous counter-clockwise ordered
///boundary of the surface.
std::vector<UniqueCurve> MGFace::outer_boundary()const{
	std::vector<MGFOuterCurve> cid=outer_curve();

	const MGSurface& srf=*(surface());
	std::vector<UniqueCurve> perim(srf.perimeter_num());
		//To save whole perimeter curve to avoid multiple process of
		//same perimeters.

	std::vector<UniqueCurve> crvs;
	//Convert id's in cid to curve expressions.
	std::vector<MGFOuterCurve>::iterator i = cid.begin(), ie = cid.end();
	for(; i!=ie; i++){
		if(i->is_loop()){
			std::vector<UniqueCurve> crvsWorld=i->loop()->curves_world();
			std::move(crvsWorld.begin(), crvsWorld.end(), std::back_inserter(crvs));
		}else{
			int id=i->perimeter_id();
			if(perim[id]==0) perim[id].reset(srf.perimeter_curve(id));
			double t0,t1; i->range(t0,t1);
			if(t0>t1){ double save=t0; t0=t1; t1=save;}
			crvs.emplace_back(perim[id]->part(t0,t1));
		}
	}
	return crvs;
}

//Obtain outer_boundary curves(parameter space coordinates representation)
//of the face.
//Let the output of outer_boundary() be wcurves and of outer_boundary_param()
//be pcurves, then wcurves[i] corresponds to pcurves[i] one by one.
std::vector<UniqueCurve> MGFace::outer_boundary_param()const{
	std::vector<MGFOuterCurve> cid=outer_curve();

	const MGSurface& srf=*(surface());
	std::vector<UniqueCurve> crvs;

	//Convert id's in cid to curve expressions.
	std::vector<MGFOuterCurve>::iterator i = cid.begin(), ie = cid.end();
	for(; i!=ie; i++){
		if(i->is_loop()){
			std::vector<UniqueCurve> crvsLp=i->loop()->curves();
			std::move(crvsLp.begin(), crvsLp.end(), std::back_inserter(crvs));
		}else{
			int id=i->perimeter_id();
			double t0,t1; i->range(t0,t1);
			MGPosition uv0=srf.perimeter_uv(id,t0), uv1=srf.perimeter_uv(id,t1);
			MGStraight* sl=new MGStraight(uv1,uv0,t1, t0);
			crvs.emplace_back(sl);
		}
	}
	return crvs;
}

//Obtain outer boundary curve expression as the combination of
//loop pointers and perimeter id's.
//******Currently this is valid only for MGSBRep, MGRSBRep, or Plane******
std::vector<MGFOuterCurve> MGFace::outer_curve() const{
	const MGSurface* srf=surface();

	MGBox prange=srf->param_range();
	double u0=(prange.ref(0)).low_point();
	double u1=(prange.ref(0)).high_point();
	double v0=(prange.ref(1)).low_point();
	double v1=(prange.ref(1)).high_point();
	double ts_corner[4]={u0,v0,u1,v1}, te_corner[4]={u1,v1,u0,v0};

	int i;
	if(no_outer_boundaries()){
		//If no outer boundaries, output surface perimeters.
		int n=srf->perimeter_num();
		std::vector<MGFOuterCurve> perim(n);
		for(i=0; i<n; i++) perim[i]=MGFOuterCurve(i,ts_corner[i],te_corner[i]);
		return perim;
	}

//////There exist perimeter or outer boundary.//////

	int id,ide, id1,id2;
	//id: surface perimeter id of uv.
	//ide: perimeter id of last perimeter;
	//i: loop id.

	MGPosition uv; double t,te,t2;
	const MGLoop* lp=loop(int(0)).get();
	if(lp->is_perimeter_boundary()){
	//Case that perimeter boundary.
		lp->both_end_on_perimeter(id1,id2);
		MGPosition uv1=lp->start_point(), uv2=lp->end_point();
		double s1=uv1[id1%2], s2=uv2[id2%2];
		if(id1>=2) s1*=-1.; if(id2>=2) s2*=-1.;
		bool reverse_on_1peri=(id1==id2 && s2<s1);

		std::vector<MGFOuterCurve> crvs;
		if(id1>id2 || reverse_on_1peri){
			crvs.push_back(MGFOuterCurve(lp));
			id=id2; uv=lp->end_point(); t=uv.ref(id2%2);
			ide=id1; uv=lp->start_point(); te=uv.ref(ide%2);
			i=1;	//To process next loop.
		}else{
			id=0; t=u0; 
			ide=3; te=v0;
			i=0;	//To process first loop.
		}

		int n=number_of_loops();
		double err=parameter_error();
		while(i<n){//Process all the perimeter boudary loops.
			lp=loop(i++).get();
			if(lp->both_end_on_perimeter(id1,id2)){
				while(id1>id){
					t2=te_corner[id];
					if(fabs(t2-t)>err) crvs.push_back(MGFOuterCurve(id,t,t2));
					t=ts_corner[++id];
				}
				uv=lp->start_point();
				t2=uv.ref(id1%2);
				if(fabs(t2-t)>err) crvs.push_back(MGFOuterCurve(id1,t,t2));
				crvs.push_back(MGFOuterCurve(lp));
				uv=lp->end_point(); id=id2; t=uv.ref(id%2);
			}else{break;}
		}

		while(ide>id){
			t2=te_corner[id];
			if(fabs(t2-t)>err) crvs.push_back(MGFOuterCurve(id,t,t2));
			t=ts_corner[++id];
		}
		t2=te;
		if(fabs(t2-t)>err) crvs.push_back(MGFOuterCurve(ide,t,t2));
		return crvs;
	}else{
	//Case that outer boundary.
		return std::vector<MGFOuterCurve>(1,MGFOuterCurve(lp));
	}
}

//Obtain parameter curves.
//In the case of surface, parameter curve is only one. However, in the case
//of face,  number of parameter curves are more than one.
std::vector<UniqueCurve> MGFace::parameter_curves(
	int is_u,	//True if x is u-value.(i.e. obtain u=const line)
	double x	//parameter value. u or v-value accordint to is_u.
)const{
	const MGSurface* srf=surface();
	int kcod= is_u ? 0:1;
	std::unique_ptr<MGCurve> whole_curve(srf->parameter_curve(is_u,x));
	std::vector<UniqueCurve> curves;

	int iperi;
	if(srf->on_a_perimeter2(is_u,x,iperi)){//if on perimeters.
		std::vector<MGInterval> pranges=perimeter_param_range(iperi);
		std::vector<MGInterval>::iterator
			ipr=pranges.begin(), iprend=pranges.end();
		double t0,t1;
		while(ipr!=iprend){
			t0=(*ipr).low_point(); t1=(*ipr++).high_point();
			curves.emplace_back(whole_curve->part(t0,t1));
		}
	}else{//When x is not on a perimeter.
		std::vector<double> pvec=isect1D_with_boundaries(x,kcod);
		int nm1=(int)(pvec.size()-1);
		double u,v;
		if(is_u) u=x; else v=x;
		for(int i=0;i<nm1; i++){
			double t0=pvec[i], t1=pvec[i+1];
			double tmid=t0+t1; tmid*=0.5;
			if(is_u) v=tmid; else u=tmid;
			if(in_range(u,v))
				curves.emplace_back(whole_curve->part(t0,t1));
		}
	}
	return curves;
}

//Obtain parent shell that this face belongs to.
MGShell* MGFace::parent_shell(){
	return static_cast<MGShell*>(parent_complex());
}
const MGShell* MGFace::parent_shell()const{
	return static_cast<const MGShell*>(parent_complex());
}

//Obtain perimeter boundadary loop's curve representation.
//Returned are curves of perimeter boundaries, do not contain perimeter
//of the surface.
std::vector<UniqueCurve> MGFace::PBloop_curves() const{
	std::vector<UniqueCurve> crvs;
	int n=number_of_perimeter_boundaries();
	for(int i=0; i<n; i++){
		std::vector<UniqueCurve> crvsi=loop(i)->curves();
		std::move(crvsi.begin(), crvsi.end(), std::back_inserter(crvs));
	}
	return crvs;
}

//Obtain perimeter's parameter range.
//Let rvec be std::vector<MGInterval> of the fucntion's output,
//then rvec[j] is j-th parameter range of iperi-th perimeter. 
std::vector<MGInterval> MGFace::perimeter_param_range(int iperi) const{
	assert(iperi<4);	//

	const MGSurface* srf=surface();
	int is_v=iperi%2;//if 0, u range perimeter. if 1, v range perimeter.
	MGBox prange=srf->param_range();
	std::vector<MGInterval> rvecw(1,prange.ref(is_v)), rvec;
	if(no_outer_boundaries()) return rvecw;

	const MGLoop* lp=loop(int(0)).get();
	int id1,id2;
	int id1s,id2s;
	if(!lp->is_perimeter_boundary()) return rvec;
	lp->both_end_on_perimeter(id1s,id2s); id1=int(id1s); id2=int(id2s);

	double u0=(prange.ref(0)).low_point(), u1=(prange.ref(0)).high_point();
	double v0=(prange.ref(1)).low_point(), v1=(prange.ref(1)).high_point();
	double ts_corner[4]={u0,v0,u1,v1}, te_corner[4]={u1,v1,u0,v0};

///////There exist perimeters or outer boundaries.

	int i;		//i: loop id.
	MGPosition uv,uve;	//uve: last point of loops.
	MGPosition uv1=lp->start_point(),uv2=lp->end_point();
	double t1=uv1[id1%2], t2=uv2[id2%2];
	if(id1>=2) t1*=-1.; if(id2>=2) t2*=-1.;
	bool reverse_on_1peri=(id1==id2 && t2<t1);

	int idend=id2;
	double ts,te;
	if((id2==iperi && id1!=iperi) || (reverse_on_1peri&&id1==iperi)) ts=uv2(is_v);
	else ts=ts_corner[iperi];
	if(id1>id2 || reverse_on_1peri){
		if(iperi>id1 || iperi<id2) return rvec;
		//Now  id2<= iperi <=id1.
		if(iperi==id1) uve=lp->start_point();
		else uve=srf->perimeter_uv(iperi,te_corner[iperi]);
		i=1;	//To process next loop.
	}else{
		if(iperi>id1 && iperi<id2) return rvec;
		uve=srf->perimeter_uv(iperi,te_corner[iperi]);
		i=0;	//To process first loop.
		idend=-1;
	}

	int n=number_of_perimeter_boundaries();
	while(idend<iperi && i<n){
		lp=loop(i++).get();
		lp->both_end_on_perimeter(id1s,id2s); id1=int(id1s); id2=int(id2s);
		if(iperi<id1) return rvecw;
		if(id1<iperi){
			if(id2>iperi) return rvec;
		}else{//id1=iperi.
			te=(lp->start_point()).ref(is_v);
			MGInterval rng(MGInterval(ts),te);
			rvec.push_back(rng);
			if(id2>iperi) return rvec;
		}
		ts=(lp->end_point()).ref(is_v);//Here id2==iperi.
		idend=id2;
	}
	if(i>=n){//loop ended while searching id2=iperi.
		if(id2<iperi) return rvecw;
		te=uve(is_v);
		MGInterval rng(MGInterval(ts),te);
		rvec.push_back(rng);
		return rvec;
	}

//Now lp is active and ( id1<iperi and id2=iperi) or, (id1=iperi, id2=iperi).
	while(id2==iperi && i<n){
		lp=loop(i++).get();
		lp->both_end_on_perimeter(id1s,id2s); id1=int(id1s); id2=int(id2s);
		if(id1>iperi){
			te=te_corner[iperi];
			MGInterval rng(MGInterval(ts),te);
			rvec.push_back(rng);
			return rvec;
		}else{//id1==iperi.
			te=(lp->start_point()).ref(is_v);
			MGInterval rng(MGInterval(ts),te);
			rvec.push_back(rng);
			ts=(lp->end_point()).ref(is_v);
		}
	}
//Now id2>iperi or i>=n.
//(i-1) indicates last processed loop.
	if((i>=n && id2==iperi)||(id1!=iperi)){
		te=uve.ref(is_v);
		MGInterval rng(MGInterval(ts),te);
		rvec.push_back(rng);
	}
	return rvec;
}

//Sort boundary occurreces in m_boundaries.
//Sorting is done according to operator< of MGComplex.
//parameter space box will be set.
void MGFace::sort_boundaries(){
	//Sort boundary.
	auto compareLoop = [](const UniqueLoop& b1, const UniqueLoop& b2){return (*b1)<(*b2); };
	std::sort(m_boundaries.begin(), m_boundaries.end(), compareLoop);
	m_box_param.set_null();
}

//Trim the face by the projection of a curve along a vector direction.
//If mgNULL_VEC as direction is specified, surface normal projection
//will be employed.
//When this face is divided into two faces, the left face of the crv
//is the face selected.
//
//When the projected curve on the face is connected to already existed
//boundary curve, no new boundary is generated, and inserted(connected)
//to the old boundary. However new projection is floating, that is, not
//connected to any old boundaries, or this boundary is the first boundary,
//new boundary is generated.
//Function's return value is error code:
//0= normal return
//(not error, this includes the case of inactive loop generation).
//1= input pcrv includes a part that is outside surface parameter range.
//2= tried to generate outer boundary loop inside perimeter boudary.
//3= tried to generate inner boundary loop that incudes active loop inside.
//4= tried to generate perimeter boudary loop that inactivates perimeter
//   boundary loops existed.
//5= tried to generate a loop outside the face.
//6= No projected curve by the projection.
//Other: Error code of surface projection.
int MGFace::trim_projection(const MGCurve& crv, const MGVector& direction){
	MGSurface* surf=surface();
	std::vector<UniqueCurve> uvs, worlds;
	int n=surf->project(crv,uvs,worlds,direction);
	if(n<0)
		return n;
	if(n==0)
		return 6;
	if(n==1)
		return trim(*uvs[0]);

	int error;
	MGFace new_face=*this;
	int i=n-1;
	for(; i>=0; i--){
		if(error=new_face.trim(*uvs[i]))
			break;
	}
	if(error==5 && i!=n-1)
		error=0;
	if(!error)
		*this=std::move(new_face);
	return error;
}

//Build a loop of this face, given a closed curve crv on this face. Although crv
//is generally a MGCompositeCurve, this may be not the case. Returned MGLoop is not
//added into this face as a boundary. User must add it after the direction is adjusted.
//That is, the output loop can be an outer or inner loop.
std::unique_ptr<MGLoop> MGFace::build_loop(
	const MGCurve& crv	//curve of world coordinates on this face.
)const{
	std::unique_ptr<MGLoop> loop(new MGLoop);
	const MGSurface* srf=surface();
	int nperi=srf->perimeter_num();

	std::vector<double> pspan[2]; int peri_num[2];
	std::vector<double>& pspan0=pspan[0];
	std::vector<double>& pspan1=pspan[1];
		//To store equal peri for the i-th curve if MGCompositeCurve.
		//(pspan[i], peri_num[i]) is one pair.

	MGCompositeCurve tempCcrv;
	const MGCompositeCurve* ccrv=dynamic_cast<const MGCompositeCurve*>(&crv);
	if(ccrv){
		int nComPeri=srf->getPerimeterCommon(ccrv->curve(0),pspan,peri_num);
		if(nComPeri>=2){
			tempCcrv=crv;
			MGCurve* firstC=tempCcrv.release_front();
			tempCcrv.connect_to_end(firstC);
			ccrv=&tempCcrv;
		}
	}else{
		const MGLBRep* lb=dynamic_cast<const MGLBRep*>(&crv);
		if(lb && lb->order()==2){
			//Devide the curve at C0 continuous point.
			std::vector<UniqueCurve> lbs;
			int ncurves=lb->divide_multi(lbs);
			for(int i=0; i<ncurves; i++){
				std::unique_ptr<MGCurve> lbi(lbs[i].release());
				if(lbi->start_point() != lbi->end_point())
					tempCcrv.connect_to_end(lbi.release());
			}
		}else{
			tempCcrv.connect_to_end(crv.clone());
		}
		ccrv=&tempCcrv;
	}

	double terror=crv.param_span()*MGTolerance::rc_zero();
	MGCompositeCurve::const_iterator icrv=ccrv->begin(), iecrv=ccrv->end();
	double tLast=(*icrv)->param_s();
	for(; icrv!=iecrv; icrv++){
		MGCurve& curvi=**icrv;

		int nComPeri=srf->getPerimeterCommon(curvi,pspan,peri_num);
		if(nComPeri>=2){
			if(loop->number_of_edges()){//If this is not the 1st curve.
				MGPosition uv00=srf->perimeter_uv(peri_num[0],pspan0[2]);
				MGPosition uv10=srf->perimeter_uv(peri_num[1],pspan1[2]);
				MGPosition uvE=loop->end_point();
				if((uvE-uv10).len()<(uvE-uv00).len()){
					pspan0=pspan1; peri_num[0]=peri_num[1];
				}

			}
		}
		loop->append_edge_from_crvWorld(*srf,curvi,tLast,terror,pspan[0],peri_num[0]);
	}

	loop->make_close();
	return loop;
}

//Trim the face giving parameter curve pcrv of this face.
//pcrv has a direction. That is, when this face is divided into two faces,
//left part face of the crv is the face selected.
//
//When the curve pcrv on the face is connected to already existed
//boundary curve, no new boundary is generated, and inserted(connected)
//to the old boundary. However, pcrv is floating, that is, not
//connected to any old boundaries, or this boundary is the first boundary,
//new boundary is generated.
//Function's return value is error code:
//0= normal return
//(not error, this includes the case of inactive loop generation).
//2= tried to generate outer boundary loop inside perimeter boudary.
//3= tried to generate inner boundary loop that incudes active loop inside.
//4= tried to generate perimeter boudary loop that inactivates perimeter
//   boundary loops existed.
//5= tried to generate a loop outside the face.
int MGFace::trim(
	const MGCurve& pcrv	//parameter(u,v) space curve of the face.
){
	MGEdge* e=new MGEdge(pcrv.clone());
	MGLoop new_loop(e);
	if((new_loop.start_point()- new_loop.end_point()).len()<= parameter_error())
		new_loop.make_close();
	return trim(std::move(new_loop));
}

//Trim the face giving a loop new_loop that does not have the parent face.
//new_loop must be parrameter representaion of this face.
//
//Function's return value is error code:
//0= normal return
//(not error, this includes the case of inactive loop generation).
//2= tried to generate outer boundary loop inside perimeter boudary.
//3= tried to generate inner boundary loop that incudes active loop inside.
//4= tried to generate perimeter boudary loop that inactivates perimeter
//   boundary loops existed.
//5= tried to generate a loop outside the face.
int MGFace::trim(
	MGLoop&& new_loop_in	//loop
){
	MGFace new_face=*this;	//Updation is done to this new_face.
	int j,j1,m;
	MGLoop* new_loop;
	if(new_loop_in.closed()){//When closed curve. Generate closed edge.
		new_loop=new MGLoop(std::move(new_loop_in));
		new_face.prepend_boundary(new_loop);
		m=new_face.number_of_boundaries();
		bool merged=false;
		for(j=m-1; j>0; j--){
			if(new_loop->merge_trim(*new_face.loop(j))){
				new_face.erase_boundary(j); merged=true;
			}
		}
		if(!merged && new_loop->is_outer_boundary()){
			if(!in_range(new_loop->start_point()))
				return 5;

			//New is outer, so erase perimeter or legacy outer boundary.
			new_face.free_boundary(new_loop);
			m=new_face.number_of_boundaries();
			for(j=m-1; j>=0; j--){
				UniqueLoop& lp2=new_face.loop(j);
				if(lp2->is_outer_boundary() || lp2->is_perimeter_boundary())
					new_face.erase_boundary(j);
			}
			new_face.prepend_boundary(new_loop);
		}
	}else{//When not closed curve.
		new_loop=0;
		m=new_face.number_of_boundaries();
		for(j=0; j<m; j++){
			UniqueLoop& lp1=new_face.loop(j);
			bool merged=lp1->merge_trim(new_loop_in);
			if(merged){
				new_loop=lp1.get();
				for(j1=m-1; j1>j; j1--){
					if(new_loop->merge_trim(*new_face.loop(j1))){
						new_face.erase_boundary(j1);
						new_face.m_box_param.set_null();
					}
				}
				break;
			}
		}
		if(!new_loop){
			MGPosition uv=new_loop_in.start_point();
			if(!in_range(uv))
				return 5;

			new_loop=new MGLoop(new_loop_in);
			new_face.append_boundary(new_loop);
		}	
	}

	//Sort the loop.
	//perimeter loop, ..., closed loop, ...., not active loop, ...
	new_face.sort_boundaries();

	if(new_loop->active()){
	//When new_loop is active, check the validity.
		if(new_loop->closed()){
			if(new_loop->is_outer_boundary()){
				if(new_face.number_of_perimeter_boundaries())
					return 2;
				//When new loop is outer boundary loop,
				//there must not exist perimeter boundaries.
			}else{
				m=new_face.number_of_loops();
				for(j=0; j<m; j++){
					UniqueLoop& lp=new_face.loop(j);
					if(lp.get()!=new_loop){
						if(!new_loop->inside(lp->start_point()))
							return 3;
						//When new loop is inner boundary loop, there must not
						//exist active boundaries inside the new loop.
					}
				}
			}
		}else{
			m=new_face.number_of_perimeter_boundaries();
			int ids,ide;
			new_loop->both_end_on_perimeter(ids,ide);
			MGPosition ps=new_loop->start_point(), pe=new_loop->end_point();
			double tps=ps.ref(ids%2), tpe=pe.ref(ide%2);
			if(ids>=2) tps=-tps; if(ide>=2) tpe=-tpe;
			MGLPoint ts(ids,tps), te(ide,tpe);
			//When new loop is perimeter boudary, all of the other
			//perimeter boundaries must lie inside the perimeter loop.
			if(te<ts){
				for(j=0; j<m; j++){
					UniqueLoop& lp=new_face.loop(j);
					if(lp.get()!=new_loop){
						MGPosition p=lp->start_point();
						lp->on_perimeter_start(ids);
						tps=p.ref(ids%2); if(ids>=2) tps=-tps;
						MGLPoint tis(ids,tps);
						if(tis<te || tis>ts)
							return 4;
					}
				}
			}else{//When ts<te.
				for(j=0; j<m; j++){
					UniqueLoop& lp=new_face.loop(j);
					if(lp.get()!=new_loop){
						MGPosition p=lp->start_point();
						lp->on_perimeter_start(ids);
						tps=p.ref(ids%2); if(ids>=2) tps=-tps;
						MGLPoint tis(ids,tps);
						if(ts<tis && tis<te)
							return 4;
					}
				}
			}
		}
	}
	*this = new_face;
	return 0;
}

//Shrink the base surface of this face to the part limitted by the parameter range of uvbx.
//New parameter range uvbx2 is so determined that uvbx2 is the smallest
//box tha includes uvbx, and all of the u or v values of uvbx2 is one of 
//the values of u or v knots of the surface knotvector.
//uvbx(0) is the parameter (us,ue) and uvbx(1) is (vs,ve).
//That is u range is from us to ue , and so on.
void MGFace::shrink_base_surface_to_knot(
	int multiple	//Indicates if start and end knot multiplicities
					//are necessary. =0:unnecessary, !=0:necessary.
){
	MGSurface* srf=surface();
	if(srf){
		const MGBox& range=box_param();
		srf->shrink_to_knot(range,multiple);
	}
}

//split this fsurface at the parameter param.
void MGFace::split(
	double param,//parameter value of this fsurface. if is_u is true, param is u-value,
				//else v-value.
	bool is_u,	//indicates if param is u or v of the surface parameter (u,v).
	std::vector<UniqueFSurface>& surfaces///<splitted surfaces will be output.
)const{
	int kcod0= is_u ? 0:1;
	const MGBox& bx=box_param();
	const MGInterval& rng=bx[kcod0];
	if(!rng.includes(param))
		return;

	if(hasOuterBoundaryLoop()){
		std::vector<double> pvec=isect1D_with_boundaries(param,kcod0);
		int nm1=(int)(pvec.size()-1);
		if(nm1<0)
			return;
		int kcod1=(kcod0+1)%2;

		MGPosition uv(2);
		uv(kcod0)=param;

		std::vector<UniqueLoop> networks;
		for(int i=0;i<nm1; i++){
			double t0=pvec[i], t1=pvec[i+1];
			double tmid=t0+t1; tmid*=0.5;
			uv(kcod1)=tmid;
			if(in_range(uv)){
				MGPosition uv0(uv), uv1(uv);
				uv0(kcod1)=t0;
				uv1(kcod1)=t1;
				networks.emplace_back(new MGLoop(new MGEdge(new MGStraight(uv1,uv0))));
			}
		}
		std::vector<UniqueFace> faces;
		split(networks,faces);
		int n=(int)faces.size();
		for(int i=0; i<n; i++){
			surfaces.emplace_back(faces[i].release());
		}
	}else{
		const MGSurface* srf=surface();
		if(!srf)
			return;
		int iperi;
		if(srf->on_a_perimeter2(is_u,param,iperi))//if on perimeters.
			return;
		MGFace f2(*this);
		f2.make_outer_boundary();
		f2.split(param,is_u,surfaces);
	}
}

///Rebuild this face. Rebuild means:
/// 1) Rebuild the base surface.
/// 2) Rebuild the parameter edges.
/// 3) Rebuild the binder edges from the rebuilt parameter edges.
/// This rebuild does not change the number of edges. All of the old edges
/// does exist after the rebuild.
///Base surface rebuild is done only when the base surface is MGSBRep, or MGRSBRep.
std::unique_ptr<MGFace> MGFace::rebuild(
	int how_rebuild,
		//intdicates how rebuild be done.
		// =0: no base surface approximation(only parameter change)
		// =1: If base surface is MGRSBRep, reconstruct it with new knot configuration
		//     again as rational spline(MGRSBRep). If MGSBRep, reconstruct with new knot
		//     configuration as MGSBRep.
		// =2: approximated by non-rational spline(MGSBRep) with new knot configuration,
		//     even base surface is MGRSBRep. If the base surface is MGSBRep, same as =1.
	int parameter_normalization,
		//Indicates how the parameter normalization be done:
		//=0: no surface parameter normalization.
		//=1: normalize to u_range=(0., 1.), and v_range=(0.,1.);
		//=2: normalize to make the average length of the 1st derivative along u and v 
		//    of the base surface is as equal to 1. as possible.
		//=3: Specify parameter range in range[4].
	double tol,	///<tolerance allowed for the approximation.
		///When tol<=0., MGTolerance::line_zero() will be employed.
	int* order,///<order of the new MGSBRep, >=4 is recomended.
		///order[0]:u-order, [1]:v-order.
		///When how_rebuild!=2, order is not used.
		///When order=0 is input, order[0]=order[1]=4 are assumed.
	double* range///valid only when parameter_normalization=3,
		/// and range[]={umin, umax, vmin, vmax}.
		///When parameter_normalization=3 and range=0, parameter_normalization=2
		///is assumed.
)const{
	assert(parameter_normalization<=3);
	if(parameter_normalization==3 && range==0)
		parameter_normalization=2;

	if(tol<0.)
		tol=MGTolerance::line_zero();

	//1. Surface rebuild.
	const MGSurface* srfOld=get_surface_pointer();
	double usOld=srfOld->param_s_u();
	double ueOld=srfOld->param_e_u();
	double vsOld=srfOld->param_s_v();
	double veOld=srfOld->param_e_v();

	UniqueSurface surfNew;
	const MGSBRep* srf=dynamic_cast<const MGSBRep*>(srfOld);
	if(srf){
		surfNew =srf->rebuild(how_rebuild,parameter_normalization,tol,order);
	}else{
		const MGRSBRep* rsrf=dynamic_cast<const MGRSBRep*>(srfOld);
		if(rsrf){
			surfNew =rsrf->rebuild(how_rebuild,parameter_normalization,tol,order);
		}else{
			parameter_normalization=0;
			surfNew.reset(srfOld->clone());
		}
	}
	if(parameter_normalization==3){
		surfNew->change_range(1,range[0], range[1]);
		surfNew->change_range(0,range[2], range[3]);
	}

	//Tolerance adjustment for parameter line rebuild.
	const MGBox& bx=box();
	double relTol=tol/bx.length();

	double ue=surfNew->param_e_u(), us=surfNew->param_s_u();
	double uspan=ue-us;
	double ve=surfNew->param_e_v(), vs=surfNew->param_s_v();
	double vspan=ve-vs;
	double lenParam=sqrt(uspan*uspan+vspan*vspan);
	double paraError=lenParam*relTol;
	mgTolSetWCZero wczeroSet(paraError);//Set&save the error.
	MGVector originOld(usOld, vsOld);
	double uscale=uspan/(ueOld-usOld), vscale=vspan/(veOld-vsOld);
	MGMatrix M(2,uscale);	M(1,1)=vscale;

	//2. Loops Rebuild.
	std::vector<UniqueLoop> newloops;
	int nloops=number_of_boundaries();
	for(int i=0; i<nloops; i++){
		const UniqueLoop& loopi=loop(i);
		int nedges=loopi->number_of_edges();

		//2.1 Rebuild loop i.
		MGLoop* newloopi=new MGLoop;
		for(int j=0; j<nedges; j++){

			//2.2 Rebuild edge j.
			MGEdge* ej=new MGEdge(*(loopi->edge(j)));//Delete if binder edges exist.
			MGCurve& cj=*(ej->base_curve());
			if(parameter_normalization){
				cj-=originOld;
				cj*=M;
			}
			MGLBRep* lbj=dynamic_cast<MGLBRep*>(&cj);
			MGRLBRep* rlbj=0;
			if(!lbj){
				rlbj=dynamic_cast<MGRLBRep*>(&cj);
			}
			if(lbj || rlbj){
				MGTrimmedCurve cjt=ej->trimmed_curve();
				std::unique_ptr<MGCurve> newCurve=
					cjt.rebuild(1,parameter_normalization,paraError,4);
				ej->set_start(newCurve->param_s());
				ej->set_end(newCurve->param_e());
				ej->set_extent(std::move(newCurve));
			}
			newloopi->append(ej);
		}
		newloopi->make_close();
		newloops.emplace_back(newloopi);
	}
	UniqueFace fnew(new MGFace(std::move(surfNew),std::move(newloops)));
	fnew->copy_appearance(*this);
	return fnew;
}
