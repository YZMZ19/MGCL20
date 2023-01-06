/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Box.h"
#include "mg/BPointSeq.h"
#include "mg/Position.h"
#include "mg/Position_list.h"
#include "mg/Transf.h"
#include "mg/Curve.h"
#include "mg/Ellipse.h"
#include "mg/LBRep.h"
#include "mg/SurfCurve.h"
#include "mg/CompositeCurve.h"
#include "mg/CParam_list.h"
#include "mg/CCisects.h"
#include "mg/CSisects.h"
#include "mg/SSisects.h"
#include "mg/Straight.h"
#include "mg/Surface.h"
#include "mg/Plane.h"
#include "mg/SBRep.h"
#include "mg/SBRepTP.h"
#include "mg/Tolerance.h"
using namespace std;

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// Implementation of MGSurface.
// MGSurface is an abstract class of 3D surface.
// Surface is represented using two parameter u and v.

//*******Intersection.************

//Intersection of Surface and curve.
MGCSisects MGSurface::intersect(const MGCurve& crv) const{
	assert(dynamic_cast<const MGPlane*>(this)==0);//This should not be plane.

	MGCSisects list(&crv,this);
	MGBox bx=box(); bx&=crv.box();
	if(bx.empty())
		return list;
	double tss=crv.param_s(), tee=crv.param_e(); 
	double t,tm,u,v; MGPosition uv(2);
	//MGPosition PS=crv.eval(tss);
	//MGPosition PE=crv.eval(tee);
	//if(on(PS,uv))
	//	list.append(PS,tss,uv);
	//if(on(PE,uv))
	//	list.append(PE,tee,uv);

	MGVector ctanS=crv.eval(tss,1);
	MGVector ctanE=crv.eval(tee,1);
	const double log2=0.69315;//log(2.);
	int nt1=crv.intersect_dnum();
	int nu1=bdim_u()/2;if(nu1<=2) nu1=3;
	int nv1=bdim_v()/2;if(nv1<=2) nv1=3;
	//Compute minimum length
	double Cspan_length=((ctanS.len()+ctanE.len())*.5*(tee-tss));
	double us=param_s_u(), ue=param_e_u();
	double um=(us+ue)*.5;
	double vs=param_s_v(), ve=param_e_v();
	double vm=(vs+ve)*.5;
	double dsdu=(eval(us,vs,1,0).len()+eval(ue,vs,1,0).len())*.5;
	dsdu+=(eval(us,vm,1,0).len()+eval(ue,vm,1,0).len())*.5;
	dsdu+=(eval(us,ve,1,0).len()+eval(ue,ve,1,0).len())*.5;
	dsdu/=3.;
	double Uspan_length=dsdu*(ue-us);

	double dsdv=(eval(us,vs,0,1).len()+eval(us,ve,0,1).len())*.5;
	dsdv+=(eval(um,vs,0,1).len()+eval(um,ve,0,1).len())*.5;
	dsdv+=(eval(ue,vs,0,1).len()+eval(ue,ve,0,1).len())*.5;
	dsdv/=3.;
	double Vspan_length=dsdv*(ve-vs);

	double one_span=Cspan_length/double(nt1);
	double U1span_length=Uspan_length/double(nu1);
	double V1span_length=Vspan_length/double(nv1);
	if(one_span>U1span_length){
		one_span=U1span_length;
		if(one_span>V1span_length)
			one_span=V1span_length;
	}else{
		if(one_span>V1span_length)
			one_span=V1span_length;
	}
	double min_span=MGTolerance::wc_zero()*500.;
	if(one_span<min_span)
		one_span=min_span;

	nt1=int(Cspan_length/one_span)+1;
	nu1=int(Uspan_length/U1span_length+.5)+1;
	nv1=int(Vspan_length/V1span_length+.5)+1;

	int ntstack=(int)(log(3.*nt1)/log2)+1;
	int nustack=(int)(log(3.*nu1)/log2)+1;
	int nvstack=(int)(log(3.*nv1)/log2)+1;
	int stack_len=ntstack+nustack+nvstack;
//Subdivision before call of isect_guess will be done until
//one intersect_dnum()-th of mean span length.
	double* lstack=new double[stack_len*2];
	std::vector<MGPosition> sstack(stack_len*2);
	int* ndivide=new int[stack_len*3];

//lstack and sstack are actually the arrays of [stack_len][2].
// lstack[i][0] - lstack[i][1] is the line parameter range(min and max).
// sstack[i][0] - sstack[i][1] is surface parameter range.
//ndivide[] indicate how deep subdivision performed, and which should be
//divided, line(0), u of surface(1), or v of surface(2). Maximum number's
//subdivision will be performed.
	int kdiv;
	MGPosition uvs(2), uvm(2), uve(2);

	lstack[0]=tss; lstack[1]=tee; 
	sstack[0]=MGPosition(us,vs); sstack[1]=MGPosition(ue,ve);
	ndivide[0]=nt1; ndivide[1]=nu1; ndivide[2]=nv1;
	int stack_id=1;	// stack_id indicates next available id for stack, i.e.
						// how many data in stack.
	MGCurve* lb1=0; MGSurface* sb1=0;
 
	int id1, id2;
	while(stack_id){
		stack_id--;	// Pop the stack data.
		id1=stack_id*2, id2=stack_id*3;
		double ts= lstack[id1], te=lstack[id1+1];
		uvs=sstack[id1]; uve=sstack[id1+1];
		int nt=ndivide[id2], nu=ndivide[id2+1], nv=ndivide[id2+2];
		delete lb1; lb1=crv.part(ts,te);
		delete sb1; sb1=part(MGBox(uvs,uve));
		//Check of box boundary override-ness.
		while(!(lb1->box()&sb1->box()).empty()){
			if((nt+nu+nv)<5){
				//Compute intersection using isect_guess.
				uvm=(uvs+uve)/2.;
				//if(sb1->isect_guess(*lb1,uvm,ts,uv,t))
				//			list.append(lb1->eval(t),t,uv);
				//if(sb1->isect_guess(*lb1,uvm,te,uv,t))
				//			list.append(lb1->eval(t),t,uv);
				tm=(ts+te)*.5;
				if(sb1->isect_guess(*lb1,uvm,tm,uv,t))
							list.append(lb1->eval(t),t,uv);
				break;			//To pop up stacked data.
			}
			//Subdivide and push to stack.
			if(nt>=nu){
				if(nt>=nv) kdiv=0;
				else       kdiv=2;
			}else if(nu>=nv) kdiv=1;
			else           kdiv=2;
			switch(kdiv){
			case(0):
			//Subdivide crv curve.
				tm=(ts+te)/2.;
				lstack[id1]=tm;  lstack[id1+1]=te;
				sstack[id1]=uvs; sstack[id1+1]=uve;
				nt/=2;
				ndivide[id2]=nt; ndivide[id2+1]=nu; ndivide[id2+2]=nv;
				te=tm;
				delete lb1; lb1=crv.part(ts,te);
				break;
			case(1):
			//Sbudivide this surface for u-direction.
				u=(uvs[0]+uve[0])/2.;
				uvm(0)=u;
				uvm(1)=uvs[1];
				lstack[id1]=ts;  lstack[id1+1]=te;
				sstack[id1]=uvm; sstack[id1+1]=uve;
				nu/=2;
				ndivide[id2]=nt; ndivide[id2+1]=nu; ndivide[id2+2]=nv;
				uve(0)=u;
				delete sb1; sb1=part(MGBox(uvs,uve));
				break;
			case(2):
			//Sbudivide this surface for v-direction.
				uvm(0)=uvs[0];
				v=(uvs[1]+uve[1])/2.;
				uvm(1)=v;
				lstack[id1]=ts;  lstack[id1+1]=te;
				sstack[id1]=uvm; sstack[id1+1]=uve;
				nv/=2;
				ndivide[id2]=nt; ndivide[id2+1]=nu; ndivide[id2+2]=nv;
				uve(1)=v;
				delete sb1; sb1=part(MGBox(uvs,uve));
				break;
			}
			stack_id++;
			id1=stack_id*2, id2=stack_id*3;
			//To check current spans' override-ness.
		}
	}		
	delete lb1; delete sb1; delete[] lstack; delete[] ndivide;
	return list;
}

//Intersection of Surface and ellipse.
MGCSisects MGSurface::intersect(const MGEllipse& el) const{
	MGCSisects list(&el, this);
	MGBox bx=box(); bx&=el.box();
	if(bx.empty()) return list;

	MGPlane el_plane(el.normal(), el.center()); //Plane that includes el.
	MGSSisects llist=isect(el_plane);	//Compute intersetion line of
											//the above plane and this.
	int nl=llist.entries();
	for(int i=0; i<nl; i++){
		std::unique_ptr<MGSSisect> is=llist.removeFirst();
		//Compute intersection point of llist line and ellipse.
		MGCCisects plist=el.isect(is->line());
		int np=plist.entries();
		for(int j=0; j<np; j++){
			std::unique_ptr<MGCCisect> ip=plist.removeFirst();
			MGPosition p=ip->point();
			list.append(MGCSisect(p,ip->param1(),param(p)));
		}
	}
	return list;
}

//Compute intersections with MGLBRep lb that does not have C0 continuity in it.
MGCSisects MGSurface::isect_withC1LB(const MGLBRep& lb)const{
	return intersect(lb);
}

//isect with SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCSisects MGSurface::isect_with_noCompoSC(const MGSurfCurve& scrv)const{
	return intersect(scrv);
}

void isectSlStack(std::stack<MGBox>& Rstack, MGBox& uvrng, int i){
	MGBox uvrng2=uvrng;		//Save the original.
	MGInterval& xrng=uvrng(i);
	double mid=xrng.mid_point();
	xrng.set_low_point(mid);
	uvrng2(i).set_high_point(mid);
	Rstack.push(uvrng2);
}

//Intersection of Surface and a straight line.
MGCSisects MGSurface::isectSl(
	const MGStraight& sl,
	const MGBox& uvbox //indicates if this surface is restrictied to the parameter
					//range of uvbox. If uvbox.is_null(), no restriction.
) const{
	MGCSisects list(&sl,this);
	const MGBox& sbx=box();
	if(!sbx.crossing(sl))
		return list;

	MGUnit_vector SLD=sl.direction();
	MGMatrix mat; mat.to_axis(SLD,2);	//Matrix to transform SLD to be z axis.

	double u0,u1,v0,v1;
	if(uvbox.is_null()){
		u0=param_s_u(); u1=param_e_u();
		v0=param_s_v(); v1=param_e_v();
	}else{
		u0=uvbox[0].low_point(); u1=uvbox[0].high_point();
		v0=uvbox[1].low_point(); v1=uvbox[1].high_point();
	}
	double u=(u0+u1)*.5, v=(v0+v1)*.5;
	MGPosition Puv=eval(u,v);
	MGVector PN=sl.eval(sl.closest(Puv));
		//PN is the nearest point of the sl from the center of this surface,
		//and is the vector to translate sl to pass through the origin.
	std::unique_ptr<MGSurface> sf2D(copy_surface());
	(*sf2D)-=PN; (*sf2D)*=mat;
	sf2D->change_dimension(2);
	//sf2D is the 2D surface that is transformed from the original this
	//surface so that the straight line sl is seen as a point of origin(0.,.0.).

//	double uspan=(u1-u0)/bdim_u()*.5, vspan=(v1-v0)/bdim_v()*.5;
	double uspan=(u1-u0)/bdim_u(), vspan=(v1-v0)/bdim_v();
	//Subdividion of the surfae sf2D will halt when the subdivided
	//parameter range becomes smaller than these spans.

	std::stack<MGBox> Rstack;
	Rstack.push(MGBox(MGInterval(u0,u1), MGInterval(v0,v1)));
	bool uDir=false;//To control that the subdivision is done u a v direction
		//alternately.
	MGBox Wbx(sdim());
	while(!Rstack.empty()){

	uDir=!uDir;
	MGBox& uvrng=Rstack.top();
	Wbx=sf2D->box_limitted(uvrng);
	if(Wbx.includes_origin()){
	//If the origin is included in the box,
	//subdivide or invoke isect_guess_straight.
		MGInterval& urng=uvrng(0);
		MGInterval& vrng=uvrng(1);
		if(uDir){
			if(urng.length()>uspan){
				isectSlStack(Rstack,uvrng,0);
				continue;
			}else if(vrng.length()>vspan){
				isectSlStack(Rstack,uvrng,1);
				continue;
			}
		}else{
			if(vrng.length()>vspan){
				isectSlStack(Rstack,uvrng,1);
				continue;
			}else if(urng.length()>uspan){
				isectSlStack(Rstack,uvrng,0);
				continue;
			}
		}
	//Now both u and v range are small enough to use isect_guess_straight.
		MGPosition uv(urng.mid_point(), vrng.mid_point());
		double t=sl.closest(eval(uv));
		if(isect_guess_straight(sl,t,uv,t,uv)) list.append(sl.eval(t),t,uv);
	}
	Rstack.pop();

	}

	return list;
}

#define MAXLOOP 6
#define MAXLOOP2 16
// "isect_guess" computes one intersection point of surface and a curve,
// given initail guess parameter values of surface and curve.
//Function's return value is 1(true) if i.p. obtained, and 0(false) if not.
int MGFSurface::isect_guess(
	const MGCurve& crv,		//Curve
	const MGPosition& uvi,	//Input initial guess parameter value
						// of the i.p. of the surface. 
	double ti,			//Input initial guess parameter value of the line.
	MGPosition& uv,		// Output parameter value obtained. 
	double& t)			// Output parameter value obtained. 
const{
	const MGStraight* sl=dynamic_cast<const MGStraight*>(&crv);
	if(sl) return isect_guess_straight(*sl,ti,uvi,t,uv);
	const MGCompositeCurve* ccrv=dynamic_cast<const MGCompositeCurve*>(&crv);
	if(ccrv) return isect_guess_composite(*ccrv,uvi,ti,uv,t);

	MGPosition A,P;
	MGVector AP,Su,Sv,SuSv,Sud,Svd,Lt,PQ;
	double du,dv,dt,Lt_sqr,AP_sqr;
	MGPlane plane; MGStraight line; MGCSisect is;
	double error_sqr=MGTolerance::wc_zero_sqr();//*.5;
	double u0=param_s_u(), u1=param_e_u();
	double v0=param_s_v(), v1=param_e_v();
	double t0=crv.param_s(), t1=crv.param_e();
	double SuSud,SvSvd;

	int loop=0, ulow=0, uhigh=0, vlow=0, vhigh=0, tlow=0, thigh=0;
	t=ti; double u=uvi.ref(0), v=uvi.ref(1);
	int ret_code=0;//Return code.
	while(loop++<MAXLOOP2 && ulow<MAXLOOP && uhigh<MAXLOOP && vlow<MAXLOOP && vhigh<MAXLOOP
                  && tlow<MAXLOOP && thigh<MAXLOOP){
		A=eval(u,v); P=crv.eval(t); AP=P-A;
	//A is guess point of the surface. P is guess point of the line.
		AP_sqr=AP%AP;
		if(AP_sqr<=error_sqr){
			uv=MGPosition(u,v);
			ret_code=1;
			break;
		}

		Su=eval(u,v,1,0); Sv=eval(u,v,0,1); SuSv=Su*Sv;
		Sud=Sv*SuSv; Svd=SuSv*Su;
		SuSud=Su%Sud; if(MGMZero(SuSud)) break;
		SvSvd=Sv%Svd; if(MGMZero(SvSvd)) break;

		plane=MGPlane(Su,Sv,MGDefault::origin());
		Lt=crv.eval(t,1);
		Lt_sqr=Lt%Lt; if(MGMZero(Lt_sqr)) break;
		line=MGStraight(MGSTRAIGHT_TYPE::MGSTRAIGHT_UNLIMIT,Lt,MGPosition(AP));
		if(line.relation(plane,is)!= MGPSRELATION::MGPSREL_ISECT) break;
		const MGPosition& Q=is.point();	//Q is the intersection point of osculating plane
	// of the surface at uv(A) and tangent straight line of the line at t(P).
		PQ=Q-AP;
		dt=(PQ%Lt)/Lt_sqr;
		du=(Sud%Q)/SuSud; dv=(Svd%Q)/SvSvd;
		
	// Update uv and t.
		u+=du; v+=dv; t+=dt;
		if(u<u0)     {ulow+=1; uhigh=0; u=u0;}
		else if(u>u1){ulow=0; uhigh+=1; u=u1;}
		else         {ulow=0; uhigh=0;}
	
		if(v<v0)     {vlow+=1; vhigh=0; v=v0;}
		else if(v>v1){vlow=0; vhigh+=1; v=v1;}
		else         {vlow=0; vhigh=0;}
	
		if(t<t0)     {tlow+=1; thigh=0; t=t0;}
		else if(t>t1){tlow=0; thigh+=1; t=t1;}
		else         {tlow=0; thigh=0;}
	}	

	//if(AP_sqr<=error_sqr*5.){uv=MGPosition(u,v); return 1;}
	if(ret_code==0){
		if(ulow>=MAXLOOP || uhigh>=MAXLOOP || vlow>=MAXLOOP || vhigh>=MAXLOOP
			|| tlow>=MAXLOOP || thigh>=MAXLOOP){
			double lzero2=MGTolerance::line_zero()*.5;
			lzero2*=lzero2;
			if(AP_sqr<=lzero2){
				uv=MGPosition(u,v);
				ret_code=1;
			}
		}
	}
	if(ret_code){
		if(!crv.in_range(t))
			ret_code=0;
	}

	return ret_code;
}

// "isect_guess" computes one intersection point of surface and a curve,
// given initail guess parameter values of surface and curve.
//Function's return value is 1(true) if i.p. obtained, and 0(false) if not.
int MGFSurface::isect_guess_composite(
	const MGCompositeCurve& ccrv,	//Curve
	const MGPosition& uvi,//Input initial guess parameter value
						//of the i.p. of the surface. 
	double ti,			//Input initial guess parameter value of the line.
	MGPosition& uv,		//Output parameter value obtained. 
	double& t)			//Output parameter value obtained. 
const{
	int i=ccrv.find(ti);
	const MGCurve& curvei=ccrv.curve(i);
	return isect_guess(curvei,uvi,ti,uv,t);
}

// "isect_guess_straight" computes one intersection point of surface and
//a straight line, given initail guess parameter values of the surface.
int MGFSurface::isect_guess_straight(
	const MGStraight& sl,	//Straight line.
	double ti,			//Initial guess parameter value of sl.
	const MGPosition& uvi,	//Input initial guess parameter value
						// of the i.p. of the surface. 
	double& t,			//Parameter of t will be returned.
	MGPosition& uv		//Surface parameter value obtained (u,v). 
	) const
//Function's return value is 1(true) if i.p. obtained, and 0(false) if not.
{
	MGPosition A,Q;
	MGVector AP,Su,Sv,SuSv,Sud,Svd,PQ,AQ;
	double du,dv,dt ,AP_sqr;
	MGCSisect is;
	double SuSud,SvSvd;

	double error_sqr=MGTolerance::wc_zero_sqr();
	MGPosition P,PS=sl.root_point(); MGPSRELATION rl;
	MGVector Lt=sl.direction(); double Lt_sqr=Lt%Lt;
	double u0=param_s_u(), u1=param_e_u();
	double v0=param_s_v(), v1=param_e_v();

	int loop=0, ulow=0, uhigh=0, vlow=0, vhigh=0, tlow=0, thigh=0;
	uv.resize(2);
	double& u=uv(0); double& v=uv(1);
	u=uvi.ref(0), v=uvi.ref(1); t=ti;
	int ret_code=0;
	while(loop++<16 && ulow<MAXLOOP && uhigh<MAXLOOP && vlow<MAXLOOP && vhigh<MAXLOOP
                  && tlow<MAXLOOP && thigh<MAXLOOP)
	{
		A=eval(u,v); P=PS+Lt*t; AP=P-A;
			//A is guess point of the surface.
		AP_sqr=AP%AP;
		if(AP_sqr<=error_sqr){
			ret_code=1;
			break;
		}

		Su=eval(u,v,1,0); Sv=eval(u,v,0,1); SuSv=Su*Sv;
		Sud=Sv*SuSv; Svd=SuSv*Su;
		SuSud=Su%Sud; if(MGMZero(SuSud)) break;
		SvSvd=Sv%Svd; if(MGMZero(SvSvd)) break;

		rl=sl.relation(MGPlane(Su,Sv,A),is);
		if(rl== MGPSRELATION::MGPSREL_COIN || rl== MGPSRELATION::MGPSREL_PARALLEL) break;
		Q=is.point();
			//Q is the intersection point of osculating plane
			// of the surface at uv(A) and the straight sl.
		PQ=Q-P; dt=(PQ%Lt)/Lt_sqr;
		AQ=Q-A; du=(Sud%AQ)/SuSud; dv=(Svd%AQ)/SvSvd;
		
	// Update uv and t.
		u+=du; v+=dv; t+=dt;
		if(u<u0)     {ulow+=1; uhigh=0; u=u0;}
		else if(u>u1){ulow=0; uhigh+=1; u=u1;}
		else         {ulow=0; uhigh=0;}
		if(v<v0)     {vlow+=1; vhigh=0; v=v0;}
		else if(v>v1){vlow=0; vhigh+=1; v=v1;}
		else         {vlow=0; vhigh=0;}
	}

	if(ret_code==0 && (ulow>=MAXLOOP || uhigh>=MAXLOOP || vlow>=MAXLOOP
		|| vhigh>=MAXLOOP || tlow>=MAXLOOP || thigh>=MAXLOOP)){
		double lzero2=MGTolerance::line_zero()*.5;
		lzero2*=lzero2;
		if(AP_sqr<=lzero2){
			ret_code=1;
		}
	}
	if(ret_code)
		if(!sl.in_range(t))
			ret_code=0;

	return ret_code;
}

//Compute average chord length along u(is_u==false) or v(is_u==true).
//Function's return value is the chord length.
double MGSurface::average_chord_length(
	int is_u,// =0, or 1. Indicates if para is u-value(is_u=1) or v(is_u=0).
	const double para[3],	//three parameter value of srf, start, middle, and end
							//of u parameter value if is_u==1, of v if is_u==0.
	const MGNDDArray& tau	//data points to evaluate on is_u parameter line.
							//tau[.] is v values if is_u==1, u values if 0.
)const{
	double span=0.;
	int other=is_u ? 0:1;
	MGPosition uv(2);
	int n=tau.length();
	for(int j=0; j<3; j++){
		uv(other)=para[j];
		uv(is_u)=tau[0];
		MGPosition Ppre=eval(uv);
		for(int i=1; i<n; i++){
			uv(is_u)=tau[i];
			MGPosition Paft=eval(uv);
			span+=Ppre.distance(Paft);
			Ppre=Paft;
		}
	}
	span/=3.;
	return span;
}

///Construct a tangent plane LBRep of a surface's perimeter.
///Output MGLBRep@tp's knot vector is the same as srf's perimeter's 
///if pcrv is MGLBRep. Let ft(t) is pcrv and g(t) be the output tp ,
///then g(t)=normal at srf(u,v) where (u,v)=f(t).
///The direction of the output is the same as perimeter_curve().
void MGSurface::TPatPerimeter(
	int perimeterNum,
	MGLBRep& tp///<Obtained tp is output.
)const{
	UniqueCurve pcrv = perimeterUVCurve(perimeterNum);
	bool uPeri = (perimeterNum % 2 == 0);
	const MGKnotVector& t = uPeri ? knot_vector_u() : knot_vector_v();
	MGNDDArray tau; tau.buildByKnotVector(t);

	int n = t.bdim();
	MGBPointSeq uvs(n, 2);
	int sd = sdim(); if(sd<3) sd = 3;
	MGBPointSeq normals(n, sd);
	for(int i = 0; i < n; i++){
		MGPosition uv=pcrv->eval(tau[i]);
		MGVector N = unit_normal(uv);
		normals.store_at(i, N);
	}
	int ordr = t.order(); if(ordr<4) ordr = 4;
	if(ordr>n) ordr = n;
	tp.setKnotVector(MGKnotVector(tau,ordr));
	tp.buildByInterpolationWithKTV(tau, normals);
}
