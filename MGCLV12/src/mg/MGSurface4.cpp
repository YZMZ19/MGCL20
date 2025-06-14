/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Box.h"
#include "mg/Position.h"
#include "mg/Position_list.h"
#include "mg/Transf.h"
#include "mg/LBRepEndC.h"
#include "mg/Straight.h"
#include "mg/LBRep.h"
#include "mg/CParam_list.h"
#include "mg/CCisects.h"
#include "mg/CSisects.h"
#include "mg/SSisects.h"
#include "mg/Surface.h"
#include "mg/Plane.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/Tolerance.h"
#include "topo/Face.h"

#include "cskernel/Bkdnp.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
using namespace std;

// Implementation of intersection of MGSurface and MGPlane.

//*******Intersection.************

//Compute intersection points of perimeters of this surface and surf, and
//of perimeters of surf and this surface.
//These intersection points are used to compute surface to surface
//intersection lines.
void MGFSurface::intersect12Boundary(
	const MGFSurface& sf2,		//The second surface.
	MGPosition_list& uvuv_list	//The intersection points will be output.
			//One member of uvuv_list is (u1,v1,u2,v2), where (u1,v1) is
			//a parameter of this surface and (u2,v2) is a parameter of 
			//surf.
)const{
	mgTolSetWCZero wczeroSet(MGTolerance::line_zero()*.5);//Set&save the error.
	int numi1=isect_boundary(sf2,uvuv_list,0);
	int numi2=sf2.isect_boundary(*this,uvuv_list,2);

	if(uvuv_list.size()>2){
		//Sort the ip's.
		int id=0;
		if(numi1>numi2) id=2;
		uvuv_list.sort_uv_space(id);
	}
	return;
}

//isect_startPlanePt compute an array of parameter value pairs of this surf and pl2
//for one intersection line of a surface and a plane,
// given starting intersetion point uvuv((u1,v1) of this and (u2,v2) of pl2)
// and direction in uvuv_startI(4-6) (optionally).
//isect_startPlanePt is a dedicated function for isect_startPlane.
//isect_startPlanePt halts the computation when intersection
//reached to a boundary of this, or reached to one of the points
//in uvuv_list.
//The function's return value is:
//  =0: Intersection was not obtained.
// !=0: Intersection was obtained as follows:
//  =1: End point is a point on a perimeter of this surfaces.
//  =3: End point is one of boundary points in uvuv_list.
//  =4: End point is the starting point.
//  =7: isect_startPlanePt halted the computation since intersection was lost
//     during the computation.
int MGFSurface::isect_startPlanePt(
	const MGPosition& uvuv_startIn, //Starting point of the intersection line.
	MGPosition_list& uvuv_list,	//isect_startPlanePt will halt when ip reached one of 
		//the point in uvuv_list. isect_startPlanePt does not change uvuv_list(actually
		//uvuv_list is const.) uvuv's space dimension is at least 4,
		//and the first 2 is (u,v) of this and the next 2 is (u,v) of pl2. 
	const MGPlane& pl2,	//2nd surface(MGPlane).
	double acuRatio,//Accurate ratio, should be decreased by multiplyng .2
		//(or a number less than 1.).
	MGBPointSeq& point,		//Surface-surface intersection parameter values
		//will be returned as:point(.,0) and point(.,1) for(u,v) of this surface
		//point(.,2) and point(.,3) for(u,v) of pl2.
		//point will have the dimension of(.,7).
	MGPosition_list::iterator& uvuv_id
			//When the end point of ip was one of the points of uvuv_list,
			//uvuv_list's iterator of the point will be returned, that is, 
			//when the function's return value was 3.
			//When was not a point of uvuv_list, end() of uvuv_list will be
			//returned.
)const{
	//Prepare work area for isect_start_boundary and isect_start_incr.
	int i, lmax=2*isect_area_length();

	MGBox pr1=param_range();
	double tol=MGTolerance::rc_zero();
	double dtmin[2]={pr1(1).length().value()*tol,pr1(0).length().value()*tol};
	//dtmin are minimum dt of below.

	tol*=(1000.*acuRatio); if(tol > .001) tol=.001;
	double tTol[4]={
		tol*knot_vector_u().param_span(),
		tol*knot_vector_v().param_span(),
		tol*pl2.u_deriv().len(), tol*pl2.v_deriv().len()
	}, tTolNow[4];	//tTol is the box tolerance to 
		//terminate the marching. if a point in uvuv_list is within this
		//tolerance of current intersection point, the marching will be
		//terminated.

	MGPosition uvuv_start(4,uvuv_startIn), uvuv_now(4);
	MGPosition uv1i(2,uvuv_start,0,0);
	MGPosition uv2i(2,uvuv_start,0,2);
	MGPosition_list::iterator uvuvi=uvuv_list.begin(), uvuve=uvuv_list.end();
	uvuv_id=uvuve;//Initialize the return value of uvuv_id.

	//Get maximum dt to not pass too far.
	double dt_max=-1.;
	for(;uvuvi!=uvuve; uvuvi++){
		double uvlen=(uv1i-MGPosition(2,*uvuvi)).len();
		if(dt_max<0. || dt_max>uvlen) dt_max=uvlen;
	}
	dt_max/=4.;

//Define initial parameter direction.
	int kdt,kdt2;;
		//When kdt=1: dt=u-value(isect of u=const parameter line)
		//     kdt=0: dt=v-value(isect of v=const parameter line)
	double du,dv, du3,dv3; 

	//////////Check too short intersection line.///////////
	if(uvuv_list.size()){
		MGPosition& uvuvend2=uvuv_list.front();
		double u0=uvuv_start[0], v0=uvuv_start[1];
		isect_dt(u0,v0,du,dv,acuRatio);
		double du2=uvuvend2[0]-u0, dv2=uvuvend2[1]-v0;
		if(fabs(du2)<du*.8 && fabs(dv2)<dv*.8){
		//Now this can be too short intersection line candidate.
			double u=u0+du2*.5, v=v0+dv2*.5;
			if(in_range(u,v)){
				point.resize(2,7);
				point.store_at(0,uvuv_start,0,0,4);
				point.store_at(0,pl2.eval(uvuv_start[2],uvuv_start[3]),4,0,3);
				point.store_at(1,uvuvend2,0,0,4);
				point.store_at(1,pl2.eval(uvuvend2[2],uvuvend2[3]),4,0,3);
				uvuv_id=uvuv_list.begin();
				return 3;
			}
		}
	}
	///////////////////************

//*****************Define initial direction.**********************
	MGPosition uvuv_start2(uvuv_startIn);
	kdt=isect_direction(pl2,0,uvuv_start2,du,dv,acuRatio);
	if(kdt==-1)
		return 0;

	double dt;
	if(kdt) dt=du; else dt=dv;
	if(dt_max>0. && dt_max<fabs(dt)){
		double save=dt; dt=dt_max; if(save<0.) dt=-dt;
	}

	point.resize(lmax,7);//(.,0-1) for this surface parameter (u,v),
							//(., 2-3) for pl2 surface parameter (u,v), and
							//(., 4-6) for (x,y,z) positional data. 

// Store starting point data.
//Recompute starting point by this method.
	uvuv_now=uvuv_start;
	isect_start_incr(pl2,uvuv_start,kdt,0.,0.,0,uvuv_now);

	point.store_at(0,uvuv_now,0,0,4);
	point.store_at(0,pl2.eval(uvuv_now[2],uvuv_now[3]),4,0,3);
	int n=1, nm1=0;
//Loop until new i.p. is on a boundary of surfaces or, i.p. reached to
// one of the points of uvuv_list.
	int obtained;
	bool kdt_changed=false;
	MGPosition uvuv_pre(uvuv_now);
	while(1){
		if(!(obtained=			//===First try.
			isect_start_incr(pl2,uvuv_pre,kdt,du,dv,0,uvuv_now))){
			du3=du*.3333; dv3=dv*.3333;
			obtained=		   //===Second try of same dt.
				isect_start_incr(pl2,uvuv_pre,kdt,du3,dv3,0,uvuv_now);
		}
		if(obtained){
			if(kdt_changed){
			double du4=uvuv_now[0]-uvuv_pre[0], dv4=uvuv_now[1]-uvuv_pre[1];
			if(kdt2){
				if(dv*dv4<=0.){ obtained=7;break;}
			}else{
				if(du*du4<=0.){  obtained=7;break;}
			}
			kdt_changed=false;
			}
		}else{
			kdt=!kdt; kdt_changed=true;
			if(kdt){if(fabs(du3)<=dtmin[kdt]) du3*=2.;}
			else {if(fabs(dv3)<=dtmin[kdt]) dv3*=2.;}
			//===Third try by changing kdt(u and v parameter line.)
			obtained=isect_start_incr(pl2,uvuv_pre,kdt,du3,dv3,0,uvuv_now);
			if(obtained){
				double du4=uvuv_now[0]-uvuv_pre[0], dv4=uvuv_now[1]-uvuv_pre[1];
				if(kdt){
					if(dv*dv4<=0.){ obtained=7;break;}
				}else{
					if(du*du4<=0.){ obtained=7;break;}
				}
				kdt_changed=false;
			}else{ obtained=7;break;};
		}

		for(int ii=0; ii<4; ii++){
			double dttemp=uvuv_pre[ii]-uvuv_now[ii];
			if(dttemp<0.) dttemp=-dttemp; tTolNow[ii]=dttemp*2.;
			if(tTolNow[ii]<tTol[ii]) tTolNow[ii]=tTol[ii];
		}
		MGBox uv_passed(uvuv_now, tTolNow);
			// uv_passed is the box of one span parameter length that passed.
			// When parameter of one of uvuv_list is included in this box, 
			// computation terminates.

		//Test if a point in uvuv_list or uvuv_start is in the last 1 span
		//parameter range(uv_passed). If so, end point of intersection point
		//is found and terminate the marching.
		if(uvuv_list.in(uv_passed,uvuv_id,2)){
			uvuv_now=MGPosition(4,*uvuv_id);
			obtained=3;
		}else if(n>3 && uv_passed.includes(uvuv_start)){
			uvuv_now=uvuv_start;
			obtained=4;
		}

	//Store obtained data in line.
		int pointSize=point.capacity();
		if(n>=pointSize){
			int len=pointSize+lmax;
			point.reshape(len);
		}
		point.store_at(n,uvuv_now,0,0,4);
		point.store_at(n,pl2.eval(uvuv_now[2],uvuv_now[3]),4,0,3);
		nm1=n++; point.set_length(n);
		if(obtained!=1) break;

	//Define new dt,kdt.
		du=uvuv_now[0]-uvuv_pre[0]; dv=uvuv_now[1]-uvuv_pre[1];
		kdt2=kdt;
		uvuv_pre=uvuv_now;
		MGPosition uv(2,uvuv_now,0,0);
		isect_inner_dt(nm1,uv,du,dv,kdt2,acuRatio);
		if(kdt2==-1){ obtained=7; break;}//kdt will be used even if kdt2=-1.
		if(kdt!=kdt2){ kdt=kdt2; kdt_changed=true;}	
	}

	if(n<=1) return 0;
	if(obtained==7){
		uvuv_now=uvuv_pre;
		if(kdt_changed) kdt=!kdt;
	}
	// obtained=1: ip found as a point on a perimeter of one of the surfaces.
	//             This case occurs at on_the_perimeter() since isect_start_incr's 
	//              return value is always 1.
	//         =3: ip passed one of boundary points in uvuv_list.
	//         =4: ip returned to the starting point.
	//         =7: ip was lost during the computation(may be tangent).

	MGPosition uvuv_end=uvuv_now;//Save the end point data.
//Normalize last points. The normalization is done by deleting last
//ndel points and adding (nadd+ndel) points in replace.
	int ndel,nadd=0;
	if(n<=isect_div_id_max()){ndel=0;}
	else if(n<isect_div_id_max()+3){ndel=n-isect_div_id_max();}
	else{ndel=3;}
	//nadd is number of points to add for the normalization.
	int naddby2=nadd+ndel;
	int nold=n; n+=nadd;
	if(n>point.capacity()) point.reshape(n);

	//Compute normalized span length of the last point (dt).
	int id_sect_div=nold-2-nadd; if(id_sect_div<0) id_sect_div=0;
	double divt,span=0., maxdivt=isect_dt_coef(id_sect_div);
	int j=naddby2;
	for(i=0; i<=naddby2; i++, j--){
		divt=isect_dt_coef(j); if(divt>maxdivt) divt=maxdivt;
		span=span+divt;
	}

	double dtSum=0.;
	int id=(kdt+1)%2;
	for(i=nold-ndel-1; i<=nold-1; i++){
		double dt2=point(i,id)-point(i-1,id);
		dtSum=dtSum+fabs(dt2);
	}
	double dtsave=dtSum/span;		//dtsave is always plus.

	//Re-compute last (nadd+ndel) points in normalized spans.
	int inext=n-naddby2-2; j=naddby2;
	uvuv_pre(0)=point(inext,0); uvuv_pre(1)=point(inext,1);
	uvuv_pre(2)=point(inext,2); uvuv_pre(3)=point(inext,3);
	nm1=n-1;
	for(i=inext+1; i<nm1; i++, j--){
		divt=isect_dt_coef(j); if(divt>maxdivt) divt=maxdivt;
		double dt=dtsave*divt;//dtsave is always plus.
		if(kdt){if(du<0.) du=-dt;else du=dt;}
		else{if(dv<0.) dv=-dt;else dv=dt;}
		isect_start_incr(pl2,uvuv_pre,kdt,du,dv,0,uvuv_now);
		point.store_at(i,uvuv_now,0,0,4);
		point.store_at(i,pl2.eval(uvuv_now[2],uvuv_now[3]),4,0,3);
		uvuv_pre=uvuv_now;
	}
	//Recompute End point by this method(uvuv_now includes the end point found).
	isect_start_incr(pl2,uvuv_end,kdt,0.,0.,0,uvuv_end);
	point.store_at(nm1,uvuv_end,0,0,4);
	point.store_at(nm1,pl2.eval(uvuv_end[2],uvuv_end[3]),4,0,3);
	point.set_length(n);
	return obtained;
}

//Compute an intersection line of this surface(a surface that is converted to
//1D SBRep surface(sf1d)) and a plane pl.
//sf1D is so converted that pl be x=0. plane. This surface is the original
// surface and pl is the original plane.
//The function's return value is:
// =0: Intersection was not obtained.
// !=0: Intersection was obtained as follows:
//    =1: End point is a point on a perimeter of one of the surfaces.
//    =3: End point is one of boundary points in uvuv_list.
//    =4: End point is the starting point.
//    =7: isect_startPlane halted the computation since intersection was lost
//     during the computation.
int MGFSurface::isect_startPlane(
	const MGPosition& uvuvS,//Starting parameter value of the intersection.
	MGPosition_list& uvuv_list,
		//isect_startPlane will halt when ip reached one of the point in
		//uvuv_list. isect_startPlane does not change uvuv_list(actually
		//uvuv_list is const.) uvuv's space dimension is at least 4,
		//and the first 2 is (u,v) of this and the next 2 is (u,v) of pl. 
		//When uvuv's space dimension is more than 4, it indicates that the uvuv
		//is used to input approcimate tangent of the intersection.
	const MGPlane& pl,	//Plane expression(2nd surface for the intersection).
	MGSSisect& ssi,		//Surface-surface intersection line will be output.
	MGPosition_list::iterator& uvuv_id
			//When the end point of ip was one of the points of uvuv_list,
			//uvuv_list's iterator of the point will be returned, that is, 
			//when the function's return value was 3 or 5.
			//When was not a point of uvuv_list, end() of uvuv_list will be
			//returned.
)const{
	ssi.set_null();
	double errwcSave=MGTolerance::wc_zero();
	double errwc=errwcSave*.5;

	//Compute intersection lines using isect_startPlane.
	int obtained;
	const double dlen[2]={eval(uvuvS,1,0).len(),eval(uvuvS,0,1).len()};
	double lineError=MGTolerance::line_zero();
	double err2=lineError*.5; err2*=err2;

	double acuRatio=1.;
	double deviOld;
	MGLBRep lineb;
	for(int numRepeat=0; numRepeat<6; numRepeat++){

	MGBPointSeq point;
	//**********Compute intersection points using isect_startPlanePt.***********
	mgTolSetWCZero wczeroSet(errwc);//Set&save the error.
	obtained=isect_startPlanePt(uvuvS,uvuv_list,pl,acuRatio,point,uvuv_id);
	if(obtained==0 || obtained==1){
		obtained=0;
		break;
	}

	int n=point.length();
	MGNDDArray tau(n);
	double tauati=0.;
	tau(0)=tauati;
	MGPosition ppre(point(0,4),point(0,5),point(0,6));
	int i;
	for(i=1; i<n; i++){
		MGPosition p(point(i,4),point(i,5),point(i,6));
		tau(i)=tauati+=(ppre-p).len();
		ppre=p;
	}
	if(tau[n-1]<=errwcSave){
		return obtained;
	}

	int IMLT=1, pointSize=point.capacity(), pointDim=7;
	double ratio=4.;
	bkdnp_(&n,tau.data(),point.data(),pointSize,pointDim,IMLT,ratio);
	tau.set_length(n); point.set_length(n);
	int k=4;//Set order to 4.
	if(k>n) k=n;
	lineb.knot_vector()=MGKnotVector(tau, k);
	lineb.buildByInterpolationWithKTV(tau,point);

//Check if the obtained intersection lies on the two surface.	
	double devi1=-1.;//get the maximum deviations in devi1,2.
	for(i=1; i<n; i++){
		double t=tau(i-1)+tau(i); t*=.5;
		MGVector iP=lineb.eval(t);
		MGVector Q(3,iP,0,4), P1(eval(iP[0],iP[1]));
		MGVector def1(P1-Q);
		double len1=def1%def1;
		if(devi1<len1) devi1=len1;
	}
	if(devi1<=err2){
		if(n<4) break;
		MGVector tan1,tan2;
		MGVector* tan[2]={&tan1, &tan2};
		int ngtan=isect_start_tan(*this,pl,lineb,tan);
		if(ngtan)
			isect_start_adjustSE(ngtan,tau,point,lineb,tan);
			//Adjust start and end tangent line by value tan.
		break;
	}

	deviOld=devi1;
	acuRatio*=.3;	//Try again by raising accuracy.
	errwc*=.5;

	}

	if(obtained){
		lineb.remove_knot(4,3);
		MGLBRep* lineuv1=new MGLBRep(2,lineb,0,0);
		MGLBRep* lineuv2=new MGLBRep(2,lineb,0,2);
		MGLBRep* linexyz=new MGLBRep(3,lineb,0,4);
		ssi=MGSSisect(linexyz,lineuv1,lineuv2);
	}
	return obtained;
}

//Compute the intersection lines of this(is not a plane) surface and a plane pl.
//sd1D that is converted to surf1D() about pl is necessary to input.
MGSSisects MGFSurface::isect_with_plane(
	MGPosition_list& uvuv_list,
	//Let a member of uvuv_list be uvuv. Then, uvuv's space dimension is
	//at least 4, and the first 2 is (u,v) of this and the next 2 is (u,v) of pl. 
	//When uvuv's space dimension is more than 4, it indicates that the uvuv
	//is used to input approximate tangent of the intersection.
	const MGPlane& pl,	//Plane expression(2nd surface for the intersection).
	const MGFSurface& fsrf2	//Original FSurface before casting to plane.
)const{
	MGSSisects lst(this,&fsrf2);
	MGSSisect ssi;
	MGPosition_list::iterator uvuv_id;
	int obtained;
	int n;
	int loop_number=0, max_loop=uvuv_list.entries();
	while(n=uvuv_list.entries()){
		if(n>=2){
			const MGPosition& Ps=uvuv_list.front();
			MGVector N1=unit_normal(Ps[0],Ps[1]), N2=pl.unit_normal(Ps[2],Ps[3]);
			double saS=N1.sangle(N2);
			if(saS<.02){
				const MGPosition& Pe=uvuv_list.back();
				N1=unit_normal(Pe[0],Pe[1]); N2=pl.unit_normal(Pe[2],Pe[3]);
				double saE=N1.sangle(N2);
				if(saS<saE) uvuv_list.reverse_order();
				//If two surfaces are more parallel at the starting point than
				//at the ending point, we change the intersection line direction.
			}
		}
		MGPosition uvuvS=uvuv_list.removeFirst();
again:	if(obtained=isect_startPlane(uvuvS,uvuv_list,pl,ssi,uvuv_id)){
			MGPosition uvuvE;
			if(obtained==3) uvuvE=uvuv_list.removeAt(uvuv_id);
			if(!ssi.is_null())
				uvuv_list.removeOn(*this,pl,ssi); //Remove uvuv that is on the ssi.

			if(uvuv_list.size() && obtained==3){
				//Check if the obtained section's end point is not a terminal point
				//but only a mid point of a section(and should be neglected).

				if(!ssi.is_null()){
				const MGCurve& iline=ssi.line();
				//Check to the ending direction.
				if(uvuvE_is_a_midpoint(fsrf2,0,ssi,uvuvE)){
					uvuvS=MGPosition(7,uvuvS);
					uvuvS.store_at(4,iline.eval(iline.param_s(),1));
					goto again;
						//goto again means uvuvE is a mid point and will be neglected.
				}

				//Check to the starting direction by reversing the ssi.
				ssi.negate();
				if(uvuvE_is_a_midpoint(fsrf2,0,ssi,uvuvS)){
					uvuvS=MGPosition(7,uvuvE);
					uvuvS.store_at(4,iline.eval(iline.param_s(),1));
					goto again;
				}
				}
			}
			if((!lst.size() || obtained!=7)&& !ssi.is_null())
				lst.append(std::move(ssi));
		}else{
			uvuv_list.append(uvuvS);
		}
		loop_number++;
		if(loop_number>=max_loop) break;
	}
	return lst;
}

//Compute common curve part of srf's perimeter and the crv.
//Function's returned value is the number of common part curve part,
//which is 2 at most.
int MGSurface::getPerimeterCommon(
	const MGCurve& crv,//crv must be a cotinuous one curve, must not be MGCompositeCurve.
	std::vector<double> pspan[2],//parameter range of crv and perimeters will be output.
	int peri_num[2]	//perimeter number of pspan[i] will be output in peri_num[i].
		//(pspan[i], peri_num[i]) is one pair.
)const{
	int nperi=perimeter_num(), nComPeri=0;
    for(int i=0; i<nperi; i++){//nperi=0, or 4.
		std::unique_ptr<MGCurve> perimeterj(perimeter_curve(i));
		int nspan=crv.common(*perimeterj,pspan[nComPeri]);
		if(nspan){
			peri_num[nComPeri]=i;
			nComPeri++;
			if(nComPeri>=2)
				break;
		}
	}
	return nComPeri;
}

//split this fsurface at the parameter param.
void MGSurface::split(
	double param,//parameter value of this fsurface. if is_u is true, param is u-value,
				//else v-value.
	bool is_u,	//indicates if param is u or v of the surface parameter (u,v).
	std::vector<UniqueFSurface>& surfaces//splitted surfaces will be output.
)const{
	surfaces.clear();
	MGBox uvrange=box_param();
	int id=is_u ? 0:1;
	MGInterval& prange=uvrange[id];
	if(prange.includes(param)){
		double error=prange.relative_error();
		MGBox uvrange2(uvrange);//save the original box.
		prange.set_high(param);
		if(prange.length()>error)
			surfaces.emplace_back(part(uvrange));

		MGInterval& prange2=uvrange2[id];
		prange2.set_low(param);
		if(prange2.length()>error)
			surfaces.emplace_back(part(uvrange2));
	}else
		surfaces.emplace_back(clone());
}

//Given MGSBRep or MGRSBRep as srf, compute normalize MGKnotVector along u and v.
void MGSurface::get_new_surface_knots(
	int parameter_normalization,
		//Indicates how the parameter normalization be done:
		//=0: no surface parameter normalization.
		//=1: normalize to u_range=(0., 1.), and v_range=(0.,1.);
		//=2: normalize to make the average length of the 1st derivative along u and v 
		//    of the base surface is as equal to 1. as possible.
	MGKnotVector& uknots,//normalized MGKnotVector will be returned.
	MGKnotVector& vknots,
	double* Oldparameter//Old parameterization will be returned. [.]=(us, ue, vs,ve).
)const{
	double prange2[4];
	if(!Oldparameter)
		Oldparameter=prange2;
	double& usOld=Oldparameter[0];
	double& ueOld=Oldparameter[1];
	double& vsOld=Oldparameter[2];
	double& veOld=Oldparameter[3];

	uknots=knot_vector_u();
	vknots=knot_vector_v();
	usOld=uknots.param_s(), ueOld=uknots.param_e();
	vsOld=vknots.param_s(), veOld=vknots.param_e();
	if(parameter_normalization==1){
		uknots.change_range(0.,1.);
		vknots.change_range(0.,1.);
	}else if(parameter_normalization==2){
		MGNDDArray utau; utau.buildByKnotVector(uknots);
		double v[3]={vsOld,(vsOld+veOld)*.5, veOld};
		double uspan=average_chord_length(0,v,utau);
		uknots.change_range(0.,uspan);

		MGNDDArray vtau; vtau.buildByKnotVector(vknots);
		double u[3]={usOld,(usOld+ueOld)*.5, ueOld};
		double vspan=average_chord_length(1,u,vtau);
		vknots.change_range(0.,vspan);
	}
}

#define INCNUM 15
///Approximate this MGSBRep or MGRSBRep by new knot configuration.
///Remove the redundant surface B-coefficients within the tolerance, which is
///performed by remove_knots.
std::unique_ptr<MGSBRep> MGSurface::approximate_as_SBRep(
	int parameter_normalization,
		//Indicates how the parameter normalization be done:
		//=0: no surface parameter normalization.
		//=1: normalize to u_range=(0., 1.), and v_range=(0.,1.);
		//=2: normalize to make the average length of the 1st derivative along u and v 
		//    of the base surface is as equal to 1. as possible.
	double tol,	///<tolerance allowed for the approximation.
		///When tol<=0., MGTolerance::line_zero() will be employed.
	int* order///<order of the new MGSBRep, >=4 is recommended.
		///order[0]:u-order, [1]:v-order.
		///When order=0 is input, the original order is unchanged.
)const{
	const MGKnotVector& ut=knot_vector_u();
	const MGKnotVector& vt=knot_vector_v();
	MGNDDArray utau,vtau;
	utau.buildByKnotVector(ut);utau.change_number(utau.length()*INCNUM);
	vtau.buildByKnotVector(vt);vtau.change_number(vtau.length()*INCNUM);
	MGSPointSeq sp;
	eval_spoint(utau,vtau,sp);
	int uorder=4, vorder=4;
	if(order){
		if(order[0])
			uorder=order[0];
		if(order[1])
			vorder=order[1];
	}
	std::unique_ptr<MGSBRep> sn(new MGSBRep);
	sn->buildByInterpolationDataPoints(utau, vtau, sp, uorder, vorder);
	mgTolSetLineZero lineZeroSet(tol);
	sn->remove_knot();
	if(!parameter_normalization)
		return sn;

	MGKnotVector uknots, vknots;
	sn->get_new_surface_knots(parameter_normalization,uknots,vknots);
	std::unique_ptr<MGSBRep> surfNew(new MGSBRep);
	surfNew->buildSBRepFromMemberData(std::move(sn->surface_bcoef()),
					std::move(uknots), std::move(vknots));
	return surfNew;
}
