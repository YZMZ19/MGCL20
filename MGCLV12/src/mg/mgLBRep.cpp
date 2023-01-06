/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"

#include "cskernel/Bluprt.h"
#include "cskernel/Blumix.h"
#include "cskernel/bkdtpg.h"
#include "cskernel/Blgl2a.h"
#include "cskernel/Blgint.h"
#include "cskernel/Bkdnp.h"
#include "cskernel/bkdtkt.h"
#include "cskernel/blg4sc.h"
#include "cskernel/blg4sp2.h"
#include "cskernel/mgblgsq.h"

#include "mg/Tolerance.h"
#include "mg/Box.h"
#include "mg/Position.h"
#include "mg/KnotArray.h"
#include "mg/CParam_list.h"
#include "mg/BPointSeq.h"
#include "mg/Ellipse.h"
#include "mg/Straight.h"
#include "mg/PPRep.h"
#include "mg/LBRepEndC.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/TrimmedCurve.h"
#include "mg/SurfCurve.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

// MGLBRep.cpp
//
// Implement MGLBRep class.
//
// This file contains all of the constructors of MGLBRep.


	//**** 1. Interpolation Constructor ****

///Construct Line B-rep by intepolation from Point data only.
///If circular is true, start and end points must be the same,
///(that is points[0]==points[n-1]) and the MGLBRep constructed is 
///smoothly connected at the point.
///If circular is true, order will be always 4, and input order is neglected.
///The knot vector is generated from data points of the input points'
///chord length starting 0.
///Function's return value is 0:successfully built, !=0:too near data is input.
int MGLBRep::buildByInterpolation(
	const MGBPointSeq& points,//Point seq data
	int order,			// Order
	bool circular		///<Circular flag
){
	assert(order>=2 && points.length()>=2);
	invalidateBox();

	m_knot_vector.size_change(order, points.length()+2);
	m_line_bcoef.resize(points.length()+2, points.sdim());

	MGBPointSeq points_temp(points);//Copy, since blg4sc_ may change input data.
	int nv=points_temp.length();	//Number of input points.
	const int iv=points_temp.capacity();
	const int ncd=points_temp.sdim();
	const int irc=m_line_bcoef.capacity();
	int n;

	int error=2;//When order=2 is input, error stays equal to 2.
	if(circular){
		if(order!=4) m_knot_vector.change_order(4);
		double* work1=new double[irc*10];
		double* work2=work1+irc;
		blg4sc_(nv, &points_temp(0,0), iv, ncd, irc, work1,work2,
				&n, &m_knot_vector(0), &m_line_bcoef(0,0), &error);
		delete[] work1;
	}else if(order>2){
		if(nv<order){
			order=nv;
			m_knot_vector.size_change(order,nv);
		}
		double* work=new double[nv*order*2];
		n=nv;
		//   GENERATE DATA POINTS IN WORK(I,1) 0<=I<=N-1 . 
		double* val=&points_temp(0,0);
		bkdtpg_(val, n, ncd, iv, work);
		//   DISCARD TOO NEAR POINTS 
		int nnew = n;
	    bkdnp_(&nnew,work,val,iv,ncd,1,MGTolerance::max_knot_ratio());
		if(nnew<order){
			order=nnew;
			m_knot_vector.change_order(order);
		}
		//   GENERATE KNOT VECTOR IN T FROM DATA POINTS. 
		double* t=&m_knot_vector(0);
		bkdtkt_(work, nnew, order, t);
		//   GENERATE B-REP. 
		error=blgint_(work,val,t,order,nnew,ncd,iv,irc,work+n,&m_line_bcoef(0,0));
		n = nnew;
		delete[] work;
	}

	if(error==1) {//Return of BLG4SC,blgint_ error=1 means normal.
		m_knot_vector.set_bdim(n);
		m_line_bcoef.set_length(n);
		error=0;
	}else{
		//Generate polyline B-Rep. of order 2.
		//When input order==2 or error detected, the polyline that connects input points is set.
		m_knot_vector=MGKnotVector(2,points.length(),0.,1.);
		m_line_bcoef=points;
	}
	return error;
}

//Construct Line B-rep of a specified order, given data point abscissa and the ordinates.
//When error is returned, B-Rep of polyline of points is constructed that has uniform
//knot vector.
///Function's return value is 0:successfully built,
//!=0:tau is illegal(too near data is input).
//When error is returned, B-Rep of polyline of points is constructed that has uniform
//knot vector.
int MGLBRep::buildByInterpolationDataPoints(
	const MGNDDArray& tau,		//Data point abscissa
	const MGBPointSeq& points,	//Point seq data(ordinates).
	int order,				//order
	double ratio			//Maximum of data point ratio of pre and after spans.
			// Let d(i)=tau[i]-tau[i-1], then if d(i)/d(i-1)>ratio or
			// d(i-1)/d(i)>ratio, either tau[i] or tau[i-1] will be removed.
			// This is done to prevent control polygon computation error.
			// When ratio<0. no data point removal will be done.
){
	invalidateBox();
	m_line_bcoef.resize(points.length(), points.sdim());
	int IMLT=1;
	MGNDDArray tau2(tau);
	MGBPointSeq points2(points);
	int n=points2.length(), pointSize=points2.capacity(), pointDim=points2.sdim();
	if(ratio>0.){
		bkdnp_((int*)&n,tau2.data(),points2.data(),pointSize,pointDim,IMLT,ratio);
		tau2.set_length(n); points2.set_length(n);
	}

	int k=order;//Order
	if(k>n) k=n;
	m_knot_vector=MGKnotVector(tau2,k);
	double* work=new double[n*(2*k-1)];
	const int irc=m_line_bcoef.capacity();
	int error=blgint_(tau2.data(), points2.data(), m_knot_vector.data(), k, n, pointDim,
			pointSize, irc, work, &m_line_bcoef(0,0));
	delete[] work;

	if(error==1){
		error=0;		//Return of BLGINT error=1 means normal.
		m_line_bcoef.set_length(n);
	}else{
		//Generate polyline B-Rep. of order 2.
		//When error detected, the polyline that connects input points is set.
		m_knot_vector=MGKnotVector(2,n,0.,1.);
		m_line_bcoef=points2;
	}
	return error;
}

//Construct Line B-rep by interpolation, given the knot vector in m_knot_vector,
//and (tau, points) as arguments.
int MGLBRep::buildByInterpolationWithKTV(
	const MGNDDArray& tau,	//Data point abscissa
	const MGBPointSeq& points	//Point seq data
){
	assert(m_knot_vector.bdim()==points.length());

	int k=m_knot_vector.order(), n=points.length(), ncd=points.sdim();
	m_line_bcoef.resize(n, ncd);
	int iv=points.capacity(), irc=m_line_bcoef.capacity();

	double* work=new double[n*(2*k-1)];
	int error=blgint_(tau.data(), points.data(), m_knot_vector.data(), k, n, ncd,
			iv, irc, work, &m_line_bcoef(0,0));
	delete[] work;

	if(error==1){
		m_line_bcoef.set_length(n);
		error=0;
	}else{
		//Generate polyline B-Rep. of order 2.
		//When error detected, the polyline that connects input points is set.
		m_knot_vector=MGKnotVector(2,n,0.,1.);
		m_line_bcoef=points;
	}
	return error;
}

void build_mgblgsp_parameter(
	int& order,				//Order of the target B-Spline.
	const MGLBRepEndC& begin,	//Begin end condition
	const MGLBRepEndC& end,		//End end conditoion
	const MGNDDArray& tau,		//Data point abscissa
	const MGBPointSeq& value,	//Data point ordinate
	MGENDCOND& beginc,
	MGENDCOND& endc,
	MGNDDArray& tau_new,	//Data point abscissa
	MGBPointSeq& value_new	//Data point ordinate.
							//value_new.size() wil be value.length().
){
	int i,j;
	int ie=1, is=0;
	const int n=value.length();
	const int nm1=n-1;
	int nnew=n;
	beginc=begin.cond(); endc=end.cond();
	if(beginc== MGENDCOND::MGENDC_1D || beginc== MGENDCOND::MGENDC_2D) {nnew+=1;is=1;}
	else if(beginc== MGENDCOND::MGENDC_12D) {nnew+=2; is=2;}
	if(endc== MGENDCOND::MGENDC_1D || endc== MGENDCOND::MGENDC_2D){nnew+=1; ie=2;}
	else if(endc== MGENDCOND::MGENDC_12D) {nnew+=2; ie=3;}
	if(int(order)>nnew) order=nnew;
	const int ncd=value.sdim();

	value_new=value;
	value_new.reshape(nnew,is);
	int nnewm1=nnew-1;
	for(j=0; j<ncd; j++) {
		value_new(0,j)=value(0,j);
		value_new(nnewm1,j)=value(nm1,j);
	}
	tau_new=tau;
	tau_new.reshape(nnew,is);
	for(int k=0; k<is; k++) tau_new(k)=tau(0);
	for(i=n+is;i<nnew; i++) tau_new(i)=tau(nm1);

	if(beginc== MGENDCOND::MGENDC_1D || beginc== MGENDCOND::MGENDC_12D)
		for(j=0; j<ncd; j++) value_new(1,j)=(begin.first()).ref(j);
	if(beginc== MGENDCOND::MGENDC_2D || beginc== MGENDCOND::MGENDC_12D)
		for(j=0; j<ncd; j++) value_new(is,j)=(begin.second()).ref(j);
	if(endc== MGENDCOND::MGENDC_1D || endc== MGENDCOND::MGENDC_12D)
		for(j=0; j<ncd; j++) value_new(nnew-2,j)=(end.first()).ref(j);
	if(endc== MGENDCOND::MGENDC_2D || endc== MGENDCOND::MGENDC_12D)
		for(j=0; j<ncd; j++) value_new(nnew-ie,j)=(end.second()).ref(j);	
}

// Construct Line B-rep of input order by interpolation from Point data
//and end condition.
// Inner point may include derivative inf.
int MGLBRep::buildByInterpolationEC(
	const MGLBRepEndC& begin,	//Begin end condition
	const MGLBRepEndC& end,		//End end conditoion
	const MGNDDArray& tau,		//Data point abscissa
	const MGBPointSeq& value,	//Data point ordinate
	int order				//Order of the target B-Spline.
){
	invalidateBox();
	MGENDCOND beginc, endc;
	MGBPointSeq valuen;	//Data point ordinate
	MGNDDArray taun;	//Data point abscissa
	build_mgblgsp_parameter(order,begin,end,tau,value,beginc,endc,taun,valuen);
	int nnew=valuen.length();
	int ncd=valuen.sdim();

	m_line_bcoef.resize(nnew,ncd);
	m_knot_vector.size_change(order,nnew);
	int irc=nnew;
	int iv=nnew;
	double* work=new double[nnew*(order*2+1)];
	int error;
	int bc = static_cast<int>(beginc), ec = static_cast<int>(endc);
	mgblgsq((int)order,bc,ec,taun.data(),valuen.data(),iv,nnew,ncd,
			irc,work,&m_knot_vector(0),&m_line_bcoef(0,0),&error);
	delete[] work;
	if(error==1){
		error=0;		//Return of mgblgsq error=1 means normal.
		m_knot_vector.set_bdim(nnew);
		m_line_bcoef.set_length(nnew);
	}else{
		//Generate polyline B-Rep. of order 2.
		//When error detected, the polyline that connects input points is set.
		m_knot_vector=MGKnotVector(2,value.length(),0.,1.);
		m_line_bcoef=value;
	}
	return error;
}

///Construct Line B-rep of any order by interpolation from Point data
///with end condition as arguments:
/// (tau(i), value(i,.)) for 0<=i<=n(the length of value),
///and the knot vector as the member.
///tau(i) and knot vector t must satisfy Shoenberg's variation diminishing
///constraint.
///Before buildByInterpolationECWithKTV(), use setKnotVector to construct this knot vector.
///For the start and end point, tau does not have multiplicity. However,
///if tau has multiplicity at inner point, this means 1st derivative data is
///provided for the associated value, i.e.:
///If tau has multiplicity 2 as tau(i)=tau(i+1), value(i,.) is 1st derivative
///at tau(i) and value(i+1,.) is positional data at tau(i)(=tau(i+1)).
///If tau has multiplicity 3 as tau(i)=tau(i+1)=tau(i+2),
///value(i,.) is 1st derivative at tau(i)- ,
///value(i+1,.) is positional data at tau(i)(=tau(i+1)), 
///value(i+2,.) is 1st derivative at tau(i)+.
///Maximum multiplicity allowed is 3.
/// =0: successfully built.
/// !=0: error occured and the MGLBRep is built as uniform BSpline of order 2.
int MGLBRep::buildByInterpolationECWithKTV(
	const MGLBRepEndC& begin,	//Begin end condition
	const MGLBRepEndC& end,		//End end conditoion
	const MGNDDArray& tau,		//Data point abscissa
	const MGBPointSeq& value	//Data point ordinate
){
	MGENDCOND beginc, endc;
	MGBPointSeq valuen;	//Data point ordinate
	MGNDDArray taun;	//Data point abscissa
	int k=order(); assert(k>=2);
	build_mgblgsp_parameter(k,begin,end,tau,value,beginc,endc,taun,valuen);
	int nnew=valuen.length();
	int ncd=valuen.sdim();

	m_line_bcoef.resize(nnew,ncd);
	double* work=new double[nnew*(k*2+1)];
	int error = 2;
	int bc = static_cast<int>(beginc), ec = static_cast<int>(endc);
	for(int i=0; i<ncd; ++i){
		blg4sp2_(k,&error,bc,ec,taun.data(),valuen.data()+i*nnew,nnew,nnew,1,
			m_knot_vector.data(),1,work,work+nnew,work+2*nnew,&m_line_bcoef(0,i));
			if (error!=1) break;
	}
	delete[] work;

	if(error==1)
		error=0;//Return of mgblgsq error=1 means normal.
	else{
		//Generate polyline B-Rep. of order 2.
		//When error detected, the polyline that connects input points is set.
		m_knot_vector=MGKnotVector(2,value.length(),0.,1.);
		m_line_bcoef=value;
	}
	return error;
}

//**** 2. Approximation Constructor ****

//Construct Curve B-Rep ep.
//This is an approximation, and the tolerance is MGTolerance::line_zero().
MGLBRep::MGLBRep(
	const MGCurve& crv,	//Original Curve.
	int order,//Order. When order=0 is input, and crv was a MGLBRep,
		//the original order will be used. Otherwise(order=0 and crv was not an MGLBRep)
		//order will be set to 4.
	int parameter_normalization,
		//Indicates how the parameter normalization be done:
		//=0: no parameter normalization.
		//=1: normalize to range=(0., 1.);
		//=2: normalize to make the average length of the 1st derivative 
		//    is as equal to 1. as possible.
	bool neglectMulti///<Indicates if multiple knots be kept.
		///< true: multiplicity is removed.
		///< false: multiplicity is kept.
):MGCurve(crv){
	crv.approximate_as_LBRep(*this,order,parameter_normalization,neglectMulti);
}

// Construct 3D B-Rep by mixing two 2D B-Rep.
//The two 2D B-Rep's directions and start and end points must be the same.
void MGLBRep::buildByMixing2DLBRep(
	int coordinate1,	//Missing oordinate kind of the brep1 
							// 0:x, 1:y, 2:z.
	const MGLBRep& brep1,	//Original 2D B-Rep1. Coordinates are
							//(y,z), (z,x), (x,y) according to coordinate1. 
	int coordinate2,	//Missing coordinate kind of the brep2.
							// 0:x, 1:y, 2:z, and 3:girth rep.
	const MGLBRep& brep2)	//Original 2D B-Rep2. Coordinates are
		//(y,z), (z,x), (x,y) and (t, g2) according to coordinate2.
		//t is parameter of brep1 and g2 is x, y, or z according to coordinate1.
//Second 2D B-Rep can be girth representaion. That is,
//Let brep1 is f(t)=(f1(t),f2(t)) and brep2 is g(s)=(g1(s),g2(s)), where
//f1,f2 are two coordinates, g1 is parameter t of f(t), and 
//g2(s) is the missing coordinate of f(t). Given s value,
// ( f1(g1(s)), f2(g1(s)), g2(s) ) is a 3D space point.
{
	assert(brep1.sdim()==2 && brep2.sdim()==2);
	assert(coordinate1<=2 && coordinate2<=3 && coordinate1!=coordinate2);
	invalidateBox();

	int i;
	int kcod1=coordinate1+1;
	const int k1=brep1.order();
	const int n1=brep1.m_line_bcoef.length();
	const int irc1=brep1.m_line_bcoef.capacity();

	int kcod2=coordinate2+1;
	const int k2=brep2.order();
	const int n2=brep2.m_line_bcoef.length();
	const int irc2=brep2.m_line_bcoef.capacity();

	int k=k1; if(k<k2) k=k2; k=4*k*k+3*k;
	int irc=n1+n2; if(irc<11) irc=11;
	if(irc*4<k) irc=int(k/4)+1;

	double error=MGTolerance::line_zero();
	int* kseq=new int[irc];
	double* wk2=new double[irc*4];
	int n; MGNDDArray tau(irc); MGBPointSeq bp(irc,3);
	int iflag; MGVector vec_s(3),vec_e(3);

	//Compute initial data by BLUMIX.
	blumix_(error,kcod1,k1,n1,brep1.knot_data(),brep1.coef_data(),irc1,
		kcod2,k2,n2,brep2.knot_data(),brep2.coef_data(),irc2,
		irc,kseq,wk2,&vec_s(0),&vec_e(0),&n,&tau(0),&bp(0,0),&iflag);
	tau.set_length(n); bp.set_length(n);
	delete[] kseq;
	delete[] wk2;
	MGLBRepEndC endc_s(MGENDCOND::MGENDC_1D,vec_s);
	MGLBRepEndC endc_e(MGENDCOND::MGENDC_1D,vec_e);

	//Check tolerance.
	int cod11,cod12,cod21;
	cod11=kcod1;if(cod11>=3) cod11=0;
	cod12=cod11+1;if(cod12>=3) cod12=0;
	int com21=1;
	if(coordinate2<3){
		cod21=kcod2;if(cod21>=3) cod21=0;
		if(cod21==coordinate1)com21=0;
	}
	//cod11,12,21 are 3D coordinate id of brep1 and 2.
	//com21 is id of 2nd line coordinate that is missing in 1st line.

	MGPosition F,G(3); MGPosition b1,b2,p1,p2; 
	int added;
	MGLBRep line;
	double tmid,s1,s2,diff;
	for(int loop=0; loop<3; loop++){
		//Maximum loop counter is 3(3 times as much as output of BLUMIX)
		added=0;
		n=bp.length();
		line.buildByInterpolationEC(endc_s,endc_e,tau,bp);
		for(i=1; i<n-1; i++){
			while((tau(i)==tau(i+1)) && i<n-1) i+=1;
			tmid=(tau(i)+tau(i+1))/2.;
			F=line.eval(tmid);
			b1=MGVector(2,F,0,cod11);
			brep1.on(b1,s1); p1=brep1.eval(s1);
			diff=(p1-b1).len();
			if(diff>error){			//brep1 differs from line, add point.
				G=brep2.closest_mix(coordinate2,coordinate1,s1,p1,F);
				bp.insert_at(i+1,G); tau.add_data(tmid);
				i+=1; n+=1; added+=1;
			} else{
				if(coordinate2<3) b2= MGVector(2,F,0,cod21);
				else			  b2= MGVector(s1,F(coordinate1));
				brep2.on(b2,s2); p2=brep2.eval(s2);
				diff=(p2-b2).len();
				if(diff>error){		//brep2 differs from line, add point.
					G.set(cod11)=p1.ref(0); G.set(cod12)=p1.ref(1);
					G.set(coordinate1)=p2.ref(com21);
					bp.insert_at(i+1,G); tau.add_data(tmid);
					i+=1; n+=1; added+=1;
				}
			}
		}
		if(!added) break;
	}
	if(added){
		n=bp.length();
		line.buildByInterpolationEC(endc_s,endc_e,tau,bp);
	}
	buildLBRepFromMemberData(
		std::move(line.m_knot_vector), std::move(line.m_line_bcoef)
	);
	copy_appearance(brep1);
}

//Function for BLUMIX constructor. Given 3D point F,
//compute correct point of F that is closest to F.
MGPosition MGLBRep::closest_mix(
	int coordinate2,	//Missing coordinate kind of this(say 2nd line).
							//coordinate2 can be 3.
	int coordinate1,	//Missing coordinate kind of P(say, of 1st line).
							//coordinate1 cannot be 3.
	double tau,				//Parameter value of P, of 1st line,
							//used only when coordinate2=3.
	const MGPosition& P,	//Point of 1st line(correct coordinates).
	const MGPosition& F		//3D point, used to coose closest point to F
							//when more than one point are found.
)const{
	MGCParam_list list; int kcod;

	int com2=0, com21=1;
	if(coordinate2<3){
		kcod=coordinate2+1; if(kcod>=3) kcod=0;
		if(kcod==coordinate1){com2=1; com21=0;}
	}
	int com1=com2+1; if(com1>=2) com1=0;
	//com2 is id of 2nd line coordinate that is common to P(1st line).
	//com21 is id of 2nd line coordinate that is missing in P(1st line).
	//com1 is id of P coordinate that is common to 2nd line.

	if(coordinate2==3) list=isect_1D(tau);
	else{
		double data=P.ref(com1); list=isect_1D(data,com2);
		double tol=MGTolerance::wc_zero();
		if(!list.entries()) list=isect_1D(data+tol,com2); 
		if(!list.entries()) list=isect_1D(data-tol,com2);
	}

	kcod=coordinate1+1; if(kcod>=3) kcod=0;
	MGPosition Q(3,P,kcod);
	int n=list.entries();
	if(!n){
		Q(coordinate1)=F(coordinate1);
		return Q;
	}
	MGPosition A=eval(list.removeFirst()); Q(coordinate1)=A(com21);
	if(n==1) return Q;

	double dist=(Q-F).len(), dist2;
	MGPosition Q2(Q);
	for(int i=1; i<n; i++){
		A=eval(list.removeFirst()); Q2(coordinate1)=A(com21);
		dist2=(Q2-F).len();
		if(dist2<dist){ dist=dist2; Q=Q2;}
	}
	return Q;
}	

// Construct Line B-rep of any order number by least square approximation
//from Point data with approximation weights and knot vector of B-Rep.
void MGLBRep::buildByL2ApproximateWithKTV(
	const MGNDDArray& tau,	//Data point abscissa.
	const MGBPointSeq& points,	//Data Point ordinates.
	const double* weight	//Weights for each points 
){
	invalidateBox();
	int ntau=tau.length(), ig=points.capacity();

	int sdim=points.sdim();
	int k=m_knot_vector.order();
	const int n=m_knot_vector.bdim();
	assert(n==ntau);

	int iw=1, irc=1, kw=1, m=1;
	double* work=new double[2*n+k*n];
	double* q=work+2*n;

	m_line_bcoef.resize(n, sdim);
	for(int i=0; i<sdim; i++)
		blgl2a_(ntau,tau.data(),points.data(0,i),weight,ig,iw,irc,
			kw,m_knot_vector.data(),n,k,m, work, q, &m_line_bcoef(0,i));
	m_line_bcoef.set_length(n);
	delete[] work;
}

//**** 3.Conversion Constructor.****

MGLBRep::MGLBRep(
	const MGPPRep& pprep	//PP-rep
//Convert PP-Rep to B-rep.
):MGCurve(){
	int i,j;

	const int ncd=pprep.sdim();
	const int k=pprep.order();	//order.
	const int l=pprep.nbreak()-1;//Number of spans.
	int n=l+k-1;					//Initial B-Rep dimension.
	MGPPRep pp(pprep); pp.normalize();
	//Generate initial knot vector.
	MGKnotVector t(k,n);
	for(i=0; i<k; i++) {
		t(i)=pp.break_point(0);
		t(i+n)=pp.break_point(l);
	}
	j=1; i=k;
	while(i<n) t(i++)=pp.break_point(j++);
	//Add multiple knot, taking continuity into account.
	MGVector v1,v2; int imk,imkp1; double brk;
	for(i=k; i<n; i++){
		imk=i-k; imkp1=imk+1;
		brk=pp.break_point(imkp1);
		v1=pp.eval_i(imk,brk,0);
		v2=pp.coef(0,imkp1);
		j=0;
		if(v1==v2){
			for(j=1; j<=k-2; j++){
				v1=pp.eval_i(imk,brk,j);
				v2=pp.coef(j,imkp1);
				if(v1!=v2) break;
			}
		}
		int mult=k-1-j;	 //// mult=k-j
		if(mult>0)		 //// mult>1
			t=MGKnotVector(t,MGKnotArray(brk,mult));
	}
	buildLBRepFromPPRep(std::move(t), pp);
}

///Gets new B-Rep by subdividing the original one into a part.
///New one is exactly the same as the original except that it is partial.
///id1 and id2 are id's of this knot_vector, and indicate the parameter range
///as from t[id1] to t[id2]. Here t=this->knot_vector().
///shrinkToKnots() employs the partial knot vector of t and this B-coefficients.
///And so, knot multiplicity of start and end of the new knot vector is not guaranteed.
///It depends on the original one.
void MGLBRep::shrinkToKnots(
	int id1,///< start id of this knot vector.
	int id2,///< End id of this knot vector.
	MGLBRep& NewBrep///<shrinked B-Rep@is output, can be this.
)const{
	if(id1>=id2){
		int idSave=id1;
		id1=id2;
		id2=idSave;
	}
	int k=order();
	int km1=k-1;
	if(id1<km1)
		id1=km1;
	int n=bdim();
	if(id2>n)
		id2=n;

	const MGKnotVector& t=knot_vector();
	const MGBPointSeq& bcoef=line_bcoef();
	while(id1<n && t[id1]==t[id1+1]){
		id1++;
	}
	while(id2>=k && t[id2]==t[id2-1]){
		id2--;
	}
	MGKnotVector t2;
	MGBPointSeq b2;
	if(id1>=id2 || t[id1]>=t[id2]){
		//The case of illegal input data.
		t2=t;
		b2=bcoef;
	}else{
		int startID=id1-km1;
		int newBdim=id2-startID;
		t2=MGKnotVector(startID, newBdim, t);
		b2=MGBPointSeq(startID, newBdim, bcoef);
	}
	NewBrep.buildLBRepFromMemberData(std::move(t2), std::move(b2));
	NewBrep.invalidateBox();
	NewBrep.copy_appearance(*this);
}

///Gets new B-Rep by computing a part of the original. New one is exactly
///the same as the original except that it is partial.
///If multiple==true(!=0), knot(i)=t1 and knot(n+i)=t2 for i=0,..., k-1
///are guaranteed. Here, n=bdim() and k=order().
///Both t1 and t2 must be inside te range of this.
void MGLBRep::shrinkToParameters(
	double t1, double t2,	//New parameter range. t1 must be less than t2.
	MGLBRep& NewBrep,	///<subdivided B-Rep is output.
	int multiple   //Indicates if start and end knot multiplicities
					//are necessary. =0:unnecessary, !=0:necessary.
//If multiple==true(!=0), knot(i)=t1 and knot(n+i)=t2 for i=0,..., k-1 will be
//guaranteed.
//Both t1 and t2 must be inside the range of old_brep.
)const{
	MGKnotVector kv2;
	MGBPointSeq b2;

	int k=order(), n1=bdim(), ncd=sdim();
	if(t1>=t2){
		double tm=(t1+t2)*.5;
		kv2.size_change(1,1);//order=1, brep dimension=1.
		b2.resize(1,ncd);
		kv2(0)=kv2(1)=tm;
		b2.store_at(0,eval(tm));
	}else{
		kv2.size_change(k, n1);
		b2.resize(n1, ncd);

		int i1=m_knot_vector.locate(t1), i2=m_knot_vector.locate(t2, 1);
		double error=m_knot_vector.param_error();
		double ti1=m_knot_vector[i1];
		if((t1-ti1)<=error) t1=ti1;

		double ti2p1=m_knot_vector[i2+1];
		if((ti2p1-t2)<=error)
			t2=ti2p1;

		const int irc1=m_line_bcoef.capacity();
		const int irc2=n1;

		int n2;
		double* work=new double[k*k];
		bluprt_(k, n1, knot_data(), coef_data(),
				irc1, ncd, t1, t2, irc2, work, &n2, &kv2(0), &b2(0, 0)
				, multiple);
		delete[] work;
		kv2.set_bdim(n2);
		b2.set_length(n2);
	}
	NewBrep.buildLBRepFromMemberData(std::move(kv2), std::move(b2));
	NewBrep.invalidateBox();
	NewBrep.copy_appearance(*this);
}

// Construct a Line B-Rep by changing space dimension and order of coordinate.
MGLBRep::MGLBRep(
		int dim,					// New space dimension.
		const MGLBRep& old_brep,	// Original Line B-rep.
		int start1, 				// Destination order of new line.
		int start2	 			// Source order of original line.
):MGCurve(old_brep)
	,m_knot_vector(old_brep.knot_vector())
	,m_line_bcoef(dim,old_brep.m_line_bcoef,start1,start2)
{
	invalidateBox();
}

///Approximate this curve as a MGLBRep curve
///within the tolerance MGTolerance::line_zero().
///When parameter_normalization=0, reparameterization will not done, and
///the evaluation at the same parameter has the same values before and after
///of approximate_as_LBRep.
void MGLBRep::approximate_as_LBRep(
	MGLBRep& lb,///<Approximated lbrep is set, can be this.
	int ordr,	///<new order. When this is MGLBRep, if ordr=0,
				///ordr=order() will be assumed, else ordr=4 is assumed.
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
	int k=order();
	if(!ordr)
		ordr=k;

	if(ordr==2 && k==ordr && !parameter_normalization){
		lb=*this;
		return;
	}

	int km1=k-1, nold=bdim();
	int noldpk=nold+k;
	double ts=param_s(), te=param_e();//Save the original parameter range.
	double mZero=MGTolerance::mach_zero();

	MGLBRep lb2(*this);
	MGKnotVector& toriginal=lb2.knot_vector();
	int start=k, index=k;
	bool modified=false;
	while(index<nold){	//Count the new B-Rep dimension.
		int multi_found=toriginal.locate_multi(start,km1,index);//Get the next multiplicity.
		if(multi_found){
			start=index+multi_found;
			double tfound=toriginal(index);
			double a=lb2.eval(tfound,1,1).len();//Get left continuous deriv at multi_found.
			double b=lb2.eval(tfound,1).len();//Get left continuous deriv at multi_found.
			if(a<=mZero || b<=mZero)
				continue;

			//Make the tangent length equal at multi_found. 
			double tanRatio=b/a;
			if(MGREqual(tanRatio, 1.))
				continue;

			int indexPmulti=index+multi_found;
			double spanOld=toriginal(indexPmulti)-tfound;
			double sdiff=spanOld*(tanRatio-1.);
			for(int i=indexPmulti; i<noldpk; i++)
				toriginal(i)+=sdiff;
			modified=true;
		}
	}

	lb2.approximate_as_LBRep2(lb,ordr,km1,nold,neglectMulti);
	if(parameter_normalization){	
		ts=lb.param_s(), te=lb.param_e();
		te= parameter_normalization==1 ? 1. : (te-ts)*lb.get_average_tangent_length();
		modified=true;
		ts=0.;
	}
	if(modified)
		lb.change_range(ts,te);			
}
