/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Box.h"
#include "mg/Vector.h"
#include "mg/Position.h"
#include "mg/KnotVector.h"
#include "mg/LBRep.h"
#include "mg/SPointSeq.h"
#include "mg/Surface.h"
#include "mg/SBRepEndC.h"
#include "mg/SBRepTP.h"
#include "mg/RSBRep.h"
#include "mg/Tolerance.h"

#include "cskernel/Blgi2d.h"
#include "cskernel/Bvstan.h"
#include "cskernel/Bsunk.h"
#include "cskernel/Bsudec.h"
#include "cskernel/Bsepl.h"
#include "cskernel/blg4sp2.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGSBRep.cpp
//
// Implements Surface B-Representation class MGSBRep.

//Add start and end data points of u and v direction to the original
//necessary for the End Condition.
void addDataPoint1direction(
	bool alongU,//=true when u direction
	const MGSBRepEndC& endc,	//End Condition of the SBRep.
	const MGNDDArray& taui,	// Original data points of u
								// without End Condition.
	MGNDDArray& tau			// u data point will be output.
){
	MGENDCOND ec[2];
	if(alongU){
		ec[0]=endc.cond(3); ec[1]=endc.cond(1);
	}else{
		ec[0]=endc.cond(0); ec[1]=endc.cond(2);
	}
	int len1,len;
	len1=len=taui.length();
	int is=0;

	//tau=min condition.
	if(ec[0]== MGENDCOND::MGENDC_1D || ec[0]== MGENDCOND::MGENDC_2D)
		is=1;
	else if(ec[0]== MGENDCOND::MGENDC_12D)
		is=2;
	len+=is;

	//tau=max condition.
	if(ec[1]== MGENDCOND::MGENDC_1D || ec[1]== MGENDCOND::MGENDC_2D)
		len+=1;
	else if(ec[1]== MGENDCOND::MGENDC_12D)
		len+=2;

	tau.resize(len);
	int i=0;
	for(; i<is; i++)//At the starting points.
		tau(i)=taui(0);
	for(int j=0; j<len1; j++)//Middle points.
		tau(i++)=taui(j);
	for(;i<len; i++)
		tau(i)=taui(len1-1);//Ending points.
	tau.set_length(len);
}

//Add start and end data points of u and v direction to the original
//necessary for the End Condition.
void add_data_point(
	const MGSBRepEndC& endc,	//End Condition of the SBRep.
	const MGNDDArray& utaui,	// Original data points of u
								// without End Condition.
	const MGNDDArray& vtaui,	// .... of v.
	MGNDDArray& utau,			// u data point will be output.
	MGNDDArray& vtau			// v data point will be output.
){
	addDataPoint1direction(true, endc, utaui, utau);
	addDataPoint1direction(false, endc, vtaui, vtau);
}

///Make m_uknot(alongU=true) or m_vknot of this from MGSBRepEndC, data points, and an order.
void MGSBRep::makeKnotVectorFromEC1Dire(
	bool alongU,
	const MGSBRepEndC & endc,
	const MGNDDArray & tau,
	int order
){
	MGNDDArray tau2;
	addDataPoint1direction(alongU,endc,tau,tau2);
	int n=tau2.length();
	if(order>n)
		order=n;
	MGKnotVector& t=alongU ? m_uknot:m_vknot;
	t = MGKnotVector(tau2, order);
}

///Make (m_uknot, m_vknot) of this from MGSBRepEndC, data points, and orders.
void MGSBRep::makeKnotVectorFromEC(
	const MGSBRepEndC & endc,//End conditions.
	const MGNDDArray & utau, const MGNDDArray & vtau,// data point seq of (u,v).
	int orderu, int orderv //order of (u,v).
){
	makeKnotVectorFromEC1Dire(true, endc, utau, orderu);
	makeKnotVectorFromEC1Dire(false, endc, vtau, orderv);
}

//Member Data
//	MGKnotVector m_uknot;			// Knot Vector of u-direction
//	MGKnotVector m_vknot;			// Knot Vector of v-direction
//	MGSPointSeq  m_surface_bcoef;	// Surface B-Coef.

//<< Constructor >>

//**** 1. Interpolation ****

// Construct Surface B-rep by intepolation from Point data only.
// Function's return value is: =0:success, !=0: failure.
//          =12: too near points included along u direction,
//          =13: too near points included along v direction
//When failure is returned(points are illegal(e.g. same points are input)),
//uniform BSpline of order 2 is built(surface B-coefficients are input points).
int MGSBRep::buildByInterpolation(
	const MGSPointSeq& points,	//Point seq data
	int orderu,		// Order of u-direction
	int orderv		// Order of v-direction
){
	assert(points.length_u()>=2 && points.length_v()>=2);
	m_surface_bcoef.resize(points.length_u(), points.length_v(), points.sdim());
	invalidateBox();

	MGNDDArray utau, vtau;
	int lenu=points.length_u(), lenv=points.length_v();
	int lenuv=lenu*lenv;
	int len=(lenu>=lenv) ? lenu:lenv;
	if(lenu<orderu)
		orderu=lenu;
	if(lenv<orderv)
		orderv=lenv;

	int order=(orderu>=orderv) ? orderu:orderv;
	int usize=points.capacity_u();
	int ncd=points.sdim();
	double* work=new double[len+lenuv*ncd+len*(2*order-1)];
	double* work2=work+len;
	double* q=work2+lenuv*ncd;

	compute_knot(points, orderu, orderv, utau, vtau);
	int error=2;
	for(int k=0; k<ncd; k++){
		blgi2d_(&error,utau.data(),points.data(0,0,k),m_uknot.data(),
			orderu,lenu,lenv,usize,lenv,work,q,work2+lenuv*k);
		if(error!=1)
			break;
	}
	if(error==1){
		error=2;
		for(int k=0; k<ncd; k++){
			blgi2d_(&error,vtau.data(),work2+lenuv*k,m_vknot.data(),
			orderv,lenv,lenu,lenv,lenu,work,q,&m_surface_bcoef(0,0,k));
			if(error!=1)
				break;
		}
		if(error!=1) error=13;	//Error detected in v-direction.
	}else
		error=12;				//Error detected in u-direction.
	if(error==1){
		error=0;
		m_surface_bcoef.set_length(lenu, lenv);
	}else{
		m_uknot=MGKnotVector(2,lenu,0.);
		m_vknot=MGKnotVector(2,lenv,0.);
		m_surface_bcoef=points;
	}
	delete[] work;
	return error;
}

// Construct Surface B-rep of any order number by interpolation 
//from data point and point data with knot vector.
// Function's return value is: =0:success, !=0: failure.
//          =12: too near data points included in utau,
//          =13: too near data points included in vtau.
//When failure is returned,
//uniform BSpline of order 2 is built(B-coefficients are input points).
int MGSBRep::buildByInterpolationWithKTV(
	const MGNDDArray& utau,		//Data point abscissa of u-direction
	const MGNDDArray& vtau,	//Data point abscissa of v-direction
	const MGSPointSeq& points	//Point seq data
){
	assert(utau.length()==points.length_u() && points.length_u()==m_uknot.bdim());
	assert(vtau.length()==points.length_v() && points.length_v()==m_vknot.bdim());
	invalidateBox();

	m_surface_bcoef.resize(utau.length(), vtau.length(), points.sdim());
	int orderu=m_uknot.order(), orderv=m_vknot.order();
	int lenu=points.length_u(), lenv=points.length_v();
	int lenuv=lenu*lenv;
	int len=lenu; if(len<lenv) len=lenv;
	int order=orderu; if(order<orderv) order=orderv;
	int usize=points.capacity_u();
	int ncd=points.sdim();
	double* work=new double[len+lenuv*ncd+len*(2*order-1)];
	double* work2=work+len;
	double* q=work2+lenuv*ncd;

	int error=2;
	for(int k=0; k<ncd; k++){
		blgi2d_(&error,utau.data(),points.data(0,0,k),m_uknot.data(),
			orderu,lenu,lenv,usize,lenv,work,q,work2+lenuv*k);
		if(error!=1)
			break;
	}
	if(error==1){
		error=2;
		for(int k=0; k<ncd; k++){
			blgi2d_(&error,vtau.data(),work2+lenuv*k,m_vknot.data(),
			orderv,lenv,lenu,lenv,lenu,work,q,&m_surface_bcoef(0,0,k));
			if(error!=1)
				break;
		}
		if(error!=1)
			error=13;	//Error detected in v-direction.
	}else
		error=12;				//Error detected in u-direction.
	if(error==1){
		error=0;
		m_surface_bcoef.set_length(lenu, lenv);
	}else{
		m_uknot=MGKnotVector(2,lenu,0.);
		m_vknot=MGKnotVector(2,lenv,0.);
		m_surface_bcoef=points;
	}
	delete[] work;
	return error;
}

// Construct Surface B-rep of order 4 by interpolation from Point data
//and tangent plane end condition.
// Inner point must be homogeneous, shoule not include first derivative inf.
// Function's return value is: =0:success, !=0: failure.
//          =12: too near data points included in utau,
//          =13: too near data points included in vtau.
//When failure is returned,
//uniform BSpline of order 2 is built(B-coefficients are input points).
int MGSBRep::buildByInterpolationTPCornerDeriv(			//Tangent Plane Version
	const MGSBRepTP& tp,		//end condition	of Tangent Plane
	const MGNDDArray& utau,		//Data point of u-direction of value
	const MGNDDArray& vtau,		//Data point of v-direction	of value
	const MGVector   uvec[4],	//Tangent vector of u-direction for 
						// v=min and v=max boundary line.
				//uvec[0], uvec[1]: start and end of v=min line
				//uvec[2], uvec[3]: end and start of v=max line.
	const MGVector   vvec[4],	//Tangent vector of v-direction for 
						// u=min and u=max boundary line.
				//vvec[0]:start of u=min line, vvec[1]:start of u=max line
				//vvec[2]: end  of u=max line, vvec[3]: end  of u=min.
		// It is allowed that uvec[i] or vvec[j] is null vector(i.e. sdim==0).
		//When this is the case, it means no inf provided
		//about the tangent vector.
	const MGSPointSeq& points	//Data point ordinate
){
	assert(utau.length()==points.length_u()
		&& vtau.length()==points.length_v());
	assert(points.sdim()==3);

	const int dim=3; int m;

	int i,k, iend;
	int lenu=utau.length(), lenv=vtau.length();
	int lenum1=lenu-1, lenvm1=lenv-1;
	double ur[2], vr[2];
	ur[0]=utau(0); ur[1]=utau(lenum1);
	vr[0]=vtau(0); vr[1]=vtau(lenvm1);
	int ktp, ntp, n;
	const double* ttp = 0; const double* rtp = 0;
	int itanse[2]; double tanse[6];
	int irtp, ip1, ip2, psizeu, psizev;
	points.capacity(psizeu, psizev);
	double work[105], tangen[3];
	MGSBRepEndC endc;

	int perimeter; double tptau;
	const int ipse[2]={1,2};
	int i1,i2;

	const int iu[2]={0,2}; const int ivvec[4]={0,1,3,2};
	n=lenv; ip1=psizeu; ip2=psizev;
	for(m=0; m<2; m++){
		//Process of perimeter num 0 and 2(v=min and max line)

	perimeter=iu[m]; i1=ivvec[2*m]; i2=ivvec[2*m+1];
	if(tp.specified(perimeter)){
		ktp=tp.TP(perimeter).order();
		ntp=tp.TP(perimeter).bdim();
		ttp=tp.TP(perimeter).knot_data();
		rtp=tp.TP(perimeter).coef_data();
		irtp=tp.TP(perimeter).line_bcoef().capacity();
	}else
		ntp=0;
	if(ntp || vvec[i1].sdim() || vvec[i2].sdim()){
	//We have to generate first derivative data for this perimeter.
		MGBPointSeq first(lenu, dim);
		i=0; itanse[0]=2; itanse[1]=2;
		if(vvec[i1].sdim()){
			for(k=0; k<dim; k++) tanse[k]=first(0,k)=vvec[i1].ref(k);
			i=1; itanse[0]=1;
		}
		iend=lenu;
		if(vvec[i2].sdim()){
			for(k=0; k<dim; k++) tanse[k+3]=first(lenum1,k)=vvec[i2].ref(k);
			iend=lenum1; itanse[1]=1;
		}
		for(; i<iend; i++){
			tptau=utau(i);
			bvstan_(ur,ktp,ntp,ttp,rtp,tptau,n,vtau.data(),
				points.data(i,0,0),ipse[m],itanse,tanse,irtp,
				ip1,ip2,work,tangen);
			for(k=0; k<dim; k++) first(i,k)=tangen[k];
		}
		endc.set_1st(perimeter,std::move(first));
	}

	}

	const int iv[2]={3,1}; const int iuvec[4]={0,3,1,2};
	n=lenu; ip1=1; ip2=psizeu*psizev;
	for(m=0; m<2; m++){
		//Process of perimeter num 3 and 1(u=min and max line)

	perimeter=iv[m]; i1=iuvec[2*m]; i2=iuvec[2*m+1];
	if(tp.specified(perimeter)){
		ktp=tp.TP(perimeter).order();
		ntp=tp.TP(perimeter).bdim();
		ttp=tp.TP(perimeter).knot_data();
		rtp=tp.TP(perimeter).coef_data();
		irtp=tp.TP(perimeter).line_bcoef().capacity();
	}else
		ntp=0;

	if(ntp || uvec[i1].sdim() || uvec[i2].sdim()){
	//We have to generate first derivative data for this perimeter.
		MGBPointSeq first(lenv, dim);
		i=0; itanse[0]=2; itanse[1]=2;
		if(uvec[i1].sdim()){
			for(k=0; k<dim; k++) tanse[k]=first(0,k)=uvec[i1].ref(k);
			i=1; itanse[0]=1;
		}
		iend=lenv;
		if(uvec[i2].sdim()){
			for(k=0; k<dim; k++) tanse[k+3]=first(lenvm1,k)=uvec[i2].ref(k);
			iend=lenvm1; itanse[1]=1;
		}
		for(; i<iend; i++){
			tptau=vtau(i);
			bvstan_(vr,ktp,ntp,ttp,rtp,tptau,n,utau.data(),
				points.data(0,i,0),ipse[m],itanse,tanse,irtp,
				ip1,ip2,work,tangen);
			for(k=0; k<dim; k++) first(i,k)=tangen[k];
		}
		endc.set_1st(perimeter,std::move(first));
	}

	}
	return buildByInterpolationEC(endc,utau,vtau,points);
}

// Construct Surface B-rep of specified order by interpolation from Point data.
///This is an easy-to-use-version of buildByInterpolationEC.
// Function's return value is: =0:success, !=0: failure.
//          =12: too near data points included in utau,
//          =13: too near data points included in vtau.
//When failure is returned,
//uniform BSpline of order 2 is built(B-coefficients are input points).
int MGSBRep::buildByInterpolationDataPoints(				//Derivative Inf.
	const MGNDDArray& utau,	//Data point of u-direction of value
	const MGNDDArray& vtau,	//Data point of v-direction	of value
	const MGSPointSeq& value,	//Data point ordinate
	int orderu, int orderv//order along u or v direction.
){
	assert(utau.length()==value.length_u()
		&& vtau.length()==value.length_v());

	if(utau.length()<orderu) orderu=utau.length();
	if(vtau.length()<orderv) orderv=vtau.length();
	m_uknot = MGKnotVector(utau, orderu), m_vknot = MGKnotVector(vtau, orderv);
	MGSBRepEndC endc;
	return buildByInterpolationECWithKTV(endc, utau, vtau, value);
}

// Construct Surface B-rep of specified order by interpolation from Point data
///and 4 corner end condition.
//Inner points must not include first derivative inf if the corresponding data
//points are multiple if the order is not 4.
// Function's return value is: =0:success, !=0: failure.
//          =12: too near data points included in utau,
//          =13: too near data points included in vtau.
//When failure is returned,
//uniform BSpline of order 2 is built(B-coefficients are input points).
int MGSBRep::buildByInterpolationEC(				//Derivative Inf.
	MGSBRepEndC& endc,		///< corner end condition
	const MGNDDArray& utaui,	///<Data point of u-direction of value
	const MGNDDArray& vtaui,	///<Data point of v-direction	of value
	const MGSPointSeq& value,	///<Data point ordinate
	int orderu,///<order along u direction,
	int orderv///<order along v direction,
){
	assert(utaui.length()==value.length_u()
		&& vtaui.length()==value.length_v());
	assert(orderu>=2 && orderv>=2);
	makeKnotVectorFromEC(endc, utaui, vtaui, orderu, orderv);
	return buildByInterpolationECWithKTV(endc, utaui, vtaui, value);
}

// Construct Surface B-rep of order 4 by interpolation from Point data
//and end condition with knot vector.
// Inner point may include derivative inf.
// Function's return value is: =0:success, !=0: failure.
//          =12: too near data points included in utau,
//          =13: too near data points included in vtau.
//When failure is returned,
//uniform BSpline of order 2 is built(B-coefficients are input points).
int MGSBRep::buildByInterpolationECWithKTV(//blg4sp2_
	MGSBRepEndC& endc,//end condition
	const MGNDDArray& utaui,	//Data point of u-direction
	const MGNDDArray& vtaui,	//Data point of v-direction
	const MGSPointSeq& value	//Data point ordinate
){
	assert(utaui.length() == value.length_u());
	assert(vtaui.length() == value.length_v());
	invalidateBox();

	int orderu = m_uknot.order(), orderv = m_vknot.order();
	const int ncd = value.sdim();
	m_surface_bcoef.resize(m_uknot.bdim(), m_vknot.bdim(), ncd);
	endc.complete_corner_deriv(utaui, vtaui);	//Twist or other data.

	int i,j,k;

	MGENDCOND ec[4];
	int ecint[4];
	for (i = 0; i < 4; i++) {
		ec[i] = endc.cond(i);
		ecint[i] = static_cast<int>(ec[i]);
	}

	int lenu1,lenu, lenv1,lenv;
	lenu1=lenu=value.length_u(); lenv1=lenv=value.length_v();
	int ius=0,ivs=0;
	int lenum1,lenvm1, lenum2,lenvm2, lenum3,lenvm3;

// Compute b-rep dimension along u and v in lenu and lenv.
	//u=min and max condition.
	if(ec[3]== MGENDCOND::MGENDC_1D || ec[3]== MGENDCOND::MGENDC_2D){
		lenu+=1; ius=1;
	}else if(ec[3]== MGENDCOND::MGENDC_12D){
		lenu+=2; ius=2;
	}
	if(ec[1]== MGENDCOND::MGENDC_1D || ec[1]== MGENDCOND::MGENDC_2D)
		lenu+=1;
	else if(ec[1]== MGENDCOND::MGENDC_12D)
		lenu+=2;

	//v=min and max condition.
	if(ec[0]== MGENDCOND::MGENDC_1D || ec[0]== MGENDCOND::MGENDC_2D){
		lenv+=1; ivs=1;
	}else if(ec[0]== MGENDCOND::MGENDC_12D){
		lenv+=2; ivs=2;
	}
	if(ec[2]== MGENDCOND::MGENDC_1D || ec[2]== MGENDCOND::MGENDC_2D)
		lenv+=1;
	else if(ec[2]== MGENDCOND::MGENDC_12D)
		lenv+=2;

	assert(lenu== m_uknot.bdim() && lenv== m_vknot.bdim());

// Prepare data point ordinate.
	// 1. Copy original data.
	int is, js;
	for(k=0; k<ncd; k++){
		for(j=0; j<lenv1; j++){
			js=ivs+j;
			for(i=0; i<lenu1; i++){
				is=ius+i;
				m_surface_bcoef(is,js,k)=value(i,j,k);
			}
		}
	}
	// 2. Copy derivative data along perimeter from endc.
	lenum1=lenu-1; lenum2=lenu-2; lenum3=lenu-3;
	lenvm1=lenv-1; lenvm2=lenv-2; lenvm3=lenv-3;
	//u=min condition.
	if(ec[3]== MGENDCOND::MGENDC_1D || ec[3]== MGENDCOND::MGENDC_12D){
		for(k=0; k<ncd; k++){
			for(j=0; j<lenv1; j++){
				js=ivs+j;
				m_surface_bcoef(0,js,k)=(endc.first(3))(j,k);
			}
		}
	}
	if(ec[3]== MGENDCOND::MGENDC_2D || ec[3]== MGENDCOND::MGENDC_12D){
		if(ec[3]== MGENDCOND::MGENDC_2D) i=0; else i=1;
		for(k=0; k<ncd; k++){
			for(j=0; j<lenv1; j++){
				js=ivs+j;
				m_surface_bcoef(i,js,k)=(endc.second(3))(j,k);
			}
		}
	}
	//u=max condition.
	if(ec[1]== MGENDCOND::MGENDC_1D || ec[1]== MGENDCOND::MGENDC_12D){
		for(k=0; k<ncd; k++){
			for(j=0; j<lenv1; j++){
				js=ivs+j;
				m_surface_bcoef(lenum1,js,k)=(endc.first(1))(j,k);
			}
		}
	}
	if(ec[1]== MGENDCOND::MGENDC_2D || ec[1]== MGENDCOND::MGENDC_12D){
		if(ec[1]== MGENDCOND::MGENDC_2D) i=lenum1; else i=lenum2;
		for(k=0; k<ncd; k++){
			for(j=0; j<lenv1; j++){
				js=ivs+j;
				m_surface_bcoef(i,js,k)=(endc.second(1))(j,k);
			}
		}
	}

	//v=min condition.
	if(ec[0]== MGENDCOND::MGENDC_1D || ec[0]== MGENDCOND::MGENDC_12D){
		for(k=0; k<ncd; k++){
			for(i=0; i<lenu1; i++){
				is=ius+i;
				m_surface_bcoef(is,0,k)=(endc.first(0))(i,k);
			}
		}
	}
	if(ec[0]== MGENDCOND::MGENDC_2D || ec[0]== MGENDCOND::MGENDC_12D){
		if(ec[0]== MGENDCOND::MGENDC_2D) j=0; else j=1;
		for(k=0; k<ncd; k++){
			for(i=0; i<lenu1; i++){
				is=ius+i;
				m_surface_bcoef(is,j,k)=(endc.second(0))(i,k);
			}
		}
	}
	//v=max condition.
	if(ec[2]== MGENDCOND::MGENDC_1D || ec[2]== MGENDCOND::MGENDC_12D){
		for(k=0; k<ncd; k++){
			for(i=0; i<lenu1; i++){
				is=ius+i;
				m_surface_bcoef(is,lenvm1,k)=(endc.first(2))(i,k);
			}
		}
	}
	if(ec[2]== MGENDCOND::MGENDC_2D || ec[2]== MGENDCOND::MGENDC_12D){
		if(ec[2]== MGENDCOND::MGENDC_2D) j=lenvm1; else j=lenvm2;
		for(k=0; k<ncd; k++){
			for(i=0; i<lenu1; i++){
				is=ius+i;
				m_surface_bcoef(is,j,k)=(endc.second(2))(i,k);
			}
		}
	}

	//3. Copy corner derivative inf.
	int i1[4]={lenum1,lenum1,0,0}, j1[4]={0,lenvm1,lenvm1,0},
		i2[4]={lenum2,lenum2,1,1}, j2[4]={1,lenvm2,lenvm2,1};

	int idu,idv;  int m,mp1;
	for(m=0; m<4; m++){
		mp1=m+1; if(mp1==4) mp1=0;
		if(m==0 || m==2) {idu=m; idv=mp1;}
		else             {idu=mp1; idv=m;}

		//Copy twist vectors( d2f/(du*dv) ).
		if((ec[idu]== MGENDCOND::MGENDC_1D || ec[idu]== MGENDCOND::MGENDC_12D) &&
		   (ec[idv]== MGENDCOND::MGENDC_1D || ec[idv]== MGENDCOND::MGENDC_12D)){
			for(k=0; k<ncd; k++)
				m_surface_bcoef(i1[m],j1[m],k)=endc.deriv11(mp1).ref(k);
		}
		//                    d3f/(du**2*dv).
		if((ec[idu]== MGENDCOND::MGENDC_1D || ec[idu]== MGENDCOND::MGENDC_12D) &&
		   (ec[idv]== MGENDCOND::MGENDC_2D || ec[idv]== MGENDCOND::MGENDC_12D)){
			if(ec[idv]== MGENDCOND::MGENDC_2D) i=i1[m]; else i=i2[m];
			for(k=0; k<ncd; k++)
				m_surface_bcoef(i,j1[m],k)=endc.deriv21(mp1).ref(k);
		}
		//                    d3f/(du*dv**2).
		if((ec[idu]== MGENDCOND::MGENDC_2D || ec[idu]== MGENDCOND::MGENDC_12D) &&
		   (ec[idv]== MGENDCOND::MGENDC_1D || ec[idv]== MGENDCOND::MGENDC_12D)){
			if(ec[idu]== MGENDCOND::MGENDC_2D) j=j1[m]; else j=j2[m];
			for(k=0; k<ncd; k++)
				m_surface_bcoef(i1[m],j,k)=endc.deriv12(mp1).ref(k);
		}
		//                    d4f/(du**2*dv**2).
		if((ec[idu]== MGENDCOND::MGENDC_2D || ec[idu]== MGENDCOND::MGENDC_12D) &&
		   (ec[idv]== MGENDCOND::MGENDC_2D || ec[idv]== MGENDCOND::MGENDC_12D)){
			if(ec[idv]== MGENDCOND::MGENDC_2D) i=i1[m]; else i=i2[m];
			if(ec[idu]== MGENDCOND::MGENDC_2D) j=j1[m]; else j=j2[m];
			for(k=0; k<ncd; k++)
				m_surface_bcoef(i,j,k)=endc.deriv22(mp1).ref(k);
		}
	}

	//Exchange positional data for blg4sp2_.
	//Perimeter 0.
	double save;
	if(ec[0]== MGENDCOND::MGENDC_1D || ec[0]== MGENDCOND::MGENDC_2D){
		for(i=0; i<lenu; i++)
			for(k=0; k<ncd; k++){
				save=m_surface_bcoef(i,1,k);
				m_surface_bcoef(i,1,k)=m_surface_bcoef(i,0,k);
				m_surface_bcoef(i,0,k)=save;
			}
	} else if(ec[0]== MGENDCOND::MGENDC_12D){
		for(i=0; i<lenu; i++)
			for(k=0; k<ncd; k++){
				save=m_surface_bcoef(i,2,k);
				m_surface_bcoef(i,2,k)=m_surface_bcoef(i,1,k);
				m_surface_bcoef(i,1,k)=m_surface_bcoef(i,0,k);
				m_surface_bcoef(i,0,k)=save;
			}
	}
	//Perimeter 2.
	if(ec[2]== MGENDCOND::MGENDC_1D || ec[2]== MGENDCOND::MGENDC_2D){
		for(i=0; i<lenu; i++)
			for(k=0; k<ncd; k++){
				save=m_surface_bcoef(i,lenvm2,k);
				m_surface_bcoef(i,lenvm2,k)=m_surface_bcoef(i,lenvm1,k);
				m_surface_bcoef(i,lenvm1,k)=save;
			}
	} else if(ec[2]== MGENDCOND::MGENDC_12D){
		for(i=0; i<lenu; i++)
			for(k=0; k<ncd; k++){
				save=m_surface_bcoef(i,lenvm3,k);
				m_surface_bcoef(i,lenvm3,k)=m_surface_bcoef(i,lenvm2,k);
				m_surface_bcoef(i,lenvm2,k)=m_surface_bcoef(i,lenvm1,k);
				m_surface_bcoef(i,lenvm1,k)=save;
			}
	}
	//Perimeter 1.
	if(ec[1]== MGENDCOND::MGENDC_1D || ec[1]== MGENDCOND::MGENDC_2D){
		for(j=0; j<lenv; j++)
			for(k=0; k<ncd; k++){
				save=m_surface_bcoef(lenum2,j,k);
				m_surface_bcoef(lenum2,j,k)=m_surface_bcoef(lenum1,j,k);
				m_surface_bcoef(lenum1,j,k)=save;
			}
	} else if(ec[1]== MGENDCOND::MGENDC_12D){
		for(j=0; j<lenv; j++)
			for(k=0; k<ncd; k++){
				save=m_surface_bcoef(lenum3,j,k);
				m_surface_bcoef(lenum3,j,k)=m_surface_bcoef(lenum2,j,k);
				m_surface_bcoef(lenum2,j,k)=m_surface_bcoef(lenum1,j,k);
				m_surface_bcoef(lenum1,j,k)=save;
			}
	}
	//Perimeter 3.
	if(ec[3]== MGENDCOND::MGENDC_1D || ec[3]== MGENDCOND::MGENDC_2D){
		for(j=0; j<lenv; j++)
			for(k=0; k<ncd; k++){
				save=m_surface_bcoef(1,j,k);
				m_surface_bcoef(1,j,k)=m_surface_bcoef(0,j,k);
				m_surface_bcoef(0,j,k)=save;
			}
	} else if(ec[3]== MGENDCOND::MGENDC_12D){
		for(j=0; j<lenv; j++)
			for(k=0; k<ncd; k++){
				save=m_surface_bcoef(2,j,k);
				m_surface_bcoef(2,j,k)=m_surface_bcoef(1,j,k);
				m_surface_bcoef(1,j,k)=m_surface_bcoef(0,j,k);
				m_surface_bcoef(0,j,k)=save;
			}
	}

	int lenuv=lenu*lenv;
	int len=lenu; if(len<lenv) len=lenv;
	int order=orderu; if(order<orderv) order=orderv;
	double* wk=new double[len*2+lenuv*ncd+len*(2*order-1)];
	double* wk2=wk+len*2;
	double* q=wk2+lenuv*ncd;

	MGNDDArray utau,vtau;
	add_data_point(endc,utaui,vtaui,utau,vtau);	//Compute data points.

	int error=2;
	for(k=0; k<ncd; k++){
		blg4sp2_(orderu,&error,ecint[3],ecint[1],utau.data(),m_surface_bcoef.data(0,0,k)
			,lenu,lenu,lenv,m_uknot.data(),lenv,wk,wk+len,q,wk2+lenuv*k);
		if(error!=1) break;
	}
	if(error==1){
		error=2;
		for(k=0; k<ncd; k++){
			blg4sp2_(orderv,&error,ecint[0],ecint[2],vtau.data(),wk2+lenuv*k,lenv,lenv,
			lenu,m_vknot.data(),lenu,wk,wk+len,q,&m_surface_bcoef(0,0,k));
			if(error!=1) break;
		}
		if(error!=1)
			error=13;	//Error detected in v-direction.
	}else
		error=12;				//Error detected in u-direction.

	if(error==1){
		error=0;
		m_surface_bcoef.set_length(lenu, lenv);
	}else{
		m_uknot=MGKnotVector(2,value.length_u(),0.);
		m_uknot.change_range(utau[0], utau[lenu1-1]);
		m_vknot=MGKnotVector(2,value.length_v(),0.);
		m_vknot.change_range(vtau[0], vtau[lenv1-1]);
		m_surface_bcoef=value;
	}
	delete[] wk;
	return error;
}

//Member Function

// 入力のパラメータ範囲の曲面部分を囲むボックスを返す。
MGBox MGSBRep::box_limitted(const MGBox& uvbox) const{
	double u1=uvbox(0).low_point(), u2=uvbox(0).high_point();
	double us=param_s_u(), ue=param_e_u();
	if(u1<us) u1=us; if(u2>ue) u2=ue;

	double v1=uvbox(1).low_point(), v2=uvbox(1).high_point();
	double vs=param_s_v(), ve=param_e_v();
	if(v1<vs) v1=vs; if(v2>ve) v2=ve;

	if(MGREqual_base(u1,u2,knot_vector_u().param_span())){
		MGLBRep line=parameter_line(1,(u1+u2)*.5);
		return line.box_limitted(uvbox(1));
	}else if(MGREqual_base(v1,v2,knot_vector_v().param_span())){
		MGLBRep line=parameter_line(0,(v1+v2)*.5);
		return line.box_limitted(uvbox(0));
	}
	MGSBRep temp;
	shrinkToParameters(uvbox,temp);
	return temp.m_surface_bcoef.box();
}

void MGSBRep::compute_box(MGBox& bx) const{
	m_surface_bcoef.compute_box(bx);
}

//Changing this object's space dimension.
void MGSBRep::change_dimension(
	int dim,		// new space dimension
	int start1, 		// Destination order of new object.
	int start2) 		// Source order of this object.
{
	m_surface_bcoef=MGSPointSeq(dim,surface_bcoef(),start1,start2);
	invalidateBox();
}

MGSBRep& MGSBRep::change_range(
	int is_u,				//if true, (t1,t2) are u-value. if not, v.
	double t1,				//Parameter value for the start of original. 
	double t2)				//Parameter value for the end of original. 
//Change parameter range, be able to change the direction by providing
//t1 greater than t2.
{
	double ts, te;
	if(t1>t2){
		ts=t2; te=t1;
		m_surface_bcoef.reverse(is_u);
		if(is_u) m_uknot.reverse();
		else	 m_vknot.reverse();
	}
	else{
		ts=t1; te=t2;
	}
	if(is_u) m_uknot.change_range(ts,te);
	else	 m_vknot.change_range(ts,te);
	invalidateBox();
	return *this;
}

//Construct new surface object by copying to newed area.
//User must delete this copied object by "delete".
MGSBRep* MGSBRep::clone() const{return new MGSBRep(*this);}

//Construct new surface object by changing
//the original object's space dimension.
//User must delete this copied object by "delete".
MGSBRep* MGSBRep::copy_change_dimension(
	int sdim,		// new space dimension
	int start1, 		// Destination order of new line.
	int start2 		// Source order of this line.
)const{
	return new MGSBRep(sdim,*this,start1,start2);
}

//Evaluate right continuous ndu'th and ndv'th derivative data.
//Function's return value is (d(ndu+ndv)f(u,v))/(du**ndu*dv**ndv).
// ndu=0 and ndv=0 means positional data evaluation.
MGVector MGSBRep::eval(
	double u, double v,	// Parameter value of the surface.
	int ndu,			// Order of Derivative along u.
	int ndv			// Order of Derivative along v.
	) const
{
	const int ku=order_u(), kv=order_v();
	double cu[10],cv[10];
	double *cup=cu, *cvp=cv;
	if(ku>10) cup=new double[ku];	//This is done to save "new" when ku<=10.
	if(kv>10) cvp=new double[kv];	//This is done to save "new" when kv<=10.
	const int num1=bdim_u()-1, nvm1=bdim_v()-1;
	int idu=m_uknot.eval_coef(u,cup,ndu);
	int idv=m_vknot.eval_coef(v,cvp,ndv);
	int i,j,m; int jj;
	const int ncd=sdim();
	MGVector S(ncd); double coefall,coefu;
	for(m=0; m<ncd; m++){
		coefall=0.;
		for(j=0;j<kv;j++){
			coefu=0.;
			jj=idv+j;
			for(i=0;i<ku;i++) coefu=coefu+coef(idu+i,jj,m)*cup[i];
			coefall=coefall+coefu*cvp[j];
		}
		S.set(m)=coefall;
	}
	if(ku>10) delete[] cup; if(kv>10) delete[] cvp;
	return S;
}

//Evaluate surface data.
MGVector MGSBRep::eval(
	const MGPosition& uv	// Parameter value of the surface.
	, int ndu			// Order of derivative along u.
	, int ndv			// Order of derivative along v.
) const{return eval(uv.ref(0),uv.ref(1),ndu,ndv);}

//Evaluate right continuous surface data.
//Evaluate all positional data, 1st and 2nd derivatives.
void MGSBRep::eval_all(
		double u, double v,	// Parameter value of the surface.
		MGPosition& f,			// Positional data.
		MGVector&   fu,			// df(u,v)/du
		MGVector&   fv,			// df/dv
		MGVector&   fuv,		// d**2f/(du*dv)
		MGVector&   fuu,		// d**2f/(du**2)
		MGVector&   fvv			// d**2f/(dv**2)
		) const
{
	int ku=order_u(), kv=order_v(); 
	int ku2=ku+ku, kv2=kv+kv;
	double* ucoef=new double[ku*3+kv*3];
	double* vcoef=ucoef+ku*3;
	int bdum1=bdim_u()-1, bdvm1=bdim_v()-1;
	int uid=m_uknot.eval_coef(u,ucoef,0);
	int vid=m_vknot.eval_coef(v,vcoef,0);
	m_uknot.eval_coef(u,ucoef+ku,1); m_uknot.eval_coef(u,ucoef+ku2,2);
	m_vknot.eval_coef(v,vcoef+kv,1); m_vknot.eval_coef(v,vcoef+kv2,2);
	double s,su,suu,c,vj,vj1; int i,j,k,dim=sdim();
	int ii,jj;
	MGPosition p(dim);
	MGVector pu(dim),pv(dim),puv(dim),puu(dim),pvv(dim);
	for(k=0; k<dim; k++){
		p(k)=0.0;pu.set(k)=0.0;pv.set(k)=0.0;puv.set(k)=0.0;
		puu.set(k)=0.0;pvv.set(k)=0.0;
		for(j=0; j<kv; j++){
			s=su=suu=0.0;
			jj=vid+j;
			for(i=0; i<ku; i++){
				ii=uid+i;
				c=coef(ii,jj,k);
				s=s+ucoef[i]*c;
				su=su+ucoef[i+ku]*c;
				suu=suu+ucoef[i+ku2]*c;
			}
			vj=vcoef[j]; vj1=vcoef[j+kv];
			p(k)=p(k)+vj*s;
			pu.set(k)=pu.ref(k)+vj*su;
			pv.set(k)=pv.ref(k)+vj1*s;
			puv.set(k)=puv.ref(k)+vj1*su;
			puu.set(k)=puu.ref(k)+vj*suu;
			pvv.set(k)=pvv.ref(k)+vcoef[j+kv2]*s;
		}
	}
	f=p;fu=pu;fv=pv;fuv=puv;fuu=puu;fvv=pvv;
	delete[] ucoef;
}

//Evaluate all of derivative data (d(i+j)f(u,v))/(du**i*dv**j),
//for 0<=i<=ndu and 0<=j<=ndv.
void MGSBRep::eval_all(
		double u, double v,	// Parameter value of the surface.
		int ndu,		//Order of Derivative along u.
		int ndv,		//Order of Derivative along v.
		double* deriv	//Output. (d(i+j)f(u,v))/(du**i*dv**j) in
						//deriv[r+j*dim+i*(ndv+1)*dim] for 0<=r<dim=sdim().
						//for 0<=i<=ndu and 0<=j<=ndv.
	//deriv is an array of deriv[ndu+1][ndv+1][r],
	//(d(i+j)f(u,v))/(du**i*dv**j) is returned in deriv[i][j][r].
	) const{
	int i,j,jj,r, dim=sdim();
	int ider,jder;
	const int ku=order_u(), kv=order_v();
	int kundu=ku*(ndu+1), kvndv=kv*(ndv+1);
	double cu[30],cv[30]; double *cup=cu, *cvp=cv;
	if(kundu>30) cup=new double[kundu];	//Done to save "new".
	if(kvndv>30) cvp=new double[kvndv];	//Done to save "new".
	int idu,idv;
	for(ider=0; ider<=ndu; ider++) idu=m_uknot.eval_coef(u,cup+ider*ku,ider);
	for(jder=0; jder<=ndv; jder++) idv=m_vknot.eval_coef(v,cvp+jder*kv,jder);

	double coefall,coefu;
	for(ider=0; ider<=ndu;ider++){
		int iderndv1dim=ider*(ndv+1)*dim;
		for(jder=0; jder<=ndv; jder++){
			int jderdim=jder*dim;
			for(r=0; r<dim; r++){
				coefall=0.;
				for(j=0;j<kv;j++){
					coefu=0.;
					jj=idv+j;
					for(i=0;i<ku;i++)
						coefu=coefu+coef(idu+i,jj,r)*cup[i+ider*ku];
					coefall=coefall+coefu*cvp[j+jder*kv];
				}
				deriv[r+jderdim+iderndv1dim]=coefall;
			}	
		}
	}

	if(kundu>30) delete[] cup;
	if(kvndv>30) delete[] cvp;
}

//Test if input parameter value is inside parameter range of the surface.
bool MGSBRep::in_range(double u, double v) const{
	return m_uknot.in_range(u) && m_vknot.in_range(v);
}

//Test if input parameter value is inside parameter range of the surface.
bool MGSBRep::in_range(const MGPosition& uv) const{
	return m_uknot.in_range(uv.ref(0)) && m_vknot.in_range(uv.ref(1));
}

//Compare two parameter values. If uv1 is less than uv2, return true.
//Comparison is done after prjected to i-th perimeter of the surface.
bool MGSBRep::less_than(
	int i,	//perimeter number.
	const MGPosition& uv1,
	const MGPosition& uv2) const
{
	assert(i<4);

	switch(i){
	case 0: return uv1.ref(0)<uv2.ref(0);
	case 1: return uv1.ref(1)<uv2.ref(1);
	case 2: return uv1.ref(0)>uv2.ref(0);
	case 3: return uv1.ref(1)>uv2.ref(1);
	}
	return true;
}

// 自身に指定したパラメータ範囲のｌｉｍｉｔをつける。
void MGSBRep::limit(
	const MGBox& uvrange
){
	shrinkToParameters(uvrange, *this);
}

// Return ending parameter value.
MGPosition MGSBRep::param_e() const{
	return MGPosition(m_uknot.param_e(),m_vknot.param_e());
}
double MGSBRep::param_e_u() const{
	return m_uknot.param_e();
}
double MGSBRep::param_e_v() const{
	return m_vknot.param_e();
}

// Compute parameter curve.
//Returned is newed area pointer, and must be freed by delete.
MGCurve* MGSBRep::parameter_curve(
	int is_u				//Indicates x is u-value if is_u is true.
	, double x				//Parameter value.
							//The value is u or v according to is_u.
)const{
	int ku=order_u(); int lud=bdim_u();
	int kv=order_v(); int lvd=bdim_v();
	int is1,is2; surface_bcoef().capacity(is1,is2);
	int ncd=sdim(), len;
	int kx;	int k;
	if(is_u){ kx=1; len=lvd; k=kv; }
	else    { kx=0; len=lud; k=ku; }
	int nderiv=0;

	MGLBRep* lb=new MGLBRep(len,k,ncd);
	MGBPointSeq& rcoef=lb->line_bcoef();
	MGKnotVector& t=lb->knot_vector();
	int n;
	bsepl_(ku,lud,knot_data_u(),kv,lvd,knot_data_v(),coef_data(),
		is1,is2,ncd,kx,x,nderiv,len,&k,&n,&t(0),&rcoef(0,0));
	return lb;
}

// Compute parameter line.
MGLBRep MGSBRep::parameter_line(
	int is_u		//Indicates x is u-value if is_u is true.
	, double x		//Parameter value.
					//The value is u or v according to is_u.
	, int nderiv		//Order of derivative.
)const{
	int ku=order_u(); int lud=bdim_u();
	int kv=order_v(); int lvd=bdim_v();
	int is1,is2; surface_bcoef().capacity(is1,is2);
	int ncd=sdim(), len;
	int kx;	int k;
	if(is_u){ kx=1; len=lvd; k=kv; }
	else    { kx=0; len=lud; k=ku; }
	MGLBRep lb(len,k,ncd);
	MGBPointSeq& rcoef=lb.line_bcoef();
	MGKnotVector& t=lb.knot_vector();
	int n;
	bsepl_(ku,lud,knot_data_u(),kv,lvd,knot_data_v(),coef_data(),
		is1,is2,ncd,kx,x,nderiv,len,&k,&n,&t(0),&rcoef(0,0));
	return lb;
}

// Return starting parameter value.
MGPosition MGSBRep::param_s() const{
	return MGPosition(m_uknot.param_s(),m_vknot.param_s());
}
double MGSBRep::param_s_u() const{
	return m_uknot.param_s();
}
double MGSBRep::param_s_v() const{
	return m_vknot.param_s();
}

//Compute part of the surface limitted by the parameter range bx.
//bx(0) is the parameter (us,ue) and bx(1) is (vs,ve).
//That is u range is from us to ue , and so on.
MGSBRep* MGSBRep::part(const MGBox& uvbx,int multiple) const{
	MGSBRep* sb=new MGSBRep;
	shrinkToParameters(uvbx, *sb,multiple);
	return sb;
}

// Compute perimeter line B-Rep.
MGLBRep MGSBRep::perimeter(int i) const
// i is perimeter number:
// =0: v=min line, =1: u=max line, =2: v=max line, =3: u=min line
{
	assert(i<4);

	int is_u; double x;
	switch(i){
		case 0:  is_u=0; x=param_s_v(); break;
		case 1:  is_u=1; x=param_e_u(); break;
		case 2:  is_u=0; x=param_e_v(); break;
		default: is_u=1; x=param_s_u(); break;
	}
	return parameter_line(is_u, x);
}

//Change the B-Rep by decreasing B-Rep dimension by ndec. This is
//an approximation of the original B-Rep. Return value is error flag.
int MGSBRep::reduce(		//BSUDEC
	int is_u,		//if true, reduce b-rep dimension of u-direction.
	int ndec		//Number of B-rep dimension to decrease 
){
	assert(ndec>=0);
	if(is_u) assert(bdim_u()-ndec>=order_u());
	else 	 assert(bdim_v()-ndec>=order_v());

	if(ndec<=0) return 0;

	int kdec=2; if(is_u) kdec=1;
	int ku=order_u(); int lud=bdim_u();
	int kv=order_v(); int lvd=bdim_v();
	int maxlen=lud; if(maxlen<lvd) maxlen=lvd;
	int is1,is2; surface_bcoef().capacity(is1,is2);
	int lud2, lvd2;
	MGKnotVector tu(ku,lud), tv(kv,lvd);
	MGSPointSeq bcoef(lud,lvd,3);
	int iflag;
	int k=ku; if(k<kv) k=kv;
	double* work2=new double[maxlen+lvd*2+maxlen*(2*k-1)];
	double* work3=work2+maxlen;
	double* work1=work3+lvd*2;
	bsudec_(
		ku,lud,knot_data_u(),kv,lvd,knot_data_v(),coef_data(),
		kdec,ndec,is1,is2,lud,lvd,work1,work2,work3,
		&lud2,&tu(0),&lvd2,&tv(0),&bcoef(0,0,0),&iflag);
	if(iflag==1){
		tu.set_bdim(lud2); tv.set_bdim(lvd2);
		m_uknot=tu; m_vknot=tv;
		bcoef.set_length(lud2,lvd2);
		m_surface_bcoef=bcoef;
		iflag=0;
	}
	delete[] work2;
	invalidateBox();
	return iflag;
}

MGSBRep& MGSBRep::refine(					//BSUNK
	const MGKnotVector& uknot,	// knot of u-direction
	const MGKnotVector& vknot)	// knot of v-direction
//Change an original B-Rep to new one with subdivided knot configuration.
//Knots t must be subdivided knots.
{
	int refine_kind;
	if(uknot!=m_uknot){
		if(vknot!=m_vknot) refine_kind=3;
		else               refine_kind=1;
	}else if(vknot!=m_vknot) refine_kind=2;
	else               return *this;
	// refine_kind=1:u-knot only, =2:v-knot only, =3:both u and v-knot. 

	int ku=order_u(); int lud=bdim_u();
	int kv=order_v(); int lvd=bdim_v();
	int is1,is2; surface_bcoef().capacity(is1,is2);
	int lud2=uknot.bdim(), lvd2=vknot.bdim();
	MGSPointSeq bcoef(lud2,lvd2,3);
	int k=ku; if(k<kv) k=kv;
	double* work1=new double[k*k+lvd2*2];
	double* work2=work1+k*k;
	bsunk_(refine_kind,ku,lud,knot_data_u(),kv,lvd,knot_data_v(),coef_data(),
		is1,is2,lud2,uknot.data(),lvd2,vknot.data(),
		lud2,lvd2,work1,work2,&bcoef(0,0,0));
	delete[] work1;
	m_uknot=uknot; m_vknot=vknot;
	m_surface_bcoef=std::move(bcoef);
	invalidateBox();
	return *this;
}

//Operator overload

//Assignment.
//When the leaf object of this and srf2 are not equal, this assignment
//does nothing.
MGSBRep& MGSBRep::operator=(const MGGel& gel2){
	const MGSBRep* gel2_is_this=dynamic_cast<const MGSBRep*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}
MGSBRep& MGSBRep::operator=(MGGel&& gel2){
	MGSBRep* gel2_is_this=dynamic_cast<MGSBRep*>(&gel2);
	if(gel2_is_this)
		operator=(std::move(*gel2_is_this));
	return *this;
}

// 曲面の平行移動を行いオブジェクトを生成する。
MGSBRep MGSBRep::operator+ ( const MGVector& vec) const{
	MGSBRep brep(*this);
	brep.m_surface_bcoef += vec;
	return brep;
}

// 与ベクトルだけ曲面を平行移動して自身とする。
MGSBRep& MGSBRep::operator+= ( const MGVector& vec){
	m_surface_bcoef += vec;
	m_box+=vec;
	return *this;
}

// 曲面の逆方向に平行移動を行いオブジェクトを生成する。
MGSBRep MGSBRep::operator- ( const MGVector& vec) const{
	MGSBRep brep(*this);
	brep.m_surface_bcoef -= vec;
	return brep;
}

// 与ベクトルだけ曲面をマイナス方向に平行移動して自身とする。
MGSBRep& MGSBRep::operator-= ( const MGVector& vec){
	m_surface_bcoef -= vec;
	m_box-=vec;
	return *this;
}

// 与えられたスケーリングで曲面の変換を行いオブジェクトを生成する。
//Scaling.
MGSBRep MGSBRep::operator* (double scale) const{
	MGSBRep sb(*this);
	sb *=scale;
	return sb;
}

// 与えられたスケーリングで曲面の変換を行いオブジェクトを生成する。
//Scaling.
MGSBRep operator* (double scale, const MGSBRep& sb){
	return sb*scale;
}

// 与えられたスケーリングで曲面の変換を行い自身の曲面とする。
//Scaling.
MGSBRep& MGSBRep::operator*= (double scale){
	m_surface_bcoef *=scale;
	invalidateBox();
	return *this;
}

// 与えられた変換で曲面の変換を行いオブジェクトを生成する。
MGSBRep MGSBRep::operator* ( const MGMatrix& mat ) const{
	MGSBRep brep(*this);
	brep.m_surface_bcoef *= mat;
	return brep;
}

// 与えられた変換で曲面の変換を行い自身の曲面とする。
MGSBRep& MGSBRep::operator*= ( const MGMatrix& mat){
	m_surface_bcoef *= mat;
	invalidateBox();
	return *this;
}

// 与えられた変換で曲面のトランスフォームを行いオブジェクトを生成する。
MGSBRep MGSBRep::operator* ( const MGTransf & tr) const{
	MGSBRep brep(*this);
	brep.m_surface_bcoef *= tr;
	return brep;
}

// 与えられた変換で曲面のトランスフォームを行い自身とする。
MGSBRep& MGSBRep::operator*= ( const MGTransf & tr){
	m_surface_bcoef *= tr;
	invalidateBox();
	return *this;
}

bool MGSBRep::operator==(const MGRSBRep& gel2)const{
	return gel2==(*this);
}

// 与曲面と自身が等しいかの比較判定を行う。
bool MGSBRep::operator== (const MGSBRep& brep2) const{
	return 
		m_uknot == brep2.m_uknot &&
		m_vknot == brep2.m_vknot &&
		m_surface_bcoef == brep2.m_surface_bcoef;
}

std::partial_ordering MGSBRep::operator<=>(const MGSBRep& gel2)const{
	int nu1=bdim_u(), nu2=gel2.bdim_u();
	if(nu1==nu2){
		int nv1=bdim_v(), nv2=gel2.bdim_v();
		if(nv1==nv2){
			int nsd1=sdim(), nsd2=gel2.sdim();
			if(nsd1==nsd2){
				MGVector v1(nsd1), v2(nsd1);
				for(int i=0; i<nsd1; i++){
					v1(i)=m_surface_bcoef.ref(0,0,i);
					v2(i)=gel2.m_surface_bcoef.ref(0,0,i);
				}
				return v1.len()<=>v2.len();
			}else
				return nsd1<=>nsd2;
		}else
			return nv1<=>nv2;
	}else
		return nu1<=>nu2;
}

bool MGSBRep::equal_test(const MGGel& g2)const {
	auto c = typeCompare(g2);
	return c == 0 ? *this == dynamic_cast<const MGSBRep&>(g2) : false;
}
std::partial_ordering MGSBRep::ordering_test(const MGGel& g2)const {
	auto c = typeCompare(g2);
	return c == 0 ? *this <=> dynamic_cast<const MGSBRep&>(g2) : c;
}

//Compute data points and m_uknot and m_vknot from point sequence.
void MGSBRep::compute_knot(
	const MGSPointSeq& points,	//Point seq data
	int orderu,			// Order of u-direction
	int orderv,			// Order of v-direction
	MGNDDArray& utau,			// u data point will be output.
	MGNDDArray& vtau)			// v data point will be output.
{
	points.make_data_point(utau,vtau);
	m_uknot=MGKnotVector(utau, orderu);
	m_vknot=MGKnotVector(vtau, orderv);
}
