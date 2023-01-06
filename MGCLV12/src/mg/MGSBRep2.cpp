/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Box.h"
#include "mg/Vector.h"
#include "mg/Position.h"
#include "mg/KnotArray.h"
#include "mg/KnotVector.h"
#include "mg/BPointSeq.h"
#include "mg/SPointSeq.h"
#include "mg/LBRep.h"
#include "mg/Surface.h"
#include "mg/SBRep.h"
#include "mg/Plane.h"
#include "mg/Tolerance.h"

#include "cskernel/Bsgsmt.h"
#include "cskernel/Bluprt.h"
#include "cskernel/Blumor.h"
#include "cskernel/Blgl2a.h"
#include "cskernel/Bludkt.h"
#include "cskernel/Bluakt.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGSBRep2.cpp
//
// Implements Surface B-Representation class MGSBRep.

//<< Constructor >>
								  
//**** 1. Approximation Constructor ****

///Approximate this, given the new knot configuration tu and tv.
///The new approximated MGSBRep is output into newSB.
///The parameter range of tu and tv must be inside the ones of this.
///newSB is an approximation of this with the parameter range
///from tu.param_s() to tu.param_e() and tv.param_s() to tv.param_e().
///Function's return value is error code:
/// =0 successfully rebuilt, !=0: input tu or tv is invalid.
///When error!=0 is returned, this is copied int newSB.
int MGSBRep::rebuildByNewKnotVector(
	const MGKnotVector& tu,//New knot configuration along u is input.
	const MGKnotVector& tv,//New knot configuration along v is input.
	MGSBRep& newSB//Rebuilt MGSBRep is output, can be this.
)const{
	assert(tu.order()==order_u() && tv.order()==order_v());

	int error=0;					//Error flag.
	int lenu=tu.bdim(), lenv=tv.bdim(), ncd=sdim();
	int ku=tu.order(), kv=tv.order();
	int kmax=ku; if(kmax<kv) kmax=kv;
	int nu, nv, sizeu, sizev;
	m_surface_bcoef.length(nu, nv);
	m_surface_bcoef.capacity(sizeu, sizev);
	int m=(kmax<=5) ? 9:(2*kmax-1);
	int n2=(lenu<lenv) ? lenv:lenu;
	double* work1=new double[n2+n2*m];
	double* work2=work1+n2;

	MGSPointSeq surf1(lenv, lenu, ncd);
	MGSPointSeq surf2(lenu, lenv, ncd);
	MGBPointSeq temp(n2,ncd);

	int irc=sizeu*sizev;
	for(int j=0; j<nv; j++){
		bludkt_(ku,nu,knot_data_u(),coef_data(0,j,0),
				irc,ncd,lenu,tu.data(),n2,
				work1,work2,&temp(0,0),&error);
		if(error!=1){
			error=2;
			break;
		}
		for(int i=0; i<lenu; i++)
			for(int k=0; k<ncd; k++)
				surf1(j,i,k)=temp(i,k);
	}
	if(error==1){
		irc=lenu*lenv;
		for(int i=0; i<lenu; i++){
			bludkt_(kv,nv,knot_data_v(),surf1.data(0,i,0),
					irc,ncd,lenv,tv.data(),n2,
					work1,work2,&temp(0,0),&error);
			if(error!=1) break;
			for(int j=0; j<lenv; j++)
				for(int k=0; k<ncd; k++)
					surf2(i,j,k)=temp(j,k);
		}
		if(error!=1)
			error=12;
	}
							   
	delete[] work1;
	if(error==1){
		error=0;		//Return of BLUDKT:error=1 means normal.
		newSB.m_uknot=tu;
		newSB.m_vknot=tv;
		newSB.m_surface_bcoef=std::move(surf2);
		newSB.invalidateBox();
		newSB.copy_appearance(*this);
	}else{
		//error=2: udirection, =12:v-direction.
		newSB=*this;
	}
	return error;
}

//**** 2.Conversion Constructor.****

//Gets new B-Rep by adding knots to an original B-Rep.
void MGSBRep::addKnots(
	const MGKnotArray& uknots,	//Knots to add for u-direction
	const MGKnotArray& vknots	//Knots to add for v-direction.
){
	invalidateBox();

	int ncd=sdim();
	int ku=order_u(), kv=order_v();
	int nu, nv, sizeu, sizev;
	m_surface_bcoef.length(nu, nv);
	m_surface_bcoef.capacity(sizeu, sizev);

	int kmax=ku; if(kmax<kv) kmax=kv;
	const int nadu=(int)uknots.size();
	const int nadv=(int)vknots.size();
	int kmax2=kmax*kmax;
	int naduPnadv=nadu+nadv;

	double* work=new double[kmax2+naduPnadv];
	double* tadu=work+kmax2;
	double* tadv=tadu+nadu;

	int* mltu=new int[naduPnadv];
	int* mltv=mltu+nadu;

	int mltu_total=0;
	for(int i=0; i<nadu; i++){
		tadu[i]=uknots[i];
		mltu[i]=uknots[i].multiplicity();
		mltu_total+=mltu[i];
	}
	int lenu=nu+mltu_total;

	int mltv_total=0;
	for(int j=0; j<nadv; j++){
		tadv[j]=vknots[j];
		mltv[j]=vknots[j].multiplicity();
		mltv_total+=mltv[j];
	}
	int lenv=nv+mltv_total;
 
	MGSPointSeq surf1(lenv, lenu, ncd);
	MGSPointSeq surf2(lenu, lenv, ncd);
	int lmax=lenu; if(lmax<lenv) lmax=lenv;
	MGBPointSeq temp(lmax,ncd);
	MGKnotVector tu(ku,lenu), tv(kv,lenv);

	int lenunew,lenvnew;
	int irc=sizeu*sizev;
	for(int j=0; j<nv; j++){
		bluakt_(ku, nu, knot_data_u(), coef_data(0,j,0),
			irc, ncd, nadu, tadu, mltu, lmax, work,
			&lenunew, &tu(0), &temp(0,0));
		for(int i=0; i<lenunew; i++)
			for(int k=0; k<ncd; k++)
				surf1(j,i,k)=temp(i,k);
	}
	tu.set_bdim(lenunew);

	irc=lenu*lenv;
	for(int i=0; i<lenunew; i++){
		bluakt_(kv, nv, knot_data_v(), surf1.data(0,i,0),
			irc, ncd, nadv, tadv, mltv, lmax, work,
			&lenvnew, &tv(0), &temp(0,0));
		for(int j=0; j<lenvnew; j++)
			for(int k=0; k<ncd; k++)
				surf2(i,j,k)=temp(j,k);
	}
	tv.set_bdim(lenvnew);
	surf2.set_length(lenunew,lenvnew);
							   
	delete[] mltu; delete[] work;

	m_uknot=std::move(tu); m_vknot=std::move(tv);
	m_surface_bcoef=std::move(surf2);
}

//**** 3.Approximation Constructor.****

void MGSBRep::buildByL2ApproximateWithKTV(			//Least square approximation.
	const MGNDDArray& utau,	//Data point abscissa of u-direction
	const MGNDDArray& vtau,	//Data point abscissa of v-direction
	const MGSPointSeq& points,	//Point seq data
	const MGSPointSeq& weight	//Weight for each point(space dimension 1).
			//weight(i,j,0) is the weight for points(i,j,.).
){
	invalidateBox();
	const int nutau=utau.length(),nvtau=vtau.length();
	const int nu=m_uknot.bdim(),nv=m_vknot.bdim();
	assert(nu==nutau && nv==nvtau);
	int sd=points.sdim();
	m_surface_bcoef.resize(nu, nv, sd);

	const int ip1=points.capacity_u(),iw=weight.capacity_u();
	int nmax=nu; if(nmax<nv) nmax=nv;
	const int ku=m_uknot.order(),kv=m_vknot.order();
	int kmax=ku; if(kmax<kv) kmax=kv;
	int nubynvtau=nu*nvtau;

	std::unique_ptr<double[]>
		workUnique(new double[2*nmax+kmax*nmax+nubynvtau+nubynvtau]);
	double* work=workUnique.get();
	double* q=work+2*nmax;
	double* swork=q+kmax*nmax;//swork[nu][nvtau].
	double* vweight=swork+nubynvtau;//Actually vweight[nu][nvtau].

	//Generate weight of u B-coef.
	for(int i=0; i<nubynvtau; i++)
		vweight[i]=0.0;
	int id;
	for(int j=0; j<nvtau; j++){
		for(int i=0; i<nutau; i++)	{
			id=m_uknot.eval_coef(utau(i),q);
			for(int m=0; m<ku; m++) vweight[j+nvtau*(id+m)]+=q[m]*weight(i,j,0);
		}
	}

	int kw=1;
	for(int m=0; m<sd; m++){
		blgl2a_(nutau,utau.data(),points.data(0,0,m),weight.data(),
			ip1,iw,nvtau,kw,m_uknot.data(),nu,ku,nvtau,work,q,swork);

		blgl2a_(nvtau,vtau.data(),swork,vweight,
			nvtau,nvtau,nu,kw,m_vknot.data(),nv,kv,nu,work,q,
			&m_surface_bcoef(0,0,m));
	}
}

//Shoenberg and Reinch's smoothing function approximation.
void MGSBRep::buildSRSmoothedSB_of_FreeEnd(
	const MGNDDArray& utau,	//Data point abscissa of u-direction
	const MGNDDArray& vtau,	//Data point abscissa of v-direction
	const MGSPointSeq& points,	//Point seq data
	const MGSPointSeq& delp,	//Error estimate for each point
		//(space dimension 1).
		// delp(i,j,0) is for points(i,j,.) at (utau(i),vtau(j)).
	double deviation		//Mean deviation of each point.
//If delp(i,j,0) becomes bigger, deviation at (utau(i),vtau(j)) becomes bigger.:M
){
	invalidateBox();

	int nu = utau.length(), nv = vtau.length();
	int nuP2 = nu + 2, nvP2 = nv + 2;
	int ncd = points.sdim();
	m_surface_bcoef.resize(nuP2, nvP2, ncd);
	m_uknot.size_change(4, nuP2), m_vknot.size_change(4, nvP2);
		
	int npmax=nu;if(npmax<nv) npmax=nv;
	int nup2bynv=nuP2*nv;
	int ip1,ip2; points.capacity(ip1,ip2);
	int ip12=ip1*ip2;

	double* v=new double[npmax*7+npmax*4+npmax*ncd+(npmax+2)*ncd+nup2bynv*ncd+nup2bynv];
	double* a4=v+npmax*7;
	double* a3=a4+npmax*4;
	double* pwork=a3+npmax*ncd;//pwork[ncd][npmax+2]
	double* swork=pwork+(npmax+2)*ncd;//swork[ncd][nv][nup2]
	double* delpu=swork+nup2bynv*ncd;//Actually delpu[nup2][nv].

	int lud;
	for(int j=0; j<nv; j++){	//Process of u-direction.
		bsgsmt_(ncd,nu,utau.data(),points.data(0,j,0),delp.data(0,j,0),
		deviation,ip12,nv,nuP2,v,a4,a3,pwork,&lud,&m_uknot(0),
		swork+j);
	}

	//Generate delp of u B-coef.
	for(int i=0; i<nup2bynv; i++) delpu[i]=0.0;
	int id;
	double coef[4];
	for(int j=0; j<nv; j++){
		for(int i=0; i<nu; i++)	{
			id=m_uknot.eval_coef(utau(i),coef);
			for(int m=0; m<4; m++) delpu[j+nv*(id+m)]+=coef[m]*delp(i,j,0);
		}
	}
							//Process of v-direction.
	int nvi,lvd;
	for(int i=0; i<nuP2; i++){
		nvi=nv*i;
		bsgsmt_(ncd,nv,vtau.data(),swork+nvi,delpu+nvi,
		deviation,nup2bynv,nuP2,nvP2,v,a4,a3,pwork,&lvd,&m_vknot(0),
		&m_surface_bcoef(i,0,0));
	}
	delete[] v;
}

//Gets new B-Rep by connecting this andbrep2 to one.
void MGSBRep::connect(
	int which1,		//which perimeter of this.
	int continuity,	//continuity.
	const MGSBRep& brep2,//B-Rep 2.
	int which2,		//which perimeter of brep2.
	int opposite	// Input if parameter direction of which2
					// is the same as which1 along common edge.
					// If opposite is true, the direction is opposite.
){
	assert(continuity>=0);
	assert(which1>=0 && which1<4 && which2>=0 && which2<4);
	assert(sdim()<=3 && brep2.sdim()<=3);
	invalidateBox();

	MGSBRep sb2(brep2);//To update the direction of brep2.
	if(opposite){
		if(which2==0 || which2==2) sb2.negate(1); //Reverse u-direction.
		else                       sb2.negate(0); //Reverse v-direction.
	}
	if(which1==0 || which1==2){
		if(which2==1 || which2==3){
			sb2.exchange_uv(); which2++; if(which2>3) which2=0;
		}
	}else{
		if(which2==0 || which2==2){
			sb2.exchange_uv(); which2--; if(which2<0) which2=3;
		}
	}

	int lwhich;
	MGLBRep lb, lb1, lb2;
	int sd1=sdim(),sd2=sb2.sdim();
	int nu,nv;
	if(which1==0 || which1==2){
		//Connecting v parameter lines(at u=start or end).
		lb1.knot_vector()=knot_vector_v();
		lb2.knot_vector()=sb2.knot_vector_v();

		int nv1=bdim_v(), nv2=sb2.bdim_v();
		MGBPointSeq bp1(nv1, sd1); MGBPointSeq bp2(nv2, sd2);
		for(int j=0; j<nv1; j++)
			for(int m=0; m<sd1; m++)
				bp1(j,m)=coef(0,j,m);
		for(int j=0; j<nv2; j++)
			for(int m=0; m<sd2; m++)
				bp2(j,m)=sb2.coef(0,j,m);
		lb1.line_bcoef()=std::move(bp1);
		lb2.line_bcoef()=std::move(bp2);

		if(which1==0){
			if(which2==0) lwhich=0;
			else          lwhich=1;
		}else{
			if(which2==0) lwhich=2;
			else          lwhich=3;
		}
		lb=lb1;
		lb.connect(continuity,lwhich,lb2); 

		m_vknot=lb.knot_vector();
		nu=m_uknot.bdim(); nv=m_vknot.bdim();
		int ncd=lb.sdim();

		m_surface_bcoef.reshape(nu,nv,ncd);
		for(int j=0; j<nv; j++)
			for(int m=0; m<ncd; m++)
				m_surface_bcoef(0,j,m)=lb.coef(j,m);

		for(int i=1; i<nu; i++){
			bp1.resize(nv1, sd1); bp2.resize(nv2, sd2);
			for(int j=0; j<nv1; j++) 
				for(int m=0; m<sd1; m++) bp1(j,m)=coef(i,j,m);
			for(int j=0; j<nv2; j++) 
				for(int m=0; m<sd2; m++) bp2(j,m)=sb2.coef(i,j,m);
			lb1.line_bcoef()=std::move(bp1);
			lb2.line_bcoef()=std::move(bp2);

			lb=lb1;
			lb.connect(continuity,lwhich,lb2); 
			for(int j=0; j<nv; j++)
				for(int m=0; m<ncd; m++)
					m_surface_bcoef(i,j,m)=lb.coef(j,m);
		}
	}else{
		//Connecting u parameter lines(at v=start or end).
		lb1.knot_vector()=knot_vector_u();
		lb2.knot_vector()=sb2.knot_vector_u();
		int i;
		int nu1=bdim_u(), nu2=sb2.bdim_u();
		MGBPointSeq bp1(nu1, sd1); MGBPointSeq bp2(nu2, sd2);
		for(i=0; i<nu1; i++)
			for(int m=0; m<sd1; m++)
				bp1(i,m)=coef(i,0,m);
		for(i=0; i<nu2; i++)
			for(int m=0; m<sd2; m++)
				bp2(i,m)=sb2.coef(i,0,m);
		lb1.line_bcoef()=std::move(bp1);
		lb2.line_bcoef()=std::move(bp2);

		if(which1==3){
			if(which2==3) lwhich=0;
			else          lwhich=1;
		}else{
			if(which2==3) lwhich=2;
			else          lwhich=3;
		}

		lb=lb1;
		lb.connect(continuity, lwhich, lb2);
		nu=lb.bdim(); nv=m_vknot.bdim();
		int ncd=lb.sdim();
		m_vknot=knot_vector_v();
		m_uknot=lb.knot_vector();

		m_surface_bcoef.reshape(nu,nv,ncd);
		for(i=0; i<nu; i++)
			for(int m=0; m<ncd; m++)
				m_surface_bcoef(i,0,m)=lb.coef(i,m);

		for(int j=1; j<nv; j++){
			bp1.resize(nu1, sd1); bp2.resize(nu2, sd2);
			for(i=0; i<nu1; i++) 
				for(int m=0; m<sd1; m++)
					bp1(i,m)=coef(i,j,m);
			for(i=0; i<nu2; i++) 
				for(int m=0; m<sd2; m++)
					bp2(i,m)=sb2.coef(i,j,m);
			lb1.line_bcoef()=std::move(bp1);
			lb2.line_bcoef()=std::move(bp2);

			lb=lb1;
			lb.connect(continuity, lwhich, lb2);
			for(i=0; i<nu; i++)
				for(int m=0; m<ncd; m++)
					m_surface_bcoef(i,j,m)=lb.coef(i,m);
		}
	}
}

///Obtain the partial Surface B-Rep restricted by sub interval of u and v parameter
///range. New one is exactly the same as the original except that it is partial.
///If multiple==true(!=0), knot_u(i)=t1 and knot_u(n+i)=t2 for i=0,..., k-1
///will be guaranteed. Here, n=bdim_u(), k=order_u(),
///t1=uvrange(0).low_point(), and t2=uvrange(0).high_point().
///About knot_v(j), the same.
///Both u-range and v-range must be inside the range of this.
void MGSBRep::shrinkToParameters(
	const MGBox& uvrange,		//u and v parameter range.
	MGSBRep& newBrep,	///<Subdivided surface is output, which can be this.
	int multiple	//Indicates if start and end knot multiplicities
		//are necessary. =0:unnecessary, !=0:necessary.
		//If multiple==true(!=0), knot_u(i)=t1 and knot_u(n+i)=t2 for i=0,..., k-1
		//are guaranteed. Here, n=bdim_u(), k=order_u(),
		//t1=uvrange(0).low_point(), and t2=uvrange(0).high_point().
		//About knot_v(j), the same.
		// Both u-range and v-range must be inside the range of this.
)const{
	newBrep.invalidateBox();
	MGBox b2=param_range();
	MGBox b3=b2&uvrange;
	if(b3 == b2){
		newBrep.m_uknot=m_uknot;
		newBrep.m_vknot=m_vknot;
		newBrep.m_surface_bcoef=m_surface_bcoef;
	}else{
		int ncd=sdim();
		int ku=order_u(), kv=order_v();
		int kmax=ku; if(kmax<kv) kmax=kv;
		int nu, nv, sizeu, sizev;
		m_surface_bcoef.length(nu, nv);
		m_surface_bcoef.capacity(sizeu, sizev);
		int nmax=nu; if(nmax<nv) nmax=nv;

		double* work=new double[kmax*kmax];

		MGSPointSeq surf1(nv, nu, ncd);
		MGSPointSeq surf2(nu, nv, ncd);
		MGBPointSeq temp(nmax, ncd);
		MGKnotVector tu(ku, nu), tv(kv, nv);

		int nunew, nvnew;
		double t1, t2;
		int irc=sizeu*sizev;
		t1=b3(0).low_point(); t2=b3(0).high_point();
		for(int j=0; j<nv; j++){
			bluprt_(ku, nu, knot_data_u(), coef_data(0, j, 0),
					irc, ncd, t1, t2, nmax, work, &nunew, &tu(0), &temp(0, 0), multiple);
			for(int i=0; i<nunew; i++)
				for(int k=0; k<ncd; k++)
					surf1(j, i, k)=temp(i, k);
		}
		tu.reshape(nunew+ku);

		irc=nu*nv;
		t1=b3(1).low_point(); t2=b3(1).high_point();
		for(int i=0; i<nunew; i++){
			bluprt_(kv, nv, knot_data_v(), surf1.data(0, i, 0),
					irc, ncd, t1, t2, nmax, work, &nvnew, &tv(0), &temp(0, 0), multiple);
			for(int j=0; j<nvnew; j++)
				for(int k=0; k<ncd; k++)
					surf2(i, j, k)=temp(j, k);
		}
		tv.reshape(nvnew+kv);
		surf2.reshape(nunew, nvnew);

		delete[] work;

		newBrep.m_uknot=std::move(tu);
		newBrep.m_vknot=std::move(tv);
		newBrep.m_surface_bcoef=std::move(surf2);
	}
	newBrep.copy_appearance(*this);
}

// Construct a Surface B-Rep by changing space dimension and
//ordering of coordinate.
MGSBRep::MGSBRep(
	int dim,				// New space dimension.
	const MGSBRep& sbrep,	// Original Surface B-rep.
	int start1, 			// Destination order of new Surface.
	int start2) 			// Source order of original Surface.
:MGSurface(sbrep)
,m_uknot(sbrep.knot_vector_u()), m_vknot(sbrep.knot_vector_v())
,m_surface_bcoef(dim,sbrep.surface_bcoef(),start1,start2){
	invalidateBox();
}

// Construct a Surface B-Rep (order = 2 ) from Plane and Parameter ranges.
void MGSBRep::buildPlaneSBRep(
		const MGPlane& plane,	// Original Plane
		const MGBox& prange)	// parameter range of new Surface.
{
	invalidateBox();
	m_surface_bcoef.resize(2, 2, plane.sdim());
	m_uknot.size_change(2, 2), m_vknot.size_change(2, 2);
	double u[2] = {prange(0).low_point(), prange(0).high_point()};
	double v[2] = {prange(1).low_point(), prange(1).high_point()};		
	
	// m_surface_bcoef‚ÉControl Polygon‚Ìƒf[ƒ^‚ð‘ã“ü
	for(int i = 0; i < 2; i++){
		m_uknot(2*i)=m_uknot(2*i+1) = u[i];
		m_vknot(2*i)=m_vknot(2*i+1) = v[i];
		m_surface_bcoef.store_at(0, i, plane.eval(u[0], v[i]));
		m_surface_bcoef.store_at(1, i, plane.eval(u[1], v[i]));
	}
}

// Exchange parameter u and v.
MGSurface& MGSBRep::exchange_uv(){
	int nu=bdim_u(), nv=bdim_v(), dim=sdim();
	// Exchange surface coefficients.
	MGSPointSeq temp(nv,nu,dim);
	for(int k=0; k<dim; k++)
		for(int j=0; j<nv; j++)
			for(int i=0; i<nu; i++) temp(j,i,k)=coef(i,j,k);
	MGKnotVector tu(std::move(m_vknot)), tv(std::move(m_uknot));
	buildSBRepFromMemberData(std::move(temp), std::move(tu), std::move(tv));
	return *this;
}

//Modify the original Surface by extrapolating the specified perimeter.
//The extrapolation is C2 continuous if the order >=4.
//The extrapolation is done so that extrapolating length is "length"
//at the position of the parameter value "param" of the perimeter.
MGSBRep& MGSBRep::extend(
	int perimeter,	//perimeter number of the Surface.
					// =0:v=min, =1:u=max, =2:v=max, =3:u=min.
	double param,	// parameter value of above perimeter.
	double length,	//chord length to extend at the parameter param of the perimeter.
	double dk      //Coefficient of how curvature should vary at
//    extrapolation start point. When dk=0, curvature keeps same, i.e.
//    dK/dS=0. When dk=1, curvature becomes zero at length extrapolated point,
//    i.e. dK/dS=-K/length at extrapolation start point.
//    (S=parameter of arc length, K=Curvature at start point)
//    That is, when dk reaches to 1 from 0, curve changes to flat.
){
	assert(sdim()<=3);
	assert(perimeter>=0 && perimeter<4);
	invalidateBox();
	
	const int ncd=surface_bcoef().sdim();
	int at_start=1;//starting perimeter
	int nu, nv;
	int order; int n,m; MGKnotVector* t;
	if(perimeter==1 || perimeter==3){	// Extrapolate to u-direction
		order=order_u();
		n=bdim_u();
		t=&(knot_vector_u());
		if(perimeter==1) at_start=0;//ending perimeter
		nu=n+1; m=nv=bdim_v();
	}else{
		// Extrapolate to v-direction
		order=order_v();
		n=bdim_v();
		t=&(knot_vector_v());
		if(perimeter==2) at_start=0;//ending perimeter
		nv=n+1; m=nu=bdim_u();
	}
	//(nu,nv) are new surface B-Rep dimensions of u and v.
	//(order,n,t) is line B-rep to extrapolate.
	//m is the number of line B-reps to extrapolate.

	double* dcoef=new double[order];
	int id;
	if(at_start)
		id=t->eval_coef(t->param_s(),dcoef,1);	//To get 1st derivative.
	else
		id=t->eval_coef(t->param_e(),dcoef,1);	//To get 1st derivative.

	MGSPointSeq surf(nu,nv,ncd);
	MGLBRep lbtemp(n,order,ncd);
	MGKnotVector& t1=lbtemp.knot_vector();
	MGBPointSeq& coeftemp=lbtemp.line_bcoef();

	MGPosition uv=perimeter_uv(perimeter,param);//Surface parameter value of param.
	int ndu=0,ndv=0;
	if(perimeter==0 || perimeter==2) ndv=1;
	else                             ndu=1;
	double slen=length/(eval(uv,ndu,ndv)).len();

	int nnew,cid; double firstd_len,data,dlen;
	for(int i=0; i<m; i++){
		if(perimeter==0 || perimeter==2){
			for(int j=0; j<n; j++)
				for(int k=0; k<ncd; k++) coeftemp(j,k)=coef(i,j,k);
		}else{
			for(int j=0; j<n; j++)
				for(int k=0; k<ncd; k++) coeftemp(j,k)=coef(j,i,k);
		}
		coeftemp.set_length(n);

	//Compute first derivative length at the end of the extrapolating line.
		firstd_len=0.;
		for(int k=0; k<ncd; k++){
			data=0.;
			for(int j=0; j<order; j++){
				cid=id+j; data=data+coeftemp(cid,k)*dcoef[j];
			}
			firstd_len+=data*data;
		}
		firstd_len=sqrt(firstd_len);

		t1=*t; dlen=firstd_len*slen;
		lbtemp.extend(at_start,dlen,dk);
		nnew=lbtemp.bdim();
		if(perimeter==0 || perimeter==2){
			for(int j=0; j<nnew; j++)
				for(int k=0; k<ncd; k++) surf(i,j,k)=coeftemp(j,k);
		}else{
			for(int j=0; j<nnew; j++)
				for(int k=0; k<ncd; k++) surf(j,i,k)=coeftemp(j,k);
		}
	}

	*t=std::move(t1);
	surf.set_length(nu,nv);
	surface_bcoef()=std::move(surf);

	delete[] dcoef;
	return *this;
}

MGSBRep& MGSBRep::move(		//BLUMOR
	int move_kind_u,		//Indicates how to move Surface for u direction.
	int move_kind_v,		//Indicates how to move Surface for v direction.
	const MGPosition& move_point_param, //indicates object point to move
							// by the (u,v) parameter value.
	const MGPosition& to_point,	//destination point of the abve source point.
	const MGPosition fix_point[2])
//Modify the original Surface by moving move_point to to_point. fix_point can be
//applied according to move_kind.
{
	int i,id;

	double u=move_point_param(0), v=move_point_param(1);
	if(u<param_s_u() || u>param_e_u() || v<param_s_v() || v>param_e_v())
		return *this;

	invalidateBox();
	MGVector delta=to_point-eval(move_point_param);
	//delta is a vector from original to destination point.

	int ku=order_u(), kv=order_v();
	int nu=bdim_u(), nv=bdim_v();
	double* ucoef=new double[ku+kv+nu+nv];
	double* vcoef=ucoef+ku;
	double* uratio=vcoef+kv;
	double* vratio=uratio+nu;

	int ui,uj,vi,vj,umovek, vmovek;
	int uki=m_uknot.eval_coef(u,ucoef);
	int vki=m_vknot.eval_coef(v,vcoef);
	//ucoef[i] and vcoef[j] are coefficients that are multiplied to
	//B-Coef's. uki and vki are the ids of first B-Coef's. 0<=i<ku, 0<=j<kv. 
	double us,ue,vs,ve;			//For fixed point parameter values.
	int uis,uie,vis,vie;	//For start and end id's of B-Coef's that
								//should be modified.
	double usum,vsum;

	//Compute u-direction ratio in uratio.
	switch(move_kind_u){
	case (1): 
		ui=ku; uj=nu; umovek=1; break;
	case (2):
		us=fix_point[0](0); us=m_uknot.range(us);
		ui=m_uknot.locate(us)+1;
		if(us<=u){umovek=2; uj=nu;}
		else     {umovek=3; uj=ui; ui=1;}
		break;
	case (3):
		us=fix_point[0](0); us=m_uknot.range(us);
		ue=fix_point[1](0); ue=m_uknot.range(ue);
		if(us>ue){ double uu=us; us=ue; ue=uu;}
		ui=m_uknot.locate(m_uknot.range(us))+1;
		uj=m_uknot.locate(m_uknot.range(ue))+1;
		if(ue<=u)    {umovek=2; ui=uj; uj=nu;}
		else if(u<us){umovek=3; uj=ui; ui=1;}
		else         {umovek=1;}
		break;
	default:
		umovek=1; uj=ui=(uki+ku); break;
	}
	blumor_(umovek,&ui,&uj,u,ku,nu,m_uknot.data(),&uis,&uie,uratio);
	uis-=1; uie-=1; 
	id=0; usum=0;
	for(i=uki; i<uki+ku; i++){
		if(i>=uis && i<=uie) usum=usum+ucoef[id]*uratio[i];
		id+=1;
	}
	if(MGMZero(usum)) goto end_p;

	//Compute v-direction ratio in vratio.
	switch(move_kind_v){
	case (1): 
		vi=kv; vj=nv; vmovek=1; break;
	case (2):
		vs=fix_point[0](1); vs=m_vknot.range(vs);
		vi=m_vknot.locate(vs)+1;
		if(vs<=v){vmovek=2; vj=nv;}
		else     {vmovek=3; vj=vi; vi=1;}
		break;
	case (3):
		vs=fix_point[0](1); vs=m_vknot.range(vs);
		ve=fix_point[1](1); ve=m_vknot.range(ve);
		if(vs>ve){ double vv=vs; vs=ve; ve=vv;}
		vi=m_vknot.locate(m_vknot.range(vs))+1;
		vj=m_vknot.locate(m_vknot.range(ve))+1;
		if(ve<=v)    {vmovek=2; vi=vj; vj=nv;}
		else if(v<vs){vmovek=3; vj=vi; vi=1;}
		else         {vmovek=1;}
		break;
	default:
		vmovek=1; vj=vi=(vki+kv); break;
	}
	blumor_(vmovek,&vi,&vj,v,kv,nv,m_vknot.data(),&vis,&vie,vratio);
	vis-=1; vie-=1;
	id=0; vsum=0.;
	for(int j=vki; j<vki+kv; j++){
		if(j>=vis && j<=vie) vsum=vsum+vcoef[id]*vratio[j];
		id+=1;
	}
	if(MGMZero(vsum)) goto end_p;

	//New B-Coefficients.
	double ratio1,ratio2;

	for(i=uis; i<=uie; i++){
		ratio1=uratio[i]/usum;
		for(int j=vis; j<=vie; j++){
			ratio2=ratio1*vratio[j]/vsum;
			for(int m=0; m<sdim(); m++)
				m_surface_bcoef(i,j,m)+=delta.ref(m)*ratio2;
		}
	}

end_p:
	delete[] ucoef;
	return *this;
}

void MGSBRep::negate(	//Change direction of the surface.
	int is_u)				// Negate along u-direction if is_u is ture,
							// else along v-direction.
{
	if(is_u){		//u-direction.
		m_uknot.reverse(); m_surface_bcoef.reverse(is_u);
	}else{
		m_vknot.reverse(); m_surface_bcoef.reverse(is_u);
	}
}

//Obtain parameter value if this surface is negated by "negate()".
// Negate along u-direction if is_u is ture,
// else along v-direction.
MGPosition MGSBRep::negate_param(const MGPosition& uv, int is_u)const{
	double u=uv(0), v=uv(1);
	if(is_u){		//u-direction.
		u=m_uknot.reverse_param(u);
	}else{
		v=m_vknot.reverse_param(v);
	}
	return MGPosition(u,v);
}

#define INCNUM 15
///Rebuild this MGSBRep. Rebuild means:
/// 1) Change the parameterization.
/// 2) Remove the redundant surface B-coefficients within the tolerance, which is
///    performed by remove_knots.
std::unique_ptr<MGSBRep> MGSBRep::rebuild(
	int how_rebuild,
		//intdicates how rebuild be done.
		// =0: no approximation(only parameter change)
		// !=0: approximated by non-rational spline(MGSBRep) with new knot configuration.
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
	std::unique_ptr<MGSBRep> srfNew;
	if(how_rebuild){
		int order2[2]={order_u(), order_v()};
		if(order){
			if(order[0])
				order2[0]=order[0];
			if(order[1])
				order2[1]=order[1];
		}
		srfNew=approximate_as_SBRep(parameter_normalization,tol,order2);
	}else{
		MGKnotVector uknots, vknots;
		get_new_surface_knots(parameter_normalization,uknots,vknots);
		srfNew=std::unique_ptr<MGSBRep>(new MGSBRep);
		srfNew->buildSBRepFromMemberData(
			surface_bcoef(), std::move(uknots), std::move(vknots));
	}
	srfNew->copy_appearance(*this);
	return srfNew;
}

//Shrink this surface to the part limitted by the parameter range of uvbx.
//New parameter range uvbx2 is so determined that uvbx2 is the smallest
//box tha includes uvbx, and all of the u or v values of uvbx2 is one of 
//the values of u or v knots of the surface knotvector.
//uvbx(0) is the parameter (us,ue) and uvbx(1) is (vs,ve).
//That is u range is from us to ue , and so on.
void MGSBRep::shrink_to_knot(
	const MGBox& uvbx,
	int multiple	//Indicates if start and end knot multiplicities
					//are necessary. =0:unnecessary, !=0:necessary.
){
	MGBox uvb2=box_param();
	uvb2&=uvbx;
	invalidateBox();

	MGInterval& uspan=uvb2[0];
	double u0=uspan.low_point(), u1=uspan.high_point();
	int idu0=m_uknot.locate(u0), idu1=m_uknot.locate(u1)+1;
	uspan.set_low(m_uknot[idu0]);
	uspan.set_high(m_uknot[idu1]);

	MGInterval& vspan=uvb2[1];
	double v0=vspan.low_point(), v1=vspan.high_point();
	int idv0=m_vknot.locate(v0), idv1=m_vknot.locate(v1)+1;
	vspan.set_low(m_vknot[idv0]);
	vspan.set_high(m_vknot[idv1]);
	shrinkToParameters(uvb2, *this, multiple);
}
