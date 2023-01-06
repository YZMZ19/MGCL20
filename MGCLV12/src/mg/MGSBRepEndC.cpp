/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/NDDArray.h"
#include "mg/BPointSeq.h"
#include "mg/Position.h"
#include "mg/Curve.h"
#include "mg/SBRepEndC.h"
#include "mg/Coons.h"
#include "mg/BSumSurf.h"

#include "cskernel/Bvltn2.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGSBRepEndC.cc
//
// Implements MGSBRepEndC class.

//Member Data
//	MGENDCOND m_cond;	//Type of end condition.
//	MGVector  m_1deriv;	//1st derivative stored 
//						//when m_cond=MGENDC_1D or MGENDC_12D
//	MGVector  m_2deriv;	//2nd derivative stored

// MGSBRepEndC_twist computes twist data of a corner, given
// df/dv along u(vderiv) and df/du along v(uderiv).
//This is a private function for complete_corner_deriv.
MGVector MGSBRepEndC_twist(
	const MGNDDArray& utau, //u-direction data point
	const MGBPointSeq& vderiv,
	const MGNDDArray& vtau,	//v-direction data point
	const MGBPointSeq& uderiv,
	int ncd,				//space dimension
	int u_dorder, int v_dorder,
	int endu, int endv);

//Constructor

//Default Constructor.
MGSBRepEndC::MGSBRepEndC():m_sdim(0), m_nu(0), m_nv(0){
	for(int i=0; i<4; i++){
		m_cond[i]= MGENDCOND::MGENDC_NO;
	}
}

//Construct from a Coon's patch.
MGSBRepEndC::MGSBRepEndC(
	const MGNDDArray& utau,
	const MGNDDArray& vtau,
	const MGCoons& coons)
:m_sdim(coons.sdim()), m_nu(utau.length()), m_nv(vtau.length()){
	const MGNDDArray* tau;
 	for(int i=0; i<4; i++){
		m_cond[i]= MGENDCOND::MGENDC_1D;
		if(i%2) tau=&vtau; else tau=&utau;
		coons.derivative(i).eval_line(*tau,m_1deriv[i]);
		m_11d[i]=coons.d2fdvu(i);
		m_2deriv[i].set_null();
	}
}

//Construct from a Boolean sum surface's patch.
MGSBRepEndC::MGSBRepEndC(
	const MGNDDArray& utau,
	const MGNDDArray& vtau,
	const MGBSumSurf& bssurf)
:m_sdim(bssurf.sdim()), m_nu(utau.length()), m_nv(vtau.length()){
	const MGNDDArray* tau;
	int m=utau.length(), n=vtau.length();
	double u0=utau[0], u1=utau[m-1];
	double v0=vtau[0], v1=vtau[n-1];
	MGPosition uvcorner[4]=
		{MGPosition(u0,v0), MGPosition(u1,v0),MGPosition(u1,v1),MGPosition(u0,v1)};
	MGPosition uv;
	int id,id2;
 	for(int i=0; i<4; i++){
		uv=uvcorner[i];
		m_11d[i]=bssurf.eval(uv,1,1);
		m_cond[i]= MGENDCOND::MGENDC_1D;
		if(id=(i%2)) tau=&vtau; else tau=&utau;
		int taulen=tau->length();
		MGBPointSeq& derivi=m_1deriv[i];
		derivi.resize(taulen,m_sdim);
		id2=(id+1)%2;
		for(int j=0; j<taulen; j++){
			uv(id)=(*tau)[j];
			derivi.store_at(j,bssurf.eval(uv,id,id2));
		}
		m_2deriv[i].set_null();
	}
}

//Copy Constructor.
MGSBRepEndC::MGSBRepEndC(const MGSBRepEndC& ec)
:m_sdim(ec.m_sdim), m_nu(ec.m_nu), m_nv(ec.m_nv){
 	for(int i=0; i<4; i++){
		m_11d[i]=ec.m_11d[i];
		m_12d[i]=ec.m_12d[i];
		m_21d[i]=ec.m_21d[i];
		m_22d[i]=ec.m_22d[i];
		m_cond[i]=ec.m_cond[i];
		if(m_cond[i]== MGENDCOND::MGENDC_1D || m_cond[i]== MGENDCOND::MGENDC_12D)
			m_1deriv[i]=ec.m_1deriv[i];
		else
			m_1deriv[i].set_null();
		if(m_cond[i]== MGENDCOND::MGENDC_2D || m_cond[i]== MGENDCOND::MGENDC_12D)
			m_2deriv[i]=ec.m_2deriv[i];
		else
			m_2deriv[i].set_null();
	}
}

//Move Constructor.
MGSBRepEndC::MGSBRepEndC(MGSBRepEndC&& ec)
:m_sdim(ec.m_sdim), m_nu(ec.m_nu), m_nv(ec.m_nv){
 	for(int i=0; i<4; i++){
		m_11d[i]=std::move(ec.m_11d[i]);
		m_12d[i]=std::move(ec.m_12d[i]);
		m_21d[i]=std::move(ec.m_21d[i]);
		m_22d[i]=std::move(ec.m_22d[i]);
		m_cond[i]=ec.m_cond[i];
		m_1deriv[i]=std::move(ec.m_1deriv[i]);
		m_2deriv[i]=std::move(ec.m_2deriv[i]);
	}
}

////////// Copy assignment //////////
MGSBRepEndC& MGSBRepEndC::operator=(const MGSBRepEndC& ec){
	m_sdim=ec.m_sdim;
	m_nu=ec.m_nu;
	m_nv=ec.m_nv;
 	for(int i=0; i<4; i++){
		m_11d[i]=ec.m_11d[i];
		m_12d[i]=ec.m_12d[i];
		m_21d[i]=ec.m_21d[i];
		m_22d[i]=ec.m_22d[i];
		m_cond[i]=ec.m_cond[i];
		if(m_cond[i]== MGENDCOND::MGENDC_1D || m_cond[i]== MGENDCOND::MGENDC_12D)
			m_1deriv[i]=ec.m_1deriv[i];
		else
			m_1deriv[i].set_null();
		if(m_cond[i]== MGENDCOND::MGENDC_2D || m_cond[i]== MGENDCOND::MGENDC_12D)
			m_2deriv[i]=ec.m_2deriv[i];
		else
			m_2deriv[i].set_null();
	}
	return *this;
}

////////// Move assignment//////////
MGSBRepEndC& MGSBRepEndC::operator=(MGSBRepEndC&& ec){
	m_sdim=ec.m_sdim; m_nu=ec.m_nu; m_nv=ec.m_nv;
 	for(int i=0; i<4; i++){
		m_11d[i]=std::move(ec.m_11d[i]);
		m_12d[i]=std::move(ec.m_12d[i]);
		m_21d[i]=std::move(ec.m_21d[i]);
		m_22d[i]=std::move(ec.m_22d[i]);
		m_cond[i]=std::move(ec.m_cond[i]);
		m_1deriv[i]=std::move(ec.m_1deriv[i]);
		m_2deriv[i]=std::move(ec.m_2deriv[i]);
	}
	return *this;
}

//Member Function

void MGSBRepEndC::complete_corner_deriv
(const MGNDDArray& utau,	//u-direction data point
 const MGNDDArray& vtau)	//v-direction data point
// Function return inf:
{
	int i,ip1, idu,idv;
	int ncd=m_sdim;
	const int ise1[4]={2,2,1,1}, ise2[4]={1,2,2,1};
	
	for(i=0; i<4; i++){
		ip1=i+1; if(ip1==4) ip1=0;
		if(i==0 || i==2) {idu=i; idv=ip1;}
		else             {idu=ip1; idv=i;}

		//Construct twist vectors( d2f/(du*dv) ).
		if((m_cond[idu]== MGENDCOND::MGENDC_1D || m_cond[idu]== MGENDCOND::MGENDC_12D) &&
		   (m_cond[idv]== MGENDCOND::MGENDC_1D || m_cond[idv]== MGENDCOND::MGENDC_12D)
		   && m_11d[ip1].sdim()==0){
			m_11d[ip1]=
				MGSBRepEndC_twist(utau,m_1deriv[idu],vtau,m_1deriv[idv],
				ncd,1,1,ise1[i],ise2[i]);
		}

		//Construct d3f/(du**2*dv).
		if((m_cond[idu]== MGENDCOND::MGENDC_1D || m_cond[idu]== MGENDCOND::MGENDC_12D) &&
		   (m_cond[idv]== MGENDCOND::MGENDC_2D || m_cond[idv]== MGENDCOND::MGENDC_12D)
		   && m_21d[ip1].sdim()==0){
			m_21d[ip1]=
				MGSBRepEndC_twist(utau,m_1deriv[idu],vtau,m_2deriv[idv],
				ncd,2,1,ise1[i],ise2[i]);
		}

		//Construct d3f/(du*dv**2).
		if((m_cond[idu]== MGENDCOND::MGENDC_2D || m_cond[idu]== MGENDCOND::MGENDC_12D) &&
		   (m_cond[idv]== MGENDCOND::MGENDC_1D || m_cond[idv]== MGENDCOND::MGENDC_12D)
		   && m_12d[ip1].sdim()==0){
			m_12d[ip1]=
				MGSBRepEndC_twist(utau,m_2deriv[idu],vtau,m_1deriv[idv],
				ncd,1,2,ise1[i],ise2[i]);
		}

		//Construct d4f/(du**2*dv**2).
		if((m_cond[idu]== MGENDCOND::MGENDC_2D || m_cond[idu]== MGENDCOND::MGENDC_12D) &&
		   (m_cond[idv]== MGENDCOND::MGENDC_2D || m_cond[idv]== MGENDCOND::MGENDC_12D)
		   && m_22d[ip1].sdim()==0){
			m_22d[ip1]=
				MGSBRepEndC_twist(utau,m_2deriv[idu],vtau,m_2deriv[idv],
				ncd,2,2,ise1[i],ise2[i]);
		}
	}
}

//Initialize the instance. Will be set to the same as constructed by the void 
//constructor.
void MGSBRepEndC::initialize(){
	m_sdim=m_nu=m_nv=0;
 	for(int i=0; i<4; i++){
		m_1deriv[i].set_null();
		m_2deriv[i].set_null();
		m_cond[i]= MGENDCOND::MGENDC_NO;
	}
}

void MGSBRepEndC::setLength(int i, const MGBPointSeq & derivs){
	assert(i<4);
	int dim=derivs.sdim();
	int n=derivs.length();
	if(!m_sdim || m_sdim<dim) m_sdim=dim;
	if(i==0 || i==2){
		if(m_nu) assert(m_nu==n);
		else m_nu=n;
	} else {
		if(m_nv) assert(m_nv==n);
		else m_nv=n;
	}

}

//Set 1st deriv and change condition type to MGENDC_1D or MGENDC_12D.
void MGSBRepEndC::set_1st(int i, const MGBPointSeq& first_deriv)
{
	setLength(i, first_deriv);
	m_1deriv[i]=first_deriv;
	if(m_cond[i]== MGENDCOND::MGENDC_2D)
		m_cond[i]= MGENDCOND::MGENDC_12D;
	else
		m_cond[i]= MGENDCOND::MGENDC_1D;
}
void MGSBRepEndC::set_1st(int i, MGBPointSeq&& first_deriv)
{
	setLength(i, first_deriv);
	m_1deriv[i]=std::move(first_deriv);
	if(m_cond[i]== MGENDCOND::MGENDC_2D)
		m_cond[i]= MGENDCOND::MGENDC_12D;
	else
		m_cond[i]= MGENDCOND::MGENDC_1D;
}

//Set 1st deriv and change condition type to MGENDC_1D or MGENDC_12D.
void MGSBRepEndC::set_1st(int i, std::unique_ptr<MGBPointSeq>&& first_derivp)
{
	setLength(i, *first_derivp);
	m_1deriv[i]=std::move(*first_derivp);
	if(m_cond[i]== MGENDCOND::MGENDC_2D)
		m_cond[i]= MGENDCOND::MGENDC_12D;
	else
		m_cond[i]= MGENDCOND::MGENDC_1D;
}

//Set 2nd deriv and change condition type to MGENDC_2D or MGENDC_12D.
void MGSBRepEndC::set_2nd(int i, const MGBPointSeq& second_deriv)
{
	setLength(i, second_deriv);
	m_2deriv[i]=second_deriv;
	if(m_cond[i]== MGENDCOND::MGENDC_1D)
		m_cond[i]= MGENDCOND::MGENDC_12D;
	else
		m_cond[i]= MGENDCOND::MGENDC_2D;
}

void MGSBRepEndC::set_2nd(int i, MGBPointSeq && second_deriv){
	setLength(i, second_deriv);
	m_2deriv[i]=std::move(second_deriv);
	if(m_cond[i]== MGENDCOND::MGENDC_1D)
		m_cond[i]= MGENDCOND::MGENDC_12D;
	else
		m_cond[i]= MGENDCOND::MGENDC_2D;
}
void MGSBRepEndC::set_2nd(int i, std::unique_ptr<MGBPointSeq>&& second_derivp)
{
	setLength(i, *second_derivp);
	m_2deriv[i]=std::move(*second_derivp);
	if(m_cond[i]== MGENDCOND::MGENDC_1D)
		m_cond[i]= MGENDCOND::MGENDC_12D;
	else
		m_cond[i]= MGENDCOND::MGENDC_2D;
}

//Set m_11d[i] inf.
void MGSBRepEndC::set_11d(int i, const MGVector& deriv)
{ 	
	assert(i<4);
	int dim=deriv.sdim();
	if(!m_sdim || m_sdim<dim) m_sdim=dim;
	m_11d[i]=deriv; }

//Set m_12d[i] inf.
void MGSBRepEndC::set_12d(int i, const MGVector& deriv)
{ 
	assert(i<4); 
	int dim=deriv.sdim();
	if(!m_sdim || m_sdim<dim) m_sdim=dim;
    m_12d[i]=deriv; }

//Set m_21d[i] inf.
void MGSBRepEndC::set_21d(int i, const MGVector& deriv)
{ 
	assert(i<4); 
	if(m_sdim) assert(m_sdim==deriv.sdim());
	else m_sdim=deriv.sdim();
	m_21d[i]=deriv; }

//Set m_22d[i] inf.
void MGSBRepEndC::set_22d(int i, const MGVector& deriv)
{ 
	assert(i<4); 
	int dim=deriv.sdim();
	if(!m_sdim || m_sdim<dim) m_sdim=dim;
	m_22d[i]=deriv; }

MGVector MGSBRepEndC_twist(
	const MGNDDArray& utau, //u-direction data point
	const MGBPointSeq& vderiv,
	const MGNDDArray& vtau,	//v-direction data point
	const MGBPointSeq& uderiv,
	int ncd,
	int u_dorder, int v_dorder,
	int endu, int endv
){
	//Construct twist vectors( d2f/(du*dv) ).
	int j, ipu=vderiv.capacity(), ipv=uderiv.capacity();
	const double* pu=vderiv.data();
	const double* pv=uderiv.data();
	int ncd1=vderiv.sdim(), ncd2=uderiv.sdim();
	if(ncd<ncd1) ncd=ncd1; if(ncd<ncd2) ncd=ncd2;
	const int nu=utau.length(), nv=vtau.length();
	int idu=0, idv=0;
	if(endu==2) idu=nu-2;
	if(endv==2) idv=nv-2;
	int error;

	double* tangen1=new double[ncd+ncd+8*ncd+81];
	double* tangen2=tangen1+ncd;
	double* work12=tangen2+ncd;
	for(j=0; j<ncd; j++)
		tangen1[j]=tangen2[j]=0.;

	bvltn2_(nu,pu,utau.data(),ipu,ncd1,endu,u_dorder,
			work12,work12+70,work12+81,tangen1,&error);
	bvltn2_(nv,pv,vtau.data(),ipv,ncd2,endv,v_dorder,
			work12,work12+70,work12+81,tangen2,&error);
	for(j=0; j<ncd; j++) { tangen1[j]=tangen1[j]+tangen2[j];}

	MGVector tan12,tan22;
	if(u_dorder==1) tan12=(vderiv(idu+1)-vderiv(idu))*
						(0.4/(utau.ref(idu+1)-utau.ref(idu)));
	else tan12=MGVector(ncd,0.0);

	if(v_dorder==1) tan12+=(uderiv(idv+1)-uderiv(idv))*
						(0.4/(vtau.ref(idv+1)-vtau.ref(idv)));
	MGVector twist=MGVector(ncd,tangen1)-tan12;

	delete[] tangen1;
	return twist;
}
