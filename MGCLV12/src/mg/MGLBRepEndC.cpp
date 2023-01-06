/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/NDDArray.h"
#include "mg/BPointSeq.h"
#include "mg/LBRepEndC.h"
#include "mg/Curve.h"

#include "cskernel/Bvltan.h"
#include "cskernel/bvutan.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGLBRepEndC.cc
//
// Implements MGLBRepEndC class.

//Member Data
//	MGENDCOND m_cond;	//Type of end condition.
//	MGVector  m_1deriv;	//1st derivative stored 
//						//when m_cond=MGENDC_1D or MGENDC_12D
//	MGVector  m_2deriv;	//2nd derivative stored

//Constructor

MGLBRepEndC::MGLBRepEndC(	//Construct of m_cond=MGENDC_1D or MGENDC_2D.
	MGENDCOND cond,			//Type of end condition
							//(MGENDC_1D or MGENDC_2D)
	const MGVector& deriv)	//Derivative inf according to cond
	: m_cond(cond) {
	switch (cond){
	case MGENDCOND::MGENDC_1D:
		m_1deriv=deriv; break;
	case MGENDCOND::MGENDC_2D:
		m_2deriv=deriv; break;
	default:
		m_cond= MGENDCOND::MGENDC_NO; break;
	}
}

MGLBRepEndC::MGLBRepEndC(	//Construct of m_cond=MGENDC_12D
	const MGVector& first_deriv,	//1st derivative
	const MGVector& second_deriv)	//2nd derivative
	: m_cond(MGENDCOND::MGENDC_12D)
	, m_1deriv(first_deriv), m_2deriv(second_deriv) {;}

MGLBRepEndC::MGLBRepEndC(		//BVLTAN
	int start,					//Indicates start(start==true) condition
								// or end.
	const MGNDDArray& tau,		//Data point abscissa
	const MGBPointSeq& points,	//Point seq data
	int &error)					//Error flag.
	:m_cond(MGENDCOND::MGENDC_1D)
{
	assert(tau.length()==points.length() && tau.length()>1);

	const int np=points.length();
	const int ip=points.capacity();
	const int ncd=points.sdim();
	int ise=2; if(start) ise=1;
	double* work=new double[70+11+8*ncd+ncd]; 
	double* work1=work+70;
	double* work2=work1+11;
	double* tangen=work2+8*ncd;
	bvltan_(np,points.data(),tau.data(),ip,ncd,ise,
			work,work1,work2,tangen,&error);
	if(error==1){
		error=0;
		m_1deriv=MGVector(ncd,tangen);
	}
	delete[] work;
}

MGLBRepEndC::MGLBRepEndC(		//BVUTAN
	int start,			//Indicates start(start==true)
						// condition or end.
	const MGBPointSeq& points)	//Point seq data
	:m_cond(MGENDCOND::MGENDC_1D)
{
	assert(points.length()>1);

	const int np=points.length();
	const int ip=points.capacity();
	const int ncd=points.sdim();
	int ise=2; if(start) ise=1;
	double* tangen=new double[ncd+8*ncd+81];
	double* work=tangen+ncd;
	int error;
	bvutan_(ise,np,points.data(),ip,ncd,work,tangen,&error);
	m_1deriv=MGVector(ncd,tangen);
	delete[] tangen;
}

// Given MGCurve, construct the curve's end condition.
MGLBRepEndC::MGLBRepEndC(
	int start,		//Indicates start(start==true) condition or end.
	MGENDCOND cond,	//Type of end condition(MGENDC_1D, MGENDC_2D, or MGENDC_12D)
	const MGCurve& curve//Curve
):m_cond(cond){
	double tau;
	if(start) tau=curve.param_s();
	else tau=curve.param_e();

	if(cond== MGENDCOND::MGENDC_1D || cond== MGENDCOND::MGENDC_12D){
		m_1deriv=curve.eval(tau,1);
	}
	if(cond== MGENDCOND::MGENDC_2D || cond== MGENDCOND::MGENDC_12D){
		m_2deriv=curve.eval(tau,2);
	}
}

//Destructor
//	~MGLBRepEndC();	  We use default destructor.

//Member Function

//Initialize the instance. Will be set to the same as constructed by the void 
//constructor.
void MGLBRepEndC::initialize(){
	m_cond= MGENDCOND::MGENDC_NO;
	m_1deriv.set_null();
	m_2deriv.set_null();
}

//Set 1st deriv and change condition type to MGENDC_1D or MGENDC_12D.
void MGLBRepEndC::set_1st(const MGVector& first_deriv)
{
	switch (m_cond){
	case MGENDCOND::MGENDC_NO:
		m_cond= MGENDCOND::MGENDC_1D; m_1deriv=first_deriv; break;
	case MGENDCOND::MGENDC_2D:
		m_cond= MGENDCOND::MGENDC_12D; m_1deriv=first_deriv; break;
	default:
		m_1deriv=first_deriv; break;
	}
}

//Set 2nd deriv and change condition type to MGENDC_2D or MGENDC_12D.
void MGLBRepEndC::set_2nd(const MGVector& second_deriv)
{
	switch (m_cond){
	case MGENDCOND::MGENDC_NO:
		m_cond= MGENDCOND::MGENDC_2D; m_2deriv=second_deriv; break;
	case MGENDCOND::MGENDC_1D:
		m_cond= MGENDCOND::MGENDC_12D; m_2deriv=second_deriv; break;
	default:
		m_2deriv=second_deriv; break;
	}
}

//Operator overload.
