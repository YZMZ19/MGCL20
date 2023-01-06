/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Vector.h"
#include "mg/PPRep.h"
#include "mg/LBRep.h"
#include "mg/Tolerance.h"

#include "cskernel/bpval2.h"
#include "cskernel/bpval.h"
#include "cskernel/Blcbpn.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
// Implements MGPPRep class.

//Constructor
MGPPRep::MGPPRep()
//Default constructor.
:m_order(0), m_sdim(0),m_nbreak(0), m_coef(nullptr){;}

//Copy constructor.
MGPPRep::MGPPRep(const MGPPRep& rhs)
:m_break_point(rhs.m_break_point), m_nbreak(rhs.m_nbreak)
,m_order(rhs.m_order), m_sdim(rhs.m_sdim){
	int len=m_order*m_nbreak*m_sdim;
	m_coef=new double[len];
	for(int i=0; i<len; i++) m_coef[i]=rhs.m_coef[i];
}

//Move constructor.
MGPPRep::MGPPRep(MGPPRep&& rhs)
:m_break_point(std::move(rhs.m_break_point)), m_nbreak(rhs.m_nbreak)
,m_order(rhs.m_order), m_sdim(rhs.m_sdim), m_coef(rhs.m_coef){
	rhs.m_order=rhs.m_nbreak=rhs.m_sdim=0;
	rhs.m_coef=nullptr;
}

MGPPRep::MGPPRep(int order, int nbreak, int sdim)
	//Constructor of dummy PPRep of specified size.
	:m_order(order),m_nbreak(nbreak),m_sdim(sdim),
	m_break_point(nbreak),m_coef(new double[order*nbreak*sdim])
{
	if(order*nbreak*sdim==0) m_coef=0;
	else{
		int n=order*nbreak*sdim;
		for(int i=0; i<n; i++) m_coef[i]=0.;
	}
	assert(nbreak>=2);
}

//Constructor of dummy PP-Rep(no data except data points) of specified size.
MGPPRep::MGPPRep(int order, int sdim, const MGNDDArray& tau)
:m_order(order),m_nbreak(tau.length()),m_sdim(sdim),
m_break_point(tau),m_coef(new double[order*tau.length()*sdim]
){
	if(order*m_nbreak*sdim==0) m_coef=0;
	else{
		int n=order*m_nbreak*sdim;
		for(int i=0; i<n; i++) m_coef[i]=0.;
	}
	assert(m_nbreak>=2);
}

MGPPRep::MGPPRep(const MGLBRep& lbrep)
	//Constructor to convert from Line B-Representation.
	:m_order(lbrep.order()),				//Order
	m_nbreak(lbrep.bdim()-m_order+2),		//Num of Break points
	m_sdim(lbrep.sdim()), 					//Space dimension
	m_break_point(m_nbreak),				//Break point area
	m_coef(new double[m_order*m_nbreak*m_sdim])	//Coef area
{
	const int irc=lbrep.line_bcoef().capacity();
	double* work=new double[m_order*m_order*m_sdim]; //Work array for BLCBP
	int nbdim= lbrep.bdim();
	const double* knotp= lbrep.knot_data();
	int nbreak;
	double* breakp=&(m_break_point(0));
	blcbpn_(m_order,nbdim,knotp,lbrep.coef_data(),irc,
		   m_sdim,m_nbreak,work,breakp,&coef(0,0,0),&nbreak);
	delete[] work;
	nbreak+=1; //Since ouput nbreak of BLCBP is interval number.
	if(nbreak<int(m_nbreak))
		reshape(nbreak);
}

MGPPRep::MGPPRep(int order, const MGPPRep& pp1)
//Constructor to change order of original PP-Representation.
//New order may be greater or less than the original one. However,
//if new one is less than the original, PP-Rep constructed may not
//be able to hold the same shape.
:m_order(order)
,m_nbreak(pp1.m_nbreak)
,m_sdim(pp1.m_sdim)
,m_break_point(pp1.m_break_point)
,m_coef(new double[m_order*m_nbreak*m_sdim])
{
	int i,j,k;
	for(k=0; k<m_sdim; k++)
		for(j=0; j<m_nbreak-1; j++)
			for(i=0; i<m_order; i++)
				coef(i,j,k)=pp1.ref(i,j,k);
}


double MGPPRep::coef(int i, int j, int k) const
    //Returns (i,j,k)-th coef value.
{
	return ref(i,j,k);
}

double& MGPPRep::coef(int i, int j, int k)
    //Returns a pointer to coef(i,j,k) to access.
{
	assert(i<m_order && j<m_nbreak && k<m_sdim);

	return m_coef[i+m_order*j+m_order*m_nbreak*k];
}

MGVector MGPPRep::coef(int i, int j) const
//Extract (i,j,k)elements for 0<=k<sdim() as a vector.
{
	MGVector v(m_sdim);
	for(int k=0; k<m_sdim; k++) v(k)=ref(i,j,k);
	return v;
}

//Returns a pointer to the PPCoef data.
const double* MGPPRep::coef_data(int i, int j, int k) const
{return &m_coef[i+m_order*j+m_order*m_nbreak*k];}

MGVector MGPPRep::eval(	//Evaluate right continuous n'th derivative(BPVAL)
		double t,		//Parameter value to evaluate.
		int n		//Dgree of derivative.
		 ) const{		//When n=0, compute positional data.
	MGVector v(m_sdim);
	for(int i=0; i<m_sdim; i++){
		v(i)=bpval_(break_point_data(), &m_coef[m_order*m_nbreak*i]
					, m_nbreak, m_order, t, n);
	}
	return v;
}

//Evaluate i-th span's n'th derivative(BPVAL)
MGVector MGPPRep::eval_i(
		int i,		//span number(from 0) of the PP-Rep.
		double t,		//Parameter value to evaluate.
		int n		//Dgree of derivative.
)const{					//When n=0, compute positional data.
	const int nbreak=1; double bp=m_break_point(i);
	MGVector v(m_sdim);
	for(int m=0; m<m_sdim; m++){
		v(m)=bpval2_(
			&bp,&m_coef[m_order*i+m_order*m_nbreak*m]
			,m_order,t,nbreak,n);
	}
	return v;
}

//"normalize" normalizes the PP-Rep, i.e. changes break point data
//and pp-coefficients so as that length of first derivatives
// from left and right at each break point are the same.
MGPPRep& MGPPRep::normalize(){
	if(m_nbreak<=2) return *this;

	MGNDDArray bp(break_point());
	int i,j,k;
	double len1,len2,len;
	MGVector v1, v2;
	v1=eval_i(0,bp(1),1); len1=v1.len();
	//len1 is length of 1st derivatives at break_point(1) of the 1st span.
	if(!MGRZero(len1)){
		//Length of 1st deriv at the end of 1st span is set to unit.
		break_point(0)=bp(0)*len1;
		break_point(1)=bp(1)*len1;
		len=1./len1;
		for(i=1; i<m_order; i++){ 
			for(k=0; k<m_sdim; k++) coef(i,0,k)=coef(i,0,k)*len;
			len /=len1;
		}
	}
	for(j=1; j<m_nbreak-1; j++){
		//Change the length of the next span(len2) so as to be the same as
		//previous span's end point(len1).
		v1=eval_i(j-1,break_point(j),1); len1=v1.len();
		v2=coef(1,j); len2=v2.len();
		if(!MGRZero(len2)){
			len1=len1/len2; len=len1;
			for(i=1; i<m_order; i++){ 
				for(k=0; k<m_sdim; k++) coef(i,j,k)=coef(i,j,k)*len;
				len *=len1;
			}
			len=bp(j+1)-bp(j);
			len=len/len1;
			break_point(j+1)=break_point(j)+len;
		}
	}
	return *this;
}

double MGPPRep::ref(int i, int j, int k) const{
	assert(j<=m_nbreak-1);

	if(i>=m_order || k>=m_sdim) return 0.0;
	return m_coef[i+m_order*j+m_order*m_nbreak*k];
}

//Change size. Change of sdim not allowed.
//Stored data so far will be guarateed to hold in the same id of coef(i,j,k).
void MGPPRep::reshape(int nbreak){
	if(nbreak==m_nbreak) return ;

	double* data=new double[m_order*nbreak*m_sdim];
	//Reshape of pp coef.
	int nb=nbreak; if(nb>m_nbreak) nb=m_nbreak;
	for(int k=0; k<m_sdim; k++){
		for(int j=0; j<nb; j++){
			for(int i=0; i<m_order ; i++){
				data[i+m_order*j+m_order*nbreak*k]=coef(i,j,k);
			}
		}
	}

	delete[] m_coef; m_coef=data;
	m_break_point.reshape(nbreak);	//Reshape of break poits.
	m_nbreak=nbreak;
}

//Resize to (order, nbreak, dim).
//Resutl will contain garbages.
void MGPPRep::resize(
	int order,	//Order number.
	int nbreak,	//number of break points, includes last points.
	int dim)		//Space dimension.
{
	if(m_coef) delete[] m_coef;
	m_order=order; m_nbreak=nbreak; m_sdim=dim;
	m_coef=new double[order*nbreak*dim];
	m_break_point.resize(nbreak);
}

//Store the vector v at coef(i,j).
void MGPPRep::store_at(int i, int j, const MGVector& v){
	int sd=sdim();
	for(int m=0; m<sd; m++) coef(i,j,m)=v[m];
}

//Operator Function

//Copy Assignment
MGPPRep& MGPPRep::operator=(const MGPPRep& rhs){
	int order=rhs.m_order, nbreak=rhs.m_nbreak, dim=rhs.m_sdim;
	int len=order*nbreak*dim;
	if(m_order*m_nbreak*m_sdim!=len) resize(order,nbreak,dim);
	else{
		m_order=order; m_nbreak=nbreak; m_sdim=dim;
	}
	for(int i=0; i<len; i++) m_coef[i]=rhs.m_coef[i];
	m_break_point=rhs.m_break_point;
	return *this;
}

//Move Assignment
MGPPRep& MGPPRep::operator=(MGPPRep&& rhs){
	m_nbreak=rhs.m_nbreak;
	m_order=rhs.m_order;
	m_sdim=rhs.m_sdim;
	m_coef=rhs.m_coef;
	m_break_point=std::move(rhs.m_break_point);

	rhs.m_order=rhs.m_nbreak=rhs.m_sdim=0;
	rhs.m_coef=nullptr;
	return *this;
}
