/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "cskernel/Blctpb.h"
#include "cskernel/Bluakt.h"
#include "cskernel/Bludkt.h"
#include "mg/Tolerance.h"
#include "mg/KnotArray.h"
#include "mg/PPRep.h"
#include "mg/LBRep.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGLBRep.cpp
//
// Implement MGLBRep class.

///This is a polar coordinate system data.
///Given polar coordinate LBRep, update this to ordinary coordinates system MGLBRep.
///This curve's (x,y) coordinates are polar coordinates system(r,theta), 
///where r is the distance from origin and theta is the angel with x coordinate.
///When return from the function (x,y) are ordinary coordinate system.
///The space dimension of this curve must be >=2;
///If this space dimension is lager than 2, the remaining coordinates are unchanged.
void MGLBRep::updatePolarCoordinates2Ordinary(){
	int n=bdim();
	assert(bdim()>=2);
	MGBPointSeq& bp=line_bcoef();
	for(int i=0; i<n; i++){
		double r=bp(i,0);
		double theta=bp(i,1);
		double x,y;
		if(r<=MGTolerance::wc_zero()){
			x=y=0.;
		}else{
			x=r*cos(theta);
			y=r*sin(theta);
		}
		bp(i,0)=x;
		bp(i,1)=y;
	}
	invalidateBox();
}

//Build MGLBRep, given PP-Rep(convert PP-Rep to MGLBRep) and
//Knot Vector in m_knot_vector of this.
//Each knot of the knot vector is break point of pprep.
//The continuities at all the break points must be C(k-2) where
//k is the order of pprep.
void MGLBRep::buildLBfromPPRep(
	const MGPPRep& pp	//PP-rep
){
	assert(m_knot_vector.order()==pp.order());
	invalidateBox();

	int n=m_knot_vector.bdim(), ncd=pp.sdim();
	m_line_bcoef.resize(n, ncd);

	int k=pp.order();	//order.
	int ipc=pp.nbreak();
	int l=ipc-1;//Number of spans.
	const int ipcw=n+k-2;
	int irc=m_line_bcoef.capacity();

	double* pcwork=new double[ipcw*k*ncd];
	blctpb_(pp.break_point_data(),pp.coef_data(),l,k,ncd,n,
			m_knot_vector.data(),ipc,ipcw,irc, pcwork, &m_line_bcoef(0,0));
	delete[] pcwork;
}

//Gets new B-Rep by adding knots to an original B-Rep.
void MGLBRep::addKnots(
	const MGKnotArray& knots ///<Knots to add.
){
	invalidateBox();

	const int k=order();//order
	const int n1=bdim();
	const int ncd=m_line_bcoef.sdim();

	MGBPointSeq bpSave(m_line_bcoef);//save the original.
	const int irc1=bpSave.capacity();

	const int nad=(int)knots.size();
	double* tad=new double[nad+k*k];
	double* work=tad+nad;

	int* mlt=new int[nad];
	int mlt_total=0;
	MGKnotVector t(m_knot_vector);//save the original.
	int km1=k-1;
	for(int i=0; i<nad; i++){
		double tadi=tad[i]=knots[i];
		int j=t.locate(tadi);
		int mlt1=0,mlt2=knots[i].multiplicity();
		if(tadi==t[j]){
			int j1=j;
			while(j1>=km1 && t[j1-1]==tadi){
				j1--;
				mlt1++;
			};
		}
		mlt2+=mlt1;
		if(mlt2>=k)
			mlt2=km1;
		mlt[i]=mlt2;
		mlt_total+=mlt[i];
	}

	int n2 = n1+mlt_total;
	m_knot_vector.reshape(n2+k);
	m_line_bcoef.reshape(n2);
	int irc2 = m_line_bcoef.capacity();

	bluakt_( k, n1, t.data(), bpSave.data(),
			irc1, ncd, nad, tad, mlt, irc2, work,
			&n2, &m_knot_vector(0), &m_line_bcoef(0,0));
	m_knot_vector.set_bdim(n2);
	m_line_bcoef.set_length(n2);

	delete[] tad; delete[] mlt;
}

///Approximate this, given by a new knot configuration newT.
///The new approximated MGLBRep is output into newLB.
///The parameter range of newT must be inside the one of this.
///newLB is an approximation of this with the parameter range
///from newT.param_s() to newT.param_e().
///Function's return value is error code:
/// =0 successfully rebuilt, !=0: input newT is invalid.
///When error!=0 is returned, this is copied int newLB.
int MGLBRep::rebuildByNewKnotVector(
	const MGKnotVector& newT,//New knot configuration is input.
	MGLBRep& newLB//Rebuilt MGLBRep is output, can be this.
)const{
	assert(newT.order()==order());

	int k=newT.order(), ncd=sdim();
	int n1=m_line_bcoef.length(), n2=newT.bdim();
	int irc1=m_line_bcoef.capacity();

	MGBPointSeq bp(n2, ncd);
	int irc2=bp.capacity();

	int error;
	int m=(k<=5) ? 9:2*k-1;
	double* work1=new double[n2+n2*m];
	double* work2=work1+n2;
	bludkt_(k,n1,knot_data(),coef_data(),
			irc1,ncd,n2, newT.data(),irc2,
			work1, work2, &bp(0,0), &error);
	delete[] work1;

	if(error==1){
		error=0;		//Return of bludkt_:error=1 means normal.
		newLB.m_knot_vector=newT;
		newLB.m_line_bcoef=std::move(bp);
		newLB.invalidateBox();
		newLB.copy_appearance(*this);
	}else if(this!=&newLB){
		newLB=*this;
	}
	return error;
}
