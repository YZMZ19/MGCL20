/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/CompositeCurve.h"
#include "mg/SSisect.h"
#include "mg/Tolerance.h"
#include "topo/HHisect.h"
#include "topo/HHisects.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

//
//Implements MGHHisect Class.
//MGHHisect is to represent one continuous intersection line of a shell
//with a shell, a face, or a surface.
//(MGCompositeCurve* iline, deque<MGFPline> uvl1, deque<MGFPline> uvl2)
//where iline is a world coordinate rep of the line, uvl1 is a deque of
//1st shell's face parameter rep, uvl2 is a deque of the 2nd shell's or
//face's parameter rep of the intersection line.
//uvl1[i] corresponds to uvl2[i] one by one for all i.
//The parameter ranges of all the uvl1[i] are continuous and the total of them is
//equal to the parameter range of iline. For uvl2, the same.
//Let uvl1[i]'s start parameter be t1, and end parameter t2, then
//uvl2[i]'s parameter range is also from t1 to t2.
//Let sf1 be MGSurfCurve(f1's surface, uvl1[i]), then sf1 is the same curve as
//the iline's part of the parameter range t1 to t2. And sf1 is also equal to
//MGSurfCurve(f2's surface, uvl2[i]).
//
//MGHHisect uses MGFPline to represent the intersection lines.

//
//MGHHisect is also used to represent a pojection curve. In this case the size of
//uvline2 is zero.
//**** Projection line rep and intersection line rep cannot be mixed. ****

MGHHisect::MGHHisect(MGHHisect&& rhs):m_iline(rhs.m_iline.release()),
m_uvlines1(std::move(rhs.m_uvlines1)), m_uvlines2(std::move(rhs.m_uvlines2)){ ; }

MGHHisect& MGHHisect::operator=(MGHHisect&& rhs){
	m_iline.reset(rhs.m_iline.release());
	m_uvlines1=std::move(rhs.m_uvlines1);
	m_uvlines2=std::move(rhs.m_uvlines2);
	return *this;
}

//Construct from MGSSisect.
MGHHisect::MGHHisect(
	const MGFSurface* face1,	//face1. This must not be null.
	const MGFSurface* face2,	//face2. This may be null
							//(e.g. for face2 that is actually a surface).
	std::unique_ptr<MGSSisect>&& ssi///<intersection line of face1 and face2 expressed as MGSSisect.
):m_iline(new MGCompositeCurve(ssi->release_line()))
,m_uvlines1(1),m_uvlines2(1){
	m_uvlines1[0]= MGFPline(face1, ssi->release_param1());
	m_uvlines2[0]= MGFPline(face2, ssi->release_param2());
}

//Construct two faces intersection lines.
//uvline1 and 2 makes uvlines whose vector length is 1.
//MGHHisect takes the ownership of iline, uvline1, and uvline2
//(These must be newed objects).
MGHHisect::MGHHisect(
	MGCurve* iline,		//Intersection line of world coordinates.
	const MGFSurface* face1,//face1. This must not be null.
	MGCurve* uvline1,	//Intersection line of face1's (u,v) coordinates.
	const MGFSurface* face2,//face2. This may be null
						//(e.g. for face2 that is actually a surface).
	MGCurve* uvline2)	//Intersection line of face2's (u,v) coordinates.
						//takes the ownership of all the curves of ssi.
	:m_iline(new MGCompositeCurve(iline)), m_uvlines1(1){
	m_uvlines1[0]= MGFPline(face1, uvline1);
	if(uvline2){
		m_uvlines2.resize(1);
		m_uvlines2[0]=MGFPline(face2, uvline2);
	}
}

///////Operator overload///////

bool MGHHisect::operator< (const MGHHisect& hhi2)const{
	return m_uvlines1.front()<hhi2.m_uvlines1.front();
}
bool MGHHisect::operator<(const MGisect & is) const{
	const MGHHisect* is2=dynamic_cast<const MGHHisect*>(&is);
	if(!is2)
		return false;
	return (*this)<(*is2);
}
bool MGHHisect::operator== (const MGHHisect& hhi2)const{
	if(m_uvlines1!=hhi2.m_uvlines1) return false;
	return m_uvlines2==hhi2.m_uvlines2;
}

bool MGHHisect::operator==(const MGisect & is) const{
	const MGHHisect* is2=dynamic_cast<const MGHHisect*>(&is);
	if(!is2)
		return false;
	return (*this)==*is2;
}

///////Member function///////

bool MGHHisect::build_one(MGHHisect& hhi){
	if(hhi.is_null())
		return false;
	if(is_null()){
		(*this)=std::move(hhi);
		return true;
	}

	double error=MGTolerance::wc_zero_sqr();
	error*=4.;

	MGPosition p0=iline().start_point(), p1=iline().end_point();
	MGCurve& il=hhi.iline();
	MGPosition q0=il.start_point();
	MGVector dif(p0-q0);
	if(dif%dif<error){
		hhi.reverse_direction();
		connect_line_to_start(std::move(hhi));
		return true;
	}
	dif=(p1-q0);
	if(dif%dif<error){
		connect_line_to_end(std::move(hhi));
		return true;
	}
	MGPosition q1=il.end_point();
	dif=(p0-q1);
	if(dif%dif<error){
		connect_line_to_start(std::move(hhi));
		return true;
	}
	dif=(p1-q1);
	if(dif%dif<error){
		hhi.reverse_direction();
		connect_line_to_end(std::move(hhi));
		return true;
	}
	return false;
}

//Extract connected lines from hhivec one by one and build one continuous
//line. Extracted lines will be released from hhivec.
bool MGHHisect::build_one(MGHHisects& hhivec){
	bool extracted=false;
	size_t n=hhivec.size();
	MGHHisects::iterator j=hhivec.begin();
	for(size_t i=0; i<n; i++, j++){
		MGHHisect* hhi=static_cast<MGHHisect*>(j->get());
		if(build_one(*hhi))
			extracted=true;
	}
	return extracted;
}

//Change parameter range, be able to change the direction by providing
//t1 greater than t2.
void MGHHisect::change_range(
	double t0,		//Parameter value for the start of original. 
	double t1		//Parameter value for the end of original. 
){
	int n=m_iline->number_of_curves();
	if(!n) return;

	double oldsp=param_e()-param_s();
	double newsp=(t1-t0);
	double ratio=newsp/oldsp;
	double ts, te;
	int nm1=n-1;
	if(t0<t1){	//Case of no direction change.
		ts=t0;
		for(int i=0; i<nm1; i++){
			MGCurve& ci=m_iline->curve(i);
			te=ts+ci.param_span()*ratio;
			ci.change_range(ts,te);
			m_uvlines1[i].change_range(ts,te);
			ts=te;
		}
		m_iline->curve(nm1).change_range(ts,t1);
	}else{		//Case of direction change.
		ts=t1;
		int nhalf=n/2, i;
		for(i=0; i<nhalf; i++){//Change the curve ordering.
			int nm1mi=nm1-i;
			MGCurve* crv=m_iline->m_composite[i];
			m_iline->m_composite[i]=m_iline->m_composite[nm1mi];
			m_iline->m_composite[nm1mi]=crv;
			MGFPline fp1i(std::move(m_uvlines1[i]));
			m_uvlines1[i]=std::move(m_uvlines1[nm1mi]);
			m_uvlines1[nm1mi]=std::move(fp1i);
			if(has_face2_data()){
				MGFPline fp2i(std::move(m_uvlines2[i]));
				m_uvlines2[i]=std::move(m_uvlines2[nm1mi]);
				m_uvlines2[nm1mi]=std::move(fp2i);
			}
		}
		for(i=0; i<nm1; i++){//Change each parameter range(and direction).
			double te=ts-ratio*(m_iline->m_composite[i]->param_span());
			m_iline->m_composite[i]->change_range(te,ts);
			m_uvlines1[i].change_range(te,ts);
			if(has_face2_data())
				m_uvlines2[i].change_range(te,ts);
			ts=te;
		}
		m_iline->m_composite[nm1]->change_range(t0,ts);
		m_uvlines1[i].change_range(t0,ts);
		if(has_face2_data())
			m_uvlines2[i].change_range(t0,ts);
	}
}

//Connect a line to this HHisect.
//When both of face2 and uvlines2 are null, it indivates (face2, uvline2) are not
//used. This case occurs when MGHHisect is used to represent projection lines.
//iline, uvline1, and uvline2 must have the same direction.
//iline's direction must be equal to this HHisect's.
//MGHHisect takes the ownership of iline, uvline1, and uvline2
//(These must be newed objects).
void MGHHisect::connect_line_to_end(
	MGCurve* iline,		//Intersection line of world coordinates.
	const MGFSurface* face1,//face1. This must not be null.
	MGCurve* uvline1,	//Intersection line of face1's (u,v) coordinates.
	const MGFSurface* face2,//When face2 is null, and uvlines2!=null, it indicates
						//face2 is actually a surface.
	MGCurve* uvline2)	//Intersection line of face2's (u,v) coordinates.
						//takes the ownership of all the curves of ssi.
{
	assert(!uvline2 || m_uvlines1.size()==m_uvlines2.size());
	assert(uvline2 || m_uvlines2.size()==0);

	MGInterval rng=m_iline->connect_to_end(iline);
	double t0=rng.low_point(), t1=rng.high_point();
	uvline1->change_range(t0,t1);
	m_uvlines1.push_back(MGFPline(face1,uvline1));
	if(uvline2){
		uvline2->change_range(t0,t1);
		m_uvlines2.push_back(MGFPline(face2,uvline2));
	}
}
void MGHHisect::connect_line_to_end(
	MGHHisect&& hhi2		//After connected, this hhi2's member data's ownership
						//will be transfered to this MGHHisect,
						//just like std::auto_ptr's assignment.
){
	if(num_of_uvline()==0){
		*this=std::move(hhi2);
		return;
	}

	int n2=hhi2.num_of_uvline();
	for(int i=0; i<n2; i++){

	MGInterval rng=m_iline->connect_to_end(&(hhi2.m_iline->curve(i)));
	hhi2.m_iline->m_composite[i]=0;
	double t0=rng.low_point(), t1=rng.high_point();
	MGFPline fpl(std::move(hhi2.m_uvlines1[i]));
	fpl.uvline().change_range(t0,t1);
	m_uvlines1.push_back(std::move(fpl));

	if(!has_face2_data() || !hhi2.has_face2_data())
		continue;

	MGFPline fp2(std::move(hhi2.m_uvlines2[i]));
	MGCurve& uvline2=fp2.uvline();
	if(&uvline2){
		uvline2.change_range(t0,t1);
		m_uvlines2.push_back(std::move(fp2));
	}

	}
	hhi2.m_iline.reset();
	hhi2.m_uvlines1.clear();
	hhi2.m_uvlines2.clear();
}

//Connect a line to this HHisect.
//iline, uvline1, and uvline2 must have the same direction.
//iline's direction must be opposite to this HHisect's.
//MGHHisect takes the ownership of iline, uvline1, and uvline2
//(These must be newed objects).
void MGHHisect::connect_line_to_start(
	MGCurve* iline,		//Intersection line of world coordinates.
	const MGFSurface* face1,//face1. This must not be null.
	MGCurve* uvline1,	//Intersection line of face1's (u,v) coordinates.
	const MGFSurface* face2,//When face2 is null, and uvlines2!=null, it indicates
						//face2 is actually a surface.
	MGCurve* uvline2)	//Intersection line of face2's (u,v) coordinates.
						//takes the ownership of all the curves of ssi.
{
	assert(!uvline2 || m_uvlines1.size()==m_uvlines2.size());
	assert(uvline2 || m_uvlines2.size()==0);

	iline->negate();
	MGInterval rng=m_iline->connect_to_start(iline);
	double t0=rng.low_point(), t1=rng.high_point();
	uvline1->change_range(t1,t0);
	m_uvlines1.push_front(MGFPline(face1,uvline1));
	if(uvline2){
		uvline2->change_range(t1,t0);
		m_uvlines2.push_front(MGFPline(face2,uvline2));
	}
}
void MGHHisect::connect_line_to_start(
	MGHHisect&& hhi2		//After connected, this hhi2's member data's ownership
						//will be transfered to this MGHHisect,
						//just like std::auto_ptr's assignment.
){
	hhi2.connect_line_to_end(std::move(*this));
	*this=std::move(hhi2);
	hhi2.m_iline=0;
	hhi2.m_uvlines1.clear();
	hhi2.m_uvlines2.clear();
}

//Release the pointer of the last curve.
//Returned will be the released MGCurve pointer.
void MGHHisect::release_back(
	MGCurve*& ilineLast,
	MGFPline& uvline1Last,
	MGFPline& uvline2Last
){
	ilineLast=m_iline->release_back();
	uvline1Last=std::move(m_uvlines1.back()); m_uvlines1.pop_back();
	if(has_face2_data()){
		uvline2Last=std::move(m_uvlines2.back()); m_uvlines2.pop_back();
	}
}

//Release the pointer of the 1st curve.
//Returned be the released MGCurve pointer.
void MGHHisect::release_front(
	MGCurve*& iline1st,
	MGFPline& uvline11st,
	MGFPline& uvline21st
){
	iline1st=m_iline->release_front();
	uvline11st=std::move(m_uvlines1.front()); m_uvlines1.pop_front();
	if(has_face2_data()){
		uvline21st=std::move(m_uvlines2.front()); m_uvlines2.pop_front();
	}
}

//Reverse the direction of this intersection line.
void MGHHisect::reverse_direction(){
	double t0=param_s(), t1=param_e();
	change_range(t1,t0);
}

//Return i-th uvline.
const MGFPline& MGHHisect::uvline1(int i)const{
	assert(size_t(i)<m_uvlines1.size());
	return m_uvlines1[i];
}
const MGFPline& MGHHisect::uvline2(int i)const{
	assert(size_t(i)<m_uvlines2.size());
	return m_uvlines2[i];
}

//Replace 1st and 2nd lines.
void MGHHisect::exchange12(){
	if(!has_face2_data())
		return;
	assert(m_uvlines1.size()==m_uvlines2.size());
	container_type work(std::move(m_uvlines1));
	m_uvlines1=std::move(m_uvlines2);
	m_uvlines2=std::move(work);
}

std::ostream & MGHHisect::toString(std::ostream & ostrm) const{
	const MGHHisect& hhi=*this;
	if(hhi.is_null())
		return ostrm;
	ostrm<<"MGHHisect::m_iline="<<*(hhi.m_iline);
	int n=hhi.num_of_uvline();
	int n2=(int)hhi.m_uvlines2.size();
	ostrm<<"num_of_uvline="<<n<<endl;
	for(int i=0; i<n; i++){
		ostrm<<"m_uvlines1["<<i<<"]="<<hhi.m_uvlines1[i];
		if(n2) ostrm<<"m_uvlines2["<<i<<"]="<<hhi.m_uvlines2[i];
	}
	return ostrm;
}