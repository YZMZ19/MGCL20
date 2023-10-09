/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Box.h"
#include "mg/BPointSeq.h"
#include "mg/CCisect.h"
#include "mg/CCisects.h"
#include "mg/CompositeCurve.h"
#include "mg/CParam_list.h"
#include "mg/CSisect.h"
#include "mg/CSisects.h"
#include "mg/Ellipse.h"
#include "mg/EReal.h"
#include "mg/Interval.h"
#include "mg/Knot.h"
#include "mg/KnotArray.h"
#include "mg/KnotVector.h"
#include "mg/LBRep.h"
#include "mg/LBRepEndC.h"
#include "mg/Matrix.h"
#include "mg/NDDArray.h"
#include "mg/Object.h"
#include "mg/OscuCircle.h"
#include "mg/OscuCircleData.h"
#include "mg/Plane.h"
#include "mg/Point.h"
#include "mg/Position.h"
#include "mg/Position_list.h"
#include "mg/PPRep.h"
#include "mg/RLBRep.h"
#include "mg/RSBRep.h"
#include "mg/SBRep.h"
#include "mg/BSumSurf.h"
#include "mg/SBRepTP.h"
#include "mg/SBRepEndC.h"
#include "mg/SPointSeq.h"
#include "mg/SSisect.h"
#include "mg/SSisects.h"
#include "mg/Straight.h"
#include "mg/SurfCurve.h"
#include "mg/Transf.h"
#include "mg/TrimmedCurve.h"
#include "mg/FSurface.h"
#include "mg/MGStl.h"
#include "mg/Tolerance.h"
#include "mgGL/Appearance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#include <iomanip>
using namespace std;

//////**********Object output*********//////////

//////////// MGObject ////////////
// Output function.
ostream& MGObject::toString(std::ostream& ostrm) const{
	ostrm << ", "<<whoami()<< "." << m_box;
	if(m_appearance)
		ostrm<<","<<(*m_appearance);

	return ostrm;
}

// Output virtual function.
std::ostream& MGCurve::toString(std::ostream& ostrm) const{
	MGGeometry::toString(ostrm);
	return ostrm;
}

// Output virtual function.
std::ostream& MGSurface::toString(std::ostream& ostrm) const{
	MGGeometry::toString(ostrm);
	return ostrm;
}

//////////// MGPoint ////////////
ostream& MGPoint::toString(ostream& ostrm)const{
	ostrm<<"MGPoint::"<< (const MGGel*)this;
	MGGeometry::toString(ostrm);

	int dim=sdim();
	ostrm <<",Q(";
	for(int i=0; i<dim; i++){
		ostrm<<m_point(i);
		if(i<dim-1) ostrm<<",";
	}
	ostrm<<")";
	return ostrm;
}

//////////// MGCompositeCurve ////////////
ostream& MGCompositeCurve::toString(ostream& ostrm) const{
	ostrm<<"MGCompositeCurve::"<< (const MGGel*)this;
	MGCurve::toString(ostrm);

	ostrm<<", number_of_curves="<<number_of_curves();
	deque<MGCurve*>::const_iterator
		i=m_composite.begin(), ie=m_composite.end();
	int j=0;
	for(; i!=ie; i++,j++)
		ostrm<<endl <<"Curve-"<<j<<"::"<<(**i);
	return ostrm;
}

//////////// MGEllipse ////////////
ostream& MGEllipse::toString(ostream& ostrm) const{
//	ostrm.setf(ios::scientific,ios::floatfield);
//	ostrm.precision(10);
	ostrm<<"MGEllipse::"<< (const MGGel*)this;
	MGCurve::toString(ostrm);

	double t0=m_prange[0], t1=m_prange[1];
	ostrm<<", circle="<<(int)m_circle
		<<", r="<<m_r
		<<", center="<<m_center<<endl
		<<", normal="<<m_normal
		<<", m="<<m_m
		<<", n="<<m_n<<endl
		<<", prange=("<<t0<<","<<t1<<")"
		<<" , start="<<eval(t0)
		<<" , end="<<eval(t1);
	if(m_gprange)
		ostrm<<endl<<", gprange=("<<m_gprange[0]<<","<<m_gprange[1]<<")";
	return ostrm;
}

//////////// MGLBRep ////////////
ostream& MGLBRep::toString(ostream& ostrm) const{
	ostrm<<"MGLBRep::"<< (const MGGel*)this;
	MGCurve::toString(ostrm);

	ostrm<<", Space Dim="<<sdim()<<endl;
	ostrm<<", "<<knot_vector();
	ostrm<<", "<<line_bcoef();
	return ostrm;
}

//////////// MGRLBRep ////////////
ostream& MGRLBRep::toString(ostream& ostrm) const{
	ostrm<<"MGRLBRep(non homogeneous form)::"<< (const MGGel*)this;
	MGCurve::toString(ostrm);

	ostrm<<", Space Dimension="<<sdim()<<endl;
	ostrm<<", "<<knot_vector();
	ostrm<<", "<<non_homogeneous_bcoef();
	return ostrm;
}

//////////// MGPlane ////////////
ostream& MGPlane::toString(ostream& ostrm) const{
//	o.setf(ios::scientific,ios::floatfield);
//	o.precision(10);
	ostrm<<"MGPlane::"<< (const MGGel*)this;
	MGSurface::toString(ostrm);

	ostrm<<endl<<", root_point="<<m_root_point
		<<", normal="<<m_normal<<", d="<<m_d
		<<", uderiv="<<m_uderiv<<", vderiv="<<m_vderiv<<endl;
	return ostrm;
}

//////////// MGRSBRep ////////////
ostream& MGRSBRep::toString(ostream& ostrm) const{
	ostrm<<"MGRSBRep(non homogeneous form)::"<< (const MGGel*)this;
	MGSurface::toString(ostrm);

	ostrm<<", Space Dim="<<sdim()<<endl;
	ostrm<<", uknot="<<knot_vector_u();
	ostrm<<", vknot="<<knot_vector_v();
	ostrm<<", surface_bcoef="<<non_homogeneous_bcoef();
	return ostrm;
}

//////////// MGSBRep ////////////
ostream& MGSBRep::toString(ostream& ostrm) const{
	ostrm<<"MGSBRep::"<< (const MGGel*)this;
	MGSurface::toString(ostrm);

	ostrm<<std::endl<<", uknot="<<m_uknot;
	ostrm<<", vknot="<<m_vknot;
	ostrm<<", surface_bcoef="<<m_surface_bcoef;
	return ostrm;
}

//////////// MGStraight ////////////
ostream& MGStraight::toString(ostream& ostrm) const{
	ostrm<<"MGStraight::"<< (const MGGel*)this;
	MGCurve::toString(ostrm);

	ostrm<<endl<<", root_point="<<m_root_point
		 <<", dir="<<m_direction
		 <<", sparam="<<m_sparam<<", eparam="<<m_endparam;
	if(m_sparam.finite())
		ostrm<<", start="<<eval(m_sparam.value());
	if(m_endparam.finite())
		ostrm<<", end="<<eval(m_endparam.value());
	//ostrm<<endl;
	return ostrm;
}

//////////// MGSurfCurve ////////////
ostream& MGSurfCurve::toString(ostream& ostrm) const{
	ostrm<<"MGSurfCurve::"<< (const MGGel*)this;
	MGCurve::toString(ostrm);

	ostrm<<", start="<<start_point()<<", end="<<end_point();
	ostrm<<", surface="<<m_surface
		<<", curve="<<m_curve << endl;
	return ostrm;
}

//////////// MGTrimmedCurve ////////////
ostream& MGTrimmedCurve::toString(ostream& ostrm) const{
	ostrm<<"MGTrimmedCurve::"<< (const MGGel*)this;
	MGCurve::toString(ostrm);

	ostrm<<", curve="<<m_curve;
	if(m_curve)
		ostrm<<", sameRange="<<m_sameRange
		<<", range="<<m_range
		<<", start="<<start_point()<<", end="<<end_point()
		<<", curve="<<endl<<(*m_curve);
	return ostrm;
}

//Debug Function
// Output virtual function.
ostream& MGBSumSurf::toString(std::ostream& ostrm) const{
	ostrm<<"MGBSumSurf::"<< (const MGGel*)this;
	MGSurface::toString(ostrm);

	ostrm<<", g1="<<m_g1;if(m_g1) ostrm<<","<<(*m_g1);
	ostrm<<", g2="<<m_g1;if(m_g2) ostrm<<","<<(*m_g2);
	ostrm<<", g12="<<m_g1;if(m_g12) ostrm<<","<<(*m_g12);
	return ostrm;
}

ostream& MGStl::toString(std::ostream& ostrm) const{
	ostrm<<"MGStl::"<< (const MGGel*)this;
	MGObject::toString(ostrm);

	// �{�b�N�X�g�̍��W���o��
	ostrm<<",box="<<m_box;

	//Print out triangle data.
	ostrm<<std::endl;
	int nIndices = int(m_indices.size());
	int nTri=nIndices/3;
	assert(nTri*3==nIndices);
	ostrm<<"Number of triangles="<<nTri<<endl;
	for(int i=0, j=0; j<nTri; i+=3, j++){
		int i0=m_indices[i], i1=m_indices[i+1], i2=m_indices[i+2];
		const MGPosition& P0=m_vecPos[i0];
		const MGPosition& P1=m_vecPos[i1];
		const MGPosition& P2=m_vecPos[i2];
		ostrm<<j<<": Normal="<<m_vecNormlTriang[j]<<endl;
		ostrm<<"0("<<P0[0]<<","<<P0[1]<<","<<P0[2]<<"), 1(";
		ostrm<<P1[0]<<","<<P1[1]<<","<<P1[2]<<"), 2(";
		ostrm<<P2[0]<<","<<P2[1]<<","<<P2[2]<<")"<<endl;
	}
	return ostrm;
}

////////************Non Object output************///////////////

//////////// Box ////////////
ostream& operator<< (ostream& out, const MGBox& box){
	int dim=box.sdim();
	if (dim) {
		out << "Box(";
		for (int i = 0; i < dim; i++) {
			out << box.m_range[i];
			if (i < dim - 1) out << ",";
		}
		out << ")";
	}
	else {
		out << "Box(Null)";
	}

	return out;
}

//////////// BPointSeq ////////////
ostream& operator<<(ostream& out, const MGBPointSeq& bps){
	int space_dim=bps.sdim();
	int length=bps.length();
	out<<"BPointSeq:: capacity="<<bps.capacity()<<", sdim="<<space_dim;
	out<<", length="<<length;
	int m;
	for(int i=0; i<length ; i++){
		m=space_dim*i;
		if((m-int(m/12)*12)==0) out<<endl;
		out<<i<<"(";
		for(int j=0; j<(space_dim-1) ; j++) out<<bps(i,j)<<",";
		out<<bps(i,space_dim-1)<<") ";
	}
	//out<<endl;
	return out;
}

/// Output virtual function.
std::ostream& MGCCisects::toString(std::ostream& ostrm)const{
	const MGCCisects& lst=*this;
	ostrm<<"MGCCisects::m_error="<<lst.m_error;
	MGisects::toString(ostrm);
	return ostrm;
}

/// Output virtual function.
std::ostream& MGSSisects::toString(std::ostream& ostrm)const{
	const MGSSisects& lst=*this;
	ostrm<<"MGSSisects::surface1="<<(lst.surface1());
	ostrm<<",surface2="<<(lst.surface2());
	return MGisects::toString(ostrm);
}

//////////// MGCParam_list ////////////
ostream& operator<< (ostream& out, const MGCParam_list& lst){
	out<<"MGCParam_list::curve="<<(lst.m_curve);
	int n=lst.entries();
	out<<", number of param="<<n<<endl;
	MGCParam_list::const_iterator itr; int i=0;
	for(itr=lst.begin(); itr!=lst.end(); itr++)
		out<<i++<<":"<<(*itr)<<" ";
	out<<endl;
	return out;
}

/// Output virtual function.
std::ostream& MGCSisects::toString(std::ostream& ostrm)const{
	const MGCSisects& lst=*this;
	ostrm<<"MGCSisects::m_errort="<<(lst.m_errort)<<",m_errou="<<lst.m_erroru;
	ostrm<<",m_errorv="<<lst.m_errorv<<std::endl;
	MGisects::toString(ostrm);
	return ostrm;
}

//////////// MGDefault ////////////
ostream& operator<<(ostream& o, const MGDefault& def){
	o<<"mgNULL_VEC="<<mgNULL_VEC;
	o<<", mgZERO_VEC="<<mgZERO_VEC;
	o<<", mgZERO_VEC_2D="<<mgZERO_VEC_2D<<","<<endl;
	o<<"mgNULL_Pos="<<mgNULL_Pos;
	o<<", mgORIGIN="<<mgORIGIN;
	o<<", mgORIGIN_2D="<<mgORIGIN_2D<<","<<endl;
	o<<"mgX_UVEC="<<mgX_UVEC;
	o<<", mgY_UVEC="<<mgY_UVEC;
	o<<", mgZ_UVEC="<<mgZ_UVEC<<","<<endl;
	o<<"mgX_UVEC_2D="<<mgX_UVEC_2D;
    o<< ", mgY_UVEC_2D="<<mgY_UVEC_2D<<","<<endl;
	o<<"mgEMP_INTERV="<<mgEMP_INTERV;
	o<<", mgNULL_BOX="<<mgNULL_BOX<<","<<endl;
	o<<"mgEMP_BOX="<<mgEMP_BOX<<","<<endl;
	o<<"mgEMP_BOX_2D="<<mgEMP_BOX_2D<<","<<endl;
	o<<"mgNULL_MATR="<<mgNULL_MATR<<","<<endl;
	o<<"mgUNIT_MATR="<<mgUNIT_MATR<<","<<endl;
	o<<"mgUNIT_MATR_2D="<<mgUNIT_MATR_2D<<","<<endl;
	o<<"mgNULL_TRANSF="<<mgNULL_TRANSF<<","<<endl;
	o<<"mgID_TRANSF="<<mgID_TRANSF<<","<<endl;
	o<<"mgID_TRANSF_2D="<<mgID_TRANSF_2D<<endl;
	o<<"mgNULL_KNOT_VECTOR="<<mgNULL_KNOT_VECTOR;
	o<<"mgINFINITE_KNOT_VECTOR="<<mgINFINITE_KNOT_VECTOR;
	return o;
}

//////////// MGEReal ////////////
ostream& operator<< (ostream& out, const MGEReal& er){
	if(er.finite()) out<<er.value();
	else if(er.minus_infinite()) out<<"INFINITE_MINUS";
	else out<<"INFINITE_PLUS";
	return out;
}

//////////// MGInterval ////////////
ostream& operator<< (ostream &o, const MGInterval& i2){
//		o.setf( ios::scientific, ios::floatfield );
//        o.precision( 10 );
	if(i2.empty()) o<<"I"<<"[EMPTY]";
	else           o<<"I"<<"["<<i2.m_low<<","<<i2.m_high<<"]";
	return o;
}

//////////// MGKnot ////////////
ostream& operator<<(ostream& out, const MGKnot& knot){
	out<< "(" << knot.value() <<" " << knot.multiplicity() <<") ";
	return out;
}

//////////// MGKnotArray ////////////
ostream& operator<<(ostream& out, const MGKnotArray& klist){
	int n=(int)klist.size();
	out<<"MGKnotArray::length="<<n<<endl;
	for(int i=0; i<n; i++) out<<klist[i]<<" ";
	out<<endl;
	return out;
}

//////////// MGKnotVector ////////////
ostream & operator << (ostream & out, const MGKnotVector & t){
	int order=t.order(); int bdim=t.bdim();
	out<<"MGKnotVector::order="<< order<<", bdim=" <<bdim;
	out<<", current=" << t.m_current <<endl;
	int n=order+bdim; int nm1=n-1; int mult=1;
	if(n<=0) return out;
	double tprev=t(0)-(t(nm1)-t(0));
	for(int i=0; i<n; i++){
		if(t(i)==tprev) mult+=1;
		else if(i){
			if((i-int(i/12)*12)==0) out<<endl;
			out<<mult<<")"<<t(i-1)<<" ";
			mult=1;
		}
		tprev=t(i);
		if(i==nm1){
			if((i-int(i/12)*12)==0) out<<endl;
			out<<mult<<")";
			out<<tprev<<" ";
		}
	}
	out<<endl;
	return out;
}

//////////// MGLBRepEndC ////////////
ostream& operator<< (ostream& out, const MGLBRepEndC& endc){
	out<<"MGLBRepEndC::cond="<<endc.cond();	
	if(endc.cond()== MGENDCOND::MGENDC_1D || endc.cond()== MGENDCOND::MGENDC_12D)
		out<<", 1st deriv="<<endc.first();
	if(endc.cond()== MGENDCOND::MGENDC_2D || endc.cond()== MGENDCOND::MGENDC_12D)
		out<<", 2nd deriv="<<endc.second();
	out<<endl;
	return out;
}

//////////// MGMatrix ////////////
ostream& operator<< (ostream& out, const MGMatrix& mat){
	int dim=mat.sdim();
    out<<endl<<"MGMatrix(";
	for(int i=0; i<dim; i++){
		if(i>0) out<<"         ";
		for(int j=0; j<dim; j++){
			out<<mat(i,j); if(j<dim-1) out<<",";
		}
		if(i<dim-1) out<<endl;
	}
	out<<")";
	return out;
}

//////////// MGNDDArray ////////////
ostream& operator<<(ostream & out, const MGNDDArray& ndda){
	int n=ndda.length();
	out<<"MGNDDArray::capacity="<<ndda.capacity()<<",length="<<n<<endl;
	for(int i=0; i<n; i++) out<<i<<":"<<ndda(i)<<" ";
	out<<endl;
	return out;
}

//////////// MGOscuCircle ////////////
ostream& operator<<(ostream& out, const MGOscuCircle& circle){
	int n=circle.length();
	out<<"MGOscuCircle::length="<<n<<endl;
	for(int i=0; i<n; i++) out<<circle(i)<<" ";
	out<<endl;
	return out;
}

//////////// MGOscuCircleData ////////////
ostream& operator<<(ostream& out, const MGOscuCircleData& circle){
	out<< "(" << circle.index() <<" " << circle.radius() <<") ";
	return out;
}

//////////// MGPosition ////////////
/*ostream& operator<< (ostream& out, const MGPosition& p){
	int dim=p.sdim();
	out <<"P(";
	for(int i=0; i<dim; i++){
		out<<p.m_element(i);
		if(i<dim-1) out<<",";
	}
	out<<")";
	return out;
}*/

//////////// MGPosition_list ////////////
ostream& operator<< (ostream& out, const MGPosition_list& lst){
	size_t n=lst.size();
	out<<"MGPosition_list::number of position="<<n<<endl;
	MGPosition_list::const_iterator i;
	int j=0;
	for(i=lst.begin(); i!=lst.end(); i++) out<<j++<<":"<<(*i)<<" ";
	out<<endl;
	return out;
}

//////////// MGPPRep ////////////
ostream& operator<< (ostream& out, const MGPPRep& pprep)
//	friend istream& operator>> (istream&, MGPPRep& );
{
	int order=pprep.order();
	int nbreak=pprep.nbreak();
	int sdim=pprep.sdim();
	out<<"MGPPRep::"<<&pprep<<"order=" <<order <<", nbreak=" <<nbreak;
	out<<", sdim=" <<sdim <<endl;
	if(!nbreak)
		return out;

	out<<"break_id break_point:: sdim:pp-coef(form order 0 to order-1)" <<endl;
	int j;
	for(j=0; j<nbreak-1; j++){
		out<<j<<" "<<pprep.break_point(j)<<"::";
		for(int k=0; k<sdim; k++){
			out<<" "<<k <<":";
			for(int i=0; i<order ; i++){
				out<<pprep.coef(i,j,k)<<" ";
			}
		}
		out<<endl;
	}
	out<<j<<" "<<pprep.break_point(j)<<"::"<<endl;
	return out;
}

//////////// MGSBRepEndC ////////////
ostream& operator<< (ostream& out, const MGSBRepEndC& endc){
	out<<"MGSBRepEndC::cond={"
		<<endc.m_cond[0]<<","<<endc.m_cond[1]<<","
		<<endc.m_cond[2]<<","<<endc.m_cond[3]<<"}"
		<<endl;	
	for(int i=0; i<4; i++){
		out<<" 11d["<<i<<"]="<<endc.m_11d[i];
		out<<" 12d["<<i<<"]="<<endc.m_12d[i];
		out<<" 21d["<<i<<"]="<<endc.m_21d[i];
		out<<" 22d["<<i<<"]="<<endc.m_22d[i];
		if(endc.cond(i)== MGENDCOND::MGENDC_1D || endc.cond(i)== MGENDCOND::MGENDC_12D)
			out<<" 1st deriv("<<i<<")="<<endc.first(i);
		if(endc.cond(i)== MGENDCOND::MGENDC_2D || endc.cond(i)== MGENDCOND::MGENDC_12D)
			out<<" 2nd deriv("<<i<<")="<<endc.second(i);
	}
	return out;
}

//////////// MGSBRepTP ////////////
ostream& operator<< (ostream& out, const MGSBRepTP& tp)
{
	out<<"MGSBRepTP::"<<endl;
	for(int i=0; i<4; i++)
		if(tp.specified(i)) out<<" TP"<<i<<":"<<tp.TP(i);
	return out;
}

//////////// MGSPointSeq ////////////
ostream& operator<<(ostream& out, const MGSPointSeq& sp){
	int w=8;
	int space_dim=sp.sdim();
	int lenu=sp.length_u(), lenv=sp.length_v();
	out<<"MGSPointSeq::capacityu="<<sp.capacity_u()<<", capacityv="<<sp.capacity_v()
	   <<", sdim="<<space_dim;
	out<<", lengthu="<<lenu<<", lengthv="<<lenv;
	int m;

	out <<setiosflags( ios::fixed )<<setprecision(3);
	for(int j=0; j<lenv ; j++){
		for(int i=0; i<lenu ; i++){
			m=space_dim*i;
			if((m-int(m/12)*12)==0) out<<endl<<i<<","<<j<<"(";
			else out<<i<<"(";
			for(int k=0; k<(space_dim-1) ; k++){
				out<<setw(w);
				out<<sp(i,j,k)<<",";
			}
			out<<setw(w);
			out<<sp(i,j,space_dim-1)<<") ";
		}
	}
	out<<setprecision(6)<<resetiosflags(ios::fixed)<<endl;
	return out;
}
//////////// MGTolerance ////////////
ostream & operator << (ostream & o, const MGTolerance & tl){
	o << "MGTolerance::mach_zero="<<tl.m_mach_zero
	  << ", wc_zero="<<tl.m_wc_zero
		<<", wc_zero_sqr="<<tl.m_wc_zero_sqr<<endl;
	o <<", rc_zero="<<tl.m_rc_zero<<", rc_zero_sqr="<<tl.m_rc_zero_sqr
		<<", angle_zero="<<tl.m_angle_zero<<endl;
	o <<", line_zero="<<tl.m_line_zero<<", max_knot_ratio="
		<<tl.m_max_knot_ratio<<", stack_count="<<tl.m_count<<endl;
	return o;
}

//////////// MGTransf ////////////
ostream& operator<< (ostream& o, const MGTransf& trn1){
    o<<"MGTransf::"<<trn1.m_affine<<", "
		<<trn1.m_translation;
	return o;
}
//////////// MGUnit_vector ////////////
/*ostream& operator<< (ostream& out, const MGUnit_vector& unit) {
	int dim=unit.sdim();
	out <<"U(";
	for(int i=0; i<dim; i++){
		out<<unit.m_element[i];
		if(i<dim-1) out<<",";
	}
	out<<")";
	return out;
}*/

//////////// MGVector ////////////
ostream& operator<< (ostream& out, const MGVector& vec){
	int dim=vec.sdim();
	out <<vec.whoami()<<"(";
	for(int i=0; i<dim; i++){
		out<<vec.m_element[i];
		if(i<dim-1) out<<",";
	}
	out<<")";
	return out;
}

// MGINTERVAL_TYPE��W���o�͂ɏo�͂���B
ostream & operator << (ostream & out, MGINTERVAL_TYPE rel){
	const char* name[5]={
     "MGINTERVAL_EMPTY",             // Interval �͋�W��
     "MGINTERVAL_FINITE",            // ����A�����Ƃ��L��
     "MGINTERVAL_FINITE_ABOVE",      // ����L���A��������
     "MGINTERVAL_FINITE_BELOW",      // ��������A�����L��
     "MGINTERVAL_INFINITE"};         // ����A�����Ƃ��ɖ��� 
	 int i=(int)rel;
	 if(i>=5) i=4; if(i<0) i=0;
	out<<name[i];
	return out;
}

// MGPSRELATION��W���o�͂ɏo�͂���B
ostream & operator << (ostream & out, MGPSRELATION rel){
	const char* name[6]={
     "MGPSREL_UNKNOWN",              // �s��
     "MGPSREL_TORSION",              // �˂���
     "MGPSREL_ISECT",               // ����
     "MGPSREL_PARALLEL",             // ���s
     "MGPSREL_COIN",                 // ���
     "MGPSREL_VIRTUAL_ISECT"};      // �w�蒼���̉�����̓_�ł̌���
	int i=(int)rel;
	if(i>=6) i=5; if(i<0) i=0;
	out<<name[i];
	return out;
}

// MGCURVE_TYPE��W���o�͂ɏo�͂���B
ostream & operator << (ostream & out, MGCURVE_TYPE rel){
	const char* name[7]={
	"MGCURVE_UNKNOWN",		// �s��
	"MGCURVE_STRAIGHT",		// ����
	"MGCURVE_ELLIPSE",		// �ȉ~
	"MGCURVE_SPLINE",		//B-Spline(LBRep). �X�v���C��
	"MGCURVE_RSPLINE",		//B-Spline(Rational)
	"MGCURVE_USER1",		//Auxiliary curve type 1
	"MGCURVE_USER2"};		//Auxiliary curve type 2
	int i=(int)rel;
	if(i>=7 || i<0) i=0;
	out<<name[i];
	return out;
}

// MGELLIPSE_TYPE��W���o�͂ɏo�͂���B
ostream & operator << (ostream & out, MGELLIPSE_TYPE rel){
	const char* name[3]={
	 "MGELLIPSE_EMPTY"		// ��
	,"MGELLIPSE_SEGMENT"	// �Z�O�����g
	,"MGELLIPSE_CLOSED"};	// �����ȉ~
	int i=(int)rel;
	if(i>=3) i=2; if(i<0) i=0;
	out<<name[i];
	return out;
}

// MGSTRAIGHT_TYPE��W���o�͂ɏo�͂���B
ostream & operator << (ostream & out, MGSTRAIGHT_TYPE rel){
	const char* name[4]={
	 "MGSTRAIGHT_EMPTY" 		// ��
	,"MGSTRAIGHT_SEGMENT"		// ����
	,"MGSTRAIGHT_HALF_LIMIT"	// ������
	,"MGSTRAIGHT_UNLIMIT"};		// ��������
	 int i=(int)rel;
	 if(i>=4) i=3; if(i<0) i=0;
	out<<name[i];
	return out;
}

// MGSURFACE_TYPE��W���o�͂ɏo�͂���B
ostream & operator << (ostream & out, MGSURFACE_TYPE rel){
	const char* name[10]={
	"MGSURFACE_UNKNOWN",	// �s��
	"MGSURFACE_PLANE",		// ����
	"MGSURFACE_CONE",		// �~��
	"MGSURFACE_SPHERE",		// ����
	"MGSURFACE_TORUS",		// �~��
	"MGSURFACE_SPLINE"		//Free form surface ���R�Ȗ�                                 
							//(Tensor product surface of LBRep)
	"MGSURFACE_RSPLINE",	//Free form surface                             
							//(Tensor product surface of Rational B-Spline) 
	"MGSURFACE_CYLINDER",	//Cylinder surface(A special case of MGSURFACE_CONE)
	"MGSURFACE_USER1",		//Auxiliary surface type 1
	"MGSURFACE_USER2"		//Auxiliary surface type 2
	};
	int i=(int)rel;
	if(i>=10 || i<0) i=0;
	out<<name[i];
	return out;
}

// MGCCRELATION��W���o�͂ɏo�͂���B
ostream & operator << (ostream & out, MGCCRELATION rel){
	const char* name[5]={
	"MGCCREL_UNKNOWN",	// ���m�i���ׂĂ��Ȃ��āA�킩��Ȃ��A�ȉ��̂����ꂩ�j
	"MGCCREL_ISECT",	// �����i�����A�ڂ���A��v�ȊO�j
	"MGCCREL_NORMAL",	// ����
	"MGCCREL_TANGENT",	// �ڂ���
	"MGCCREL_COIN"};	// ��v����
	 int i=(int)rel;
	 if(i>=5) i=4; if(i<0) i=0;
	out<<name[i];
	return out;
}

// MGCSRELATION��W���o�͂ɏo�͂���B
ostream & operator << (ostream & out, MGCSRELATION rel){
	const char* name[6]={
	"MGCSREL_UNKNOWN",		// ���m
	"MGCSREL_IN",			// �Ȗʂ̓�������̌�_
	"MGCSREL_OUT",			// �Ȗʂ̊O������̌�_
	"MGCSREL_IN_TAN",		// �Ȗʂ̓�������ڂ��Ă���
	"MGCSREL_OUT_TAN",		// �Ȗʂ̊O������ڂ��Ă���
	"MGCSREL_COIN"};		// �Ȑ����ȖʂɊ܂܂�Ă���
	 int i=(int)rel;
	 if(i>=6) i=5; if(i<0) i=0;
	out<<name[i];
	return out;
}

// MGSSRELATION��W���o�͂ɏo�͂���B
ostream & operator << (ostream & out, MGSSRELATION rel){
	const char* name[4]={
	"MGSSREL_UNKNOWN",	// ���m�i���ׂĂ��Ȃ��āA�킩��Ȃ��A�ȉ��̂����ꂩ�j
	"MGSSREL_ISECT",	// �����i�ڂ���A��v�ȊO�j
	"MGSSREL_TANGENT",	// �ڂ���
	"MGSSREL_COIN"};	// ��v����
	 int i=(int)rel;
	 if(i>=4) i=3; if(i<0) i=0;
	out<<name[i];
	return out;
}

// MGENDCOND��W���o�͂ɏo�͂���B
ostream & operator << (ostream & out, MGENDCOND rel){
	const char* name[5]={
	"MGENDC_UNKNOWN",	// ���m
	"MGENDC_1D",		// 1st deravative provided.
	"MGENDC_2D",		// 2nd deravative provided.
	"MGENDC_NO",		// no end cond(only positional data)
	"MGENDC_12D"};		// both 1st and 2nd deravatives provided.
	 int i=(int)rel;
	 if(i>=5) i=4; if(i<0) i=0;
	out<<name[i];
	return out;
}
