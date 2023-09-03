/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Box.h"
#include "mg/Position.h"
#include "mg/BPointSeq.h"
#include "mg/SPointSeq.h"
#include "mg/OscuCircle.h"
#include "mg/KnotArray.h"
#include "mg/PPRep.h"
#include "mg/Straight.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/Plane.h"
#include "mg/Tolerance.h"

#include "cskernel/blurev.h"
#include "cskernel/Blunk.h"
#include "cskernel/Blumov.h"
#include "cskernel/Bludec.h"
#include "cskernel/Blucpr.h"
#include "cskernel/Blucon.h"
#include "cskernel/Blqbox.h"
#include "cskernel/blel.h"
#include "cskernel/blelin.h"

using namespace std;

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGLBRep.cpp
//
// Implement MGLBRep class.

#define local_basis_size 20
//This is the size of local(automatic) variable to store B-Spline
//basis functinon value of order k.
//When order exceeds local_basis_size, the area will be obtained by new.

//Member Function

// ���͂̃p�����[�^�͈͂̋Ȑ��������͂ރ{�b�N�X��Ԃ��B
MGBox MGLBRep::box_limitted(
	const MGInterval& intrvl// Parameter Range of the curve.
)const{
	if(intrvl.infinite()) return box_unlimit();
	else{
		double t1=intrvl.low_point(), t2=intrvl.high_point();
		if(MGREqual_base(t1,t2,knot_vector().param_span())){
			return MGBox(eval(t1),eval(t2));
		}
		const int k=order();
		const int n=bdim();
		const int irc=m_line_bcoef.capacity();
		const int ncd=sdim();
		MGInterval i2=param_range()&intrvl;
		double ts=i2.low_point(), te=i2.high_point();
		double* tw=new double[(n+k)+k*k+ncd*2];
		double* wk2=tw+(n+k);
		double* fbox=wk2+k*k;
		blqbox_(k,n,knot_data(),coef_data(),irc,ncd,ts,te,tw,wk2,fbox);
		MGBox rbox(ncd); MGInterval i1;
		for(int i=0; i<ncd; i++){
			i1.set_low_point(fbox[i*2]); i1.set_high_point(fbox[i*2+1]);
			rbox(i)=i1;
		}
		delete[] tw;
		return rbox;
	}
}

//Return box of whole of the curve.
void MGLBRep::compute_box(MGBox& bx) const{
	m_line_bcoef.compute_box(bx);
}

//Connect brep2 to this brep to make one B-Representation.
//This parameter range will not be changed, instead brep2's range
//will be so changed that brep2 has the same 1st derivative magnitude
//as the original this brep's at the connecting point
//(start or end point of this).
//continuity and which can be obtained using the fucntion continuity().
void MGLBRep::connect(
	int continuity,	//continuity.
	int which,		//which point of this to which of brep2.
				// =0: start of this and start of brep2.
				// =1: start of this and end of brep2.
				// =2: end of this and start of brep2.
				// =3: end of this and end of brep2.
	const MGLBRep& brep2	//B-Rep 2.
){
	assert(0<=which && which<=3);
	assert(continuity>=-1);
	invalidateBox();

	MGLBRep b2t;
	if(which<=1){//When connecting to the start of this.
		b2t=std::move(*this); *this=brep2;
		if(which==0) negate();
		double ts1=param_s(), te1=param_e();
		double a=eval(te1,1).len();
		double ts2=b2t.param_s();
		double b=b2t.eval(ts2,1).len();
		if(MGREqual_base(a,b,1.)){
			double new_ts=ts2-(te1-ts1);
			change_range(new_ts, ts2);
		}else if(b<=MGTolerance::mach_zero()){
			if(continuity>=1) continuity=0;
		}else{
			double new_ts=ts2-(te1-ts1)*a/b;
			change_range(new_ts, ts2);
			//This change_range() guarantees that b2t's(original this)
			//parameter range will not be modified.
		}
	}else{//When connecting to the end of this.
		b2t=brep2;
		if(which==3) b2t.negate();
	}
	//Now this end is connected to the start of the b2t.

	//Make order even for brep1 and 2.
	int k=order(), k2=b2t.order();
	if(k<k2){
		k=k2; change_order(k);
	}else if(k>k2){
		b2t.change_order(k);
	}
	int km2=k-2;
	if(continuity>km2)
		continuity=km2;

	//Make space dimension even for brep1 and 2.
	int dim=sdim(), dim2=b2t.sdim();
	if(dim<dim2){
		dim=dim2;
		m_line_bcoef=MGBPointSeq(dim,m_line_bcoef);
	}else if(dim>dim2){
		b2t.m_line_bcoef=MGBPointSeq(dim,b2t.m_line_bcoef);
	}

	//This will be the whole curve container, and its knot vector and
	//bcoef area should be large enough to store the both curves.
	int n1=bdim(), n2=b2t.bdim();
	m_knot_vector.reshape(k+n1+n2);
	m_line_bcoef.reshape(n1+n2);
	int irc1=m_line_bcoef.capacity(), irc2=b2t.m_line_bcoef.capacity();

	int kk=k+2;
	if(kk<dim)
		kk=dim;
	double* work=new double[k*k+k*kk];

	double ratio;
	int it2s;
	blucon_(k, &n1, &knot(0), &coef(0,0), irc1, dim,
			n2, b2t.knot_data(), b2t.coef_data(), irc2,
			continuity, work, work+k*k, &ratio, &it2s);
	delete[] work;
	m_knot_vector.set_bdim(n1);
	m_line_bcoef.set_length(n1);
}

//Exchange ordering of the coordinates.
//Exchange coordinates (j1) and (j2).
void MGLBRep::coordinate_exchange(int j1, int j2){
	assert(j1<sdim() && j2<sdim());
	m_line_bcoef.coordinate_exchange(j1,j2);
	invalidateBox();
}

//Construct new curve object by copying to newed area.
//User must delete this copied object by "delete".
MGLBRep* MGLBRep::clone() const{return new MGLBRep(*this);}

void convert_to_Bezier_knot(
	MGKnotVector& t
){
	int k=t.order(), n=t.bdim();
	int km1=k-1;
	for(int i=0; i<k; i++)
		t(i)=0.;
	int idSpan=1;
	for(int i=k; i<n;){
		for(int j=0; j<km1; j++)
			t(i++)=double(idSpan);
		idSpan++;
	}
	for(int i=n; i<n+k; i++)
		t(i)=double(idSpan);
}

///Convert this curve to Bezier curve.
///If this is MGLBRep or MGStraight, the shape is exactly the same
///as the original. Otherwise, this is apporoximated by MGLBRep.
///The result MGLBRep is of order 2 if original order is 2 
///and is of order4 otherwise.
///Output bezier can be this.
void MGLBRep::convert_to_Bezier(MGLBRep& bezier)const{
	int k=order();
	if(k==2){
		bezier.m_knot_vector=m_knot_vector;
		convert_to_Bezier_knot(bezier.m_knot_vector);
	}else if(k==4){
		int km1=k-1;
		MGKnotArray knots;
		int n=bdim();
		const MGKnotVector& told=knot_vector();
		knots.push_back(MGKnot(told(km1),k));//Starting point parameter.
		for(int i=k; i<n;){
			knots.push_back(MGKnot(told(i++),km1));
			while(i<n && (told(i)==told(i+1)))
				i++;
		}
		knots.push_back(MGKnot(told(n),k));//Ending point parameter.

		MGLBRep lbtemp(*this);
		lbtemp.addKnots(knots);

		MGBPointSeq& bpnew=bezier.line_bcoef();
		bpnew=std::move(lbtemp.line_bcoef());
		int nnew=bpnew.length();
		MGKnotVector& tnew=bezier.knot_vector();
		tnew.size_change(k,nnew);
		convert_to_Bezier_knot(bezier.m_knot_vector);
	}else{
		MGLBRep lb2;
		approximate_as_LBRep(lb2,4);
		lb2.convert_to_Bezier(bezier);
	}
	bezier.copy_appearance(*this);
}

//Construct new curve object by changing
//the original object's space dimension.
//User must delete this copied object by "delete".
MGLBRep* MGLBRep::copy_change_dimension(
	int sdim,		// new space dimension
	int start1, 		// Destination order of new line.
	int start2) 		// Source order of this line.
	const{
	return new MGLBRep(sdim,*this,start1,start2);
}

//Construct new curve object by copying to newed area,
//and limitting the parameter range to prange.
//Returned is newed object and must be deleted.
MGCurve* MGLBRep::copy_limitted(const MGInterval& prange) const{
	MGInterval pr=param_range()&prange;
	MGLBRep* lb=new MGLBRep;
	shrinkToParameters(pr.low_point(), pr.high_point(), *lb);
	return lb;
}

//Compute curvilinear integral of the 1st two coordinates.
//(���ϕ��j�����߂�B
//This integral can be used to compute area sorrounded by the curve.
double MGLBRep::curvilinear_integral(double t1, double t2) const{
	if(sdim()<2) return 0.;
	return MGCurve::curvilinear_integral(t1,t2);
}

//Divide this curve at the designated knot multiplicity point.
//Function's return value is the number of the curves after divided.
int MGLBRep::divide_multi(
	std::vector<UniqueCurve>& crv_list,	//divided curves are appended.
	int multiplicity	//designates the multiplicity of the knot to divide at.
						//When multiplicity<=0, order()-1 is assumed.
						//When multiplicity>=order(), order() is assumed.
)const{
	int k=order();assert(k>=2);
	if(k<=2){
		int sd=sdim();
		const MGKnotVector& t=knot_vector();
		const MGBPointSeq& bp=line_bcoef();
		int n=bdim();
		for(int i=1; i<n; i++){
			MGBPointSeq bpi(2,sd);
			bpi.store_at(0,bp(i-1));
			bpi.store_at(1,bp(i));
			double ti[4];
			ti[0]=ti[1]=t[i];
			ti[2]=ti[3]=t[i+1];
			MGLBRep* lb=new MGLBRep;
			lb->buildLBRepFromMemberData(MGKnotVector(2, 2, ti),std::move(bpi));
			crv_list.emplace_back(lb);
		}
		return n-1;
	}

	int start_index=k-1;
	int bd=bdim();
	if(multiplicity<=0)
		multiplicity=k-1;
	else if(multiplicity>int(k))
		multiplicity=k;

	int nold=(int)crv_list.size();
	const MGKnotVector& t=knot_vector();
	int index, multi;
	do{
		multi = t.locate_multi(start_index, multiplicity, index);
		MGCurve* curvei = part(t[start_index],t[index]);
		crv_list.emplace_back(curvei);
		start_index = index + multi - 1;
	}while(index!=bd);	//���d�x��������Ȃ�������I���

	return (int)(crv_list.size()-nold);
}

MGVector MGLBRep::eval(		//Evaluate right continuous n'th derivative
	double x,		//Parameter value to evaluate.
	int nderiv,	//degree of derivative.
					//When 0, compute positional data.
	int leftcon		//Left continuous(leftcon=true)
					//or right continuous(leftcon=false).
)const{
	const int k=order();
	double c[local_basis_size]; double* cp=c;
	if(k>local_basis_size) cp=new double[k];
	//This is done to save "new" when k<=local_basis_size.
	int id=m_knot_vector.eval_coef(x,cp,nderiv,leftcon);
	int i, j; int ii;
	const int ncd=sdim();
	MGVector v(ncd); double data;
	for(j=0; j<ncd; j++){
		data=0.;
		for(i=0; i<k; i++){
			ii=id+i;
			data=data+coef(ii,j)*cp[i];
		}
		v.set(j)=data;
	}
	if(k>local_basis_size) delete[] cp;
	return v;
}

//Compute position, 1st and 2nd derivatives.
// �p�����[�^�l��^���Ĉʒu�A�ꎟ�����l�A�񎟔����l�����Ƃ߂�B
void MGLBRep::eval_all(
	double tau,			//Input parameter value(�p�����[�^�l)
	MGPosition& P,		//Position(�ʒu)
	MGVector& V1,		//1st derivative(1�������l)
	MGVector& V2		//2nd derivative(2�������l)
)const{
	int sd=sdim();
	double* data=new double[sd*3];
	eval_all(tau,2,data);
	P=MGPosition(sd,data);
	V1=MGVector(sd,data+sd);
	V2=MGVector(sd,data+2*sd);
	delete[] data;
}

//Evaluate all of i'th derivative data for 0<=i<=nderiv.
//Output will be put on deriv[j+i*sdim()]
//for 0<=i<=nderiv and 0<=j<sdim(), i.e. 
//deriv[j+i*sdim()] is i-th derivative data for 0<=j<sdim(). 
void MGLBRep::eval_all(
	double tau,		// Parameter value to evaluate.
	int nderiv,	// Order of Derivative.
	double* deriv,	// Output area of size (nderiv+1)*sdim().
	int leftcon		//Left continuous(leftcon=true)
					//or right continuous(leftcon=false).
)const{
	int i,j,m;
	int k=order(); int km1=k-1;
	int sd=sdim();
	int nd=nderiv;
	if(nderiv>km1) nd=km1;
    double	biatx_buf[local_basis_size],
			deltal_buf[local_basis_size],
			deltar_buf[local_basis_size];
	double* deltal=deltal_buf; double* deltar=deltar_buf;
	double* biatx=biatx_buf;
	if(k>local_basis_size){
		biatx=new double[k+k+k]; deltal=biatx+k; deltar=deltal+k;
	}
	const double* t=knot_data();
	int left=knot_vector().locate(tau,leftcon);
	MGSPointSeq alpha(k,nd+1,sd);	//Work area.
// STORE THE K B-SPLINE COEFF.S RELEVANT TO CURRENT KNOT INTERVAL
// IN alpha(.,0,.) .
	int leftmkp1=left-k+1;
	for(i=0; i<k; ++i)
	    for(m=0; m<sd; ++m) alpha(i,0,m) = coef(leftmkp1+i,m);

// FOR J=1,...,nd, COMPUTE THE  K-J  B-SPLINE COEFF.S RELEVANT TO
// CURRENT KNOT INTERVAL FOR THE J-TH DERIVATIVE BY DIFFERENCING
// THOSE FOR THE (J-1)ST DERIVATIVE, AND STORE IN alpha(.,J,.).
	double diff; int kmj;
	for (j=1; j<=nd; ++j) {
	    kmj = k - j;
		double fkmj=double(kmj);
	    for (i=1; i<=kmj; ++i) {
			diff = t[left+i] - t[left+i-kmj];
			for(m=0; m<sd; ++m)
				alpha(i-1,j,m)=(alpha(i,j-1,m)-alpha(i-1,j-1,m))/diff*fkmj;
	    }
	}

	//Actually deriv's array is deriv[nd][sd], deriv[i] will contain
	//i-th derivative for 0<=i<=nderiv.
	if(nderiv>km1){
		for(j=k; j<=nderiv; ++j)
			for(m=0; m<sd; ++m) deriv[m+j*sd]=0.;
	}

//     FOR  J = 0, ..., nd, FIND THE VALUES AT  T(LEFT)  OF THE  J+1
//     B-SPLINES OF ORDER  J+1  WHOSE SUPPORT CONTAINS THE CURRENT
//     KNOT INTERVAL FROM THOSE OF ORDER  J  (IN  BIATX ), THEN COMB-
//    INE WITH THE B-SPLINE COEFF.S (IN alpha(.,nd-J,.) ) FOUND EARLIER
//     TO COMPUTE THE (K-J-1)ST DERIVATIVE AT  T(LEFT)  OF THE GIVEN
//     SPLINE.
	biatx[0] = 1.;
	double saved,sum,term;
	for(m=0; m<sd; ++m) deriv[m+nd*sd] = alpha(0,nd,m);
	for (j=1; j<k; ++j) {
		int jm1=j-1;
	    deltar[jm1] = t[left+j] - tau;
	    deltal[jm1] = tau - t[left-jm1];
	    saved = 0.;
	    for(i=0; i< j; ++i) {
			term = biatx[i] /(deltar[i]+deltal[jm1-i]);
			biatx[i] = saved+deltar[i] * term;
			saved = deltal[jm1-i] * term;
	    }
	    biatx[j] = saved;
		int kmjm1=k-j-1;
		if(kmjm1<=nd){
			for(m=0; m<sd; ++m){
				sum=0.;
				for (i=0; i<=j; ++i) sum = biatx[i]*alpha(i,kmjm1,m)+sum;
				deriv[m+kmjm1*sd] = sum;
			}
		}
	}

	if(k>local_basis_size)
		delete[] biatx;
}
									
void MGLBRep::eval_line(			//Evaluate data for data point seq.(BLELIN)
			MGENDCOND beginc,		//Begin end condition
			MGENDCOND endc,			//End end conditoion 
			const MGNDDArray& tau,	//Data points.
			MGBPointSeq& value		//Values evaluated.
) const{
	const int k=order();
	const int n=bdim();
	const int irc=m_line_bcoef.capacity();
	const int ncd=sdim();
	const int ntau=tau.length();
	int iv2=value.capacity();
	if(iv2<ntau || value.sdim()!=ncd){
		value=MGBPointSeq(ntau,ncd);
		iv2=ntau;
	}
	const int iv1=1;

	assert( sdim()<=3);
	assert(beginc != MGENDCOND::MGENDC_12D &&	beginc != MGENDCOND::MGENDC_UNKNOWN);
	assert(endc != MGENDCOND::MGENDC_12D && endc != MGENDCOND::MGENDC_UNKNOWN);
	int bc = static_cast<int>(beginc), ec = static_cast<int>(endc);
	blelin_(k,n,knot_data(),coef_data(),irc,ncd,
			bc, ec, ntau, tau.data(), iv1, iv2, &value(0,0));
	value.set_length(ntau);
}

//Extrapolate the curve by the chord length.
void MGLBRep::extend(
		int start,			//Flag of start or end poit of the line.
							//If start is true extend on the start point.
		double length,		//chord length to extend. 
		double dk          //Coefficient of how curvature should vary at
//    extrapolation start point. When dk=0, curvature keeps same, i.e.
//    dK/dS=0. When dk=1, curvature becomes zero at length extrapolated point,
//    i.e. dK/dS=-K/length at extrapolation start point.
//    (S=parameter of arc length, K=Curvature at start point)
//    That is, when dk reaches to 1 from 0, curve changes to flat.
//The extrapolation is C2 continuous if the order >=4.
){
	assert(sdim()<=3);
	double t=start ? param_s():param_e();
	MGVector V1=eval(t,1);
	double tan_ratio=V1.len();
	length/=tan_ratio;
	if(start)
		length*=-1.;
	double tau=t+length;
	MGPPRep extention;
	extrapolated_pp(tau,dk,extention);

// ***** NOW PP-REP OBTAINED IN extention, CONVERT PP-REP TO B-REP and connect to this.
	MGLBRep extentionLB(extention);
	if(start)
		connect(2,0,extentionLB);
	else
		connect(2,2,extentionLB);
}

//Extrapolate this curve by an (approximate) chord length.
//The extrapolation is C2 continuous.
void MGLBRep::extend(
	double length,	//approximate chord length to extend. 
	bool start		//Flag of which point to extend, start or end point of the line.
					//If start is true extend on the start point.
){
	int se=start ? 1:0;
	extend(se,length,0.);
}

//Extrapolate the curve by the parameter value.
void MGLBRep::extend_with_parameter(
	double tau,		//The parameter value at the end of extended point.
					//When tau<param_s(), extension will be done at the starting point.
					//When tau>param_e(), extension will be done at the end point.
	double dk     //Coefficient of how curvature should vary at the connecting point.
){
	MGPPRep pp;
	extrapolated_pp(tau,dk,pp);

	MGLBRep lbext(pp);
	if(tau<param_s()) connect(2,0,lbext);
	else connect(2,2,lbext);
}

//Extracts control points.
//Fucntion's return value is 
//true if control points was obtained, false if not.
bool MGLBRep::get_control_points(
	MGBPointSeq& cpoints	//Control points will be output.
)const{
	cpoints=line_bcoef();
	return true;
}

///Test if this is a Bezier Curve.
///Functions's return value is MGLBRep* if Bezier, null if not.
///If input ordr>=2, order is also tested if this Bezier's order is the same as input order.
///If input ordr<=1, any ordr>=2 is allowed for Bezier curve.
///Bezier curve is defined as follows. Here t=knot_vector(), k is this LBRep's order,
///n=bdim(), and m=(n-k)/(k-1).
///(1) n=k+(k-1)*m.
///(2) t(0)=t(1)=,...,=t(k-1)=0
///(3) t(i)=t(i+1)=,...,=t(i+k-2)=j+1
///         for i=k, k+(k-1),...,k+j*(k-1) and j=0,...,m-1.
///(4) t(n)=t(n+1)=,...,=t(n+k-1)=m+1
const MGLBRep* MGLBRep::is_Bezier(int ordr)const{
	int k=order();
	if(k<=1)
		return 0;
	if(ordr>=2){
		if(ordr!=k)
			return 0;
	}

	int n=bdim(), km1=k-1;
	int m=(n-k)/km1;
	if(n!=(k+km1*m))
		return 0;
	
	const MGKnotVector& t=knot_vector();
	double ts=0.;//Starting parameter.
	int i=0;//Index of t.
	for(; i<k; i++){//Starting parameter multiplicity.
		double ti=t(i);
		if(!MGREqual(ts,ti))
			return 0;
	}

	for(int j=0; j<m; j++){
		double ti=double(j+1);
		for(int p=0; p<km1; p++)//Middle knot parameter multiplicity.
			if(!MGREqual(ti,t(i++)))
				return 0;
	}

	double te=double(m+1);//Ending parameter.
	for(int p=0; p<k; p++){//Ending parameter multiplicity.
		double ti=t(i++);
		if(!MGREqual(te,ti))
			return 0;
	}

	return this;
}

//Test if this cure is co-planar with the 2nd curve curve2.
//MGPlane expression will be out to plane if this is co-planar.
//Function's return value is true if co-planar.
bool MGLBRep::is_coplanar(
	const MGCurve& curve2, MGPlane& plane
)const{
	const MGStraight* sl11=dynamic_cast<const MGStraight*>(&curve2);
	if(sl11)
		return sl11->is_coplanar(*this,plane);

	MGPosition point;
	MGStraight sl1;
	int plkind=planar(plane,sl1,point);
	if(plkind==0)
		return false;

	MGPosition point2;
	MGStraight sl2;
	int plkind2=-1;
	MGPlane plane2;

	const MGLBRep* lb2=dynamic_cast<const MGLBRep*>(&curve2);
	if(lb2)
		plkind2=lb2->planar(plane2,sl2,point2);
	else{
		const MGRLBRep* rlb=dynamic_cast<const MGRLBRep*>(&curve2);
		if(rlb)
			plkind2=rlb->planar(plane2,sl2,point2);
	}
	if(plkind2==0) return false;

	MGPosition uv;
	if(plkind2==1){
		if(plkind==1){
			plane=MGPlane(point,point2,point);
			return true;
		}else if(plkind==2){
			plane=MGPlane(sl1,point2);
			return true;
		}else{
			return plane.on(point2,uv);
		}
	}else if(plkind2==2){
		if(plkind==1){
			plane=MGPlane(sl2,point);
			return true;
		}else if(plkind==2){
			return sl1.is_coplanar(sl2,plane);
		}else{
			return plane.on(sl2);
		}
	}else if(plkind2==3){
		if(plkind==1){
			plane=plane2;
			return plane.on(point,uv);
		}else if(plkind==2){
			plane=plane2;
			return plane.on(sl1);
		}else{
			return plane==plane2;
		}
	}

	//When curve2 is neither MGLBRep nor MGRLBRep.
	if(!curve2.is_planar(plane2)) return false;
	if(plkind==1){
		plane=plane2;
		return plane.on(point,uv);
	}else if(plkind==2){
		plane=plane2;
		return plane.on(sl1);
	}
    return plane2==plane;
}

//Test if this cure is planar or not.
//MGPlane expression will be out to plane if this is planar.
//Function's return value is true if planar.
bool MGLBRep::is_planar(MGPlane& plane)const{
	MGStraight line;
	MGPosition point;
	int isp=planar(plane,line,point);
	if(isp==1){//IF this is within a point.
		plane=MGPlane(mgZ_UVEC,point);
	}else if(isp==2){//IF this is within a straight.
		return line.is_planar(plane);
	}

	return isp>0;
}

// ���g�Ɏw�肵���p�����[�^�͈͂̂���������������B
void MGLBRep::limit(const MGInterval& i1){
	MGInterval i2=param_range();
	MGInterval i3=i2 & i1;
	if(i3 != i2)
		shrinkToParameters(i3.low_point(), i3.high_point(), *this);
}

//Modify the original line by moving move_point to to_point. fix_point can be
//applied according to move_kind.
void MGLBRep::move(
		int move_kind,				//Indicates how to move line.
		double move_point_param,	//indicate object point to move by the
									//parameter value.
		const MGPosition& to_point,	//destination point of the abve source
									//point.
		const double fix_point[2])
{
	if(move_kind<1 || move_kind>4)
		move_kind=4;

	const int k=order();
	const int n=bdim();
	const int irc1=m_line_bcoef.capacity();
	const int ncd=sdim();
	const int irc2=n;
	double p[3]; p[0]=to_point(0); p[1]=to_point(1); p[2]=to_point(2);
	MGBPointSeq bc(n,ncd);

	double tfo[2];
	double* work=new double[n];
	blumov_( k,n,knot_data(),coef_data(),irc1,ncd,
			move_point_param, p, move_kind, fix_point,
			irc2, work, &bc(0,0), tfo);
	delete[] work;
	bc.set_length(n);

	m_line_bcoef=std::move(bc);
	invalidateBox();
}

void MGLBRep::negate()	//BLUREV
							//Change direction of the line.
{
	const int k=order();
	int n=bdim();
	const int irc=m_line_bcoef.capacity();
	const int ncd=sdim();
	blurev_(k, n, knot_data(), coef_data(), irc, ncd,
			irc, &n, &knot(0), &coef(0,0));
}

//Obtain parameter value if this curve is negated by "negate()".
double MGLBRep::negate_param(double t)const{
	double tspte=param_s()+param_e();
	return tspte-t;
}

//Changing this object's space dimension.
void MGLBRep::change_dimension(
	int sdim,		// new space dimension
	int start1, 		// Destination order of new object.
	int start2) 		// Source order of this object.
{
	m_line_bcoef=MGBPointSeq(sdim,m_line_bcoef,start1,start2);
	invalidateBox();
}

//Change order of the B-Rep. When new order is greater than the original,
//new B-rep is guaranteed to be the same line as the original. However,
//if new order is less than the original one, new line is not the same
//in general.
void MGLBRep::change_order(
	int knew)		//New order number. 
{
	int kold=order(); int kdif=knew-kold;
	if(kdif==0)
		return;

	int nold=bdim();
//1.Generate new knot vector.
	int i,m,n; int j,mold,mnew;
	int nspan=nold-kold+1; int n1=nspan*kdif+nold;
	int n2=nspan+knew-1; if(n1<n2) n1=n2;
	MGKnotVector t(knew,n1);
	//1.1 First knew knots.
	for(i=1; i<=knew; i++){
		j=kold-i; if(j<0) j=0;
		t(knew-i)=knot(j);
	}
	//1.2 Internal break point.
	n=knew; m=kold;
	while(m<nold){
		//Count old knots' multiplicity. (mold=multiplicity)
		mold=1; while(knot(m)==knot(m+1)){mold++; m++;}
		mnew=knew-kold+mold; if(mnew<=0) mnew=1;
		//mnew is the new knots' multiplicity.
		while(mnew--) t(n++)=knot(m);
		m+=1;
	}
	//1.3 Last knew knots.
	for(i=0; i<knew; i++){
		m=nold+i; if(m>=nold+kold) m=nold+kold-1;
		t(n++)=knot(m);
	}
	t.set_bdim(n-knew);

//2. Convert to PP Rep, then change order.
	MGPPRep pp(knew,MGPPRep(*this));

//3. Construct new B-Rep.
	buildLBRepFromPPRep(std::move(t), pp);
}

//Change order of the B-Rep by approximation.
void MGLBRep::change_order_by_approximation(
	int ordr		//New order number. 
){
	assert(ordr>=order());
	int k=order(), n=bdim(), new_n;
	if(ordr<=order())
		return;
	int dif=ordr-k;
	new_n=n-k+ordr;
	MGLBRep temp;
	MGKnotVector& new_t=temp.knot_vector();
	new_t=MGKnotVector(ordr, new_n);
	int i,j;
	for(i=0; i<k; i++) new_t(i)=m_knot_vector[i];
	for(j=0; j<dif; j++) new_t(i++)=m_knot_vector[k-1];
	for(j=k; j<n; j++) new_t(i++)=m_knot_vector[j];
	for(j=0; j<dif; j++) new_t(i++)=m_knot_vector[n];
	for(j=n; j<n+k; j++) new_t(i++)=m_knot_vector[j];

	MGNDDArray tau; tau.buildByKnotVector(new_t);
	MGBPointSeq bp(new_n,sdim());
	eval_line(MGENDCOND::MGENDC_NO, MGENDCOND::MGENDC_NO,tau,bp);

	temp.buildByInterpolationWithKTV(tau, bp);
	buildLBRepFromMemberData(std::move(temp.m_knot_vector), std::move(m_line_bcoef));
	//This is to avoid destroying the other data of this(e.g. appearance).
}

//Change parameter range, be able to change the direction by providing
//t1 greater than t2.
void MGLBRep::change_range(	//BLUCPR
		double t1,			//Parameter value for the start of original. 
		double t2)			//Parameter value for the end of original. 
{
	const int k=order();
	int n=bdim();
	const int irc=m_line_bcoef.capacity();
	const int ncd=sdim();
	blucpr_(t1,t2,k,n, knot_data(), coef_data(),irc,ncd,
			irc, &n, &knot(0), &coef(0,0));
	m_knot_vector.set_bdim(n);
	m_line_bcoef.set_length(n);
	invalidateBox();
}

// Return ending parameter value.
double MGLBRep::param_e() const{
	return m_knot_vector.param_e();
}

//Normalize parameter value t to the nearest knot if their distance is
//within tolerance.
double MGLBRep::param_normalize(double t) const{
	return m_knot_vector.param_normalize(t);
}

// Return starting parameter value.
double MGLBRep::param_s() const{
	return m_knot_vector.param_s();
}
	
//Compute part of this curve from parameter t1 to t2.
//Returned is the pointer to newed object, and so should be deleted
//by calling program, or memory leaked.
MGCurve* MGLBRep::part(double t1, double t2,int multiple) const{
	double ts=param_s();
	if(t1<ts)
		t1=ts;
	if(t2-t1>param_error()){
		auto lb=new MGLBRep;
		shrinkToParameters(t1, t2, *lb, multiple);
		return lb;
	}else{
		MGPosition P1=eval(t1), P2=eval(t2);
		return new MGStraight(t2,t1,P2-P1,P1);;
	}
}

//Change an original B-Rep to new one with subdivided knot configuration.
//Knots t must be refined knots.
void MGLBRep::refine(
	const MGKnotVector& t	//refined knot vector.
){
	assert(t.order()==order());

	const int k=order();
	const int n1=bdim();
	const int ncd=sdim();
	MGKnotVector knot1(m_knot_vector);//save the original
	MGBPointSeq bcoef1(m_line_bcoef);//save the original
	const int irc=bcoef1.capacity();

	m_knot_vector=t;
	int n2=t.bdim();
	m_line_bcoef.reshape(n2);
	int irc2 = m_line_bcoef.capacity();

	double* work=new double[k*k];
	blunk_(k,n1,knot1.data(),bcoef1.data(),irc,ncd,n2,t.data(), irc2,work,&coef(0,0));
	delete[] work;
}

//Change the  B-Rep by decreasing B-Rep dimension.	This is an approximation
//of the origimal B-Rep.
int MGLBRep::reduce(			//BLUDEC
			int ndec)			//Number of B-rep dimension to decrease 
{
	assert(bdim()-ndec >= order() && ndec>0);

	int error;
	const int ism=1;
	const int k=order();
	const int n1=bdim();
	const MGKnotVector& oldt=knot_vector();
	MGBPointSeq oldbcoef(line_bcoef());
	const int irc1=oldbcoef.capacity();
	const int ncd=sdim();
	const int irc2=m_line_bcoef.capacity();
	double* work2=new double[n1+n1*(2*k-1)];
	double* work1=work2+n1;
	int n2;

	bludec_(ism,k,n1, oldt.data(), oldbcoef.data(),
			irc1, ncd, ndec, irc2, work1, work2,
			&n2, &m_knot_vector(0), &m_line_bcoef(0,0), &error);
							   
	delete[] work2;
	if(error==1){
		error=0;		//Return of BLUDEC:error=1 means normal.
		m_knot_vector.set_bdim(n2);
		m_line_bcoef.set_length(n2);
	}else{
		m_knot_vector=oldt;
		m_line_bcoef=oldbcoef;
	}
	invalidateBox();
	return error;
}

// Compute box of unlimitted.
const MGBox& MGLBRep::box_unlimit() const{
	return box();
}

//Operator overload
MGLBRep& MGLBRep::operator=(const MGGel& gel2){
	const MGLBRep* gel2_is_this=dynamic_cast<const MGLBRep*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}
MGLBRep& MGLBRep::operator=(MGGel&& gel2){
	MGLBRep* gel2_is_this=dynamic_cast<MGLBRep*>(&gel2);
	if(gel2_is_this)
		operator=(std::move(*gel2_is_this));
	return *this;
}


// �Ȑ��̕��s�ړ����s���Ȑ��𐶐�����B
MGLBRep MGLBRep::operator+ (const MGVector& vec) const{
	MGLBRep brep(*this);
	brep += vec;
	return brep;
}
MGLBRep operator+ (const MGVector& v, const MGLBRep& lb){
	return lb+v;
}

// �^�x�N�g�������Ȑ��𕽍s�ړ����Ď��g�Ƃ���B
MGLBRep& MGLBRep::operator+= (const MGVector& vec){
	m_line_bcoef += vec;
	m_box+=vec;
	return *this;
}

// �Ȑ��̋t�����ɕ��s�ړ����s���Ȑ��𐶐�����B
MGLBRep MGLBRep::operator- (const MGVector& vec) const{
	MGLBRep brep(*this);
	brep-= vec;
	return brep;
}

// �^�x�N�g�������Ȑ����}�C�i�X�����ɕ��s�ړ����Ď��g�Ƃ���B
MGLBRep& MGLBRep::operator-= (const MGVector& vec){
	m_line_bcoef -= vec;
	m_box-=vec;
	return *this;
}

// �^����ꂽ�X�P�[���������I�u�W�F�N�g�𐶐�����B
//generate line by scaling.
MGLBRep MGLBRep::operator* (double scale) const{
	MGLBRep lb(*this);
	lb*=scale;
	return lb;
}

// �^����ꂽ�X�P�[���������I�u�W�F�N�g�𐶐�����B
//generate line by scaling.
MGLBRep operator* (double scale, const MGLBRep& lb){
	return lb*scale;
}

// ���g�̋Ȑ��ɗ^����ꂽ�X�P�[����������B
//Scale the curve.
MGLBRep& MGLBRep::operator*= (double scale){
	m_line_bcoef *= scale;
	m_knot_vector *= scale;
	invalidateBox();
	return *this;
}

// �^����ꂽ�ϊ��ŋȐ��̕ϊ����s���Ȑ��𐶐�����B
MGLBRep MGLBRep::operator* (const MGMatrix& mat) const{
	MGLBRep brep(*this);
	brep *= mat;
	return brep;
}

// �^����ꂽ�ϊ��ŋȐ��̕ϊ����s�����g�̋Ȑ��Ƃ���B
MGLBRep& MGLBRep::operator*= (const MGMatrix&  mat){
	m_line_bcoef *= mat;
	invalidateBox();
	return *this;
}

// �^����ꂽ�ϊ��ŋȐ��̃g�����X�t�H�[�����s���Ȑ��𐶐�����B
MGLBRep MGLBRep::operator* (const MGTransf& tr) const{
	MGLBRep brep(*this);
	brep *= tr;
	return brep;
}

// �^����ꂽ�ϊ��ŋȐ��̃g�����X�t�H�[�����s�����g�Ƃ���B
MGLBRep& MGLBRep::operator*= (const MGTransf& tr){
	m_line_bcoef *= tr;
	invalidateBox();
	return *this;
}

bool MGLBRep::operator==(const MGRLBRep& gel2)const{
	return gel2.MGRLBRep::operator==(*this);
}
bool MGLBRep::operator==(const MGLBRep& gel2)const {
	return m_knot_vector == gel2.m_knot_vector &&
		m_line_bcoef == gel2.m_line_bcoef;
}

std::partial_ordering MGLBRep::operator<=>(const MGLBRep& gel2)const{
	int n1=bdim(), n2=gel2.bdim();
	if(n1==n2){
		int sd1=sdim(), sd2=gel2.sdim();
		if(sd1==sd2){
			MGPosition v1(sd1), v2(sd1);
			m_line_bcoef.point(0,0,sd1,v1);
			gel2.m_line_bcoef.point(0,0,sd1,v2);
			return v1.len()<=>v2.len();
		}else
			return sd1<=>sd2;
	}else
		return n1<=>n2;
}
bool MGLBRep::equal_test(const MGGel& g2)const {
	auto c = typeCompare(g2);
	return c == 0 ? *this == dynamic_cast<const MGLBRep&>(g2) : false;
}
std::partial_ordering MGLBRep::ordering_test(const MGGel& g2)const {
	auto c = typeCompare(g2);
	return c == 0 ? *this <=> dynamic_cast<const MGLBRep&>(g2) : c;
}
