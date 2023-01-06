/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/NDDArray.h"
#include "mg/KnotVector.h"
#include "mg/KnotArray.h"
#include "mg/Tolerance.h"

#include "cskernel/B1nk.h"
#include "cskernel/bkcrng.h"
#include "cskernel/bkdtkt.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGKnotVector
// MGKnotVector Implementation.

// Private member function.
void MGKnotVector::new_knot_coef(
	const MGKnotVector& old_knot, // old knot vector.
	int j_of_t,	// Index of new Knot vector to get coef.
	double* coef	// coef's to mutiply are returned.
){
	int mu= old_knot.m_current+1;
	b1nk_(m_order, old_knot.data(), mu, data(), j_of_t, coef);
}

//Constructor

//
MGKnotVector::MGKnotVector(
	int k,		// order
	int brdim,	// B-Rep dimension
	const double* data	// all of the knot data sequence if data!=NULL.
	// Construct garbage knot vector of spedified size.
):MGNDDArray(k+brdim,data),m_order(k){
	if(k)
		m_current=k-1;
	else
		m_current=0;
}

// Construct uniform knot vector with inital value init_value
// and unit incremental.
MGKnotVector::MGKnotVector(
	int k,			// order
	int brdim,		// B-Rep dimension.
	double init_value,		// Initial value
	double increment		//Incremental.
):MGNDDArray(k+brdim),m_order(k){
	assert (k>0 && brdim>=k && increment>0.);

	int km1=k-1;
	double end_value=init_value+double(brdim-km1)*increment;
	int i;
	for(i=0; i<k; i++){
		(*this)(i)=init_value;
		(*this)(brdim+i)=end_value;
	}
	for (i=k; i<brdim; i++){
		(*this)(i)=init_value+double(i-km1)*increment;
	}
	m_current=km1;
}

MGKnotVector::MGKnotVector(const MGNDDArray& dtp, int k)
// From Data Point, obtain knot vector.
	:MGNDDArray(dtp.length()+k),m_order(k){
	assert (k>0 && dtp.length()>=k);

	int brdim=dtp.length();
	bkdtkt_(dtp.data(), brdim, m_order, data());
	m_current=k-1;
}

//Add knots into original knot vector.
MGKnotVector::MGKnotVector(
	const MGKnotVector& vec2,	//Original Knot Vector
	const MGKnotArray& knots	//Knots to add with multiplicity.
):MGNDDArray(vec2),m_order(vec2.m_order){
	const int nad=(int)knots.size();
	int mlt_total=0;
	int i;
	for(i=0; i<nad; i++)
		mlt_total+=knots[i].multiplicity();

	int n2 = vec2.length();
	int nNew = n2+mlt_total;
	reshape(nNew);
	set_length(n2);

	double tau;
	double ts=param_s(); double te=param_e();
	mlt_total=0;
	for(i=0; i<nad; i++){
		tau=knots[i];
		if(ts<tau && tau<te){
			mlt_total += add_data(knots[i], m_order);
		}
	}
	m_current=m_order-1;
}

// Construct by extracting sub interval of vec2.
MGKnotVector::MGKnotVector(
	int start_id,			//Start id of vec2(from 0).
	int num,					//new B-Rep dimension.
	const MGKnotVector& vec2	//Original knot vector.
):MGNDDArray(start_id,vec2.order()+num,vec2),m_order(vec2.m_order){
	assert((start_id+num)<=vec2.bdim() && num>=vec2.order());
	m_current=m_order-1;
}

//Construct new order knot vector. 
MGKnotVector::MGKnotVector(
	const MGKnotVector& knotv,	// input knot vector
	int k				// new order
):MGNDDArray(0),m_order(k){
	assert(k>=1);
	int kold=knotv.order(), nold=knotv.bdim();
	if(kold == k){*this = knotv; return;}

	resize(nold-kold+2*k);
		//nold-kold+2*k is the maximum length of the knot vector.
	int i;
	double sparam=knotv[kold-1];
	for(i=0; i<k; i++) (*this)(i)=sparam;
		//Starting knot has multiplicity of order.
	int mult=1;
	double prev=sparam;
	int n=k, km1=k-1;
	for(i=kold; i<nold; i++){
		double data=knotv[i];
		if(data!=prev){ mult=1; prev=data;}
		else if(mult<km1) mult+=1;
		else continue;
		(*this)(n++)=data;
	}
	double eparam=knotv[nold];
	for(i=0; i<k; i++) (*this)(n++)=eparam;
		//Ending knot has multiplicity of order.
	set_length(n);
	m_current=m_order-1;
}

//Generate knot vector as a part of original knotv.
//New knot vector's parame range is from t1 to t2.
//Should hold t1>=knotv.param_s() && t2<=knotv.param_e().
//Knots between t1< <t2 are copied from knotv.
MGKnotVector::MGKnotVector(
	const MGKnotVector& knotv,double t1, double t2
){
	assert(t1>=knotv.param_s() && t2<=knotv.param_e() && t1<t2);

	t1=knotv.param_normalize(t1); int id1=knotv.locate(t1);
	t2=knotv.param_normalize(t2); int id2=knotv.locate(t2,1);
	int k=knotv.order();
	int n=id2-id1+k;
	if(t2==knotv(id2))
		n-=1;
	*this=MGKnotVector(k,n);
	int i,j, id1mkm1=id1-(k-1);

	if(t1==knotv(id1))
		for(i=0; i<k; i++)
			(*this)(i)=knotv(id1mkm1+i);
	else
		for(i=0; i<k; i++)
			(*this)(i)=t1;
	for(j=id1+1; j<=id2; j++)
		(*this)(i++)=knotv(j);
	if(t2==knotv(id2))
		for(j=1; j<k; j++)
			(*this)(i++)=knotv(id2+j);
	else
		for(j=0; j<k; j++)
			(*this)(i++)=t2;
	m_current=m_order-1;
}

//
// Operator overload /////////////////////////////
//

//Addition and subtraction of real number.
//All of the elements will be added or subtracted.
MGKnotVector MGKnotVector::operator+ (double a) const{
	MGKnotVector t(*this); t.MGNDDArray::operator+= (a); return t;
}
MGKnotVector& MGKnotVector::operator+= (double a){
	MGNDDArray::operator+= (a); return *this;
}
MGKnotVector MGKnotVector::operator- (double a) const{
	MGKnotVector t(*this); t.MGNDDArray::operator-= (a); return t;
}
MGKnotVector& MGKnotVector::operator-= (double a){
	MGNDDArray::operator-= (a); return *this;
}

// 単項マイナス。
//Unary minus. Reverse the ordering of elements by changing all of the
//signs.
MGKnotVector MGKnotVector::operator- () const{
	MGKnotVector t(*this); t.MGNDDArray::operator- (); return t;
}

//Scaling.
MGKnotVector MGKnotVector::operator* (double scale) const{
	MGKnotVector t(*this); t.MGNDDArray::operator*= (scale); return t;
}
MGKnotVector operator* (double scale, const MGKnotVector& t){
	return t*scale;
}
MGKnotVector& MGKnotVector::operator*= (double scale){
	MGNDDArray::operator*= (scale); return *this;
}

//Compare two KnotVector if they are equal.
//Return true if equal.
bool MGKnotVector::operator== (const MGNDDArray& t2) const{
	auto kv = dynamic_cast<const MGKnotVector*>(&t2);
	if(!kv)
		return false;
	return operator== (*kv);
}
bool MGKnotVector::operator== (const MGKnotVector& t2) const {
	if (m_order != t2.m_order) return false;
	return MGNDDArray::operator== (t2);
}

//////// Public member function //////

MGKnotVector& MGKnotVector::change_knot_number(int nnew){
	assert(nnew>=order());

	int k = order();
	int km1 = k - 1;
	int brdim=bdim();
	MGNDDArray tDDArray(km1, brdim-km1+1, *this);
	tDDArray.change_number(nnew);
	MGKnotVector knotVector(k, k);
	int i = 0;
	for(i=0; i<k; i++){
		knotVector(i) = (*this)(i);
		knotVector(k+i) = (*this)(brdim+i);
	}
	int nnewm1=nnew-1;
	for(i=1; i<nnewm1; i++){
		knotVector.add_data(tDDArray(i), k);
	}
	*this = std::move(knotVector);
	return *this;
}

//Change order for area adjustment. Allowed only for garbage knot vector.
MGKnotVector& MGKnotVector::change_order(int k){
	assert(k>=1);

	reshape(k+bdim());
	m_order=k;
	m_current=k-1;
	return *this;
}

// Change value range. //BKCRNG
void MGKnotVector::change_range(double ts, double te){
	assert(ts<te);
	int brdim=bdim();
	bkcrng_(m_order, brdim, data(), ts, te, data()); 
}

//Divide every spans. Every spans are subdivided into num equal spans.
//Result knot vector's bdim() becomes approximately bdim()*num.
void MGKnotVector::divide_span(int num){
	int k = order();
	int km1 = k - 1;
	MGKnotVector copyKnotVector(*this);
	int brdim=copyKnotVector.bdim();
	for(int i = km1; i < brdim; i++){
		double	spara = copyKnotVector.ref(i),      //スパンの始点
				epara = copyKnotVector.ref(i + 1);  //スパンの終点
		if(epara - spara == 0.0)
			continue;  //マルチノットのときの処理
		double shortspan = (epara-spara)/num, tmpParam = spara;
		tmpParam += shortspan;
		for(int j=1; j<num; j++){
			add_data(tmpParam); tmpParam += shortspan;
		}
	}
}

#define MAX_ORDER 20
// Function's return value id is the index of B-coefficients that
// should be multiplied to.
// coef[a] is for (id+a)-th B-Coefficients, 0<= a <=order-1.
// Multiplication with coef should be done like:
//
// data=0.;
// for(int a=0; a<order; a++) data+=coef[a]*B_coef[id+a].
//
//left indicates whether left continuous(left=true), or right continuous
//(left=false) evaluation.
int MGKnotVector::eval_coef(
	double x,		// Parameter value to evaluate
	double* coef,	// coef's to mutiply are returned.
	int nderiv,//order of derivative, =0 means position
	int left		//Left continuous(left=true)
					//or right continuous(left=false).
)const{
	int k=m_order;
	int mu=locate(x,left);// k-1<= mu <=bdim()-1 is guaranteed.
	int mup1=mu+1;
	int id=mup1-k;
	if(nderiv>=int(k)){
		for(int a=0; a<k; ++a)
			coef[a]=0.f;
		return id;
	}

    double deltal_buf[MAX_ORDER], deltar_buf[MAX_ORDER];
	double* deltal=deltal_buf; double* deltar=deltar_buf;
	if(k>MAX_ORDER){
		deltal=new double[k+k]; deltar=deltal+k;
	}

	//deltar[j] will be t[mu+j+1]-x
	//deltal[j] will be x-t[mu-j]            for j=0, ... ,k-1-nderiv
	const MGKnotVector& t=*this;
    int kmnd = k-nderiv;
    coef[0] = 1.f;//=Bmu,1(x)
    int j;
	for(j=1; j<kmnd; j++){
		int jm1=j-1;
		int mup1mj = mup1-j;
		deltar[jm1] = t[mu+j]-x;
		deltal[jm1] = x-t[mup1-j];
		double saved = 0.f;
		for(int a=0; a<j; ++a) {
			double term = coef[a]/(t[mup1+a]-t[mup1mj+a]);
			coef[a] = saved+deltar[a]*term;
			saved = deltal[jm1-a]*term;
		}
		coef[j] = saved;
	//Here coef[a] contains Bmu-j+a,j+1(x) 
	//(B-coefficient of (mu-j+a)th and of order j+1) for a=0, ... , j.
	}
	if(k>MAX_ORDER){
		delete[] deltal;
	}

	// NOW B-SPLINE coef(a) 0<=a<=K-JDERIV(ORDER OF K-JDERIV)-1 ARE OBTAINED
	// GET DERIVATIVE PART OF THE COEFFICIENTS.
	// Here j == kmnd, coef[a] contains Bmu-k+nderiv+1+a,k-nderiv(x) 
	//(B-coefficient of (mu-k+nderiv+1+a)th and of order k-nderiv) for a=0, ... , k-nderiv-1.
	for(;j<k;j++){
	    double fj = double(j);
		int mup1mj = mup1-j;
		double saved = 0.f;
		for(int a=0; a<j; ++a){
			double term=coef[a]*fj/(t[mup1+a]-t[mup1mj+a]);
			coef[a] = saved-term;
			saved =term;
	    }
		coef[j] = saved;
	}

	return id;
}

//Test if input parameter value is inside parameter range of the line.
bool MGKnotVector::in_range(double t) const{
	const double t1=param_s(), t2=param_e();
	double error=(t2-t1)*MGTolerance::rc_zero();
	return (t>=t1-error && t<=t2+error);
}

// finds where tau is located in MGKnotVector.
// Returned id is :
//      order()-1<= id < bdim().
//left indicates whether left continuous(left=true), or right continuous
//(left=false) location.
int MGKnotVector::locate(double tau, int left) const{
	const MGKnotVector& t=*this;
	double tc=t[m_current], tcp1=t[m_current+1];
	if(left){
		if(tc<tau && tau<=tcp1)
			return m_current;
	}else{
		if(tc<=tau && tau<tcp1)
			return m_current;
	}

	int nm1=length()-1;
	m_current=MGNDDArray::locate(tau);
    if(m_current<0){
		//FIND SMALLEST LEFT SUCH THAT t(m_current)<t(m_current+1)
		m_current = 0;
		while(m_current<nm1){
			int cp1 = m_current + 1;
			if(t[m_current]<t[cp1])
				break;
			m_current = cp1;
		}
    }else if(m_current>=nm1){
		//FIND LARGEST LEFT SUCH THAT t(m_current) < t(m_current+1)
		m_current = nm1 - 1;
		int cp1=nm1;
		while(0<m_current){
			if(t[m_current]<t[cp1])
				break;
			cp1=m_current--;
		}
    }

	if(left || tau>=param_e()){
		while(tau==t[m_current] && m_current>=int(m_order))
			m_current--;
	}

	int km1=m_order-1,brdim=bdim();
	if(m_current<km1) m_current=km1;
	if(m_current>=brdim) m_current=brdim-1;
	return m_current ;
}

// finds where tau is located in MGKnotVector.
// Returned id is :
//      0<= id <= bdim()-1.
int MGKnotVector::locate(double tau) const{
	return locate(tau,0);
}

//Locate where data of multiplicity of multi is after start and before
//bdim(). index is the starting point index of this found first
//after start, index>=start.
//Function's return value locate_multi is actual multiplicity at the
//index, i.e. locate_multi>=multi if found.
//If position of the multiplicity is not found before bdim(),
//index=bdim() (index of the param_e()) and locate_multi=0 will be returned.
//multi must be >=1 and start must be <=bdim().
int MGKnotVector::locate_multi(
	int start, int multi, int& index
)const{
	assert(start<=bdim());

	int brdim=bdim(),multi_found;
	if(start>=brdim)
		index=start;
	else
		multi_found=MGNDDArray::locate_multi(start, multi, index);
	if(index>=brdim){
		multi_found=0;
		index=brdim;
	}
	return multi_found;
}

//ノットベクトルをショートスパンのできないように足しあわせる
//ノットベクトルのパラメータ範囲は等しいものとする。
//エラーのとき元のノットベクトルを返す。
MGKnotVector& MGKnotVector::mix_knot_vector(const MGKnotVector& knot){
	double span1 = param_span(), span2 = knot.param_span();
	if(!MGRZero((span1 - span2) / span1))return *this;
	int k1 = order(), k2 = knot.order();
	int p1 = k1 - 1, p2 = k2 - 1;
	//引数のノットベクトルに多重度があるときには、thisノットベクトルに加えておく
	int i = 0;
	for(i = k2; i < knot.bdim(); i++){
		int index = 0, multi = 0;
		multi = knot.locate_multi(i, 2, index);
		if(multi > 1){	//多重度があったときの処理
			add_data(MGKnot(knot(index), multi), multi);
		}
		i = index;
	}
	MGNDDArray tDDArray1(p1, bdim() - (p1 - 1), *this);
	MGNDDArray tDDArray2(p2, knot.bdim() - (p2 - 1), knot);
	MGNDDArray mixDDArray(0, tDDArray1.length(), tDDArray1,	0, tDDArray2.length(), tDDArray2);

	MGKnotVector knotVector(k1, k1);	//ミックス後のオーダはthisノットベクトルのを用いる
	for(i = 0;i < k1; i++){
		knotVector(i) = mixDDArray(0);
		knotVector(knotVector.bdim() + i) = mixDDArray(mixDDArray.length() - 1);
	}
	for(i = 1; i < mixDDArray.length() - 1; i++){
		knotVector.add_data(mixDDArray(i), k1 - 1);
	}
	*this = std::move(knotVector);
	return *this;
}

//Return tolerance allowed in knot vector parameter space.
double MGKnotVector::param_error()const{
	return param_span()*MGTolerance::rc_zero();
}

//Normalize parameter value t to the nearest knot if their distance is
//within tolerance.
double MGKnotVector::param_normalize(double t) const{
	double plen=param_span();
	int i=locate(t);
	double tnew1=(*this)(i), tnew2=(*this)(i+1);
	double dif1=t-tnew1, dif2=tnew2-t;
	if(dif1<=dif2){
		if(MGRZero2(dif1,plen)) return tnew1;
	}else{
		if(MGRZero2(dif2,plen)) return tnew2;
	}
	return t;
}

double MGKnotVector::param_span() const{
	return (*this)(bdim())-(*this)(order()-1);
}

// 入力パラメータをパラメータ範囲でまるめて返却する。
double MGKnotVector::range(double t) const{
	const double t1=(*this)(order()-1);
	double s;
	if(t<t1)	s=t1;
	else{
		const double t2=(*this)(bdim());
		if(t>t2) s=t2;
		else s=t;
	}
	return s;
}

//Resize this so that this order is order and this b-rep dimension is brdim.
//Result of size_change will contain garbages.
void MGKnotVector::size_change(int k, int brdim){
	MGNDDArray::resize(k+brdim);
	m_order=k;
	m_current=k-1;
}

//Reverse the ordering of knots.
void MGKnotVector::reverse(){
	double add=param_s()+param_e();
	int len=bdim()+order();
	int half=len/2, i2=len-1;
	double save;
	int i;
	for(i=0; i<half; i++){
		save=(*this)(i);
		(*this)(i)=add-(*this)(i2);
		(*this)(i2--)=add-save;
	}
	if(i2>=i)
		(*this)(i)=add-(*this)(i);
}

//Obtain parameter value if this knot vector is reversed by reverse().
double MGKnotVector::reverse_param(double t) const{
	return param_s()+param_e()-t;
}

//Set B-rep dimension.
//Only set the B-Rep Dimension.
MGKnotVector& MGKnotVector::set_bdim(int brdim){
	assert(m_order+brdim<=capacity());

	set_length(m_order+brdim);
	if(m_current>=int(brdim))
		m_current=brdim-1;
	return *this;
}
