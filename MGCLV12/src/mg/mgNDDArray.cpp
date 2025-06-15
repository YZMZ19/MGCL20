/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/NDDArray.h"
#include "mg/KnotVector.h"
#include "mg/Knot.h"
#include "mg/BPointSeq.h"
#include "mg/Tolerance.h"

#include "cskernel/bkdtpg.h"
#include "cskernel/Bkdtwe.h"
#include "cskernel/bkktdp.h"
#include "cskernel/Bkcdtn.h"
#include "cskernel/Bkcrng.h"
#include "cskernel/Bkdmix.h"
#include "cskernel/Bkdnp.h"
#include "cskernel/bkmlt.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGNDDArray
// Implementation of Class MGNDDArray


//Copy constructor.
MGNDDArray::MGNDDArray(const MGNDDArray& nd)
:m_capacity(nd.m_length), m_length(nd.m_length),
m_element(new double[nd.m_length]),m_current(nd.m_current){
	for(int i=0; i<m_length; i++) m_element[i]=nd.m_element[i];
}

//Copy Assignment
MGNDDArray& MGNDDArray::operator=(const MGNDDArray& vec2){
	if(this!=&vec2){
		int len=vec2.length();
		if(capacity()<len) resize(len);
		else m_length=len;
		for(int i=0; i<len; i++) m_element[i]=vec2.m_element[i];
		m_current=vec2.m_current;
	}
	return *this;
}

//Move constructor.
MGNDDArray::MGNDDArray(MGNDDArray&& nd)
:m_capacity(nd.m_capacity), m_length(nd.m_length),
m_element(nd.m_element),m_current(nd.m_current){
	nd.m_element=nullptr;
	nd.m_capacity=nd.m_length=nd.m_current=0;
}

//Move Assignment
MGNDDArray& MGNDDArray::operator=(MGNDDArray&& nd){ 
	m_capacity=nd.m_capacity;
	m_length=nd.m_length;
	m_current=nd.m_current;
	delete[] m_element;
	m_element=nd.m_element;

	nd.m_element=nullptr;
	nd.m_capacity=nd.m_length=nd.m_current=0;
	return *this;
}

/// Constructor MGNDDArray of size n and lenght=n.
/// When data is given, data is an array of length n, and
/// data construct the element data of this MGNDDArray.
/// data[i] must be <= data[i+1].
MGNDDArray::MGNDDArray(
	int n,		//size of this.
	const double* data	//data array of length n if data!=NULL.
):m_capacity(n),m_length(n),m_current(0),m_element(nullptr){
	if(n){
		m_element=new double[n];
		if(data)
			for(int i=0; i<n; i++) m_element[i]=data[i];
	}
}

// Construct MGNDDArray	with initial init and incrementation
// of increment.
MGNDDArray::MGNDDArray(int n, double init, double increment)
:m_length(n),m_capacity(n),m_current(0),m_element(nullptr){
	if(n){
		m_element=new double[n];
		double value=init;
		for(int i=0; i<n; i++) {
			(*this)(i)=value;
			value=value+increment;
		}
	}
}

// From Data Point ordinate, obtain data point seq abscissa
// by adding chord length to starting value zero.
MGNDDArray::MGNDDArray(const MGBPointSeq& bps)
:m_element(new double[bps.length()]),m_length(bps.length())
,m_capacity(bps.length()),m_current(0){
	const int ncd=bps.sdim();
	const int ip=bps.capacity();
	bkdtpg_(bps.data(), m_length, ncd, ip, data()); 
}

// From Data Point ordinate with End condition, obtain data point seq
// abscissa by adding chord length to starting value zero.
MGNDDArray::MGNDDArray(MGENDCOND beginc, MGENDCOND endc, const MGBPointSeq& bps)
:m_element(new double[bps.length()]),m_length(bps.length())
,m_capacity(bps.length()),m_current(0){
	assert(beginc!= MGENDCOND::MGENDC_12D && endc!= MGENDCOND::MGENDC_12D);
	const int ncd=bps.sdim();
	const int ip=bps.capacity();
	int bc = static_cast<int>(beginc), ec = static_cast<int>(endc);
	bkdtwe_(bc,ec,bps.data(),m_length,ncd,ip, data());
}

//From knot vector with End conditions, obtain data point
//seq abscissa. This data point can be input of MGSBRep constructor
//of the following form:
//MGSBRep::MGSBRep(
//MGSBRepEndC& endc,//end condition
//const MGNDDArray& utaui,	//Data point of u-direction
//const MGNDDArray& vtaui,	//Data point of v-direction
//const MGSPointSeq& value,	//Data point ordinate
//const MGKnotVector& tu,	//knot vector of u-direction, of order 4
//const MGKnotVector& tv,	//knot vector of v-direction, of order 4
//int &error)				//Error flag.
//That is, generate utaui or vtaui from tu or tv each, taking endc 
//into account.
//The order is assumed to be 4.
MGNDDArray::MGNDDArray(
	MGENDCOND begin, MGENDCOND end, const MGKnotVector& t
):m_current(0){
	assert(t.order()>=4 && t.bdim()>=4);

	int km1=t.order()-1, n=t.bdim();
	int nstart=0;
	if(begin== MGENDCOND::MGENDC_1D || begin== MGENDCOND::MGENDC_2D) nstart=1;
	else if(begin== MGENDCOND::MGENDC_12D) nstart=2;
	int nend=0;
	if(end== MGENDCOND::MGENDC_1D || end== MGENDCOND::MGENDC_2D) nend=1;
	else if(end== MGENDCOND::MGENDC_12D) nend=2;

	int newn=n-nstart-nend;
	m_element=new double[newn];
	m_capacity=m_length=newn;
	int newnm1=newn-1;

	int id,id1;
    double fkm1 = double(km1);
    m_element[0]= t(km1);
	for(int i=1;i<newnm1; i++){
		double r=0.;
		id1=i+nstart;
		for (int j = 1; j <= km1; ++j) {
		    id=id1+j;
			if (id < km1) id=km1;
			if (id > n)   id=n;
		    r+= t(id);
		}
		m_element[i]=r/fkm1;
	}
	m_element[newnm1]=t(n);
}

// From data point, obtain data point of updated number.
MGNDDArray::MGNDDArray(const MGNDDArray& tau1, int nnew)
:m_element(0),m_length(0),m_capacity(0),m_current(0){
	*this=tau1;
	change_number(nnew);
}

// From data point, obtain data point of updated value range.
MGNDDArray::MGNDDArray(const MGNDDArray& tau1, double ts, double te) //BKCRNG
:m_element(new double[tau1.length()]),m_length(tau1.length())
,m_capacity(tau1.length()),m_current(0){
	assert(ts<te);

	const int k=0;
	bkcrng_(k, m_length, tau1.data(), ts, te, data()); 
}

// Construct by extracting sub interval of array2.
MGNDDArray::MGNDDArray(
	int start_id,			//Start id of array2(from 0).
	int num,					//new array length.
	const MGNDDArray& array2	//Original NDDArray.
):m_element(new double[num]),m_length(num),m_capacity(num),m_current(0){
	assert((start_id+num)<=array2.length());
	int is=start_id;
	for(int i=0; i<num; i++, is++) (*this)(i)=array2(is);
}

//Construct by mixing two arrays.
//Mixing is so done as that no too close points are not included.
//Although data point multiplicities of array1 are preserved,
//multiplicities of array2 are not.
//DATA POINT MULTIPLICITY IS ALLOWED ONLY IN array1.
MGNDDArray::MGNDDArray(
	int id1,			//Start id of array1(from 0).
	int num1,				//new array length to use of array1.
	const MGNDDArray& array1,	//Original NDDArray1.
	int id2,			//Start id of array2(from 0).
	int num2,				//new array length to use of array2..
	const MGNDDArray& array2	//Original NDDArray2.
):m_element(new double[num1+num2]),m_capacity(num1+num2),m_current(0){
	assert((id1+num1)<=array1.length());
	assert((id2+num2)<=array2.length());

	bkdmix_(num1,array1.data(id1),num2,array2.data(id2),&m_length,data());
}

//Set the length of effective data.
void MGNDDArray::set_length(int len){
	assert(len<=capacity());
	m_length=len;
	m_current=0;
}

//Add data of multiplicity 1 into data points. mult_max is the maximum
//multiplicity allowed for NDDArray.
//Return value is number of data actually added.
int MGNDDArray::add_data(double value, int mult_max){
	return add_data(MGKnot(value, 1),mult_max);
}

//Add data with multiplicity into data points.
int MGNDDArray::add_data(
	const MGKnot& knot,	//data with multiplicity.
	int mult_max	//maximum multiplicity allowed.
){
	assert(mult_max>0 && knot.multiplicity()>0);

	double tau=knot; int imult=knot.multiplicity();
	int mult = mult_max>imult ? imult : mult_max;
	int i=MGNDDArray::locate(tau);
	if(i>=0){
		int current_mult=i;
		while(i>=0 && (*this)(i)==tau) i-=1;
		current_mult=current_mult-i; //current_mult is the multiplicity of tau.
		int able_to_add = mult_max - current_mult; //number of multiplicity  
		if(able_to_add < mult) mult=able_to_add;
	}
	if(mult<0) mult=0;
	int oldlength = m_length;
	int newlength=m_length+mult;
	reserve(newlength);//guarantee capacity()>=newlength
	set_length(newlength);
	int rval=mult;
	if(mult){
		int endid=newlength-1; int j= oldlength-1;
		while(j>i) (*this)(endid--)=(*this)(j--);
		while(mult--) (*this)(++i)=tau;
	}
	return rval;	//Return number of data added.
}

//Delete one data at index.
int MGNDDArray::del_data(int index){
	assert(length()>0);
	if(index<m_length){
		int i=index; int j=i+1;
		while(j<m_length) (*this)(i++)=(*this)(j++);
		m_length -=1;
	}
	m_current=0;
	return m_length; //Return new length of the array.
}

// finds index where tau is located in MGNDDArray.
// tau < (*this)(0) : index=-1.
// (*this)(0) <= tau< (*this)(n-1) : 0<= index <n-1, such that
//                         (*this)(index) <= tau < (*this)(index+1)
// (*this)(n-1) <= tau : index=n-1.
//Here n=lenght().
int MGNDDArray::locate(double tau) const{
	int nm1=m_length-1;
	if(tau>=m_element[nm1])
		return nm1;
	if(tau<m_element[0])
		return -1;

    int ihi = m_current+1;
	if(ihi>nm1){
	    m_current = nm1-1;
		ihi = nm1;
	}

    int istep=1;
    if(tau<m_element[ihi]){
	    if(m_element[m_current]<=tau)
			return m_current;

		//**** NOW tau<m_element[m_current]. DECREASE m_current TO CAPTURE tau. 
		do{
		    ihi = m_current;
			m_current = ihi - istep;
		    istep <<= 1;//Multiply by 2.
		}while(m_current>0 && tau<m_element[m_current]);
	}else{
		//**** NOW m_element[ihi]<=tau. INCREASE ihi TO CAPTURE  tau . 
		do{
		    m_current = ihi;
			ihi = m_current + istep;
		    istep <<= 1;//Multiply by 2.
		}while(ihi<nm1 && m_element[ihi]<=tau);
	}
	if(m_current<0)
		m_current=0;
	if(ihi>nm1)
		ihi=nm1;

// **** NOW m_element[m_current] <= tau < m_element[ihi] . NARROW THE INTERVAL. 
	while(1){
		int middle=(m_current+ihi)/2;
		if(middle==m_current)//IT IS ASSUMED THAT MIDDLE=m_current IN CASE ihi=m_current+1. 
			break;
		if(tau<m_element[middle])
			ihi=middle;
		else
			m_current=middle;
	};

	return m_current;
}

//Locate where data of multiplicity of multi is after start.
int MGNDDArray::locate_multi(int start, int multi, int& index)
const{
	assert( multi>=1 );

	int is=start+1;
	int multi_found; int isf;
	bkmlt_(m_length, data(), is, multi, &isf, &multi_found);
	if(isf>0)
		isf=isf-1;
	else
		isf=length()-1;
	index=isf;
	return multi_found;
}

// Change data point range.
void MGNDDArray::change_range(double ts, double te){
	assert(ts<te);

	const int k=0;
	bkcrng_(k, m_length, data(), ts, te, data());
}

//Update number of data points.
MGNDDArray& MGNDDArray::change_number(int nnew){
	assert(nnew>0 && length()>0);
	if(m_length == nnew) return *this;

	MGNDDArray temp(nnew);
	const int lenm1=m_length-1;
	double dlenm1=ref(lenm1), d0=ref(0);
	if(dlenm1<=d0)
		for(int i=0; i<nnew; i++)
			temp(i)=d0;
	else{

	if(nnew<m_length)
		bkcdtn_(m_length, data(), nnew, temp.data());
	else{
		const int inc_total=nnew-m_length;
		const double ratio=double(inc_total)/(dlenm1-d0);
		int iend;
		int mult=locate_multi(0, 1, iend);
		int i1=0, j1=0, inc_now=0;
		int inext=i1+mult;
		while(inext<=lenm1){
			mult=locate_multi(inext, 2, iend); if(!mult) mult++;
			int n1=iend-i1+1;
			int inc;
			if(iend+mult>lenm1)
				inc=inc_total-inc_now;
			else{
				inc=int((ref(iend)-ref(i1))*ratio+.5);
				if(inc+inc_now>inc_total)
					inc=inc_total-inc_now;
			}
			int n2=n1+inc;
			bkcdtn_(n1, data()+i1, n2, temp.data()+j1);
			i1+=n1-1; j1+=n2-1; inc_now+=inc;
			inext=i1+mult;
		}
		while(j1<nnew)
			temp(j1++)=ref(i1);
	}

	}
	*this=temp;
	return *this;
}

//Copy data points from array by removing the mltiple knot.
void MGNDDArray::copy_removing_multi(
	int start_id,			//Start id of array2(from 0).
	int num,					//extracting length from start_id.
	const MGNDDArray& array2	//Original NDDArray.
){
	assert((start_id+num)<=array2.length());
	resize(num);
	int j=start_id, i=0;
	double prev=array2(j++);
	(*this)(i++)=prev;
	int ie=start_id+num;
	for(; j<ie; j++){
		double current=array2(j);
		if(current>prev){
			(*this)(i++)=current;
			prev=current;
		}
	}
	set_length(i);
}

//Remove too near data points. Removal will be done with the ordinates.
void MGNDDArray::remove_too_near(
	bool allow_multi,	//indicates if multiple data point is allowed or not.
							//when allow_multi=false, multiple data points will be removed.
							//when allow_multi=true, will not be removed.
	double ratio	//maximum ratio allowed for neighboring span.
			//let ti=(*this)[i], then
			//if (t(i+1)-ti)/(ti-t(i-1))>ratio or (t(i+1)-ti)/(ti-t(i-1))<1/ratio,
			//a data point will be removed(along with the ordinates).
){
	int n=length(), ordinates_size=0, sd=0;
	int IMLT=1;
	if(allow_multi)
		IMLT=0;
	bkdnp_(&n,data(),0,ordinates_size,sd,IMLT,ratio);		
	set_length(n);
}

//Remove too near data points. Removal will be done with the ordinates.
void MGNDDArray::remove_too_near(
	MGBPointSeq& ordinates,	//ordinate.
							//ordinates.length() must be equal to this->length().
	bool allow_multi,	//indicates if multiple data point is allowed or not.
							//when allow_multi=false, multiple data points will be removed.
							//when allow_multi=true, will not be removed.
	double ratio	//maximum ratio allowed for neighboring span.
			//let ti=(*this)[i], then
			//if (t(i+1)-ti)/(ti-t(i-1))>ratio or (t(i+1)-ti)/(ti-t(i-1))<1/ratio,
			//a data point will be removed(along with the ordinates).
){
	assert(ordinates.length()==length());

	int n=length(), ordinates_size=ordinates.capacity(), sd=ordinates.sdim();
	int IMLT=1;
	if(allow_multi) IMLT=0;
	bkdnp_(&n,data(),ordinates.data(),ordinates_size,sd,IMLT,ratio);		
	set_length(n); ordinates.set_length(n);
}

//Change the size. 
void MGNDDArray::reshape(int size, int start){
	int new_length=m_length+start;
	if(size<new_length) new_length=size;
	int n=new_length-start;	//n is number of data to move.
 
	double* data=(size>m_capacity) ? new double[size] :m_element;
	if(start || m_element!=data)
		for(int i=n-1; i>=0 ; i--)
			data[start+i]=m_element[i];

	if(m_element!=data){
		delete[] m_element;
		m_element=data;
		m_capacity=size;
	}
	m_length= size;
	m_current=0;
}

//Resize the array. Result will contain garbages.
void MGNDDArray::resize(int nsize){
	if(m_element) delete[] m_element;
	m_element=new double[nsize];
	m_capacity=m_length=nsize;
	m_current=0;
}

//Set this as a null NDDArray.
void MGNDDArray::set_null(){
	if(m_element) delete[] m_element;
	m_element=0;
	m_capacity=m_length=0;
	m_current=0;
}

#define INITIAL_LENGTH 20
#define INCREMENTAL 1.5
/// Reserve capacity.
/// When newCapacity>=1 is input, adequate new capacity() is made to guarantee capacity() >= newCapacity.
/// When newCapacity<=0 is input, the most adequate new capacity() is made .
/// The length() and the data stored are unchanged.
void MGNDDArray::reserve(int newCapacity) {
	if (1<=newCapacity && newCapacity <= m_capacity)
		return;

	int new_capa =
<<<<<<< HEAD
		m_capacity >= INITIAL_LENGTH ? int(m_capacity*INCREMENTAL) : INITIAL_LENGTH;
	while(new_capa<=newCapacity)
		new_capa = int(new_capa*INCREMENTAL);
=======
		m_capacity >= INITIAL_LENGTH ? m_capacity * INCREMENTAL : INITIAL_LENGTH;
	while(new_capa<=newCapacity)
		new_capa *= INCREMENTAL;
>>>>>>> 43f8b42bf874d0b96a44d789ec156c4bbb7cabd9
	int lenSave = length();
	reshape(new_capa);
	set_length(lenSave);
}

/// Reserve length n area before the index id.
/// On return, length() will be increased by n.
/// Area from [id] to [id+n-1] will contain garbage.
/// The old data from [id] to [lenOld-1] will be transfered to [id+n] ...
void MGNDDArray::reserveInsertArea(int n, int id) {
	assert(id <= m_length );
	int lenNew = m_length + n;
	double* dataNew;
	if (m_capacity >= lenNew) {
		dataNew = m_element;
	}else{
		dataNew = new double[lenNew];
		for (int i = 0; i<id; ++i)
			dataNew[i] = m_element[i];
	}
	int lenOldm1 = m_length - 1;
	int lenNewm1 = lenNew - 1;
	for (int i = lenOldm1, j = 0; id <= i; --i, ++j)
		dataNew[lenNewm1 - j] = m_element[lenOldm1 - j];
	if (dataNew != m_element) {
		delete[] m_element;
		m_element = dataNew;
		m_capacity = lenNew;
	}
	m_length = lenNew;
}

///store data at the position i of this array.
///tau(i) must be => tau(i-1).
///When i>=capacity(), reshape will take place.
///The data validity is not checked.
///The length is set same as the capacity. User must set the length by set_length().
void MGNDDArray::store_with_capacityCheck(
	int i,	///the postion to store at.
	double data///the data to store.
){
	assert(i>=0);
	if(i>=m_capacity)
		reserve();
	m_element[i]=data;
	if (m_length <= i)
		set_length(i + 1);
}

//From Knot vector, obtain data point seq abscissa.
MGNDDArray& MGNDDArray::buildByKnotVector(const MGKnotVector& t){
	m_length=t.bdim();
	if(capacity()<m_length) resize(m_length);
	const int k=t.order();
	bkktdp_(m_length, k, t.data(), data()); 
	m_current=0;
	return *this;
}

//Addition and subtraction of real number.
//All of the elements will be added or subtracted.
MGNDDArray MGNDDArray::operator+ (double a) const{
	MGNDDArray t(*this); return t+=a;
}
MGNDDArray& MGNDDArray::operator+= (double a){
	for(int i=0; i<m_length; i++) m_element[i]+=a;
	return *this;
}
MGNDDArray MGNDDArray::operator- (double a) const{
	MGNDDArray t(*this); return t-=a;
}
MGNDDArray& MGNDDArray::operator-= (double a){
	for(int i=0; i<m_length; i++)
		m_element[i]-=a;
	return *this;
}

// 単項マイナス。
//Unary minus. Reverse the ordering of elements by changing all of the
//signs.
MGNDDArray MGNDDArray::operator- ()const{
	MGNDDArray t(m_length);
	for(int i=0, j=m_length-1; i<m_length; i++, j--)
		t.m_element[i]=-m_element[j];
	return t;
}

//Scaling
MGNDDArray MGNDDArray::operator* (double scale) const{
	MGNDDArray t(*this); return t*=scale;
}
MGNDDArray operator* (double scale, const MGNDDArray& nd){
	return nd*scale;
}
MGNDDArray& MGNDDArray::operator*= (double scale){
	if(scale>=0.){
		for(int i=0; i<m_length; i++)
			m_element[i]*=scale;
	}else{
		int half=m_length/2;
		for(int i=0, j=m_length-1; i<half; i++, j--){
			double x=m_element[i];
			m_element[i]=scale*m_element[j];
			m_element[j]=scale*x;
		}
		if(half*2!=m_length)
			m_element[half]*=scale;
	}
	return *this;
}

//Copmarison operator.
bool MGNDDArray::operator== (const MGNDDArray& t2) const{
	if(m_length != t2.m_length) return false;
	if(m_length<=0) return true;
	if(m_length==1) return MGREqual((*this)(0),t2(0));
	int n=m_length-1;
	for(int i=0; i<n; i++){
		if(!MGREqual2(m_element[i+1]-m_element[i],
			t2.m_element[i+1]-t2.m_element[i]))
			return false;
	}
	return true;
}
