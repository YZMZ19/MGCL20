/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "cskernel/Bvavpl.h"
#include "mg/Tolerance.h"
#include "mg/Vector.h"
#include "mg/Matrix.h"
#include "mg/Transf.h"
#include "mg/Box.h"
#include "mg/BPointSeq.h"
#include "mg/SPointSeq.h"
#include "mg/Straight.h"
#include "mg/Plane.h"


#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
// Implementation of Class MGBpointSeq

//Constructor
MGBPointSeq::MGBPointSeq(int capacity, int dim)
:m_capacity(capacity),m_sdim(dim),m_length(capacity){
	int len=capacity*dim;
	if(len)
		m_bpoint=new double[len];
	else{
		m_capacity=m_length=0;
		m_bpoint=nullptr;
	}
}

//Construct a MGBPointSeq by copying original MGBPointSeq.
//Can change the order of coordinates 
MGBPointSeq::MGBPointSeq(
		int dim,				// New Space Dimension.
		const MGBPointSeq& old_bpoint,	// Origianl BPointSeq.
		int start1,			// Destination start order to store.
		int start2			// Source start order to retrieve.
):m_capacity(old_bpoint.m_length),m_sdim(dim)
, m_length(old_bpoint.m_length){
	assert(start1<dim && start2<old_bpoint.m_sdim);

	int len=m_length*dim;
	if(len){
		m_bpoint=new double[len];
		int dim2=old_bpoint.m_sdim;
		int dimmin= dim<dim2 ? dim:dim2; 
		int k,j2=start2,j1=start1;
		for(k=0; k<dimmin; k++){
			for(int i=0; i<m_length; i++)
				(*this)(i,j1)=old_bpoint(i,j2);
			j1 +=1;
			if(j1>=dim)
				j1=0; 
			j2 +=1;
			if(j2>=dim2)
				j2=0;
		}
		while(k++<dim){
			for(int i=0; i<m_length; i++)
				(*this)(i,j1)=0.0;
			j1 +=1;
			if(j1>=dim)
				j1=0; 
		}
	}else
		m_bpoint=0;
}

// Construct by extracting sub interval of array2.
MGBPointSeq::MGBPointSeq(
	int start_id,		// Start id of bp_old
	int num,						// New length(of new BPoint)
	const MGBPointSeq& bp_old		// Origianl BPoint.
): m_capacity(num),m_sdim(bp_old.sdim()), m_length(num)
, m_bpoint(new double[num*bp_old.sdim()]){
	assert(start_id<bp_old.length() && start_id+num<=bp_old.length());

	if(m_length==0){
		m_bpoint=0;
		return;
	}
	for(int m=0; m<m_sdim; m++)
		for(int i=0; i<num; i++)
			(*this)(i,m)=bp_old(start_id+i,m);
}

//Conversion constructor.
MGBPointSeq::MGBPointSeq(const std::vector<MGPosition>& poses)
:m_capacity(0),m_sdim(0),m_length(0),m_bpoint(0){
	int i, n=(int)poses.size();
	if(!n)
		return;

	int sdmax=poses[0].sdim();
	for(i=1; i<n; i++)
		if(poses[i].sdim()>sdmax)
			sdmax=poses[i].sdim();
	resize(n,sdmax);
	for(i=0; i<n; i++)
		for(int j=0; j<sdmax; j++)
			operator()(i,j)=poses[i][j];
}

//Copy constructor.
MGBPointSeq::MGBPointSeq(const MGBPointSeq& bp)
:m_sdim(bp.m_sdim),m_length(bp.m_length),m_capacity(bp.m_capacity) {
	if(m_capacity==0 || m_sdim==0){
		m_sdim=m_length=m_capacity=0;
		m_bpoint=nullptr;
		return;
	}
	m_bpoint=new double[m_sdim*m_capacity];
	for(int m=0; m<m_sdim; m++)
		for(int i=0; i<m_length; i++)
			(*this)(i,m)=bp(i,m);
}

//Move constructor.
MGBPointSeq::MGBPointSeq(MGBPointSeq&& bp)
:m_sdim(bp.m_sdim),m_length(bp.m_length),m_capacity(bp.m_capacity),
m_bpoint(bp.m_bpoint){
	bp.m_bpoint=nullptr;
	bp.m_sdim=bp.m_length=bp.m_capacity=0;
}

// Construct by extracting one line data of sp along u or v direction.
MGBPointSeq::MGBPointSeq(
	bool along_u,	//indicates which direction make a line out of sp.
		//if true, along u direction:sp(i,m,.) for i=0, ..., nu-1 makes the BPointSeq.
		//if false, along v direction:sp(m,j,.) for j=0, ..., nv-1 makes the BPointSeq.
	int m,	//index of u or v as above according to along_u.
	const MGSPointSeq& sp		// Origianl SPoint.
):m_sdim(0),m_length(0),m_capacity(0),m_bpoint(0){
	int n;
	if(along_u)
		n=sp.length_u();
	else
		n=sp.length_v();
	resize(n,sp.sdim());
	if(along_u)
		for(int i=0; i<n; i++)
			store_at(i,sp(i,m));
	else
		for(int j=0; j<n; j++)
			store_at(j,sp(m,j));
}

//////////Destructor//////////

//Member Function

//compute an average plane of the point sequence.
//Function's return value is:
// 0: Number of points included is zero.
// 1: Point seq is a point.		2: Point seq is on a line.
// 3: Plane is output. 
int MGBPointSeq::average_plane(
	MGPosition& center		// center of point seq will be output.
	, MGPlane& plane		// Plane will be output, when average_plane=3.
	, MGStraight& line		// Straight line will be output            =2.
	, double& deviation)	// Maximum deviation from point, line
	 const{					// , or plane will be output.
	const double error=MGTolerance::wc_zero();	//error allowed.
	int n=length(), ncd=sdim(), ip=capacity();
	int kplane; double g[4], fcenter[3]; const int one=1;
	bvavpl_(error,n,one,data(),ncd,ip,one,&kplane,g,fcenter,&deviation);
	MGVector direction(3,g);
	center=MGPosition(ncd, fcenter);
	switch (kplane){
	case 1:	break;	//When Point.
	case 2:			// Straight line
		line=MGStraight(MGSTRAIGHT_TYPE::MGSTRAIGHT_UNLIMIT, direction, center);
		break;
	case 3:			//Plane will be output.
		plane=MGPlane(direction, g[3]);
		break;
	default:
		break;
	}
	return kplane;
}

//Compute minimum box sorrounding the all the points.
MGBox MGBPointSeq::box()const{
	if(is_null()) return mgNULL_BOX;
	int dim=sdim();
	MGInterval* intrvl=new MGInterval[dim];
	int len=length();
	for(int j=0; j<dim; j++){
		double dmax,dmin; dmax=dmin=ref(0,j);
		for(int i=1; i<len; i++){
			double dataij=ref(i,j);
			dmax= dmax>=dataij ? dmax:dataij;
			dmin= dmin<=dataij ? dmin:dataij;
		}
		intrvl[j]=MGInterval(dmin,dmax);
	}
	MGBox ubox(dim, intrvl);
	delete[] intrvl;
	return ubox;
}

//Compute minimum box sorrounding the all the points.
void MGBPointSeq::compute_box(MGBox& bx)const{
	bx.set_null();
	if(is_null()) return;

	int dim=sdim();
	bx = MGBox(dim);
	int len=length();
	for(int j=0; j<dim; j++){
		double dmax,dmin; dmax=dmin=ref(0,j);
		for(int i=1; i<len; i++){
			double dataij=ref(i,j);
			dmax= dmax>=dataij ? dmax:dataij;
			dmin= dmin<=dataij ? dmin:dataij;
		}
		bx[j]=MGInterval(dmin,dmax);
	}
}

//Exchange ordering of the coordinates.
//Exchange coordinates (j1) and (j2).
void MGBPointSeq::coordinate_exchange(int j1, int j2){
	assert(j1<sdim() && j2<sdim());
	double w;
	for(int i=0; i<m_length; i++){
		w=ref(i,j1); (*this)(i,j1)=ref(i,j2); (*this)(i,j2)=w;
	}
}

//Returns a pointer to the area.
double* MGBPointSeq::data(int i, int j){
	assert(i<capacity()&&j<m_sdim);
	return &m_bpoint[i+m_capacity*j];
}

//Returns a pointer to the area.
const double* MGBPointSeq::data(int i, int j)const{
	assert(i<capacity()&&j<m_sdim);
	return &m_bpoint[i+m_capacity*j];
}

//Insert vector data to this(i,j)  0<=j<sdim().
void MGBPointSeq::insert_at(
	int i,	//id of this which indicates the placement.
	const MGVector& vector)	//Input vector.
{
	assert(i<=length());
	int oldlen=length();
	int newlen=oldlen+1;
	reshape(newlen);

	int j, m, ii = i;
	int sd=sdim();
	for(j=oldlen-1; j>=ii; j--)
		for(m=0; m<sd; m++)
			m_bpoint[j+1+m_capacity*m]=m_bpoint[j+m_capacity*m];
	for(m=0; m<sd; m++)
		m_bpoint[i+m_capacity*m]=vector.ref(m);
}


//Insert array data[j] to this(i,j)  0<=j<sdim().
void MGBPointSeq::insert_at(
	int i,				//id of this which indicates the placement.
	const double* data	//Input data array.
){
	assert(i<=length());
	int oldlen = length();
	int newlen = oldlen+1;
	reshape(newlen);

	int j, m, ii = i;
	int sd = sdim();
	for(j = oldlen-1; j>=ii; j--)
		for(m = 0; m<sd; m++)
			m_bpoint[j+1+m_capacity*m] = m_bpoint[j+m_capacity*m];
	for(m = 0; m<sd; m++)
		m_bpoint[i+m_capacity*m] = data[m];
}

//Transformation for rational(MGRLBRep) Control Polygon.
//When rational control polygon, coordinates are of homogeneous form,
//i.e., the last space dimension id is for weights,
//and other elements include weight multiplied.
MGBPointSeq& MGBPointSeq::homogeneous_transform(double scale){
//1. Scaling.
	--m_sdim; operator*=(scale); ++m_sdim;
	return *this;
}

//2. Vector addition.
MGBPointSeq& MGBPointSeq::homogeneous_transform(const MGVector& vec){
	int i,j;
	int dim1=sdim()-1, dim2=vec.sdim();
	int dim= dim1>=dim2 ? dim1:dim2;
	double* v=new double[dim];
	for(j=0; j<dim; j++) v[j]=vec.ref(j);
	int len=length();
	double weight;
	if(dim>dim1){
		MGBPointSeq bnew(len, dim+1);
		for(i=0; i<len; i++){
			weight=ref(i,dim1);
			for(j=0; j<dim1; j++)
				bnew(i,j)=ref(i,j)+v[j]*weight;
			for(j=dim1; j<dim; j++)
				bnew(i,j)=v[j]*weight;
			bnew(i,dim)=weight;
		}
		*this=std::move(bnew);
	}else{
		for(i=0; i<len; i++){
			weight=ref(i,dim1);
			for(j=0; j<dim2; j++)
				(*this)(i,j)+=v[j]*weight;
		}
	}
	delete[] v;
	return *this;
}

//3. Matrix multiplication.
MGBPointSeq& MGBPointSeq::homogeneous_transform(const MGMatrix& mat){
	int dim1=sdim()-1; int dim2=mat.sdim();
	int dim= dim1>=dim2 ? dim1:dim2;
	int i,j;	double a;
	int len=length();
	if(dim>dim1){
		MGBPointSeq bnew(len, dim+1);
		for(int k=0; k<len; k++){
			for(i=0; i<dim; i++){
				a=0.;
				for(j=0; j<dim1; j++) a+=ref(k,j)*mat.ref(j,i);
				bnew(k,i)=a;
			}
			bnew(k,dim)=ref(k,dim1);
		}
		*this=std::move(bnew);
	}else{
		double* v=new double[dim];
		for(int k=0; k<len; k++){
			for(i=0; i<dim; i++){
				a=0.;
				for(j=0; j<dim; j++)
					a+=ref(k,j)*mat.ref(j,i);
				v[i]=a;
			}
			for(i=0; i<dim; i++)
				(*this)(k,i)=v[i];
		}
		delete[] v;
	}

	return *this;
}

//4. Transformation multiplication.
MGBPointSeq& MGBPointSeq::homogeneous_transform(const MGTransf& tr){
	homogeneous_transform(tr.affine());
	homogeneous_transform(tr.translation());
	return *this;
}

///Negate this point sequence.
///That is, change the ordering of this(i,.) for i.
void MGBPointSeq::negate(){
	int sd=sdim(), n=length();
	int nHalf=n/2, nm1=n-1;
	for(int j=0; j<sd; j++){
		int id=m_capacity*j;
		int nm1pid=nm1+id;
		for(int i=0; i<nHalf; i++){
			double& m_bpointipid=m_bpoint[i+id];
			double& m_bpointnm1pidmi=m_bpoint[nm1pid-i];
			double save=m_bpointipid;
			m_bpointipid=m_bpointnm1pidmi;
			m_bpointnm1pidmi=save;
		}
	}
}

//Compute non_homogeneous coordonate data without w coordinate element,
//assumed that this is homogeneous coordinate data,
//i.e., maximum space dimension element is w(weight) coordinate.
MGBPointSeq MGBPointSeq::non_homogeneous() const{
	int i,j, sd=sdim()-1, n=length();
	MGBPointSeq cp(n,sd);
	for(i=0; i<n; i++){
		double weight=(*this)(i,sd);
		for(j=0; j<sd; j++)
			cp(i,j)=(*this)(i,j)/weight;
	}
	return cp;
}

//Check if the line B-rep is planar.
//Funtion's return value is;
// 0: Not planar, nor a point, nor straight line.
// 1: B-Rep is a point.		2: B-Rep is a straight line.
// 3: B-Rep is planar.
int MGBPointSeq::planar(
	MGPlane& plane			//When Brep is not straight line nor a point,
							// plane is returned.
							//Even when not planar, plane nearest is returned.
	, MGStraight& line		//When Brep is a line, line is returned.
	, MGPosition& center	//Center of the B-Rep is always returned.
)const{
	double deviation; int kplane;
	kplane=average_plane(center, plane, line, deviation);
	if(kplane==3 && deviation>MGTolerance::wc_zero())
		kplane=0;
	return kplane;
}

//Retrieve sub data of i-th point of the BPointSeq.
//That is, P(k)=(*this)(i,j+k) for k=0, ..., sd-1.
void MGBPointSeq::point(int i, int j, int sd, MGPosition& P)const{
	P.resize(sd);
	int nd=sdim();
	for(int k=0; k<sd; k++){
		int id=j+k; if(id>=nd) id-=nd;
		P(k)=ref(i,id);
	}
}

double MGBPointSeq::ref(int i, int j) const{
//Return (i,j)-th element data.
// When j>=sdim(), return 0.0  .
	assert(i< length());
	if(j>=sdim())
		return 0.;
	else
		return m_bpoint[i+m_capacity*j];
}

///Change the capacity to at least sz and update the effective length to sz,
///does not update sdim. After invoking reshape(), at least sz point data are guaranteeed 
///to store.
///reshape() guarantees the original data BPoint(i,.) before invoking reshape, which
///will be stored in the new BPoint(start+i,.).
void MGBPointSeq::reshape(int sz, int start){
	int new_length = m_length+start;
	m_length = sz;
	int capa = capacity();
	if(sz==capa && start==0) return;

	if(sz<new_length) new_length = sz;
	int n = new_length-start; //n is the number of data to move.

	//Reshape of m_bpoint.
	double* data = capa<sz ? new double[sz*m_sdim] : m_bpoint;
	int i, j; int jstart_new, jstart_old;
	if(start || data!=m_bpoint){
		for(j = 0; j<m_sdim; j++){
			jstart_new = start+sz*j;
			jstart_old = m_capacity*j;
			for(i = 0; i<n; i++)
				data[i+jstart_new] = m_bpoint[i+jstart_old];
		}
	}
	if(data!=m_bpoint){
		delete[] m_bpoint;
		m_bpoint = data;
		m_capacity = sz;
	}
}
											  
//Resize the array. The result will contain garbages.
//m_capacity,  m_sdim, and m_length will be set as
//m_capacity=sz, m_sdim=dim, m_length=sz.
//dim=0 means space dimension be set same as the old one.
void MGBPointSeq::resize(int sz, int dim){
	if(dim<=0)
		dim=m_sdim;
	int len=sz*dim;
	if(m_capacity*m_sdim != len){
		if(m_bpoint)
			delete[] m_bpoint;
		m_bpoint=new double[len];
	}
	m_capacity=m_length=sz; m_sdim=dim; 
}

//Set the length of effective data.
void MGBPointSeq::set_length(int length){
	assert(length<=m_capacity);
	m_length=length;
}

//Set this BPointSeq as a null.
void MGBPointSeq::set_null(){
	if(m_bpoint)
		delete[] m_bpoint;
	m_bpoint=0;
	m_capacity=m_length=m_sdim=0;
}

//Store vector data vector(from+j) to this(i,to+j) for 0<=j<sdim().
//When (form+j) or (to+j) reached to maximum space dimension id, next id
//becomes 0(form the start).
void MGBPointSeq::store_at(
	int i,				//id of this which indicates the placement.
	const MGVector& vector,	//Input vector.
	int to,		//Indicates to where of this in the space dimension id.
	int from	//Indicates from where of vector in the space dimension.
){
	assert(i<capacity());
	int sd1=sdim(), sd2=vector.sdim();
	int len=sd1;
	if(len>sd2)
		len=sd2;
	for(int j=0; j<len; j++){
		if(to>=sd1)
			to=0;
		if(from>=sd2)
			from=0;
		m_bpoint[i+m_capacity*(to++)]=vector.ref(from++);
	}
}

//Store vector data vector(from+j) to this(i,to+j) for 0<=j<len.
//When (form+j) or (to+j) reached to maximum space dimension id, next id
//becomes 0(form the start).
void MGBPointSeq::store_at(
	int i,	//id of this which indicates the placement.
	const MGVector& vector,	//Input vector.
	int to,	//Indicates to where of this in the space dimension id.
	int from,	//Indicates from where of vector in the space dimension.
	int len		//length of the data to store.
){
	assert(i<capacity());
	int sd1=sdim(), sd2=vector.sdim();
	if(len>sd1)
		len=sd1;
	for(int j=0; j<len; j++){
		if(to>=sd1)
			to=0; 
		if(from>=sd2)
			from=0;
		m_bpoint[i+m_capacity*(to++)]=vector.ref(from++);
	}
}

//Store data[j] to this(i,to+j) for 0<=j<sdim().
//When (to+j) reached to maximum space dimension id, next id
//becomes 0(form the start).
void MGBPointSeq::store_at(
	int i,				//id of this which indicates the placement.
	const double* data,		//Input data array.
	int to			//Indicates to where of this in the space dimension id.
){
	assert(i<capacity());
	int sd1=sdim();
	for(int j=0; j<sd1; j++){
		if(to>=sd1)
			to=0;
		m_bpoint[i+m_capacity*(to++)]=data[j];
	}
}

//Operator Definition

//Extract (i,j) coordinate values for 0<=j<sdim().
MGVector MGBPointSeq::operator()(int i)const{
	assert(i<length());

	MGVector v(m_sdim);
	for(int j=0; j<m_sdim; j++)
		v(j)=ref(i,j);
	return v;
}

// 曲線の平行移動を行いオブジェクトを生成する。
MGBPointSeq MGBPointSeq::operator+(const MGVector& v) const{
	MGBPointSeq bnew(*this);
	return bnew+=v;
}
MGBPointSeq operator+(const MGVector& v, const MGBPointSeq& b){
	return b+v;
}

// 与ベクトルだけ曲線を平行移動して自身とする。
MGBPointSeq& MGBPointSeq::operator+= (const MGVector& vec){
	int i,j; int dim1=sdim(), dim2=vec.sdim();
	int dim= dim1>=dim2 ? dim1:dim2;
	double* v=new double[dim];
	for(j=0; j<dim; j++) v[j]=vec.ref(j);
	int len=length();
	if(dim>dim1){
		MGBPointSeq bnew(len, dim);
		for(j=0; j<dim; j++)
			for(i=0; i<len; i++)
				bnew(i,j)=ref(i,j)+v[j];
		*this=std::move(bnew);
	}else{
		for(j=0; j<dim; j++){
			for(i=0; i<len; i++)
				(*this)(i,j)+=v[j];
		}
	}
	delete[] v;
	return *this;
}

// 曲線の逆方向に平行移動を行いオブジェクトを生成する。
MGBPointSeq MGBPointSeq::operator- (const MGVector& vec) const{
	MGVector v=-vec;
	return (*this)+v;
}

// 与ベクトルだけ曲線をマイナス方向に平行移動して自身とする。
MGBPointSeq& MGBPointSeq::operator-= (const MGVector& vec){
	MGVector v=-vec;
	return *this += v;
}

//Add and subtract operation of two BPointSeq.
MGBPointSeq MGBPointSeq::operator+ (const MGBPointSeq& bp2) const{
	MGBPointSeq bp(*this);
	return bp+=bp2;
}

MGBPointSeq& MGBPointSeq::operator+= (const MGBPointSeq& bp2){
	int sd1=sdim(), sd2=bp2.sdim();
	int sd = (sd1>=sd2) ? sd1:sd2;
	int n1=length(), n2=bp2.length();
	int n = n2>n1 ?  n2:n1;
	if(sd>sd1){
		MGBPointSeq bp;
		if(n>n1){
			bp=MGBPointSeq(sd, bp2);
			for(int i=0; i<n1; i++)
				bp.store_at(i,bp(i)+(*this)(i));
		}else{
			bp=MGBPointSeq(sd, *this);
			for(int i=0; i<n2; i++)
				bp.store_at(i,bp(i)+bp2(i));
		}
		(*this)=std::move(bp);
	}else{
		reshape(n);
		int i;
		for(i=n1; i<n; i++)
			store_at(i,MGVector(0.,0.));//Clear
		for(i=0; i<n2; i++)
			store_at(i,(*this)(i)+bp2(i));
	}
	return *this;
}

MGBPointSeq MGBPointSeq::operator- (const MGBPointSeq& bp2) const{
	return (*this)+(bp2*(-1.));
}
MGBPointSeq& MGBPointSeq::operator-= (const MGBPointSeq& bp2){
	return (*this)+=(bp2*-1.);
}

// 与えられたスケールで曲線の変換を行いオブジェクトを生成する。
//Generates an object by multiplying scale to the original.
MGBPointSeq MGBPointSeq::operator* (double scale) const{
	MGBPointSeq bp(m_length,m_sdim);
	for(int j=0; j<m_sdim; j++){
		for(int i=0; i<m_length; i++) bp(i,j)=(*this)(i,j)*scale;
	}
	return bp;
}

//BPointをスケーリングしてできるオブジェクトを生成する。
//Generates a object by scaling.
MGBPointSeq operator* (double scale, const MGBPointSeq& bp){
	return bp*scale;
}

// 与えられたスケールで曲線の変換を行い自身の曲線とする。
//Updates the object by multiplying scale.
MGBPointSeq& MGBPointSeq::operator*= (double scale){
	for(int j=0; j<m_sdim; j++){
		for(int i=0; i<m_length; i++)
			(*this)(i,j)*=scale;
	}
	return *this;
}

// 与えられたスケールで曲線の変換を行いオブジェクトを生成する。
//Generates an object by multiplying scale to the original.
MGBPointSeq MGBPointSeq::operator/ (double scale) const{
	return (*this)*(1./scale);
}

// 与えられたスケールで曲線の変換を行い自身の曲線とする。
//Updates the object by multiplying scale.
MGBPointSeq& MGBPointSeq::operator/= (double scale){
	return (*this)*=(1./scale);
}

// 与えられた変換で曲線の変換を行いオブジェクトを生成する。
MGBPointSeq MGBPointSeq::operator* (const MGMatrix& mat) const{
	MGBPointSeq bnew(*this);
	return bnew*=mat;
}

// 与えられた変換で曲線の変換を行い自身の曲線とする。
MGBPointSeq& MGBPointSeq::operator*= (const MGMatrix& mat){
	int dim1=sdim(); int dim2=mat.sdim();
	int dim= dim1>=dim2 ? dim1:dim2;
	int i,j;	double a;
	int len=length();
	if(dim>dim1){
		MGBPointSeq bnew(len, dim);
		for(int k=0; k<len; k++){
			for(i=0; i<dim; i++){
				a=0.;
				for(j=0; j<dim; j++)
					a=a+ref(k,j)*mat.ref(j,i);
				bnew(k,i)=a;
			}
		}
		bnew.set_length(len);
		*this=std::move(bnew);
	}else{
		double* v=new double[dim];
		for(int k=0; k<len; k++){
			for(i=0; i<dim; i++){
				a=0.;
				for(j=0; j<dim; j++)
					a=a+ref(k,j)*mat.ref(j,i);
				v[i]=a;
			}
			for(i=0; i<dim; i++)
				(*this)(k,i)=v[i];
		}
		delete[] v;
	}

	return *this;
}

// 与えられた変換で曲線のトランスフォームを行いオブジェクトを生成する。
MGBPointSeq MGBPointSeq::operator* (const MGTransf& tr) const{
	return (*this)*tr.affine()+tr.translation();
}

// 与えられた変換で曲線のトランスフォームを行い自身とする。
MGBPointSeq& MGBPointSeq::operator*= (const MGTransf& tr){
	return ((*this)*=tr.affine())+=tr.translation();
}
										  
//Assignment
MGBPointSeq& MGBPointSeq::operator= (const MGBPointSeq& bp2){
	int len=bp2.length(), dim=bp2.sdim();
	int total_len=capacity()*sdim();
	if(dim*len>total_len) resize(len, dim);
	else{
		m_sdim=dim;
		m_capacity=total_len/dim;
	}
	m_length=len;
	for(int j=0; j<m_sdim; j++)
		for(int i=0; i<m_length; i++)
			(*this)(i,j)=bp2(i,j);
	return *this;
}
										  
//Moove Assignment
MGBPointSeq& MGBPointSeq::operator= (MGBPointSeq&& bp2){
	if(this!=&bp2){
		if(m_bpoint)
			delete[] m_bpoint;
		m_bpoint=bp2.m_bpoint;
		m_capacity=bp2.m_capacity;
		m_sdim=bp2.m_sdim;
		m_length=bp2.m_length;

		bp2.m_bpoint=nullptr;
		bp2.m_sdim=bp2.m_length=bp2.m_capacity=0;
	}
	return *this;
}

//Compare two BPointSeq if they are equal.
bool MGBPointSeq::operator== (const MGBPointSeq& bpoint) const{
	int len=length();
	if(len!=bpoint.length())
		return false;
	if(len<=0)
		return true;

	int dim=sdim(), dim2=bpoint.sdim();
	if(dim<dim2)
		dim=dim2;
	int i,j; double a,b; 
	double error=MGTolerance::wc_zero(); error=error*error;
	for(i=0; i<len; i++){
		a=0.;
		for(j=0; j<dim; j++){
			b=ref(i,j)-bpoint.ref(i,j); a += b*b;
		}
		if(a>error)
			return false;
	}
	return true;
}
