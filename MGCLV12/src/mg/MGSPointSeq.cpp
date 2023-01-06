/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/NDDArray.h"
#include "mg/Box.h"
#include "mg/BPointSeq.h"
#include "mg/SPointSeq.h"
#include "mg/Vector.h"
#include "mg/Matrix.h"
#include "mg/Transf.h"
#include "mg/Straight.h"
#include "mg/Plane.h"
#include "mg/Tolerance.h"

#include "cskernel/bkdtpg.h"
#include "cskernel/Bvavpl.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
// Implementation of Class MGSPointSeq

//Constructor
MGSPointSeq::MGSPointSeq(int sizeu, int sizev, int dim)
:m_capacityu(sizeu),m_capacityv(sizev),m_sdim(dim)
,m_lengthu(sizeu),m_lengthv(sizev){
	int len=sizeu*sizev*dim;
	if(len)
		m_spoint=new double[len];
	else{
		m_sdim=m_capacityu=m_capacityv=m_lengthu=m_lengthv=0;
		m_spoint=nullptr;
	}
}

//Copy constructor.
MGSPointSeq::MGSPointSeq(const MGSPointSeq& rhs)
: m_capacityu(rhs.m_capacityu),m_capacityv(rhs.m_capacityv)
, m_sdim(rhs.m_sdim)
, m_lengthu(rhs.m_lengthu), m_lengthv(rhs.m_lengthv)
,m_spoint(nullptr){
	int len=m_capacityu*m_capacityv*m_sdim;
	if(len){
		MGSPointSeq& SP=*this;
		m_spoint=new double[len];
		for(int k=0; k<m_sdim; k++)
			for(int i=0; i<m_lengthu; i++)
				for(int j=0; j<m_lengthv; j++)
					SP(i, j, k)=rhs.ref(i, j, k);
	}else
		m_sdim=m_capacityu=m_capacityv=m_lengthu=m_lengthv=0;

}

//Move constructor.
MGSPointSeq::MGSPointSeq(MGSPointSeq&& rhs)
: m_capacityu(rhs.m_capacityu),m_capacityv(rhs.m_capacityv)
, m_sdim(rhs.m_sdim)
, m_lengthu(rhs.m_lengthu), m_lengthv(rhs.m_lengthv)
,m_spoint(rhs.m_spoint){
	rhs.m_spoint=nullptr;
	rhs.m_capacityu=rhs.m_capacityv=rhs.m_sdim=rhs.m_lengthu=rhs.m_lengthv=0;
}

//Construct a MGSPointSeq by copying original MGSPointSeq.
//Can change the order of coordinates 
MGSPointSeq::MGSPointSeq(
	int dim,				// New Space Dimension.
	const MGSPointSeq& old,	// Origianl SpointSeq.
	int start1,			// Destination start order to store.
	int start2			// Source start order to retrieve.
): m_capacityu(old.m_lengthu),m_capacityv(old.m_lengthv)
, m_sdim(dim)
, m_lengthu(old.m_lengthu), m_lengthv(old.m_lengthv){
	assert(start1<dim && start2<old.m_sdim);

	int len=m_lengthu*m_lengthv*m_sdim;
	if(len){
		MGSPointSeq& SP=*this;
		m_spoint=new double[len];
		int dim2=old.m_sdim; 
		int dimmin= dim<dim2 ? dim:dim2; 
		int k, k2=start2, k1=start1;
		for(k=0; k<dimmin; k++) {
			for(int i=0; i<m_lengthu; i++)
				for(int j=0; j<m_lengthv; j++) SP(i,j,k1)=old(i,j,k2);
			k1 +=1; if(k1>=dim) k1=0; 
			k2 +=1; if(k2>=dim2) k2=0;
		}
		while(k++<dim){
			for(int i=0; i<m_lengthu; i++)
				for(int j=0; j<m_lengthv; j++) SP(i,j,k1)=0.0;
			k1 +=1; if(k1>=dim) k1=0; 
		}
	}else
		m_spoint=0;
}

//Member Function

//compute an average plane of the point sequence.
//Function's return value is:
// 0: Number of points included is zero.
// 1: Point seq is a point.		2: Point seq is on a line.
// 3: Plane is output.
int MGSPointSeq::average_plane(
	MGPosition& center		// center of point seq will be output.
	, MGPlane& plane		// Plane will be output, when average_plane=3.
	, MGStraight& line		// Straight line will be output            =2.
	, double& deviation)	// Maximum deviation from point, line
	 const{					// , or plane will be output.
	const double error=MGTolerance::wc_zero();	//error allowed.
	int nu=length_u(), nv=length_v(), ncd=sdim(), ipu=capacity_u(), ipv=capacity_v();
	int kplane=0; double g[4], fcenter[3];
	if(ncd){
		bvavpl_(error,nu,nv,data(),ncd,ipu,ipv,&kplane,g,fcenter,&deviation);
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
	}
	return kplane;
}

//Compute minimum box sorrounding the all the points.
MGBox MGSPointSeq::box()const{
	int dim=sdim();
	MGBox ubox(dim);
	if(dim){
		int lenu=length_u(), lenv=length_v();
		for(int k=0; k<dim; k++){
			MGInterval& Ik=ubox(k);
			double max, min; max=min=ref(0, 0, k);
			for(int j=0; j<lenv; j++){
				for(int i=0; i<lenu; i++){
					double dijk=ref(i, j, k);
					max= max>=dijk ? max:dijk;
					min= min<=dijk ? min:dijk;
				}
			}
			Ik=MGInterval(min, max);
		}
	}
	return ubox;
}

//Returns a pointer to the area.
const double* MGSPointSeq::data(int i, int j, int k) const
{return &m_spoint[i+m_capacityu*j+m_capacityu*m_capacityv*k];}

//Returns a pointer to the area.
double* MGSPointSeq::data(int i, int j, int k)
{return &m_spoint[i+m_capacityu*j+m_capacityu*m_capacityv*k];}

//Transformation for rational(MGRLBRep) Control Polygon.
//When rational control polygon, coordinates are of homogeneous form,
//i.e., the last space dimension id is for weights,
//and other elements include weight multiplied.
MGSPointSeq& MGSPointSeq::homogeneous_transform(double scale)
//1. Scaling.
{
	--m_sdim; operator*=(scale); ++m_sdim;
	return *this;
}
//2. Vector addition.
MGSPointSeq& MGSPointSeq::homogeneous_transform(const MGVector& vec)
{
	int i,j,r;
	int dim1=sdim()-1, dim2=vec.sdim();
	int dim= dim1>=dim2 ? dim1:dim2;
	double* v=new double[dim];
	for(j=0; j<dim; j++) v[j]=vec.ref(j);
	int lenu=length_u(),lenv=length_v();
	double weight;
	if(dim>dim1){
		MGSPointSeq bnew(lenu,lenv, dim+1);
		for(i=0; i<lenu; i++){
		for(j=0; j<lenv; j++){
			weight=ref(i,j,dim1);
			for(r=0; r<dim1; r++) bnew(i,j,r)=ref(i,j,r)+v[r]*weight;
			for(r=dim1; r<dim; r++) bnew(i,j,r)=v[r]*weight;
			bnew(i,j,dim)=weight;
		}
		}
		*this=std::move(bnew);
	}else{
		MGSPointSeq& SP=*this;
		for(i=0; i<lenu; i++){
		for(j=0; j<lenv; j++){
			weight=ref(i,j,dim1);
			for(r=0; r<dim2; r++) SP(i,j,r)+=v[r]*weight;
		}
		}
	}
	delete[] v;
	return *this;
}
//3. Matrix multiplication.
MGSPointSeq& MGSPointSeq::homogeneous_transform(const MGMatrix& mat){
	int dim1=sdim()-1; int dim2=mat.sdim();
	int dim= dim1>=dim2 ? dim1:dim2;
	int i,j,r1,r2;	double a;
	int lenu=length_u(),lenv=length_v();
	if(dim>dim1){
		MGSPointSeq bnew(lenu,lenv, dim+1);
		for(i=0; i<lenu; i++){
		for(j=0; j<lenv; j++){
			for(r1=0; r1<dim; r1++){
				a=0.;
				for(r2=0; r2<dim1; r2++) a+=ref(i,j,r2)*mat.ref(r2,r1);
				bnew(i,j,r1)=a;
			}
			bnew(i,j,dim)=ref(i,j,dim1);
		}
		}
		*this=std::move(bnew);
	}else{
		MGSPointSeq& SP=*this;
		double* v=new double[dim];
		for(i=0; i<lenu; i++){
		for(j=0; j<lenv; j++){
			for(r1=0; r1<dim; r1++){
				a=0.;
				for(r2=0; r2<dim; r2++) a+=ref(i,j,r2)*mat.ref(r2,r1);
				v[r1]=a;
			}
			for(r1=0; r1<dim; r1++) SP(i,j,r1)=v[r1];
		}
		}
		delete[] v;
	}

	return *this;
}
//4. Transformation multiplication.
MGSPointSeq& MGSPointSeq::homogeneous_transform(const MGTransf& tr)
{
	homogeneous_transform(tr.affine());
	homogeneous_transform(tr.translation());
	return *this;
}
//Generate data point abscissa from data point ordinates SPointSeq.
//SPointSeq(*this) must be homogeneous, i.e. if SPointSeq is positional
//data, all of them must be positional data, should not include
//derivative data.
void MGSPointSeq::make_data_point(MGNDDArray& utau, MGNDDArray& vtau) const{
	int i,j,k, lenu,lenv, ipu,ipv,ipuv, mid; 
	length(lenu, lenv); const int ncd=sdim();
	if((ncd*lenu*lenv)==0){
		utau=vtau=MGNDDArray();
		return;
	}
	double d1,d2, d;
 	capacity(ipu, ipv); ipuv=ipu*ipv;
	MGNDDArray tau1(lenu), tau2(lenu), tau3(lenu);

	//Compute u data points.
	mid=lenv/2;
	bkdtpg_(data(0,0,0), lenu, ncd, ipuv, &tau1(0));
	bkdtpg_(data(0,mid,0), lenu, ncd, ipuv, &tau2(0));
	bkdtpg_(data(0,lenv-1,0), lenu, ncd, ipuv, &tau3(0));
	utau.reshape(lenu);
	for(i=0; i<lenu; i++)
		utau(i)=(tau1(i)+tau2(i)+tau3(i))/3.;

	//Compute v data points and knot vector.
	vtau.reshape(lenv);
	mid=lenu/2;
	int id[3]={0, mid, lenu-1};
	vtau(0)=0.;
	for(j=1; j<lenv; j++){
		d2=0.0;
		for(int m=0; m<3; m++){
			d1=0.;
			for(k=0; k<ncd; k++){
				d = ref(id[m],j,k)-ref(id[m],j-1,k);
				d1+=d*d;
			}
			d2 += sqrt(d1);
		}
		vtau(j)=vtau(j-1)+d2/3.;
	}
}

//Compute non_homogeneous coordonate data without w coordinate element,
//assumed that this is homogeneous coordinate data,
//i.e., maximum space dimension element is w(weight) coordinate.
//Result data does not include weight elements.
MGSPointSeq MGSPointSeq::non_homogeneous() const{
	int i,j,r, sd=sdim()-1, m=length_u(), n=length_v();
	MGSPointSeq cp(m,n,sd);
	for(i=0; i<m; i++){
	for(j=0; j<n; j++){
		double weight=ref(i,j,sd);
		for(r=0; r<sd; r++) cp(i,j,r)=ref(i,j,r)/weight;
	}
	}
	return cp;
}

double MGSPointSeq::ref(int i, int j,int k) const{
//Return (i,j,k)-th element data.
// When k>=sdim(), return 0.0  .
	assert(i<length_u() && j<length_v());
	if(k>=sdim())
		return 0.;
	else
		return m_spoint[i+m_capacityu*j+m_capacityu*m_capacityv*k];
}

///Change size. Change of sdim not allowed.
///reshape does update the effective length to sizeu and sizev.
void MGSPointSeq::reshape(
	int sizeu, int sizev
	, int startu, int startv
){
	assert(startu<=sizeu && startv<=sizev);

	int nlenu=m_lengthu+startu; int nlenv=m_lengthv+startv;
	if(sizeu<nlenu) nlenu=sizeu; if(sizev<nlenv) nlenv=sizev;
	int nu=nlenu-startu, nv=nlenv-startv;
	//nu,nv are the numbers of data to move.
	
	//Reshape of m_spoint.
	int i,j,k;
	int sizeuv=sizeu*sizev;
	double* data=new double[sizeuv*m_sdim];
	for(k=0; k<m_sdim; k++){
		int sizeuvk=sizeuv*k;
		for(j=0; j<nv; j++){
			int js=startu+(j+startv)*sizeu+sizeuvk;
			for(i=0; i<nu; i++) data[i+js]=ref(i,j,k);
		}
	}
	delete[] m_spoint;
	m_spoint=data;
	m_capacityu= m_lengthu = sizeu; m_capacityv=m_lengthv=sizev;
}
	
//Change the size of the array to
//m_lengthu=m_capacityu=lenu, m_lengthv=m_capacityv=lenv, m_sdim=dim.
void MGSPointSeq::resize(int lenu, int lenv,  int dim){
	if(m_spoint) delete[] m_spoint;
	m_spoint=new double[lenu*lenv*dim];
	m_capacityu=m_lengthu=lenu;
	m_capacityv=m_lengthv=lenv;
	m_sdim=dim;
}

//Reverse the ordering of points.
void MGSPointSeq::reverse(
	int is_u)				//if true, u-drection. if not, v.
{
	int i,i2,j,j2,k; double save;
	int ncd=sdim(), lenu=length_u(), lenv=length_v();
	int half;
	MGSPointSeq& SP=*this;
	if(is_u){
		half=lenu/2;
		for(k=0; k<ncd; k++){
			for(j=0; j<lenv; j++){
				i2=lenu-1;
				for(i=0; i<half; i++){
					save=ref(i2,j,k);
					SP(i2--,j,k)=ref(i,j,k);
					SP(i,j,k)=save;
				}
			}
		}
	}else{
		half=lenv/2;
		for(k=0; k<ncd; k++){
			for(i=0; i<lenu; i++){
				j2=lenv-1;
				for(j=0; j<half; j++){
					save=ref(i,j2,k);
					SP(i,j2--,k)=ref(i,j,k);
					SP(i,j,k)=save;
				}
			}
		}
	}
}

//Set the length of effective data.
void MGSPointSeq::set_length(
	int lengthu,
	int lengthv
){
	assert(lengthu<=m_capacityu && lengthv<=m_capacityv);
	m_lengthu=lengthu; m_lengthv=lengthv;
}

//Set this as a null.
void MGSPointSeq::set_null(){
	if(m_spoint) delete[] m_spoint;
	m_spoint=0;
	m_lengthu=m_lengthv=m_sdim=m_capacityu=m_capacityv=0;
}

//Store vector data vector(from+r) to this(i,j,to+r) for 0<=r<sdim().
//When (form+r) or (to+r) reached to maximum space dimension id, next id
//becomes 0(form the start).
void MGSPointSeq::store_at(
	int i, int j,	//id of this which indicates the placement.
	const MGVector& vector,	//Input vector.
	int to,		//Indicates to where of this in the space dimension id.
	int from)	//Indicates from where of vector in the space dimension.
{
	assert(i<m_lengthu && j<m_lengthv);
	int sd1=sdim(), sd2=vector.sdim();
	int len=sd1; if(len>sd2) len=sd2;
	for(int r=0; r<len; r++){
		if(to>=sd1) to=0; if(from>=sd2) from=0;
		m_spoint[i+m_capacityu*j+m_capacityu*m_capacityv*(to++)]=vector.ref(from++);
	}
}

//Store vector data vector(from+r) to this(i,j,to+r) for 0<=r<sdim().
//When (form+r) or (to+r) reached to maximum space dimension id, next id
//becomes 0(form the start).
void MGSPointSeq::store_at(
	int i, int j,	//id of this which indicates the placement.
	const MGVector& vector,	//Input vector.
	int to,		//Indicates to where of this in the space dimension id.
	int from,	//Indicates from where of vector in the space dimension.
	int len)		//Length of data to store.
{
	assert(i<m_lengthu && j<m_lengthv);
	int sd1=sdim(), sd2=vector.sdim();
	if(len>sd1) len=sd1;
	for(int r=0; r<len; r++){
		if(to>=sd1) to=0; if(from>=sd2) from=0;
		m_spoint[i+m_capacityu*j+m_capacityu*m_capacityv*(to++)]=vector.ref(from++);
	}
}

//Store data[r] to this(i,j,to+r)  0<=r<sdim().
//When (to+r) reached to maximum space dimension id, next id
//becomes 0(form the start).
void MGSPointSeq::store_at(
	int i, int j,		//id of this which indicates the placement.
	const double* data,		//Input data.
	int to)		//Indicates to where of this in the space dimension id.
{
	assert(i<m_lengthu && j<m_lengthv);
	int sd1=sdim();
	for(int r=0; r<m_sdim; r++){
		if(to>=sd1) to=0;
		m_spoint[i+m_capacityu*j+m_capacityu*m_capacityv*(to++)]=data[r];
	}
}

//Store BPointSeq along u at v's id j. That is, 
//data bp(i,from+r) to this(i,j,to+r) for i=0,..., length_u()-1, and  0<=r<sdim().
//When (form+r) or (to+r) reached to maximum space dimension id, next id
//becomes 0(form the start).
void MGSPointSeq::store_BP_along_u_at(
	int j,	//id of this which indicates the placement.
	const MGBPointSeq& bp,	//Input vector.
	int to,	//Indicates to where of this in the space dimension id.
	int from)	//Indicates from where of vector in the space dimension.
{
	assert(j<m_lengthv);
	int sd1=sdim(), sd2=bp.sdim();
	int len=sd1; if(len>sd2) len=sd2;
	int szubyj=m_capacityu*j, szubyszv=m_capacityu*m_capacityv;
	for(int r=0; r<len; r++, to++, from++){
		if(to>=sd1) to=0; if(from>=sd2) from=0;
		int juvto=szubyj+szubyszv*to;
		for(int i=0; i<m_lengthu; i++) m_spoint[i+juvto]=bp(i,from);
	}
}

//Store BPointSeq along v at u's id i. That is, 
//data bp(j,from+r) to this(i,j,to+r) for j=0,..., length_v()-1, and  0<=r<sdim().
//When (form+r) or (to+r) reached to maximum space dimension id, next id
//becomes 0(form the start).
void MGSPointSeq::store_BP_along_v_at(
	int i,	//id of this which indicates the placement.
	const MGBPointSeq& bp,	//Input vector.
	int to,	//Indicates to where of this in the space dimension id.
	int from)	//Indicates from where of vector in the space dimension.
{
	assert(i<m_lengthu);
	int sd1=sdim(), sd2=bp.sdim();
	int len=sd1; if(len>sd2) len=sd2;
	int szubyszv=m_capacityu*m_capacityv;
	for(int r=0; r<len; r++, to++, from++){
		if(to>=sd1) to=0; if(from>=sd2) from=0;
		int iuvto=i+szubyszv*to;
		for(int j=0; j<m_lengthv; j++)
			m_spoint[iuvto+m_capacityu*j]=bp(j,from);
	}
}

//Operator Definition

//Copy Assignment
MGSPointSeq& MGSPointSeq::operator=(const MGSPointSeq& sp2){
	if(this!=&sp2){
		int lenu=sp2.length_u(), lenv=sp2.length_v(), dim=sp2.sdim();
		int len=capacity_u()*capacity_v()*sdim();
		if(len<lenu*lenv*dim) resize(lenu, lenv, dim);
		else{
			m_sdim=dim;
			m_lengthu=m_capacityu=lenu; m_lengthv=lenv;
			m_capacityv=len/(dim*lenu);
		}
		MGSPointSeq& SP=*this;
		for(int m=0; m<m_sdim; m++)
			for(int i=0; i<lenu; i++)
				for(int j=0; j<lenv; j++)
					SP(i, j, m)=sp2(i, j, m);
	}
	return *this;
}
//Move Assignment
MGSPointSeq& MGSPointSeq::operator=(MGSPointSeq&& rhs){
	if(m_spoint)
		delete m_spoint;
	m_spoint=rhs.m_spoint;
	rhs.m_spoint=nullptr;

	m_capacityu=rhs.m_capacityu;
	m_capacityv=rhs.m_capacityv;
	m_sdim=rhs.m_sdim;
	m_lengthu=rhs.m_lengthu;
	m_lengthv=rhs.m_lengthv;
	rhs.m_capacityu=rhs.m_capacityv=rhs.m_sdim=rhs.m_lengthu=rhs.m_lengthv=0;
	return *this;
}

double& MGSPointSeq::operator()(int i, int j, int k) 
//Access to (i,j)th element
{
	assert(i<m_capacityu && j<m_capacityv && k<m_sdim);
	return m_spoint[i+m_capacityu*j+m_capacityu*m_capacityv*k];
}

MGVector MGSPointSeq::operator()(int i, int j) const
//Extract (i,j,k) elements for 0<=k<sdim() as a vector.
{
	MGVector v(m_sdim);
	for(int k=0; k<m_sdim; k++) v(k)=ref(i,j,k);
	return v;
}

// 曲線の平行移動を行いオブジェクトを生成する。
MGSPointSeq MGSPointSeq::operator + ( const MGVector& vec) const{
	MGSPointSeq bnew(*this);
	return bnew+=vec;
}

// 与ベクトルだけ曲線を平行移動して自身とする。
MGSPointSeq& MGSPointSeq::operator+= (const MGVector& vec){
	int i,j,k; int dim1=sdim(), dim2=vec.sdim();
	int dim= dim1>=dim2 ? dim1:dim2;
	int lenu=length_u(), lenv=length_v();
	if(dim>dim1){
		MGSPointSeq bnew(lenu, lenv, dim);
		for(k=0; k<dim; k++){
			double vk=vec[k];
			for(i=0; i<lenu; i++)
				for(j=0; j<lenv; j++) bnew(i, j, k)=ref(i, j, k)+vk;
		}	
		bnew.set_length(lenu,lenv);
		*this=std::move(bnew);
	}else{
		MGSPointSeq& SP=*this;
		for(k=0; k<dim; k++)
			for(i=0; i<lenu; i++)  
				for(j=0; j<lenv; j++) SP(i,j,k)+=vec[k];
	}
	return *this;
}

// 曲線の逆方向に平行移動を行いオブジェクトを生成する。
MGSPointSeq MGSPointSeq::operator- (const MGVector& vec) const{
	MGVector v=-vec;
	return (*this)+v;
}

// 与ベクトルだけ曲線をマイナス方向に平行移動して自身とする。
MGSPointSeq& MGSPointSeq::operator-= (const MGVector& vec){
	MGVector v=-vec;
	return *this += v;
}

// 与えられたスケーリングで曲線の変換を行いオブジェクトを生成する。
//Scaling.
MGSPointSeq MGSPointSeq::operator* (double scale) const{
	MGSPointSeq sp(*this);
	sp*=scale;
	return sp;
}

// 与えられたスケーリングで曲線の変換を行いオブジェクトを生成する。
//Scaling.
MGSPointSeq operator* (double scale, const MGSPointSeq& sp){
	return sp*scale;
}

// 与えられたスケーリングで曲線の変換を行い自身の曲線とする。
//Scaling.
MGSPointSeq& MGSPointSeq::operator*= (double scale){
	MGSPointSeq& SP=*this;
	for(int k=0; k<sdim(); k++)
		for(int i=0; i<length_u(); i++)  
			for(int j=0; j<length_v(); j++) SP(i,j,k)*=scale;
	return *this;
}

// 与えられた変換で曲線の変換を行いオブジェクトを生成する。
MGSPointSeq MGSPointSeq::operator* (const MGMatrix& mat) const{
	MGSPointSeq bnew(*this);
	return bnew*=mat;
}

// 与えられた変換で曲線の変換を行い自身の曲線とする。
MGSPointSeq& MGSPointSeq::operator*= (const MGMatrix& mat){
	int dim1=sdim(); int dim2=mat.sdim();
	int dim= dim1>=dim2 ? dim1:dim2;
	int i,j,k,k2;	double a;
	int lenu=length_u(), lenv=length_v();
	if(dim>dim1){
		MGSPointSeq bnew(lenu, lenv, dim);
		for(i=0; i<lenu; i++){
		for(j=0; j<lenv; j++){
			for(k=0; k<dim; k++){
				a=0.;
				for(k2=0; k2<dim; k2++) a=a+ref(i,j,k2)*mat.ref(k2,k);
				bnew(i,j,k)=a;
			}
		}
		}
		bnew.set_length(lenu,lenv);
		*this=std::move(bnew);
	}else{
		MGVector v(dim);
		MGSPointSeq& SP=*this;
		for(i=0; i<lenu; i++){
		for(j=0; j<lenv; j++){
			for(k=0; k<dim; k++){
				a=0.;
				for(k2=0; k2<dim; k2++) a=a+ref(i,j,k2)*mat.ref(k2,k);
				v(k)=a;
			}
			for(k=0; k<dim; k++) SP(i,j,k)=v[k];
		}
		}
	}

	return *this;
}

// 与えられた変換で曲線のトランスフォームを行いオブジェクトを生成する。
MGSPointSeq MGSPointSeq::operator* (const MGTransf& tr) const{
	return (*this)*tr.affine()+tr.translation();
}

// 与えられた変換で曲線のトランスフォームを行い自身とする。
MGSPointSeq& MGSPointSeq::operator*= (const MGTransf& tr){
	return ((*this)*=tr.affine())+=tr.translation();
}

//Compare two SPointSeq if they are equal.
bool MGSPointSeq::operator== (const MGSPointSeq& spoint) const{
	int lenu=length_u(), lenv=length_v();
	if(lenu!=spoint.length_u() || lenv!=spoint.length_v()) return 0;
	if(lenu<=0 && lenv<=0) return 1;

	int dim=sdim(), dim2=spoint.sdim();
	if(dim<dim2) dim=dim2;
	int i,j,k; double a,b; 
	double error=MGTolerance::wc_zero(); error=error*error;
	for(j=0; j<lenv; j++){
	for(i=0; i<lenu; i++){
		a=0.;
		for(k=0; k<dim; k++){
			b=ref(i,j,k)-spoint.ref(i,j,k); a += b*b;
		}
		if(a>error) return 0;
	}
	}
	return 1;
}

//Add and subtract operation of two MGSPointSeq.
MGSPointSeq& MGSPointSeq::operator+= (const MGSPointSeq& sp2){
	int sd=sdim(), sd2=sp2.sdim();
	if(sd2>sd) sd=sd2;
	int nu1=length_u(), nu2=sp2.length_u();
	int nu=nu1;	if(nu2>nu) nu=nu2;
	int nv1=length_v(), nv2=sp2.length_v();
	int nv=nv1;	if(nv2>nv) nv=nv2;
	int i,j;
	const MGVector zero(0.,0.,0.);
	MGSPointSeq& SP=*this;
	if(sd>sdim()){
		MGSPointSeq sp(nu,nv,sd);
		for(i=0; i<nu1; i++)
			for(j=0; j<nv1; j++)
				sp.store_at(i,j,SP(i,j));//copy
		for(i=nu1; i<nu; i++)
			for(j=nv1; j<nv; j++)
				sp.store_at(i,j,zero);//Clear
		for(i=0; i<nu2; i++)
			for(j=0; j<nv2; j++)
				sp.store_at(i,j,sp(i,j)+sp2(i,j));
		SP=std::move(sp);
	}else{
		reshape(nu,nv);
		for(i=nu1; i<nu; i++)
			for(j=nv1; j<nv; j++)
				store_at(i,j,zero);//Clear
		for(i=0; i<nu2; i++)
			for(j=0; j<nv2; j++)
				store_at(i,j,SP(i,j)+sp2(i,j));
	}
	return *this;
}
MGSPointSeq MGSPointSeq::operator+ (const MGSPointSeq& sp2)const{
	MGSPointSeq sp(*this);
	return sp+=sp2;
}
MGSPointSeq MGSPointSeq::operator- (const MGSPointSeq& sp2) const{
	MGSPointSeq sp(*this);
	return sp-=sp2;
}
MGSPointSeq& MGSPointSeq::operator-= (const MGSPointSeq& sp2){
	int sd=sdim(), sd2=sp2.sdim();
	if(sd2>sd) sd=sd2;
	int nu1=length_u(), nu2=sp2.length_u();
	int nu=nu1;	if(nu2>nu) nu=nu2;
	int nv1=length_v(), nv2=sp2.length_v();
	int nv=nv1;	if(nv2>nv) nv=nv2;
	int i,j;
	const MGVector zero(0.,0.,0.);
	MGSPointSeq& SP=*this;
	if(sd>sdim()){
		MGSPointSeq sp(nu,nv,sd);
		for(i=0; i<nu1; i++) for(j=0; j<nv1; j++) sp.store_at(i,j,SP(i,j));//copy
		for(i=nu1; i<nu; i++) for(j=nv1; j<nv; j++) sp.store_at(i,j,zero);//Clear
		for(i=0; i<nu2; i++) for(j=0; j<nv2; j++) sp.store_at(i,j,sp(i,j)-sp2(i,j));
		SP=std::move(sp);
	}else{
		reshape(nu,nv);
		for(i=nu1; i<nu; i++) for(j=nv1; j<nv; j++) store_at(i,j,zero);//Clear
		for(i=0; i<nu2; i++) for(j=0; j<nv2; j++) store_at(i,j,SP(i,j)-sp2(i,j));
	}
	return *this;
}
