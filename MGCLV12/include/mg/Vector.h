/********************************************************************/
/* Copyright (c) 2021 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#pragma once

/** @file */
/** @addtogroup BASE
 *  @{
 */
#include <stddef.h>
#include <vector>
#include "mg/MGCL.h"

// MGVector.h
// Header for MGVector.

//Forward Declaration
class MGUnit_vector;
//class MGPosition;
class MGIfstream;
class MGOfstream;
class MGIgesOfstream;

///Vector of a general n space dimension.
class MG_DLL_DECLR MGVector {
protected:
	/// Protected data member
	int m_sdim;
	double* m_element;///< data pointe, m_data when m_sdim<=3.
	double m_data[3]={0.,0.,0.};///<For vector data of space dimension less or equal to 3.
	mutable double m_length;///<To hold vector length if computed.
						///<When length not computed, negative value will be set.

public:

///Addition of two vectors.
MG_DLL_DECLR friend MGVector operator+(const MGVector& vec1,const MGVector& vec2);

///Subtraction of two vectors.
MG_DLL_DECLR friend MGVector operator-(const MGVector& vec1,const MGVector& vec2);

///Inner product of two vectors.
MG_DLL_DECLR friend double operator%(const MGVector& vec1,const MGVector& vec2);

///vector product of two vectors.
MG_DLL_DECLR friend MGVector operator*(const MGVector& vec1,const MGVector& vec2);

///Scalar multiplication.
MG_DLL_DECLR friend MGVector operator*(const MGVector& vec1,double scale);

///Scalar multiplication.
MG_DLL_DECLR friend MGVector operator* (double, const MGVector&);

///Scalar division.
MG_DLL_DECLR friend MGVector operator/(const MGVector& vec1,double scale);

///Test if this vector is less than v2.
///Comparison depends on two vectors' length.
inline friend
bool operator<(const MGVector& v1,const MGVector& v2){return v1.len()<v2.len();};
inline friend
bool operator<=(const MGVector& v1,const MGVector& v2){return v1.len()<=v2.len();};
inline friend
bool operator>(const MGVector& v1,const MGVector& v2){return v1.len()>v2.len();};
inline friend
bool operator>=(const MGVector& v1,const MGVector& v2){return v1.len()>=v2.len();};

///Test if two vectors are equal.
MG_DLL_DECLR friend bool operator==(const MGVector& v1,const MGVector& v2);

///Test if two vectors are equal.
inline friend bool operator!=(const MGVector& v1,const MGVector& v2){return !(v1==v2);}

///String stream function
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream&, const MGVector&);

///Determinant of 3 by 3 matrix of 3 vectors.
MG_DLL_DECLR friend double MGDeterminant(
		const MGVector& v1, const MGVector& v2, const MGVector& v3
);


////////Special member functions/////////
explicit MGVector(int sdim=0);
~MGVector(){if(m_sdim>3) delete[] m_element;}
MGVector(const MGVector&);///Copy constructor.
MGVector& operator= (const MGVector&);///Copy assignment.
MGVector(MGVector&&);		///Move constructor.
MGVector& operator= (MGVector&&);///Move assignment.

///Construct 2D vector by providing each element data.
MGVector(double x, double y);

///Construct 3D vector by providing each element data.
MGVector(double x, double y , double z);

///Construct 4D vector by providing each element data.
MGVector(double x, double y, double z, double w);

///Vector of same value for each coordinate element.
MGVector(int sdim, double v);

///Vector from array of double v[sdim].
///***** This is the fundamental constructor.*****
MGVector(int sdim, const double* v);

///Construct a vector from a difference of two vectors.
MGVector(
	const MGVector& dvec,///<Destination point
	const MGVector& svec///<Source point
);

///Construct Vector by copying old Vector, changing space dimension and
///ordering of old coordinates.
/// (*this)(start1+i)=vec2(start2+i).
MGVector(
	int sdim,			///<Space dimension.
	const MGVector& vec2,///<Original vector.
	int start1=0,		///<id of constructing vector that indicates
						///<from where to store the elements of vec2.
	int start2=0		///<id of vec2.
);

///Construct from std::vector<double>
MGVector(const std::vector<double>& darrays);

///Return i-th element of the vector.
double operator() (int i) const{return ref(i);}  

///Return i-th element of the vector.
double operator[] (int i) const{return ref(i);}  

///Access to i-th Inteval.
///This is left hand side value operator. If only regference is needed,
/// operator[] should be used. 
double& operator()(int i);

///Update vector data by array of double.
///This space dimension's length data are copied from data*.
MGVector& operator=(const double*);

///Addition of two vectors.
MGVector & operator+= (const MGVector&);

///Unary minus. Negate all the elements of the vector.
MGVector operator- () const;

///Subtraction of two vectors.
MGVector & operator -= ( const MGVector & );

///Scalar multiplication.
MGVector& operator*= (double scale);

///Update own vector by vector product output, changes to 3D vector.
MGVector& operator*= (const MGVector& vec2);

///Scalar division.
MGVector& operator/= (double scale);

//////////// Member Function ////////////

///Compute angle in radian of two vectors.
/// 0<= angle <pai.
double angle(const MGVector&) const;

///Compute angle in radian of two vectors.
/// 0<= angle <pai.
double anglepai(const MGVector& v2)const{return angle(v2);};

/// Compute the angle in radian measured from this to v2 around the normal N.
/// The angle's range is 0<= angle <2*pai.
/// Although N is assumed to be parallel to N2=(*this)*v2, N may not be perpendicular
/// to v1 and v2, in which case, the projected N to N2 is used to measure the angle.
/// N is used only to decide the direction of N2.
/// v1.angle2pai(v2,N)+v2.angle2pai(v1,N)=2*pai always holds.
double angle2pai(const MGVector& v2, const MGVector& N)const;

/// 自身のベクトルと与えられたベクトルのなす角度を cosΘ で返却する

/// 自身か与えられたベクトルが零ベクトルの時は、cosΘは 1.0 とする
///Compute angle in cosine of two vectors.
double cangle(const MGVector&) const;

/// 自身のベクトルと与えられたベクトルのなす角度の sin値を返却

///Compute angle in sine of two vectors.
/// sanlge>=0.
double sangle(const MGVector& ) const;

///Compute signed sangle for 2D vectors,  (*this , v2).
double sangleSigned2D(const MGVector& v2)const;

///Clear all the element by the value init.
MGVector& clear(double init=0.0);

/// <summary>
/// Get the concavity of this and V2 around N.
/// </summary>
/// <returns>concavity from 2. to -2:
/// 2: is concave, where almost closed,
/// 1: is concave, where 90 degree open(this and v2 makes 90 degree),
/// 0: is flat, where 180 degree open(=180 degree closed),
/// -1: is convex, where 270 degree open(this and v2 makes 90 degree),
/// -2: is convex, where 360 degree open.
/// </returns>
double concavity(const MGVector& V2, const MGVector& N)const;
double concavity(const MGVector& V2)const;

///Return the 1st address of the array of the vector double data.
const double* data()const{return m_element;};
double* data(){return m_element;};

/// Generate a vector by interpolating two vectors.

///Input scalar is a ratio and when zero, output vector is a copy of *this.
/// New vector vnew=(1-t)*(*this)+t*vec2.
MGVector interpolate(double t, const MGVector& vec2) const;

/// Generate a vector by interpolating two vectors by rotation.

/// Input scalar t is a ratio and when t is zero,
/// output vector is a copy of *this and when t=1., 	output is vec2.
/// New vector vnew=a*(*this)+b*vec2, where
/// a=sin(theta2)/sin(theta), b=sin(theta1)/sin(theta). Here,
/// theta=angle of *this and vec2. theta1=t*theta, theta2=theta-theta1.
/// theta may be zero.
///When ratio is not null, ratio[0]=a and ratio[1]=b will be returned.
MGVector interpolate_by_rotate(
	double t, const MGVector& vec2,
	double* ratio=0
) const;

///Test if this and v2 are on a single straight line.

///Function's return value is true if the three points are on a straight,
///false if not.
bool is_collinear(
	const MGVector& v2
)const{	return (*this).parallel(v2);}

///Test if this, v2, and v3 are on a single straight line.

///Function's return value is true if the three points are on a straight,
///false if not.
bool is_collinear(
	const MGVector& v2,
	const MGVector& v3
)const;

///Test if this is null.
bool is_null()const{return m_sdim==0;}

///Test if the vector is unit.
bool is_unit_vector() const;

///Return true when the vector is a zero vector.
bool is_zero_vector() const;

///Return vector length.
double len() const;

///Negate the vector.
void negate(){operator*=(-1.);};

///Generate unit vector from the vector.
MGUnit_vector normalize() const;

///Test if two vectors are orthogonal, i.e. cross at right angle.
bool orthogonal(const MGVector& ) const;

///Update this to unit vector, then compute orthonormal system.

///(*this, v1, v2) organizes orthonormal system of 3D, that is
///this, v1, and v2 are all unit.
///If sv.orthogonal(*this), v1=sv.normalize().
///This is supposed not to be parallel to sv.
void orthonormalize(const MGVector& sv
	, MGVector& v1, MGVector& v2);

///Iges output. PD123=Direction.
///Function's return value is the directory entry id created.
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

///Compute the vector that is orthogonal to vec2 and is closest to this.

///"closest" means that the angle of the two vectors is minimum and
///the two vector length are equal.
MGVector orthogonize(const MGVector& vec2)const;

///Test if two vectors are parallel.
bool parallel(const MGVector& ) const;

///Get the parallelism of two vectors(this and v2). 
///Let para be the function's output, then 0<= para <=2.
///The smaller para is, dir1 and dir2 are more parallel.
///when para=0., both are completely parallel.
///         =1., they are perpendicular,
///         =2., their directions are completely opposite.
double parallelism(const MGVector& v2) const;

/// 自身のベクトルをベクトル(v2)に射影したベクトルを求める。

/// v2 が 零ベクトルのとき(*this)が返る。
MGVector project(const MGVector& v2) const;

///Reference to i-th element.
double ref(int i) const{ 
	if(i<sdim()) return m_element[i];
	else         return 0.;
}

///Resize the vector, that is , change space dimension.

///When this is enlarged, the extra space will contain garbages.
void resize(int new_sdim);

///Get the space dimension
int sdim() const { return m_sdim; };

///Set this as a null vector.
void set_null();

///Change this to a unit vector.
void set_unit();

///Store vec2 data into *this.

///Store length is vec2.len().
///Storing will be done rap-around. That is, if id i or j reached to
///each sdim(), the id will be changed to 0.
void store_at(
	int i,		///<Displacement of *this.
	const MGVector& vec2,///<Vector 2.
	int j=0		///<Displacement of vec2.
);

///Store vec2 data into *this.

///Storing will be done rap-around. That is, if id i or j reached to
///each sdim(), the id will be changed to 0.
void store_at(
	int i,		///<Displacement of *this.
	const MGVector& vec2,///<Vector 2.
	int j,		///<Displacement of vec2.
	int len		///<Length to store 
);

///swap two coordinates, (i) and (j).
void swap(int i, int j);

///Calculate dump size
virtual int dump_size() const;

///Dump Function
virtual int dump(MGOfstream& ) const;

///Restore Function
virtual int restore(MGIfstream& );

///Get the name of the class.
virtual std::string whoami()const { return "V"; };

	///Set data at element i.
protected:
	///This should be used with care, since m__length will not be set.
	///Maintenance of m_length should be done by the user, that is, m_length=-1
	///must set if updated.
	double& set(int i) {return m_element[i];};

///Friend Function
friend class MGLBRep;
friend class MGSBRep;

};

namespace MGCL{

///Compute the angel around the normal N in radian range[0., 2*pia).

///angle(v1,v2,N)+angle(v2,v1,N)=2*pai always holds.
inline double angle(const MGVector& V1, const MGVector& V2, const MGVector& N) {
	return V1.angle2pai(V2, N);
};

inline double concavity(const MGVector& V1, const MGVector& V2, const MGVector& N) {
	return V1.concavity(V2, N);
};

/// <summary>
/// N is assummed to be Z axis.
/// </summary>
inline double concavity(const MGVector& V1, const MGVector& V2) {
	return V1.concavity(V2);
};

///Get the parallelism of two vectors. 
///Let para be the function's output, then 0<= para <=2.
///The smaller para is, dir1 and dir2 are more parallel.
///when para=0., both are completely parallel.
///         =1., they are perpendicular,
///         =2., their directions are completely opposite.
inline double parallelism(const MGVector& v1, const MGVector& v2) {
	return v1.parallelism(v2);
};

/// V1をベクトル(v2)に射影したベクトルを求める。

/// v2 が 零ベクトルのときV1が返る。
MG_DLL_DECLR MGVector project(const MGVector& V1, const MGVector& V2);

};

/** @} */ // end of BASE group
