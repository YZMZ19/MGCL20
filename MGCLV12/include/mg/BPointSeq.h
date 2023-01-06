/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
/** @file */ 
/** @addtogroup BASE
 *  @{
 */
#pragma once
#include <cassert>
#include <vector>
#include "mg/MGCL.h"

// MGBPointSeq.h
//

//Forward Declaration
class MGIfstream;
class MGOfstream;
class MGBox;
class MGVector;
class MGPosition;
class MGSPointSeq;
class MGMatrix;
class MGTransf;
class MGPlane;
class MGStraight;

///Defines BPoint seq of a space dimension and of a capacity.

///MGBPointSeq is an array of world coordinates of a space dimension(any number),
///desined especially for MGLBRep class's control polygon handling.
///Let n be a B-rep dimension(i.e., number of points of the control polygon), and 
///sd be a space dimension, then MGBPointSeq bp has the minimum area of n*sd.
///The subscription of MGBPointSeq bp is 2 dimension as bp(i,j), 0<=i<n, and 0<=j<sd.
///If MGBPointSeq can be considered as a Matrix of n by sd. 
class MG_DLL_DECLR MGBPointSeq {

public:

///translation by a vector.
MG_DLL_DECLR friend MGBPointSeq operator+(const MGVector& v, const MGBPointSeq& b);

///BPointをスケーリングしてできるオブジェクトを生成する。
///Generates a object by scaling.
MG_DLL_DECLR friend MGBPointSeq operator* (double scale, const MGBPointSeq&);

///String stream Function
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream&, const MGBPointSeq&);

////////Special member functions/////////

///Dimensioning of MGBPointSeq is MGBPointSeq(capacity,dim).
///effective length will be set capacity(i.e. length() will be capacity).
explicit MGBPointSeq(int capacity=0, int dim=0);

~MGBPointSeq(){if(m_bpoint) delete[] m_bpoint;};///<Dtor.
MGBPointSeq(const MGBPointSeq& rhs);///<Copy constructor.
MGBPointSeq& operator=(const MGBPointSeq& rhs);///<Copy assignment.
MGBPointSeq(MGBPointSeq&& rhs);///<Move constructor.
MGBPointSeq& operator=(MGBPointSeq&& rhs);///<Move assignment.


///Construct a MGBPointSeq by copying original MGBPointSeq,
///Can change the order of coordinates.
/// (*this)(start1,.) will be old_brep(start2,.).
MGBPointSeq(int dim,		///< New Space Dimension.
	const MGBPointSeq& old_brep,///< Origianl B-Rep.
	int start1=0,		///< Destination start order to store.
	int start2=0			///< Source start order to retrieve.
);

/// Construct by extracting one line data of sp along u or v direction.
MGBPointSeq(
	bool along_u,	///<indicates which direction make a line out of sp.
		///<if true, along u direction:sp(i,m,.) for i=0, ..., nu-1 makes the BPointSeq.
		///<if false, along v direction:sp(m,j,.) for j=0, ..., nv-1 makes the BPointSeq.
	int m,	///<index of u or v as above, according to along_u.
	const MGSPointSeq& sp	///< Origianl SPoint.
);

/// Construct by extracting sub interval of bp_old.
///(*this)(0+i,.) will be bp_old(start_id+i,.) for i=0,...,num-1.
MGBPointSeq(
	int start_id,	///< Start id of bp_old
	int num,			///< New length(of new BPoint)
	const MGBPointSeq& bp_old///< Origianl BPoint.
);

///Conversion constructor.
MGBPointSeq(const std::vector<MGPosition>& poses);
/////////Destructor////////

////////Operator Oveload////////

///Access to (i,j)th element(LHS version).
double& operator()(int i, int j){	
	assert(i<capacity()&&j<m_sdim);
	return m_bpoint[i+m_capacity*j];
}

///Access to (i,j)th element(RHS version).
double operator()(int i, int j) const{return ref(i,j);};

///Extract (i,j) coordinate values for 0<=j<sdim() as a vector.
MGVector operator() (int i) const;

/// 曲線の平行移動を行いオブジェクトを生成する。
///Generates an object by translation.
MGBPointSeq operator+(const MGVector& v) const;

/// 与ベクトルだけ曲線を平行移動して自身とする。
///Updates an object by translation.
MGBPointSeq& operator+= (const MGVector&);

/// 曲線の逆方向に平行移動を行いオブジェクトを生成する。
///Generates an object by translation.
MGBPointSeq operator- (const MGVector&) const;

/// 与ベクトルだけ曲線をマイナス方向に平行移動して自身とする。
///Updates the object by translation.
MGBPointSeq& operator-= (const MGVector&);

///Add operation of two BPointSeq.
MGBPointSeq operator+ (const MGBPointSeq& bp2) const;

///Add operation of two BPointSeq.
MGBPointSeq& operator+= (const MGBPointSeq& bp2);

///Subtract operation of two BPointSeq.
MGBPointSeq operator- (const MGBPointSeq& bp2) const;

///Subtract operation of two BPointSeq.
MGBPointSeq& operator-= (const MGBPointSeq& bp2);

/// 与えられたスケールで曲線の変換を行いオブジェクトを生成する。
///Generates an object by multiplying scale to the original.
MGBPointSeq operator* (double scale) const;

/// 与えられたスケールで曲線の変換を行い自身の曲線とする。
///Updates the object by multiplying scale.
MGBPointSeq& operator*= (double scale);

/// 与えられたスケールで曲線の変換を行いオブジェクトを生成する。
///Generates an object by multiplying scale to the original.
MGBPointSeq operator/ (double scale) const;

/// 与えられたスケールで曲線の変換を行い自身の曲線とする。
///Updates the object by multiplying scale.
MGBPointSeq& operator/= (double scale);

/// 与えられた変換で曲線の変換を行いオブジェクトを生成する。
///Generates an object by multiplying matrix to the original.
MGBPointSeq operator* (const MGMatrix&) const;

/// 与えられた変換で曲線の変換を行い自身の曲線とする。
///Updates the object by multiplying matrix.
MGBPointSeq& operator*= (const MGMatrix&);

/// 与えられた変換で曲線のトランスフォームを行いオブジェクトを生成する。
///Generates an object by multiplying transformation to the original.
MGBPointSeq operator* (const MGTransf&) const;

/// 与えられた変換で曲線のトランスフォームを行い自身とする。
///Updates the object by multiplying transformation.
MGBPointSeq& operator*= (const MGTransf&);							  

///Compare two BPointSeq if they are equal.
///Return true if equal.
bool operator== (const MGBPointSeq& ) const;

///Compare two BPointSeq if they are equal.
/// ///Return true if not equal.
bool operator!= (const MGBPointSeq& brep) const
{return !(operator== (brep));}

////////Member Function////////

///compute an average plane of the point sequence.
///Function's return value is:
/// 1: Point seq is a point.		2: Point seq is on a line.
/// 3: Plane is output. 
int average_plane(
	MGPosition& center	///< center of point seq will be output.
	, MGPlane& plane	///< Plane will be output, when average_plane=3.
	, MGStraight& line	///< Straight line will be output            =2.
	, double& deviation///< Maximum deviation from point, line, or plane will be output.
)const;

///Compute minimum box sorrounding the all the points.
MGBox box()const;

///Compute minimum box sorrounding the all the points.
void compute_box(MGBox& bx) const;

///Exchange ordering of the coordinates.
///Exchange coordinates (j1) and (j2).
void coordinate_exchange(int j1, int j2);

///Returns a pointer to the data area.
double* data(int i=0, int j=0);

///Returns a pointer to the area.
const double* data(int i=0, int j=0)const;

///Insert vector data to this(i,j)  0<=j<sdim().
void insert_at(
	int i,	///<id of this which indicates the placement.
	const MGVector& vctr///<Input vector.
);

///Insert array data[j] to this(i,j)  0<=j<sdim().
void insert_at(int i,	///<id of this which indicates the placement.
	const double* data	///<Input data array.
);

///Test if this is null.
bool is_null()const{return m_length ==0;};

///Returns the actual size of Bpoint seq.
int length()const{return m_length;};

///Negate this point sequence.
///That is, change the ordering of this(i,.) for i.
void negate();

///Compute non_homogeneous coordonate data without w coordinate element,
///assumed that this is homogeneous coordinate data,
///i.e., maximum space dimension element is w(weight) coordinate.
MGBPointSeq non_homogeneous()const;

///Check if coefficients are planar.
///Funtion's return value is;
/// 0: Not planar, nor a point, nor straight line.
/// 1: coefficients are within a point.  2: coefficients are a straight line.
/// 3: coefficients are planar.
int planar(
	MGPlane& plane	///<When coefficients are not straight line
					///<nor a point, plane is returned,
					///<Even when not planar, plane nearest is returned.
	, MGStraight& line	///<When coefficients are a line, the line is returned.
	, MGPosition& center///<Center of the coefficients is always returned.
)const;

///Retrieve sub data of i-th point of the BPointSeq.
///That is, P(k)=(*this)(i,j+k) for k=0, ..., sd-1.
void point(int i, int j, int sd, MGPosition& P)const;

///Return (i,j)-th element data, return 0.0 when j>=sdim().
double ref(int i, int j) const;

///Change the capacity to at least sz and update the effective length to sz,
///does not update sdim. After invoking reshape(), at least sz point data are guaranteeed 
///to store.
///reshape() guarantees the original data BPoint(i,.) before invoking reshape, which
///will be stored in the new BPoint(start+i,.).
void reshape(
	int capacity///<New capacity
	, int start=0///<To which place to store the original data.
);								  

///Resize the array. The result will contain garbages.
///m_capacity,  m_sdim, and m_length will be set as
///m_capacity=sz, m_sdim=dim, m_length=sz.
//dim=0 means space dimension be set same as the old one.
void resize(int sz, int dim=0);

///Returns the space dimension.
int sdim() const {return m_sdim;}

///Set the length of effective data.
void set_length(int length);

///Set this BPointSeq as a null.
void set_null();

///Returns the capacity.
int capacity() const {return m_capacity;}

///Store vector data vector(from+j) to this(i,to+j) for 0<=j<sdim().
///When (form+j) or (to+j) reached to maximum space dimension id, next id
///becomes 0(form the start).
void store_at(int i,	///<id of this which indicates the placement.
	const MGVector& vctr,///<Input vector.
	int to=0,	///<Indicates to where of this in the space dimension id.
	int from=0	///<Indicates from where of vector in the space dimension.
);

///Store vector data vector(from+j) to this(i,to+j) for 0<=j<len.
///When (form+j) or (to+j) reached to maximum space dimension id, next id
///becomes 0(form the start).
void store_at(int i,	///<id of this which indicates the placement.
	const MGVector& vctr,///<Input vector.
	int to,	///<Indicates to where of this in the space dimension id.
	int from,///<Indicates from where of vector in the space dimension.
	int len///<length of the data to store.
);

///Store data[j] to this(i,to+j) for 0<=j<sdim().
///When (to+j) reached to maximum space dimension id, next id
///becomes 0(form the start).
void store_at(int i,	///<id of this which indicates the placement.
	const double* data,	///<Input data array.
	int to=0///<Indicates to where of this in the space dimension id.
);

///Transformation of own for rational(MGRLBRep) Control Polygon.
///Scaling.
MGBPointSeq& homogeneous_transform(double);

///Add the vector.
MGBPointSeq& homogeneous_transform(const MGVector&);

///Multiply the matrix.
MGBPointSeq& homogeneous_transform(const MGMatrix&);

///Multiply the transform.
MGBPointSeq& homogeneous_transform(const MGTransf&);

///Dump Functions, Calculate dump size.
int dump_size() const;

///Dump Function
int dump(MGOfstream& ) const;

///Restore Function
int restore(MGIfstream& );
	
////////Member Data////////
private:
	int m_capacity;	///< maximum capacity to store BPoint seq of space dimension m_sdim.
	int m_sdim;		///< Space Dimension.
	int m_length;	///< length of BPoint actual data stored.
						///< This BPointSeq is null when m_length==0;
	double* m_bpoint;	///<Area to store Bpoint.

friend class MGLBRep;
friend class MGRLBRep;
};

/** @} */ // end of BASE group
