/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGPoint_HH_
#define _MGPoint_HH_

#include "mg/Geometry.h"
#include "mg/Position.h"
#include "mg/isects.h"

class MGIfstream;
class MGOfstream;
class MGInterval;
class MGBox;
class MGVector;
class MGPosition_list;
class MGMatrix;
class MGTransf;

/** @addtogroup GEO
 *  @{
 */

///MGPoint represents one dimensional manifold, a point in a space.

///See MGPosition class, which is recommended to use usually.
class MG_DLL_DECLR MGPoint:public MGGeometry{

public:

//////////// Constructor ////////////
MGPoint()=default;
~MGPoint()=default;
MGPoint(const MGPoint&)=default;
MGPoint(MGPoint&&)=default;
MGPoint& operator=(const MGPoint&)=default;
MGPoint& operator=(MGPoint&&)=default;

///Conversion constructor from a position.
MGPoint(const MGPosition& P);

///Construct a point by changing the space dimension or ordering
///the space dimension element.
MGPoint(
	int sdim,		///<new space dimension.
	const MGPoint& P		///<original point.
	, int start1=0		///<start position coordinate of new point.
	, int start2=0		///<start position coordinate of the original.
);

///Assignment.
///When the leaf object of this and obj2 are not equal, this assignment
///does nothing.
MGPoint& operator=(const MGGel& gel2);//Copy Assignment.
MGPoint& operator=(MGGel&& gel2);//Move Assignment.

///Return i-th element of the position.
double operator[] (int i) const{return m_point.ref(i);}
double operator() (int i) const{return m_point.ref(i);}

///Access to i-th element.
double& operator()(int i){return m_point(i);};

///Object transformation.
MGPoint& operator+=(const MGVector& v);
MGPoint& operator-=(const MGVector& v);
MGPoint& operator*=(double scale);
MGPoint& operator*=(const MGMatrix& mat);
MGPoint& operator*=(const MGTransf& tr);

////////////Logical operator overload/////////

///comparison
bool operator==(const MGPoint& p2)const;
std::partial_ordering operator<=>(const MGPoint& p2)const;

//gel2 must be the same class as this.
bool equal_test(const MGGel& gel2)const override;

//gel2 must be the same class as this.
std::partial_ordering ordering_test(const MGGel& gel2)const override;

////////// Member Function ////////////

///Obtain ceter coordinate of the geometry.
MGPosition center() const{return m_point;};

///Obtain ceter parameter value of the geometry.
MGPosition center_param() const{ return MGPosition();};

///Changing this object's space dimension.
void change_dimension(
	int sdim,			///< new space dimension
	int start1=0, 		///< Destination order of new object.
	int start2=0 		///< Source order of this object.
);

///Construct new geometry object by copying to newed area.
///User must delete this copied object by "delete".
MGPoint* clone() const;

///Return minimum box that includes whole of the geometry.
///曲線部分を囲むボックスを返す。
///Compute the box from the scratch.
void compute_box(MGBox& bx) const;

///Construct new geometry object by changing
///the original object's space dimension.
///User must delete this copied object by "delete".
MGPoint* copy_change_dimension(
	int sdim,			///< new space dimension
	int start1=0, 		///< Destination order of new line.
	int start2=0 		///< Source order of this line.
)const;

///Compute direction unit vector of the geometry.
///For the point, this is undefined.
MGUnit_vector direction(const MGPosition& param) const;

///Draw 3D curve in world coordinates.
///The object is converted to curve(s) and is drawn.
void drawWire(
	mgVBO& vbo,///<The target graphic object.
	int line_density=1	///<line density to draw a surface in wire mode.
)const;

///Draw 3D point(vertex) in world coordinates.

///The object is converted to point(s) and is drawn.
///This is valid only for topology objects or MGPoint.
void drawVertex(
	mgVBO& vbo///<The target graphic object.
)const;

/// Evaluate n'th derivative data. n=0 means positional data evaluation.
MGVector evaluate(
	const MGPosition& t,	///< Parameter value,
				///<t's space dimension is geometry's manifold dimension.
	const int* nderiv=0	///<Order of derivative of i-th parameter
				///<in nderiv[i],
				///<When nderiv=null, nderiv[i]=0 is assumed for all i.
)const{return m_point;}

/// Return This object's typeID
long identify_type()const;

///Test if input parameter value is inside parameter range of the line.
bool in_range(const MGPosition& t)const;

///Intersections are obtained from two objects, which are known using
///the MGisects::object1() and object2().
///****NOTE****
///When two objects' manifold dimension are the same, object1 is this object
///at the invocation of MGObject::intersection(), and object2 is the argument
///object.
///However, their manifold dimension are not the same, object1 is always
///the lower dimension's object and object2 is the higer dimension's object.
MGisects isect(const MGObject& obj2)const{return MGisects();};

///Return manifold dimension, i.e. 0:point, 1:curve, 2:surface.
int manifold_dimension() const{return 0;};

///Negate direction of this geometry.
void negate(){;};

///Transform the coordinates of boundary of this geometry so that
///new coordinate of boundary is the same coordinate as the new one of
///this geometry after negate() of this geometry is done.
///That is, boundary coordinates are parameter world of this geometry.
void negate_transform(MGGeometry& boundary)const{;};

///Test if given point is on the geometry or not. If yes, return parameter
///value of the geometry. Even if not, return nearest point's parameter.
/// 指定点が自身上にあるかを調べる。曲線上にあれば，そのパラメーター値を，
/// なくても最近傍点のパラメータ値を返す。
/// Function's return value is >0 if the point is on the geometry,
/// and 0 if the point is not on the geometry.
bool on(
	const MGPosition& P,///<Point(指定点)
	MGPosition&	param	///<Parameter of the geometry(パラメータ)
) const;

///Compute parameter value of given point.
/// 自身の上の指定点を表すパラメータ値を返す。
/// If input point is not on the geometry, return the nearest point on the
/// geometry.
MGPosition parameter(
	const MGPosition& P	///<Point(指定点)
) const;

///Return parameter range of the geometry(パラメータ範囲を返す)
MGBox parameter_range() const;

const MGPosition& position()const{return m_point;}
MGPosition& position(){return m_point;}

///Round t into geometry's parameter range.
/// 入力パラメータをパラメータ範囲でまるめて返却する。
MGPosition range(const MGPosition& t) const;

///Return i-th element of the point.
double ref(int i) const{return m_point.ref(i);};

///Return space dimension
int sdim() const;

///IGES output function.
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

/// Output function.
std::ostream& toString(std::ostream&) const;

///Get the name of the class.
std::string whoami()const{return "Point";};

protected:

///メンバデータを読み出す関数
void ReadMembers(MGIfstream& buf);

///メンバデータを書き込む関数
void WriteMembers(MGOfstream& buf) const;

private:
	MGPosition m_point;

};

/** @} */ // end of GEO group
#endif
