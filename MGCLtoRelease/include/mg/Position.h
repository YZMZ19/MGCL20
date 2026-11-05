/********************************************************************/
/* Copyright (c) 2021 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#pragma once

/** @file */
/** @addtogroup BASE
 *  @{
 */
#include "mg/Vector.h"

// MGPosition.h
// Header for MGPosition.

// Forward Declaration.
class MGMatrix;
class MGTransf;
class MGPoint;
class MGCurve;
class MGSurface;
class MGCParam_list;
class MGPosition_list;
class MGIfstream;
class MGOfstream;

///Represent a positional data.
class MG_DLL_DECLR MGPosition{
	MGVector m_vector;///<MGPosition is a MGVector, but it is not a MGVector. It is a position.
	friend class MGVector;

public:

///Translation of the position.
inline friend
MGPosition operator+(const MGPosition& p1, const MGVector& vec) {
	return operator+(p1.m_vector, vec);
};

inline friend
MGPosition operator+(const MGPosition& p1, const MGPosition& p2) {
	return p1.m_vector + p2.m_vector;
};

inline friend
MGPosition operator+ (const MGVector& v, const MGPosition& p) {
	return v+p.m_vector; 
};

///A position translated by v.
inline friend
MGVector operator-(const MGPosition& p1, const MGPosition& p2) {
	return p1.m_vector - p2.m_vector;
};

///A position translated by v.
inline friend
MGPosition operator-(const MGPosition& p1, const MGVector& v) {
		return p1.m_vector - v;
};

///A position translated by v.
inline friend
MGVector operator-(const MGVector& v, const MGPosition& p1) {
	return v - p1.m_vector;
};

///Scaling of the position.
inline friend
MGPosition operator*(double s, const MGPosition& p) {
	return p.m_vector * s; };

///Scaling of the position.
inline friend
MGPosition operator*(const MGPosition& p1, double s) {
	return p1.m_vector * s;
};

///Scaling of the position.
inline friend
MGPosition operator/(const MGPosition& p1, double s) {
	return p1.m_vector / s;
};

///Test if P1, P2, and P3 are on a single straight line.
///Function's return value is true if the three points are on a straight,
///false if not.
MG_DLL_DECLR friend bool is_collinear(
	const MGPosition& P1,
	const MGPosition& P2,
	const MGPosition& P3
);

////////// Constructor コンストラクタ ////////////

///Void constructor.
explicit MGPosition(int sdim=0):m_vector(sdim,0.0){;};

/// <summary>
/// initialize MGPosition by initializer list.
/// </summary>
/// <param name="data"></param>
MGPosition(std::initializer_list<double> data);

/// Conversion Constructor from a point
MGPosition(const MGPoint& point);

/// Conversion Constructor from a vector
MGPosition(const MGVector& vec): m_vector(vec){;};

MGPosition(MGVector&& v):m_vector(std::move(v)){;};

///Construct 2D position by providing x,y coordinate data.
MGPosition(double x, double y):m_vector(x,y){;};

///Construct 3D position by providing x,y,z coordinate data.
MGPosition(double x, double y, double z):m_vector(x,y,z){;};

///Construct 4D position by providing x,y,z,w coordinate data.
MGPosition(double x, double y, double z, double w):m_vector(x,y,z,w){;};

///Construct, given coordinate data through double array.
///***** This is the fundamental constructor.*****
MGPosition(int sdim, const double* v ):m_vector(sdim,v){;};

///Construct by copying a position, changing space dimension and
///ordering of coordinates.
MGPosition(int sdim, const MGPosition& p, int start1=0, int start2=0)
:m_vector(sdim, p.m_vector, start1, start2){;};

///Construct from std::vector<double>.
MGPosition(const std::vector<double>& darrays): m_vector(darrays){;};

///Return i-th element of the vector.
double operator() (int i) const { return m_vector.ref(i); }

///Access to i-th Inteval.
double& operator()(int i) { return m_vector(i); }

///Return i-th element of the vector.
double operator[] (int i) const { return m_vector.ref(i); }

///Assignment
///Update position data by array of double.
///This space dimension's length data are copied from data*.
MGPosition& operator=(const double* a);

///Translation of the position.
MGPosition& operator+= (const MGVector& vec);

///Unary minus. Negate all of the elements of the position.
MGPosition operator- () const;

///Translation of the position.
MGPosition& operator-= (const MGVector& vec);

///Scaling of the position.
MGPosition& operator*=(double scale);

///Scaling of the position.
MGPosition& operator/= (double);

///comparison
bool operator==(const MGPosition& p2)const;
bool operator==(const MGVector& p2)const;//p2 is treated as a position.

///Comparison depends on two vectors' length.
std::partial_ordering operator<=>(const MGVector& v2)const {
	return m_vector <=> v2;
}
std::partial_ordering operator<=>(const MGPosition& p2)const {
	return m_vector <=> p2.vector();
}

///String stream function
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream&, const MGPosition&);

///Compute an angle around a normal.

///Let this be the center of the rotation, then compute the angle rotated 
///around the normal from start to end.
///angle(start,end,normal)+angle(end,start,normal)=2*pai always holds.
double angle(
	const MGPosition& start,
	const MGPosition& end,
	const MGVector& normal
)const;

///Clear all the element by the value init.
MGPosition& clear(double init=0.0);

///Compute the closest point parameter value of the curve from this point.
///Function's return value is the parameter value of the curve.
double closest(const MGCurve& curve) const;

///Compute the closest point parameter value (u,v)of the surface
///from this point.
MGPosition closest(const MGSurface& surf) const;

///Construct new surface object by copying to newed area.
///User must delete this copied object by "delete".
MGPosition* clone() const{return new MGPosition(*this);};

///Return the 1st address of the array of the point double data.
const double* data()const{return m_vector.data();};
double* data(){return m_vector.data();};

///Return the distance of this and P2.
double distance(const MGPosition& P2)const;

///Calculate dump size
int dump_size() const { return m_vector.dump_size(); };

/// Generate a Position by interpolating two Position.
///Input scalar is:
/// a ratio t2. When t2 is zero, output position is a copy of the own position.
///Output=(*this)*(1-t2)+vec2*t2.
MGPosition interpolate(double t2, const MGPosition& vec2) const;

///Test if this, v2, and v3 are on a single straight line.
///Function's return value is true if the three points are on a straight,
///false if not.
bool is_collinear(
	const MGPosition& p2,
	const MGPosition& p3
)const{
	return m_vector.is_collinear(p2.vector(), p3.vector());
};

///Return the length from the origin to this.
double len() const{return m_vector.len();};

///Test if this position is on a curve.
///If on, return the parameter value.
///Even if not on, return the nearest point of the curve.
/// Function's return value is >0 if the point is on the curve,
/// and 0 if the point is not on the curve.
bool on(
	const MGCurve& curve,	///< Curve
	double& t	///< Parameter value of the nearest point on the curve.
) const;

///Test if this position is on a surface.

///If on, return the parameter value.
///Even if not on, return the nearest point of the surface.
/// Function's return value is >0 if the point is on the curve,
/// and 0 if the point is not on the curve.
bool on(
	const MGSurface& surf,	///< Surface pointer
	MGPosition& uv	///<Parameter value of the nearest point on surface.
) const;

///Out to iges as PD116=Point.

///Function's return value is the directory entry id created.
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

/// Return curve's parameter value of this point.

/// If this point is not on the curve, return the nearest point's parameter
/// value on the curve.
double param(const MGCurve& crv) const;

/// Return surface's parameter value of this point.

/// If this point is not on the surface, return the nearest point's parameter
/// value on the surface.
MGPosition param(const MGSurface& srf) const;

///Compute all foot points of the perpendicular line from this point to a curve.
MGCParam_list perps(
	const MGCurve& crv		///<Curve
)const;

///Compute all foot points of the perpendicular line from this point to a surface.
MGPosition_list perps(
	const MGSurface& srf	///<Surface
) const;

///Reference to i-th element.
double ref(int i) const { return m_vector.ref(i); };

///Resize the position, that is , change space dimension.
///When this is enlarged, the extra space will contain garbages.
inline
void resize(int new_sdim) {m_vector.resize(new_sdim);};

///Dump Function
int dump(MGOfstream& ) const;

///Restore Function
int restore(MGIfstream& );

///Get the space dimension
inline
int sdim() const { return m_vector.sdim(); };

///Set this as a null position.
void set_null() { m_vector.set_null(); };
///Test if this is null.
bool is_null()const { return sdim() == 0; }

///Change this to length=1.
void set_unit() { m_vector.set_unit(); };

///Store vec2 data into *this.
///Store length is vec2.len().
///Storing will be done rap-around. That is, if id i or j reached to
///each sdim(), the id will be changed to 0.
void store_at(
	int i,		///<Displacement of *this.
	const MGPosition& vec2,///<Vector 2.
	int j = 0		///<Displacement of vec2.
){m_vector.store_at(i, vec2.vector(), j);};

///Store vec2 data into *this.
///Storing will be done rap-around. That is, if id i or j reached to
///each sdim(), the id will be changed to 0.
void store_at(
	int i,		///<Displacement of *this.
	const MGPosition& vec2,///<Vector 2.
	int j,		///<Displacement of vec2.
	int len		///<Length to store 
) {	m_vector.store_at(i, vec2.vector(), j, len);};

///swap two coordinates, (i) and (j).
void swap(int i, int j){m_vector.swap(i,j);};	

//Get the vector of the position.
const MGVector& vector()const { return m_vector; };
MGVector& vector(){ return m_vector; };

///Get the name of the class.
std::string whoami()const{ return "P"; };

};

namespace MGCL{

///Compute the angel around the Normal N in radian range[0., 2*pia).

///Here V1=P1-origin, V2=P2-origin.
///Although N is assumed to be parallel to N2=V1*v2, N may not perpendicular
///to v1 and v2, in which case, the projected normal to N2 is used to measure the angle.
///angle(origin,P1,P2,N)+angle(origin,P2,P1,N)=2*pai always holds.
inline double angle(
	const MGPosition& origin,
	const MGPosition& P1,
	const MGPosition& P2,
	const MGVector& N
){
	return origin.angle(P1,P2,N);
};

};

/** @} */ // end of BASE group
