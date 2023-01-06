/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGUnit_vector_HH_
#define _MGUnit_vector_HH_
/** @file */
/** @addtogroup BASE
 *  @{
 */
class MGPosition;
#include "mg/Vector.h"

// MGUnit_vector.h 
// Header for class MGUnit_vector

/// Define a unit vector, is a MGVector.
class MG_DLL_DECLR MGUnit_vector: public MGVector {

public:

/// Void constructor.
MGUnit_vector(
	int sdim=3	///< space dimension
);

/// Unit vector from general vector.
MGUnit_vector(const MGVector& v);
MGUnit_vector(MGVector&& v);

/// Assignment From Vector.
MGUnit_vector& operator= (const MGVector& vec2) ;
MGUnit_vector& operator= (MGVector&& vec2) ;

///@cond
// スカラーの乗算を行い自身のベクトルとする.
MGUnit_vector& operator *= ( double )=delete;

// スカラー除算を行い自身のベクトルとする.
MGUnit_vector& operator /= ( double )=delete;
///@endcond

///
///Update vector data by array of double.
///The result is unit of the updated vector.
///
MGUnit_vector& operator= (const double*);

///Update the unit vector by adding vec2.
//The result is unit of the vector of two vector addition.
MGUnit_vector& operator+= (const MGVector& vec2);

///Unary minus. Negate the vector.
MGUnit_vector operator- () const;

///Update the unit vector by subtractiong vec2.
///The result is unit of the vector of two vector subtraction.
MGUnit_vector& operator-= (const MGVector& vec2);

///Update own vector by vector product output.
//Cchanges to 3D vector. The result is a unit one of two vector product.
MGUnit_vector& operator*= (const MGVector& vec2);

//////// Member function. ////////

///Compute orthonormal system, given sub(sv) vectors.
///(*this, v1, v2) organizes orthonormal system of 3D, that is
///this, v1, and v2 are all unit, and this=v1*v2, v1=v2*this, v2=this*v1.
///If sv.orthogonal(*this), v1=sv.normalize().
///This is supposed to be not parallel to sv.
void orthonormal(
	const MGVector& sv,// to be nearly equal to v1.
	MGVector& v1,
	MGVector& v2
) const;

///Get the name of the class.
std::string whoami()const override{ return "U"; };

};

///Get the unit normal of the triangle (P0, P1, P2).
///UnitNormal(P0,P1,P2)=-UnitNormal(P0,P2,P1).
///(V1,V2,UnitNOrmal) organizes orthonormal system, wher
///V1=P1-P0, V2=P2-P1.
MGUnit_vector UnitNormal(
	const MGPosition& P0, // 三角形の頂点の座標
	const MGPosition& P1,
	const MGPosition& P2
);

/** @} */ // end of BASE group
#endif
