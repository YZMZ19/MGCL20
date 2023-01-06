/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Transf.h"
#include "mg/Position.h"
#include "mg/CParam_list.h"
#include "mg/Position_list.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGPosition.cpp
// Implemetation of class MGPosition
//
// Represent a positional data. MGPosition is the same class as MGposition,
// only is an alias name.


//Let this be the center of the rotation, then compute the angle rotated 
//around the normal from start to end.
//Although normal is assumed to be parallel to N2=V1*v2, normal may not perpendicular
//to v1 and v2, in which case, the projected normal to N2 is used to measure the angle.
//Here v1=start-*this, and v2=end-*this.
//this->angle(start,end,normal)+this->angle(end,start,normal)=2*pai always holds.
double MGPosition::angle(
	const MGPosition& start,
	const MGPosition& end,
	const MGVector& normal
)const{
	MGVector V1(start-*this), V2(end-*this);
	return V1.angle2pai(V2,normal);
}

//Clear all the elements by the value init.
MGPosition& MGPosition::clear(double init){
	static_cast<MGVector&>(*this).clear(init);
	return *this;
}

//Return the distance of this and P2.
double MGPosition::distance(const MGPosition& P2)const{
	return (static_cast<const MGVector&>(*this)
			- static_cast<const MGVector&>(P2)).len();
}

// Generate a Position by interpolating two Position. Input scalar is
// a ratio and when zero, output position is a copy of the own vector.
MGPosition MGPosition::interpolate(double t2, const MGPosition& vec2) const{
	return static_cast<const MGVector&>(*this).interpolate(t2,vec2);
}

//Operator Oveload
	
//Update position data by array of double.
//This space dimension's length data are copied from data*.
MGPosition& MGPosition::operator=(const double* a){
	MGVector::operator=(a);
	return *this;
}

// 自身のPositionに与えられたVectorを加算して自身のPositionとする 
MGPosition& MGPosition::operator+= (const MGVector& vec){
	static_cast<MGVector&>(*this) +=vec;
	return *this;
}

// 単項マイナス。自身のPositionを反転し、Positionを生成
MGPosition MGPosition::operator- () const{
	return -static_cast<const MGVector&>(*this);
}

// 自身のPositionから与えられたVectorを減算し自身のPositionとする
MGPosition& MGPosition::operator-= (const MGVector& vec){
	static_cast<MGVector&>(*this) -= vec;
	return *this;
}

// Scalarの乗算を行い自身のPositionとする
MGPosition& MGPosition::operator*= (double a){
	static_cast<MGVector&>(*this)*=a;
	return *this;
}

// MatrixによるPositionの変換を行い自身のPositionとする
MGPosition& MGPosition::operator*= (const MGMatrix& mat){
	static_cast<MGVector&>(*this) *= mat;
	return *this;
}

// PositionのTransformを行いVectorを生成
MGPosition operator*(const MGPosition& p1,const MGTransf& tr){
	return p1*tr.affine()+tr.translation();
}

// PositionのTransformを行いPositionを生成して，
// 自身のPositionとする
MGPosition& MGPosition::operator*= (const MGTransf& tran){	
	MGVector& V=static_cast<MGVector&>(*this);
	V*= tran.affine();
	V+= tran.translation();
	return *this;
}

// Scalar除算を行い自身のPositionとする
MGPosition& MGPosition::operator/= (double a){
	static_cast<MGVector&>(*this) /= a;
	return *this;
}

//Friend Function

//Test if P1, P2, and P3 are on a single straight line.
//Function's return value is true if the three points are on a straight,
//false if not.
bool is_collinear(
	const MGPosition& P1,
	const MGPosition& P2,
	const MGPosition& P3
){return P1.is_collinear(P2,P3);}

