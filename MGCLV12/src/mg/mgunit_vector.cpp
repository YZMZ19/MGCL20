/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Unit_vector.h"
#include "mg/Tolerance.h"
#include "mg/Default.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGUnit_vector.cc
// Implementation of MGUnit_vector

//
//Constructor
// 
// void コンストラクタ
MGUnit_vector::MGUnit_vector(int sdim):MGVector(sdim){
	assert(sdim>0);
	int sdimm1=sdim-1;
	for(int i=0; i<sdimm1; i++) m_element[i]=0.;
	m_element[sdimm1]=1.;
	m_length=1.;
}

// ベクトルを指定してその単位ベクトルを生成する
MGUnit_vector::MGUnit_vector(const MGVector& vec) :MGVector(vec.sdim()){
	int dim=sdim();
	if(!dim) (*this)=MGUnit_vector();
	else{
		int i;
		double length = vec.len();

		// 自身のベクトルが零ベクトルの時はデフォルトベクトルを生成
		if(MGMZero(length)){
			for(i=0; i<dim-1; i++) m_element[i]=0.;
			m_element[dim-1]=1.;
		}
		// 零ベクトル以外は与えられたベクトルの成分を長さで割る
		else {
			for(i=0; i<dim; i++) m_element[i]=vec.ref(i)/length;
		}
		m_length=1.;
	}
}

// ベクトルを指定してその単位ベクトルを生成する
MGUnit_vector::MGUnit_vector(MGVector&& vec):MGVector(std::move(vec)){
	MGVector::set_unit();
}

///Compute orthonormal system, given sub(sv) vectors.
///(*this, v1, v2) organizes orthonormal system of 3D, that is
///this, v1, and v2 are all unit, and this=v1*v2, v1=v2*this, v2=this*v1.
///If sv.orthogonal(*this), v1=sv.normalize().
///This is supposed to be not parallel to sv.
void MGUnit_vector::orthonormal(
	const MGVector& sv,//nearly equal to v1.
	MGVector& v1,
	MGVector& v2
)const{
	v2=(*this)*sv;
	if(MGMZero(v2.len())){
		if(*this==mgZ_UVEC) v1=MGVector(1.,0.);
		else if(*this==-mgZ_UVEC) v1=MGVector(0.,1.);
		else{
			double dx=fabs(ref(0)), dy=fabs(ref(1)), dz=fabs(ref(2));
			v1 = MGVector(0., 0., 1.);
			if(MGRZero(dz))
				v1 = MGVector(-ref(1), ref(0)); 
			else if(dx>dz){
				if(dz>dy)
					v1 = MGVector(0., 1., 0.);//dy is min.
			}else if(dy>dx)
				v1 = MGVector(1., 0., 0.);//dx is min.
			else
				v1 = MGVector(0., 1., 0.);//dy is min.
		}
		// Normalize m and n.
		v2 = (*this)*v1;
	}
	v2.set_unit();
	v1 = v2*(*this);
	v1.set_unit();
}


MGUnit_vector& MGUnit_vector::operator= (const MGVector& vec2){
	if(vec2.is_unit_vector()) MGVector::operator= (vec2);
	else MGVector::operator= (vec2.normalize());
	return *this;
}

//Move Assignment From Vector.
MGUnit_vector& MGUnit_vector::operator= (MGVector&& vec2){
	MGVector::operator= (std::move(vec2));
	MGVector::set_unit();
	return *this;
}
	
//Update vector data by array of double.
//The result is unit of the updated vector.
MGUnit_vector& MGUnit_vector::operator=(const double* array){
	MGVector vec(sdim(),array);
	return *this=std::move(vec);
}

//Update the unit vector by adding vec2. The result is unit of the vector
//of two vector addition.
MGUnit_vector& MGUnit_vector::operator+= (const MGVector& vec2){
	return *this=(*this+vec2);
}

// 単項マイナス。自身の単位ベクトルを反転したオブジェクトを生成する
MGUnit_vector MGUnit_vector::operator- () const{
	int dim=sdim();
	MGUnit_vector temp(*this);
	for(int i=0; i<dim; i++) temp.m_element[i]=-m_element[i];
	return temp;
}

//Update the unit vector by subtractiong vec2. The result is unit of
//the vector of two vector subtraction.
MGUnit_vector& MGUnit_vector::operator-= (const MGVector& vec2){
	return *this=(*this-vec2);
}
 
//Update own vector by vector product output, changes to 3D vector.
//The result is unit of two vector product.
MGUnit_vector& MGUnit_vector::operator*= (const MGVector& vec2){
	return *this=(*this*vec2);
}

///Get the unit normal of the triangle (P0, P1, P2).
///UnitNormal(P0,P1,P2)=-UnitNormal(P0,P2,P1).
///(V1,V2,UnitNOrmal) organizes orthonormal system, wher
///V1=P1-P0, V2=P2-P1.
MGUnit_vector UnitNormal(
	const MGPosition& P0, // 三角形の頂点の座標
	const MGPosition& P1,
	const MGPosition& P2
){
	MGVector vector1(P1 - P0);
	MGVector vector2(P2 - P0);
	return vector1*vector2;
}
