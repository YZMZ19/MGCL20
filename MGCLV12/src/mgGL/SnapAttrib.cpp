/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
// SnapAttrib.cpp: MGSnapAttrib クラスのインプリメンテーション
//
//////////////////////////////////////////////////////////////////////
#include "StdAfx.h"
#include "mgGL/SnapAttrib.h"
#include "mg/Ifstream.h"
#include "mg/Ofstream.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//////////////////////////////////////////////////////////////////////
// 構築/消滅
//////////////////////////////////////////////////////////////////////

/*	enum{
		BIT_END=0,
		BIT_KNOT=1,
		BIT_NEAR=2,
		BIT_VERTEX=3,
		BIT_CENTER=4,
		BIT_GRID=5,
		//BIT_ANGLE=6
	};
*/

MGSnapAttrib::MGSnapAttrib()
:m_bitset(){
	m_dSnapAperture[0]=7.;
	m_dSnapAperture[1]=7.;
}

MGSnapAttrib::MGSnapAttrib(
	float apx, float apy,
	bool bEnd, bool bKnot, bool bNear, bool bVertex,
	bool bCenter, bool bGrid//, bool bAngle
):m_bitset(){
	m_dSnapAperture[0]=apx;
	m_dSnapAperture[1]=apy;
	m_bitset[BIT_END] = bEnd;
	m_bitset[BIT_KNOT] = bKnot;
	m_bitset[BIT_NEAR] = bNear;
	m_bitset[BIT_VERTEX] = bVertex;
	m_bitset[BIT_CENTER] = bCenter;
	m_bitset[BIT_GRID] = bGrid;
	//m_bitset[BIT_ANGLE] = bAngle;
}
MGSnapAttrib::MGSnapAttrib(
	float apx, float apy,
	const std::bitset<32>& bits
):m_bitset(bits){
	m_dSnapAperture[0]=apx;
	m_dSnapAperture[1]=apy;
}

//MGSnapAttrib::~MGSnapAttrib(){;}

//Text output to stream.
std::ostream& operator<< (std::ostream& ostrm, const MGSnapAttrib& atr){
	ostrm<<"SnapAttrib::"
		<< "End=" << atr.m_bitset[0] ;
	ostrm<< ", Knot=" << atr.m_bitset[1] << ", Near=" << atr.m_bitset[2];
	ostrm<< ", Vertex=" << atr.m_bitset[3] << ", Center=" << atr.m_bitset[4];
	ostrm<< ", Grid=" << atr.m_bitset[5] ;//<< ", Ortho=" << atr.m_bitset[6];
	ostrm<<", Aperture=("<< atr.m_dSnapAperture[0]<<","<< atr.m_dSnapAperture[1]<<")";
	return ostrm;
}

// Serialization.
MGOfstream& operator<< (MGOfstream& buf, const MGSnapAttrib& atr) {
	buf << atr.m_bitset.to_ulong();
	buf << atr.m_dSnapAperture[0];
	buf << atr.m_dSnapAperture[1];
	return buf;
}
MGIfstream& operator>> (MGIfstream& buf, MGSnapAttrib& atr) {
	long lbit;
	buf >> lbit;
	atr.m_bitset = std::bitset<32>(lbit);
	buf >> atr.m_dSnapAperture[0];
	buf >> atr.m_dSnapAperture[1];
	return buf;
}
