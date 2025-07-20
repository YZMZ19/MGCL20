/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
/// SnapAttrib.h: MGSnapAttrib クラスのインターフェイス
#if !defined(AFX_SNAPATTRIB_H__3704A1F4_AAA3_4021_98F4_087C5A16006B__INCLUDED_)
#define AFX_SNAPATTRIB_H__3704A1F4_AAA3_4021_98F4_087C5A16006B__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif /// _MSC_VER > 1000

#include <bitset>
#include <iosfwd>
#include "mg/MGCL.h"
#include "mgGL/VBO.h"
#include "mgGL/Appearance.h"

class MGIfstream;
class MGOfstream;
/** @file */

/** @addtogroup GLAttrib
 *  @{
 */

///Defines Snap attributes.

///Snap means snap to designated location when inputting positional data
///from Cursor(generally mouse position).
///Currently MGCL supports the following enum's snap kind.
class MG_DLL_DECLR MGSnapAttrib{
public:

	enum{
		BIT_END=0,
		BIT_KNOT=1,
		BIT_NEAR=2,
		BIT_VERTEX=3,
		BIT_CENTER=4,
		BIT_GRID=5
		//BIT_ANGLE=6
	};

///Text output to stream.
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream& ostrm, const MGSnapAttrib& atr);

/// Serialization.
MG_DLL_DECLR friend MGOfstream& operator<< (MGOfstream& buf, const MGSnapAttrib& atr);
MG_DLL_DECLR friend MGIfstream& operator>> (MGIfstream& buf, MGSnapAttrib& atr);

	MGSnapAttrib();
	MGSnapAttrib(float apx, float apy,
		bool bEnd=false, bool bKnot=false, bool bNear=false, bool bVertex=false,
		bool bCenter=false, bool bGrid=false//, bool bOrtho=false
	);
	MGSnapAttrib(
		float apx, float apy,
		const std::bitset<32>& bits
	);

///	virtual ~MGSnapAttrib();

public:

	/// Aperture reference
	const float* getSnapAperture()const{return m_dSnapAperture;};
	void setSnapAperture(const float Aperture[2]){
		m_dSnapAperture[0] = Aperture[0];
		m_dSnapAperture[1] = Aperture[1];
	};

	///Set aperture.
	void setSnapAperture(double x,double y){
		m_dSnapAperture[0]=float(x);m_dSnapAperture[1]=float(y);};

	///Set/get end snap status.
	bool getEnd()  const {return m_bitset[BIT_END];};
	void setEnd(bool bEnd) {m_bitset[BIT_END] = bEnd;};

	///Set/get knot snap status.
	bool getKnot() const {return m_bitset[BIT_KNOT];};
	void setKnot(bool bKnot) {m_bitset[BIT_KNOT] = bKnot;};

	///Set/get near snap status.
	bool getNear() const {return m_bitset[BIT_NEAR];};
	void setNear(bool bNear) {m_bitset[BIT_NEAR] = bNear;};

	///Set/get vertex snap status.
	bool getVertex() const {return m_bitset[BIT_VERTEX];};
	void setVertex(bool bVertex) {m_bitset[BIT_VERTEX] = bVertex;};

	///Set/get center snap status.
	bool getCenter() const { return m_bitset[BIT_CENTER];};
	void setCenter(bool bCenter) {m_bitset[BIT_CENTER] = bCenter;};

	///Set/get grid snap status.
	bool getGrid() const { return m_bitset[BIT_GRID];};
	void setGrid(bool bGrid) {m_bitset[BIT_GRID] = bGrid;};

	///////////// Proxy interfaces.//////////

	///Test if any of the attributes are set on. If on, return true.
	bool any() const{ return m_bitset.any();}

	///Test if none of the attributes is set on. If none is on, return true.
	bool none() const{ return m_bitset.none();}

	///Get the number of attributes that are on.
	int count() const{ return int(m_bitset.count());}

	///clear all the attrib to false;
	void clear_attrib(){m_bitset.reset();};

///member data
private:
	std::bitset<32> m_bitset;
	float m_dSnapAperture[2];
};

/** @} */ // end of GLAttrib group
#endif // !defined(AFX_SNAPATTRIB_H__3704A1F4_AAA3_4021_98F4_087C5A16006B__INCLUDED_)
