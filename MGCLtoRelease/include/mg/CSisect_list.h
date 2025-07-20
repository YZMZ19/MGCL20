/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGCSisect_list_HH_
#define _MGCSisect_list_HH_

/** @addtogroup IsectContainer
 *  @{
 */

#include <list>
#include "mg/CSisect.h"
#include "mg/isects.h"

//Forward class declaration.
class MGCurve;
class MGFSurface;

/// MGCSisects defines linked list of MGCSisect.

/// Used to represent Intersection points of a curve and a surface.
class MG_DLL_DECLR MGCSisects:public MGisects{
private:
	double m_errort;///<error to regard same curve point in parameter space.
	double m_erroru;///<error to regard same surface point in u-parameter space.
	double m_errorv;///<error to regard same surface point in v-parameter space.

public:

////////////// Constructor////////////
explicit MGCSisects(const MGCurve *crv=nullptr, const MGFSurface *srf=nullptr);

////////////// Member Function. ////////////

// Adds the MGCSisect to the end of the list.
//End points will be prefered.
void append(MGCSisect&& isect);

/// 全てのコンポーネントを指定して交点を生成する
void append(
		const MGPosition& point,		///<intersection point.
		double t,				///<Curve's parameter value.
        const MGPosition& uv,	///<Surface's parameter values.
		const MGCSRELATION rl=MGCSREL_UNKNOWN
								///<Curve and Surface relation
	);

void append(MGCSisects&& list);

///Return the pointer to the curve.
const MGCurve* curve() const;

///Return the pointer to the surface.
const MGFSurface* surface() const;

/// Return the number of items that are in the list.
int entries()const{	return (int)size();};

/// Return(but does not remove) first element in the list.
/// If list is empty, behavior is undefined.
MGCSisect& first(){return isectCast<MGCSisect>(begin());};

/// Return(but does not remove) last element in the list.
/// If list is empty, behavior is undefined.
const MGCSisect& last() const{return *static_cast<const MGCSisect*>(back().get());};

///Remove the parameter and return the parameter. If i is not valid, 
/// behavior is undefined.
std::unique_ptr<MGCSisect> removeAt(iterator i);

///Remove the first MGCSisect int the list and return the MGCSisect.
///If i is not valid, behavior is undefined.
std::unique_ptr<MGCSisect> removeFirst();

///Remove the last MGCSisect in the list and return the MGCSisect.
///If i is not valid, behavior is undefined.
std::unique_ptr<MGCSisect> removeLast();

/// Output virtual function.
std::ostream& toString(std::ostream& ostrm)const override;

};

/** @} */ // end of IsectContainer group
#endif
