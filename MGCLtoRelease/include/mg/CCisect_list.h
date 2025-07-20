/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGCCisect_list_HH_
#define _MGCCisect_list_HH_
/** @file */
/** @addtogroup IsectContainer
 *  @{
 */

#include <list>
#include "mg/CCisect.h"
#include "mg/isects.h"

//Forward class declaration.
class MGCurve;

/// Defines a list of MGCCisect(curve to curve intersection).

///Used to represent Intersection points of two curves.
class MG_DLL_DECLR MGCCisects:public MGisects{
//We cannot use inheritance of std::list<MYELM> to make DLL.
private:
	double m_error;			///< Square of Tolerance in parameter space.

public:

/// Constructor
explicit MGCCisects(const MGCurve *c1=NULL, const MGCurve *c2=NULL);

////////////Member Function.////////////

/// Adds the MGCCisect to the end of the list.
//void append(
//    const MGCCisect& isect///<isect to append.
//);

/// 交点の全てのコンポーネントを指定して，交点リストに追加
///Add one intersection point to the list.
void append(
	const MGPosition& point,	///<Intesection point(x,y,z)
	double t1,					///<parameter value of curve 1.
	double t2,					///<parameter value of curve 2.
	const MGCCRELATION r1=MGCCREL_UNKNOWN///<Input relation
);

/// Adds the MGCCisects to the end of the list.
void append(MGCCisects&& lst);

///Return the pointer to curve1.
const MGCurve* curve1() const;

///Return the pointer to curve2.
const MGCurve* curve2() const;

/// Return the number of items that are in the list.
int entries() const{return int(size());};

/// Return(but does not remove) first element in the list.
/// If list is empty, behavior is undefined.
MGCCisect& first(){return isectCast<MGCCisect>(begin());};

/// Adds the MGCCisect to the beginning of the list.
void prepend(MGCCisect&& isect);

///Remove the first MGCCisect int the list and return the MGCCisect.
///If i is not valid, behavior is undefined.
std::unique_ptr<MGCCisect> removeFirst();

///Remove the last MGCCisect in the list and return the MGCCisect.
///If i is not valid, behavior is undefined.
std::unique_ptr<MGCCisect> removeLast();

/// Output virtual function.
std::ostream& toString(std::ostream& ostrm)const override;

};

/** @} */ // end of IsectContainer group
#endif
