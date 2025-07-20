/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGSSisect_list_HH_
#define _MGSSisect_list_HH_
/** @file */
/** @addtogroup IsectContainer
 *  @{
 */

#include <list>
#include "mg/SSisect.h"

//Forward class declaration.
class MGFSurface;

/// MGSSisects defines linked list of MGSSisect.

/// This is value based list.
/// Used to represent intersection lines of two surfaces.
class MG_DLL_DECLR MGSSisects:public MGisects{
//We cannot use inheritance of std::list<MYELM> to make DLL.

public:

////////Special member functions/////////
~MGSSisects() = default;
MGSSisects(const MGSSisects&) = delete;
MGSSisects& operator=(const MGSSisects&) = delete;
MGSSisects(MGSSisects&&) = default;
MGSSisects& operator=(MGSSisects&&) = default;

//////////// Constructor ////////////
explicit MGSSisects(const MGFSurface *s1=NULL, const MGFSurface *s2=NULL);

/// Adds the MGSSisect to the end of the list.

///isect transfer the ownership of the curves in isect to this list.
void append(MGSSisect&& isect);
void append(MGSSisects&& isectlist);

///Add one intersection line to the list.

///iline, param1, and param2 must be newed objects, and their ownership
///are transfered to MGSSisects.
void append(
	MGCurve* iline,
	MGCurve* param1,
	MGCurve* param2,
	const MGSSRELATION r1=MGSSREL_UNKNOWN);

///Add one intersection line to the list.

///this append copies the three curves.
void append(
	const MGCurve& iline,
	const MGCurve& param1,
	const MGCurve& param2,
	const MGSSRELATION r1=MGSSREL_UNKNOWN);

///Find where in this ssi2  have common parts (in line_zero()) in 
///their world representation.

///Fucntion's return value is the iterator of this that had the common.
///		!=end():have common part. 
///		==end():no common part(except a point) found.
iterator find_common(const MGSSisect& ssi2);

///Return the pointer to surface1.
const MGFSurface* surface1() const;

///Return the pointer to surface2.
const MGFSurface* surface2() const;

/// Return the number of items that are in the list.
int entries() const{return int(size());};

///Insert MGSSisect at the index position i.

///This position must be between zero and the number of items in the list,
/// or behavior is undefined.
void insertAt(iterator i, MGSSisect&& isect){
	insert(i, std::unique_ptr<MGSSisect>(new MGSSisect(std::move(isect))));
};

///Return true (1) if there are no items in the list,
/// false(0) otherwise.
bool isEmpty() const{return empty();};

/// Return(but does not remove) lst element in the list.
/// If list is empty, behavior is undefined.
const MGSSisect& last() const{return *static_cast<const MGSSisect*>(back().get());};

/// Adds the MGSSisect to the beginning of the list.
///isect transfer the ownership of the curves in isect to this list.
void prepend(MGSSisect&& isect){
	emplace_front(new MGSSisect(std::move(isect)));
};

///Remove the MGSSisect and return the MGSSisect. If i is not valid, 
/// behavior is undefined.
std::unique_ptr<MGSSisect> removeAt(iterator i);

///Remove the first MGSSisect int the list and return the MGSSisect.
///If i is not valid, behavior is undefined.
std::unique_ptr<MGSSisect> removeFirst();

///Remove the last MGSSisect in the list and return the MGSSisect.
///If i is not valid, behavior is undefined.
std::unique_ptr<MGSSisect> removeLast();

/// Output virtual function.
std::ostream& toString(std::ostream& ostrm)const override;

};

/** @} */ // end of IsectContainer group
#endif
