/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGCFisect_vector_HH_
#define _MGCFisect_vector_HH_
/** @file */
/** @addtogroup IsectContainer
 *  @{
 */

#include <vector>
#include "mg/CSisect_list.h"
#include "topo/CFisect.h"

///MGCFisects defines a vector of MGCFisect.

///The vector is implemeted using STL's vector.
///All the methods to handle the vector are available from the STL's vector class,
///and public member m_CFivector. Refer to STL vector class.
///MGCFisects is used to represent intersection points of a curve with a shell.
///The behavior of MGCFisect is like an auto_ptr. Copy or assignment
///of MGCFisect means transfer of the ownership of all the included curves
///to copied or assigned MGCFisect and original MGCFisect does not have the
///ownership more. Users should be aware of this fact.
class MG_DLL_DECLR MGCFisects:public MGCSisects{

public:

//////////// Constructor ////////////

///Void constructor.
MGCFisects(const MGCurve* curve=0):MGCSisects(curve){;};

///Constructor of 1 MGCFisect.
MGCFisects(const MGCurve* curve,MGCFisect&& cfi);


//////////// Member Function. ////////////

///Adds one MGCFisect to the end of the vector.
///Transfers the ownership of the curves in isect to this vector.
void push_back(MGCFisect&& isect);
void emplace_back(const MGCSisect& is, const MGFSurface& f);

///Return(but does not remove) last element in the vector.
///If vector is empty, behavior is undefined.
const MGCFisect& backCF()const {
	return *static_cast<const MGCFisect*>(back().get());};
MGCFisect& backCF() {
	return *static_cast<MGCFisect*>(back().get());};

/// Return(but does not remove) first element in the vector.
/// If vector is empty, behavior is undefined.
const MGCFisect& frontCF() const{
	return *static_cast<const MGCFisect*>(front().get());};
MGCFisect& frontCF(){
	return *static_cast<MGCFisect*>(back().get());};

///Insert MGCFisect at the index position i.
void insertAt(iterator i, MGCFisect&& isect);

/// Output virtual function.
std::ostream& toString(std::ostream& ostrm)const override;

};

/** @} */ // end of IsectContainer group
#endif
