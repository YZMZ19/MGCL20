/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGHHisect_vector_HH_
#define _MGHHisect_vector_HH_
/** @file */
/** @addtogroup IsectContainer
 *  @{
 */

#include <vector>
#include "mg/isects.h"
#include "topo/HHisect.h"
class MGShell;

//Forward class declaration.

///MGHHisects defines a vector of MGHHisect.

///The vector is implemeted using STL's vector.
///All the methods to handle the vector are available from the STL's vector class,
///and public member m_HHivector. Refer to STL vector class.
///MGHHisects is used to represent intersection lines of a shell with
///another shell, a face, or a surface.
class MG_DLL_DECLR MGHHisects:public MGisects{
public:

////////Special member functions/////////
~MGHHisects() = default;
MGHHisects(const MGHHisects&) = delete;
MGHHisects& operator=(const MGHHisects&) = delete;
MGHHisects(MGHHisects&&) = default;
MGHHisects& operator=(MGHHisects&&) = default;

///shel1,2 are the target shell for the intersection.
MGHHisects(const MGShell* shel1=nullptr, const MGShell* shel2=nullptr);

///Constructor of 1 MGHHisect.
explicit MGHHisects(MGHHisect&& hhi);

const MGShell* shell1()const;
const MGShell* shell2()const;

///Adds one MGHHisect to the end of the vector.
///Transfers the ownership of the curves in isect to this vector.
//void push_back(MGHHisect& isect){ std::vector<MGHHisect>::push_back(isect); };
void append(
	const MGFSurface* face1,	//face1. This must not be null.
	const MGFSurface* face2,	//face2. This may be null
								//(e.g. for face2 that is actually a surface).
	std::unique_ptr<MGSSisect>&& ssi);			//intersection line of face1 and face2 expressed as

private:

};

/** @} */ // end of IsectContainer group
#endif
