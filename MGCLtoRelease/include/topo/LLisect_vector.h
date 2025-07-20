/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGLLisect_vector_HH_
#define _MGLLisect_vector_HH_
/** @file */
/** @addtogroup IsectContainer
 *  @{
 */

#include <vector>
#include "topo/LLisect.h"

//Forward class declaration.
class MGCurve;
class MGLoop;
class MGInterval;

/// MGLLisect_vector defines a vector of MGLLisect.

/// Used to represent Intersection points of two loops.
class MG_DLL_DECLR MGLLisect_vector:public std::vector<MGLLisect>{

public:

typedef std::vector<MGLLisect>::iterator LLiterator;
typedef std::vector<MGLLisect>::const_iterator const_LLiterator;
typedef std::vector<MGLLisect>::iterator iterator;
typedef std::vector<MGLLisect>::const_iterator const_iterator;

public:

///String stream Function
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream&, const MGLLisect_vector& );

////////////Constructor/////////////
MGLLisect_vector();
MGLLisect_vector(const MGLoop& loop);		///Loop

///Add one intersection point to the list.
void append(const MGLLisect& lli);
void append(
	const MGPosition& uv,	///<Parameter (u,v) of the parent face.
	const MGLEPoint& lp1,	///<First loop's point data.
	const MGLEPoint& lp2	///<Second loop's point data.
);

/// Adds the MGLLisect_vector to the end of the list.
void append(const MGLLisect_vector& list);

/// Return the number of items that are in the list.
int entries() const{return (int)size();};

/// Return first element in the list.
/// If list is empty, behavior is undefined.
const MGLLisect& first() const{return front();};

///Insert MGLLisect at the position i.
void insertAt(LLiterator i, const MGLLisect& llisect){insert(i, llisect);};

/// Return(but does not remove) last element in the list.
/// If list is empty, behavior is undefined.
const MGLLisect& last() const{return back();};

private:
	double m_error_square;		///<Error allowed to compute isect and 
								///<end point coincidence.

};

/** @} */ // end of IsectContainer group
#endif
