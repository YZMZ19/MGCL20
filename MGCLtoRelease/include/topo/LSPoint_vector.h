/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGLSPoint_vector_HH_
#define _MGLSPoint_vector_HH_
/** @file */
/** @addtogroup IsectContainer
 *  @{
 */

#include <vector>
#include "topo/LSPoint.h"

/// MGLSPoint_vector defines a vector of MGLSPoint.

/// Used to represent Intersection points of a loop and a surface.
class MG_DLL_DECLR MGLSPoint_vector{

public:

typedef std::vector<MGLSPoint>::iterator LSiterator;
typedef std::vector<MGLSPoint>::const_iterator const_LSiterator;

public:

///String stream Function
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream&, const MGLSPoint_vector& );

////////////Constructor/////////////
MGLSPoint_vector(){;};

/// Member Function.

MGLSPoint& operator[](size_t i){return m_lspoints[i];};
const MGLSPoint& operator[](size_t i)const{return m_lspoints[i];};

///Add one intersection point to the list.
bool append(const MGLSPoint& lsp);

/// Adds the MGLSPoint_vector to the end of the list.
void append(const MGLSPoint_vector& vec);

/// Return the number of items that are in the list.
int entries() const{return int(m_lspoints.size());};
int size() const{return int(m_lspoints.size());};

/// Return first element in the list.
/// If list is empty, behavior is undefined.
const MGLSPoint& front() const{return m_lspoints.front();};

///Insert MGLSPoint at the position i.
void insertAt(LSiterator i, const MGLSPoint& llisect)
{m_lspoints.insert(i, llisect);};

///Return true if there are no items in the list, false otherwise.
bool empty() const{return m_lspoints.empty();};

/// Return(but does not remove) last element in the list.
/// If list is empty, behavior is undefined.
const MGLSPoint& back() const{return m_lspoints.back();};
private:
	std::vector<MGLSPoint> m_lspoints;///<Vector of MGLSPoint.

};

/** @} */ // end of IsectContainer group
#endif
