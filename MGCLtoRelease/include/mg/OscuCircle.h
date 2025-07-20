/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGOscuCircle_HH_
#define _MGOscuCircle_HH_
/** @file */
/** @addtogroup GEORelated
 *  @{
 */

#include "mg/OscuCircleData.h"
#include <vector>

// MGOscuCircle.h
//

//Forward Declaration
class MGIfstream;
class MGOfstream;

/// Defines Array of OscuCircle data.

///This class is used for MGLBRep constructor:
///	MGLBRep(						///BLGCS
///			const MGLBRepEndC& begin,	///Begin end condition
///			const MGLBRepEndC& end,		///End end conditoion
///			const MGBPointSeq& points,	///Point seq data
///			const int* point_kind,		///Point kind of above point.
///			const MGOscuCircle& circle,	///Provides osculating circle data.
///			int &error);				///Error flag.
///Defines array of MGOscuCircleData which provides index of points and
///radius of the osculating circles.
class MG_DLL_DECLR MGOscuCircle {

public:

//////////// Constructor ////////////

/// Dummy constructor, setts m_n=0
MGOscuCircle():m_n(0){;};

///Reference i-th osculating circle data.
const MGOscuCircleData& operator ()(int i) const;

//////////// Member Function ////////////

///Add to the end of list.
MGOscuCircle& add(const MGOscuCircleData&);	

///Add to the end of list.
MGOscuCircle& add(int index, double radious);

///Remove i-th OscuCircleData.
MGOscuCircleData remove(int i);

///Get the number of data.
int length() const {return m_n;};

///Debug Function.
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream&, const MGOscuCircle& );

///Dump Functions
int dump_size() const;

///Dump Function
int dump(MGOfstream& ) const;

///Restore Function
int restore(MGIfstream& );

//////////// Member Data ////////////

private:
	int m_n;	///<Number of stored OscuCircle data.
	std::vector<MGOscuCircleData> m_circle;	///<OscuCircle data list.

};

/** @} */ // end of GEORelated group
#endif
