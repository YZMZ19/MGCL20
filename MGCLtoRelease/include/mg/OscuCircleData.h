/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGOscuCircleData_HH_
#define _MGOscuCircleData_HH_
/** @file */
/** @addtogroup GEORelated
 *  @{
 */
#include "mg/MGCL.h"

// MGOscuCircleData.h
//

class MGIfstream;
class MGOfstream;

///The class for MGLBRep constructor of osculating circles.

///This class is used for MGLBRep constructor:
///	MGLBRep(						///BLGCS
///			const MGLBRepEndC& begin,	///Begin end condition
///			const MGLBRepEndC& end,		///End end conditoion
///			const MGBPointSeq& points,	///Point seq data
///			const int* point_kind,		///Point kind of above point.
///			const MGOscuCircle& circle,	///Provides osculating circle data.
///			int &error);				///Error flag.
/// Defines OscuCircle data, index of BPointSeq points and circle radius.
class MG_DLL_DECLR MGOscuCircleData {

public:

///String stream Function.
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream&, const MGOscuCircleData& );

////////////Constructor////////////

///Default Constructor
MGOscuCircleData(){;};
MGOscuCircleData(
	int index,	///<Index of BPointSeq indicating where
				///<circle be inserted.
	double radius	///<Radius of circle.
);

bool operator< (const MGOscuCircleData& ocd2)const;
bool operator> (const MGOscuCircleData& ocd2)const{return ocd2<(*this);};
bool operator<= (const MGOscuCircleData& ocd2)const{return !(ocd2<(*this));};
bool operator>= (const MGOscuCircleData& ocd2)const{return !((*this)<ocd2);};
bool operator== (const MGOscuCircleData& ocd2)const;
bool operator!= (const MGOscuCircleData& ocd2)const{return !operator==(ocd2);};
	
///Get the index.
int index() const {return m_index;}

///Get the radius.
double radius() const {return m_radius;}

///Dump Functions
int dump_size() const;

///Dump Function
int dump(MGOfstream& ) const;

///Restore Function
int restore(MGIfstream& );

////////////Member Data//////////
private:
	int m_index;			///<Index of points.
	double	m_radius;		///<radius of oscurating circle.

};

/** @} */ // end of GEORelated group
#endif
