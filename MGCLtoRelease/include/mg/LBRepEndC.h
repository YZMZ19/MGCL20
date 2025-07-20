/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGLBRepEndC_HH_
#define _MGLBRepEndC_HH_
#include "mg/MGCL.h"
#include "mg/Vector.h"

// MGLBRepEndC.h
//

class MGNDDArray;
class MGBPointSeq;
class MGCurve;
class MGIfstream;
class MGOfstream;

/** @file */
/** @addtogroup GEORelated
 *  @{
 */

/// Defines End Condition of Line B-Representation.

///End condition to get spline by interpolation.
///enum MGENDCOND {
///	MGENDC_UNKNOWN=0,	/// Unknown(usually not used).–¢’m
///	MGENDC_1D  =1,		/// 1st deravative provided.
///	MGENDC_2D  =2,		/// 2nd deravative provided.
///	MGENDC_NO  =3,		/// no end cond(only positional data)
///	MGENDC_12D =4		/// both 1st and 2nd deravatives provided.
///};
class MG_DLL_DECLR MGLBRepEndC {

public:

////////////Constructor////////////

///Default Constructor.
MGLBRepEndC():m_cond(MGENDCOND::MGENDC_NO){;};

///Construct of m_cond=MGENDC_1D or MGENDC_2D.
MGLBRepEndC(
	MGENDCOND cond,	///<Type of end condition
					///<(MGENDC_1D or MGENDC_2D).
	const MGVector& deriv	///<Derivative inf according to cond
);

///Construct of m_cond=MGENDC_12D
MGLBRepEndC(
	const MGVector& first_deriv,	///<1st derivative
	const MGVector& second_deriv	///<2nd derivative
);

/// Given data points ordinates and abscissa, compute approximate
/// 1st derivative.
MGLBRepEndC(
	int start,			///<Indicates start(start==true)
						///< condition or end.
	const MGNDDArray& tau,		///<Data point abscissa
	const MGBPointSeq& points,	///<Point seq data
	int &error				///<Error flag.
);

/// Given Positional data sequence, compute approximate
/// 1st derivative. The 1st derivative is unit vector.
MGLBRepEndC(
	int start,			///<Indicates start(start==true)
						///< condition or end.
	const MGBPointSeq& points	///<Point seq data
);

/// Given MGCurve, construct the curve's end condition.
MGLBRepEndC(
	int start,		///<Indicates start(start==true) condition or end.
	MGENDCOND cond,	///<Type of end condition(MGENDC_1D, MGENDC_2D, or MGENDC_12D)
	const MGCurve& curve///<Curve
);

///Debug Function.
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream&, const MGLBRepEndC& );

////////////Member Function////////////

///Get end condition.
MGENDCOND cond() const {return m_cond;}

///Get the 1st derivative.
const MGVector& first() const {return m_1deriv;}
MGVector& first() {return m_1deriv;}

///Get the 2nd derivative.
const MGVector& second() const {return m_2deriv;}
MGVector& second() {return m_2deriv;}

///Initialize the instance. Will be set to the same as constructed by the void 
///constructor.
void initialize();

///Set 1st deriv and change condition type to MGENDC_1D or MGENDC_12D.
void set_1st(const MGVector& first_deriv);

///Set 2nd deriv and change condition type to MGENDC_2D or MGENDC_12D.
void set_2nd(const MGVector& second_deriv);

///Dump Functions
int dump_size() const;

///Dump Function
int dump(MGOfstream& ) const;

///Restore Function
int restore(MGIfstream& );

///Member Data
private:

	MGENDCOND m_cond;	///<Type of end condition.
	MGVector  m_1deriv;	///<1st derivative stored 
						///<when m_cond=MGENDC_1D or MGENDC_12D
	MGVector  m_2deriv;	///<2nd derivative stored
						///<when m_cond=MGENDC_2D	or MGENDC_12D

};

/** @} */ // end of GEORelated group

#endif
