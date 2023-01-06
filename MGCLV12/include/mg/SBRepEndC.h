/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGSBRepEndC_HH_
#define _MGSBRepEndC_HH_
#include <assert.h>
#include <memory>
#include "mg/MGCL.h"
#include "mg/Vector.h"
#include "mg/BPointSeq.h"

// MGSBRepEndC.h
//

class MGNDDArray;
class MGCoons;
class MGBSumSurf;

/** @file */
/** @addtogroup GEORelated
 *  @{
 */

/// Defines End Condition of Surface B-Representation.

///Class for MGSBRep constructor like buildByInterpolationEC().
///Provides four perimeter end conditions and four corners end conditions.
///For a perimeter, four types of conditions can be applied:
///	MGENDC_NO:   no end cond(only positional data)
///	MGENDC_1D:   1st deravative provided.
///	MGENDC_2D:   2nd deravative provided.
///	MGENDC_12D:  both 1st and 2nd deravatives provided.
///When m_cond[0]=MGENDC_1D and m_cond[1]=MGENDC_1D, twist vector d2f/(du*dv)
///is necessary at the corner of u=max and v=min. Although users can provide
///the data, complete_corner_deriv computes the data from given data
///if not provided.
///At other corners or for other conditions, the same.
///
///To construct MGSBRepEndC, first construct through void constructor.
///Then use set_1st(), set_2nd, set_11d(), set_12d(), set_21d(), or
///set_22d() to complete. Any of set_11d(), set_12d(), set_21d(), or
///set_22d() can be omitted. complete_corner_deriv will compute instead.
///This is done in SBRep constructor. 
class MG_DLL_DECLR MGSBRepEndC {
private:
//////////// Member Data /////////
	MGENDCOND m_cond[4];	///< Type of end conditions of 4 perimeters.
	///< m_cond[0]: v=min boundary line, m_cond[1]: u=max boundary line
	///< m_cond[2]: v=max boundary line, m_cond[3]: u=min boundary line
	MGBPointSeq m_1deriv[4]; ///<1st derivative stored in m_1deriv[i]
							 ///<when m_cond[i]=MGENDC_1D or MGENDC_12D.
	MGBPointSeq m_2deriv[4]; ///<2nd derivative stored in m_2deriv[i]
							 ///<when m_cond[i]=MGENDC_2D or MGENDC_12D.

	MGVector m_11d[4];		///<Twist vector of four corners. d2f/(du*dv)
	///< [0]:u-min and v-min, [1]:u-max and v-min, 
	///< [2]:u-max and v-max, [3]:u-min and v-max.
	MGVector m_12d[4];		///< d3f/(du**2*dv) of four corners.
	MGVector m_21d[4];		///< d3f/(du*dv**2) of four corners.
	MGVector m_22d[4];		///< d4f/(du**2*dv**2) of four corners.
	int   m_sdim, m_nu, m_nv;	///<space dimension,
							///< length of u and v direction

public:

//////////// String stream Function. ////////////
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream&, const MGSBRepEndC& );


////////Special member functions/////////
MGSBRepEndC();
~MGSBRepEndC()=default;
MGSBRepEndC(const MGSBRepEndC&);///Copy constructor.
MGSBRepEndC& operator= (const MGSBRepEndC&);///Copy assignment.
MGSBRepEndC(MGSBRepEndC&&);		///Move constructor.
MGSBRepEndC& operator= (MGSBRepEndC&&);///Move assignment.


///Construct from a Coon's patch.
MGSBRepEndC(
	const MGNDDArray& utau,
	const MGNDDArray& vtau,
	const MGCoons& coons
);

///Construct from a Coon's patch.
MGSBRepEndC(
	const MGNDDArray& utau,
	const MGNDDArray& vtau,
	const MGBSumSurf& bssurf
);


//////////// Member Function ////////////

/// Complete corner derivatives(twist vectors and higher order derivatives)
void complete_corner_deriv
	(const MGNDDArray& utau, ///<u-direction data point
	const MGNDDArray& vtau ///<v-direction data point
);

///Return i-th perimeter's condition.

///i=0, 2 are v=min and max u-parameter line.
///i=1, 3 are u=max and min v-parameter line.
MGENDCOND cond(int i) const {assert(i<4); return m_cond[i];}

///Return i-th perimeter's first and second derivatives.
const MGBPointSeq& first(int i) const{assert(i<4);return m_1deriv[i];}
const MGBPointSeq& second(int i) const{assert(i<4);return m_2deriv[i];}

///Return i-th corner's derivative inf.

///i=0:(u-min, v-min), =1:(u-max, v-min),
/// =2:(u-max, v-max), =3:(u-min, v-max)
const MGVector& deriv11(int i) const{assert(i<4);return m_11d[i];}
const MGVector& deriv12(int i) const{assert(i<4);return m_12d[i];}
const MGVector& deriv21(int i) const{assert(i<4);return m_21d[i];}
const MGVector& deriv22(int i) const{assert(i<4);return m_22d[i];}

///Initialize the instance. 

///Will be set to the same as constructed by the void constructor.
void initialize();

///Set 1st deriv of i-th perimeter and change condition type
///to MGENDC_1D or MGENDC_12D.

///For i=0 and 2, first_deriv=df/dv,
///for i=1 and 3  first_deriv=df/du
///1st form is to copy first_deriv, and 2nd form is to transfer the ownership
///of the first_derivp to SBRepEndC.
///2nd form is recommended since no copy takes place.
void set_1st(int i, const MGBPointSeq& first_deriv);
void set_1st(int i, MGBPointSeq&& first_deriv);
void set_1st(int i, std::unique_ptr<MGBPointSeq>&& first_derivp);

///Set 2nd deriv of i-th perimeter and change condition type
///to MGENDC_2D or MGENDC_12D.

///For i=0 and 2, first_deriv=d2f/dv**2,
///for i=1 and 3  first_deriv=d2f/du**2
///1st form is to copy second_deriv, and 2nd form is to transfer the ownership
///of the second_derivp to SBRepEndC.
///2nd form is recommended since no copy takes place.
void set_2nd(int i, const MGBPointSeq& second_deriv);
void set_2nd(int i, MGBPointSeq&& second_deriv);
void set_2nd(int i, std::unique_ptr<MGBPointSeq>&& second_derivp);

///Set m_11d[i] inf(d2f/(du*dv)) of i-th corner.
void set_11d(int i, const MGVector& deriv);

///Set m_12d[i] inf(d3f/(du*dv**2)) of i-th corner.
void set_12d(int i, const MGVector& deriv);

///Set m_21d[i] inf(d3f/(du**2*dv)) of i-th corner.
void set_21d(int i, const MGVector& deriv);

///Set m_22d[i] inf(d4f/du**2*dv**2) of i-th corner.
void set_22d(int i, const MGVector& deriv);

private:
	void setLength(int i, const MGBPointSeq& derivs);

};

/** @} */ // end of GEORelated group

#endif
