/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGCoons_HH_
#define _MGCoons_HH_

#include <iosfwd>
#include <memory>
#include "mgGL/VBO.h"

// MGCoons.h

// Forward Declaration
class MGCurve;
class MGLBRep;
class MGSPointSeq;

/** @file */
/** @addtogroup GEORelated
 *  @{
 */

/// Defines Coons Patch surface.
class MG_DLL_DECLR MGCoons {
 
public:

///String stream Function
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream&, const MGCoons& );

//////////// Constructor ////////////

///The ownership of perimeters, and derivatives are transfered to this.
MGCoons(
	std::unique_ptr<MGLBRep> perimeters[4],
	std::unique_ptr<MGLBRep> derivatives[4]
);
MGCoons(
	std::unique_ptr<MGCurve> perimeters[4],
	std::unique_ptr<MGCurve> derivatives[4]
);

//////////// Member Function ////////////

///get the derivative data d2f/((dv)(du)) at corner i.
const MGVector& d2fdvu(int i)const;

///Return perimeter i's derivative data.
const MGCurve& derivative(int i)const{return *(m_derivatives[i]);};

///Evaluate surface data.
///Currently ndu=ndv=0 is assumed.
MGVector eval(
	double u,///< Parameter value (u,v) of the surface,
	double v///< must be 0<=u,v<=1.
	, int ndu=0	///< Order of derivative along u.
	, int ndv=0	///< Order of derivative along v.
	) const;

///Evaluate at (ui, vj) data for (i, j).
///Here ui=utai(i) 0<=i<utau.length(), and for j, the same.
void eval(
	const MGNDDArray&	utau,		///<u方向のデータポイント
	const MGNDDArray&	vtau,		///<v方向のデータポイント
	MGSPointSeq&		spoint///<evaluated data will be output to spoint.
)const;

///Get space dimension.
int sdim()const;

//////////// Member Data ////////////

private:

	std::unique_ptr<MGCurve> m_perimeters[4];///<perimeters.
		///<m_perimeters[0]:v-min, [1]:u-max, [2]:v-max, [3]:u-min.
		///<the parameter ranges of all the curves must be [0,1].
	std::unique_ptr<MGCurve> m_derivatives[4];///<Derivatives along perimeters.
		///<m_derivatives[0]:along perimeters[0](df/dv(u,0))
		///<[1]:along perimeters[1](df/du(1,v))
		///<[2]:along perimeters[2](df/dv(u,1))
		///<[3]:along perimeters[3](df/du(0,v))

	///Data at each corner
	MGVector m_f00, m_f01, m_f10, m_f11;
	MGVector m_dfdu00, m_dfdu01, m_dfdu10, m_dfdu11;
	MGVector m_dfdv00, m_dfdv01, m_dfdv10, m_dfdv11;
	MGVector m_df2duv00, m_df2dvu00;
	MGVector m_df2duv01, m_df2dvu01;
	MGVector m_df2duv10, m_df2dvu10;
	MGVector m_df2duv11, m_df2dvu11;

void eval_corner();
MGVector eval_inner(
	double u, double v	///< Parameter value of the surface.
						///< must be 0<=u,v<=1.
	) const;
};

/** @} */ // end of GEORelated group
#endif
