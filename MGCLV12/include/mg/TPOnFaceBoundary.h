/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

/** @file TPOnFaceBoundary.h
 *  mgTPOnFaceBoundary is an internal class for mgTPmaker to construct an MGSBRepT object.
 *
 *  Declaration for class %mgTPOnFaceBoundary.
 */
#if !defined(__mgTPOnFaceBoundary_H__)
#define __mgTPOnFaceBoundary_H__

#include <iostream>
#include <vector>
#include "mg/TrimmedCurve.h"
#include "mg/SBRepTP.h"

 /** \addtogroup GEORelated
  *  @{
  */

//! @brief A class which represents a tangent plane line brep(unit normal vector line brep).

//! mgTPOnFaceBoundary constitutes mgPerimTP.
//!
//! %mgTPOnFaceBoundary provides the following two facilities.
//! - To hold the parameter ranges of both a piece of the boundary curve
//! of a tangent MGFSurface and the associated perimeter curve to build MGSBRep.
//! - To build tangent plane MGLBRep curve. The parameter range is same as
//! a common-part of the perimeter curve(m_periPart).
class mgTPOnFaceBoundary{
	friend std::ostream& operator<<(std::ostream&, const mgTPOnFaceBoundary&);

// Member data
private:
	MGTrimmedCurve m_periPart;//!< The common part of the original perimeter to construct MGSBRep.
	MGLBRep m_tp;//!<tangent plane line lbrep is set whose parameter range is the same as m_periPart

	const MGFSurface& m_f; //!< tangent MGFSurface of constructing MGSBRep at m_periPart.
	double m_boundaryS, m_boundaryE;//!< The parameter range of m_boundaryParam.
	    //m_boundaryS corresponds to m_periPart.param_s() and m_boundaryE to m_periPart.param_e().
	    //m_boundaryS>m_boundaryE takes place when m_boundaryParam and m_periPart have opposite directions.
	const MGCurve& m_boundaryParam;//!< The parameter (u,v) rep of the common boundary of m_f.
	bool m_necessaryToNegate;   //!< indicates if the normal of m_f is necessary to negate.

public:
////////Special member functions/////////
	mgTPOnFaceBoundary()=delete;///Default Constructor.
	~mgTPOnFaceBoundary()=default;
	mgTPOnFaceBoundary(const mgTPOnFaceBoundary&)=delete;///Copy constructor.
	mgTPOnFaceBoundary& operator= (const mgTPOnFaceBoundary&)=delete;///Copy assignment.
	mgTPOnFaceBoundary(mgTPOnFaceBoundary&&) = default;		///Move constructor.
	mgTPOnFaceBoundary& operator= (mgTPOnFaceBoundary&&) = default;///Move assignment.

// The Constructor
	mgTPOnFaceBoundary(
		const MGInterval& perimRange,//!< parameter range of perimcrv,
		   //! the part of the perimeter curve(makes the paramter range of this mgTPOnFaceBoundary).
		const MGLBRep& perimcrv,//!< perimeter curve of the target MGSBRep to build.

		const MGFSurface& face, //!< Tangent MGSurface or MGFace at the range(ts,te) of bndcrv,
		double ts, double te,   //!< parameter range of bndcrv(a boundary of face).
			  //ts>te when bndcrv and perimcrv have the opposite directions.
		const MGCurve&    bndcrv,//!< face's boundary curve that has common part with perimcrv. 
		const MGCurve&    paramcrv,//!< The parameter (u,v) rep of bndcrv.

		bool face_negated = false  //!< whether face normal is necessary to negate for the tp data.
	);

// Member functions

    //! Build tangent plane line lbrep(m_tp).
	void build();

	/// Outputs the state of this object.
	void toString(std::ostream& os = std::cout) const;

	/// Returns the tangent plane lbrep.
	MGLBRep& tp(){ return m_tp;}
	const MGLBRep& tp() const{ return m_tp; }

	/// Returns the unit normal vector at the paramete uv of m_f.
	/// Note that the direction may be negated according to m_necessaryToNegate.
	MGVector unit_normal_at(const MGPosition& uv) const;

	///Get the start parameter of m_periPart.
	double param_sPeri(){return m_periPart.param_s();}

	///Get the end parameter of m_periPart.
	double param_ePeri(){ return m_periPart.param_e(); }

private:
	void getStoreTPdata(
		double tPeri, //parameter of m_periPart to get normal of m_f.
		double tBoundary,//parameter of m_boundaryParam to get normal of m_f.
		int i,//position of tpd to store
		MGBPointSeq& tpd
	);
};
/** @}*/ //GEORelated

#endif // __mgTPOnFaceBoundary_H__
