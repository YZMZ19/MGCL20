/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGSBRepVecTP_HH_
#define _MGSBRepVecTP_HH_

#include <assert.h>
#include "mg/LBRep.h"
#include "mg/SBRepTP.h"

// MGSBRepVecTP.h
//

class MGSurface;
class MGLBRep;
class MGOfstream;
class MGIfstream;

/** @file */
/** @addtogroup GEORelated
 *  @{
 */

///Defines Tangent Plane Line B-Representation Class.

///Tangent plane is a line b-representation of (unit)normal vector of
///tangent plane along a surface perimeter.
///MGSBRepVecTP has a vector of 4 newed object pointer of MGLBRep representing the normal vetor
///B-rep along the 4 perimeter of surface.
class MG_DLL_DECLR MGSBRepVecTP{
	using MYELM = std::unique_ptr<MGLBRep>;
	using MYVEC = std::vector<MYELM>;

//////////// Member Data ////////////
	MYVEC m_TP[4];	///<Tangent Planes are stored,
			///<Tangent plane m_TP[i] is a vector of line b-representations of
			///<(unit)normal vector along the i-th perimeter,
			///<Parameter range of the TP is the same as u or v parameter
			///<range of the corresponding surface representation,
	///< m_TP[0]: v=min boundary line, m_TP[1]: u=max boundary line,
	///< m_TP[2]: v=max boundary line, m_TP[3]: u=min boundary line.

	MGInterval m_prange[4];///m_prange[j] is the parameter range of perimeter j

	bool m_to_SE[8];///m_to_SE[i*2] and [i*2+1] indicate if m_TP[i]'s parameter range covers
		///from the start point(i*2), or to the end point(i*2+1),
		///when m_to_SE[i*2] is true, m_TP[i] starts from the start point of perimeter i,
		///when m_to_SE[i*2+1] is true, m_TP[i] ends at the end point of the perimeter i,
		///when false, m_TP[i] starts from the inner point, or ends at the inner point. 

public:

///String stream Function.
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream& ostrm, const MGSBRepVecTP& vectp);


////////Special member functions/////////
MGSBRepVecTP()=default;//so set that no TPs' are specified.
~MGSBRepVecTP()=default;
MGSBRepVecTP(const MGSBRepVecTP&);///Copy constructor.
MGSBRepVecTP& operator= (const MGSBRepVecTP&)=delete;///Copy assignment.
MGSBRepVecTP(MGSBRepVecTP&&);	///Move constructor.
MGSBRepVecTP& operator= (MGSBRepVecTP&&);///Move assignment.


///Conversion constructor from ordinary MGSBRepTP.
MGSBRepVecTP(const MGSBRepTP& tp2);
MGSBRepVecTP(MGSBRepTP&& tp2);

///change the parameter range to (t0,t1).
void change_range(
	bool along_u,	///<the objective range is along u parameter or v.
	double t0,	///<target range, from t0,
	double t1	///< to t1.
);

///change the parameter range to prange.
void change_range(
	bool along_u,	///<the objective range is along u parameter or v.
	const MGInterval& prange///< target range.
);

///evaluate TP at the perimeter i's parameter t.
///Function's return value is:
///true if t was inside the parameter range of a tangent plane of m_TP[i].
///false if t was outside the parameter range of the m_TP[i] for all i.
bool eval(
	int i,		///<perimeter numeber.
	double t,		///<parameter vaule of perimeter i.
	MGVector& normal///<evaluated normal will be returned.
)const;

/// Compute the maximum (absolute) cos value of between vector deris[i](t) 
/// and vector this->TP(i)(t) for i=0,1,2,3, where t is a common
/// parameter of the data point obtained from deris[i]'s knot vector.
///Function's return value is the max out of cosmax[.].
double get_perimeters_max_cos(
	const UniqueLBRep deris[4],
	double taumax[4], ///< parameter on which the maximum value attains be stored
	double cosmax[4]  ///< the maximum value be stored
)const;

/// Compute the maximum (absolute) sin value of between vector srf.normal(uv(t))
/// and vector this->TP(i)(t) for i=0,1,2,3, where perim[i] is 
/// the same as srf.perimeter_curve(i), and t is a common parameter
/// of deris[i] and TP(i).
///Function's return value is the max out of sinmax[.].
double get_perimeters_max_sin(
	const MGSurface& srf,       ///< surface which must corresponds to this object
	double         taumax[4], ///< parameters on which the maximum value attains will be stored.
	double         sinmax[4],  ///< the maximum value will be stored.
	bool*          eval=0	///<indicates perimeters to evalate if eval!=null
			///<When eval[i] is true, perimeter i is evaluated for 0<=i<=3.
)const;

///Return true if at least one TP is specified at i-th perimeter.
///i=0, 2 are v=min and max u-parameter line.
///i=1, 3 are u=max and min v-parameter line.
bool specified(int i) const{assert(i<4); return m_TP[i].size()>0;};

///Set i-th perimeter's TP as a null, as an unspecified one.
void set_TP_null(int i);

///Set i-th perimeter's TP(std::vector version).
///vectp[i] must be newed objects, and all of the ownership will be transferer to
///this instance.
void set_TP(
	int i,		///<perimeter number.
	std::vector<UniqueLBRep>&& vectp,///<TP data.
	const MGInterval& prange		///<Whole perimeter's parameter range.
);

///Return i-th perimeter's TP.
const MYVEC& vecTP(int i) const{assert(i<4);return m_TP[i];}
MYVEC& vecTP(int i) {assert(i<4);return m_TP[i];}

MYVEC* vecTP(){return m_TP;};

};

/** @} */ // end of GEORelated group
#endif
