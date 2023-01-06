/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGSBRepTP_HH_
#define _MGSBRepTP_HH_

#include <assert.h>
#include <memory>
#include <iosfwd>
#include "mg/LBRep.h"

class MGFSurface;
class MGSurface;

// MGSBRepTP.h
//

/** @file */
/** @addtogroup GEORelated
 *  @{
 */

///Defines Tangent Plane Line B-Representation Class.

///MGSBRepTP holds 4 tangent plane MGLBRep, each one of 
///that is the tangent plane of MGSBRep's perimeter.
///m_TP[i] is the one of the perimeter i of MGSBRep.
///Tangent plane is a line b-representation of (unit)normal vector
///along surface perimeter.
class MG_DLL_DECLR MGSBRepTP{

//////////// Member Data ////////////
	MGLBRep* m_TP[4];	///<Tangent Plane will be stored.
	///<Tangent plane m_TP[i] is a line b-representation of
	//<(unit)normal vector along the i-th perimeter.
	///<Parameter range of the TP is the same as u or v parameter
	///<range of the corresponding surface representation.
	///< m_TP[0]: v=min boundary line, m_TP[1]: u=max boundary line
	///< m_TP[2]: v=max boundary line, m_TP[3]: u=min boundary line
	///< Note
	///< Let MGVector f(t) is m_TP[i] and N=f(t), then tangent T at the
	///< constructing MGSBRep' perimeter is obtained by T=N*Td*N. Here Td is
	///< a given approximate tangent.

public:

///String stream Function.
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream&, const MGSBRepTP& );


////////Special member functions/////////
MGSBRepTP();///Default Constructor, so set that no TPs' are specified.
~MGSBRepTP();
MGSBRepTP(const MGSBRepTP&);///Copy constructor.
MGSBRepTP& operator= (const MGSBRepTP&);///Copy assignment.
MGSBRepTP(MGSBRepTP&&);		///Move constructor.
MGSBRepTP& operator= (MGSBRepTP&&);///Move assignment.

///Default Constructor, will be set as no TPs' are specified.

///Compute TP of four boundaries of a Surface B-Rep.
MGSBRepTP(const MGSurface& srf);

/// Build this from surrounding 4 perimeters of a surface.
void build(
	const MGLBRep* perims[4],//4 perimeters.
	const std::vector<const MGFSurface*> faces[4]//Tangent faces.
	      //faces[i] is faces touching the perim[i]. They may have gaps on the perim[i].
);
void build(
	const UniqueLBRep perims[4],//4 perimeters.
	const std::vector<const MGFSurface*> faces[4]//Tangent faces.
		  //faces[i] is faces touching the perim[i]. They may have gaps on the perim[i].
);

///Compute the maximum (absolute) cos value of between vector deris[i](t)
///and vector this->TP(i)(t) for i=0,1,2,3, where t is a common
///parameter of the data point obtained from deris[i]'s knot vector.
///Function's return value is the max out of cosmax[.].
double get_perimeters_max_cos(
	const UniqueLBRep deris[4],
	double taumax[4], ///< parameter on which the maximum value attains be stored.
	double cosmax[4]  ///< the maximum value be stored.
)const;

///Compute the maximum (absolute) sin value of between vector srf.normal(uv(t))
///and vector this->TP(i)(t) for i=0,1,2,3, where perim[i] is 
///the same as srf.perimeter_curve(i), and t is a common parameter
///of deris[i] and TP(i).
///Function's return value is the max out of sinmax[.].
double get_perimeters_max_sin(
	const MGSurface& srf,    ///< surface which must corresponds to this object.
	double         taumax[4],///< parameters on which the maximum value attains will be stored.
	double         sinmax[4],///< the maximum value will be stored.
	bool*          eval=0	///<indicates perimeters to evalate if eval!=null,
			///<When eval[i] is true, perimeter i is evaluated for 0<=i<=3.
)const;

///Compute maximun abs(cons(theta)), where theta=angle of TP(i) and  corresponding
///edge_crvl[i]'s start and end points' tangent vector.
double max_cos(
	const MGCurve*	perimeter[4]///<‹«ŠEüƒŠƒXƒg(vmin,umax,vmax,umin‚Ì‡,•Ó”Ô†0,1,2,3‚Ì‡)
)const;

///Compute maximun abs(cons(theta)), where theta=angle of TP(i) and  corresponding
///edge_crvl[i]'s start and end points' tangent vector.
double max_cos(
	const UniqueLBRep perimeters[4]///<‹«ŠEüƒŠƒXƒg(vmin,umax,vmax,umin‚Ì‡,•Ó”Ô†0,1,2,3‚Ì‡)
)const;

///Return if i-th perimeter's TP specified(true) or not.

///i=0, 2 are v=min and max u-parameter line.
///i=1, 3 are u=max and min v-parameter line.
bool specified(int i) const{assert(i<4); return m_TP[i]!=0;};

///Set i-th perimeter's TP as a null, as an unspecified one.
void set_TP_null(int i);

///Set i-th perimeter's TP(copy version).
void set_TP(int i, const MGLBRep& tp);

///Set i-th perimeter's TP(unique_ptr version).
void set_TP(int i, std::unique_ptr<MGLBRep>&& tp);

///Return i-th perimeter's TP.
const MGLBRep& TP(int i) const{assert(i<4);return *(m_TP[i]);}
MGLBRep& TP(int i) {assert(i<4);return *(m_TP[i]);}

///Get the TP.
MGLBRep** TP(){return m_TP;};

};

/** @} */ // end of GEORelated group

#endif