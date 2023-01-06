/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

/** @file TPmaker.h
 *  Class mgTPmaker is internal for class MGCBRepTP to build.
 *
 *  Declaration for class %mgTPmaker.
 */
#if !defined(__mgTPmaker_H__)
#define __mgTPmaker_H__
#include <vector>
#include "mg/Vector.h"
#include "mg/LBRep.h"
#include "mg/PerimTP.h"

// Forward declaration.
class MGObject;
class MGCurve;
class MGFSurface;
class MGSBRepTP;
class MGSBRepVecTP;

/** \addtogroup GEORelated
 *  @{
 */

//! @brief A class which constructs an object of class %MGSBRepTP.
//!
//! Function @c make_tp() creates an object of this class. 
class MG_DLL_DECLR mgTPmaker{
	friend class mgPerimTP;

	// Data members.
private:
	const MGLBRep* m_perimeters[4]; // size is 4.
	const std::vector<const MGFSurface*>* m_faces[4];//four vectors of osculating faces.
	mgPerimTP  m_perimTP[4];
	MGVector m_cornerNormal[4];//Unit normal vector at the vertex i in [i].

public:

	//! Note that perim should be so processed in advance that is_valid_perim() returns true.
	mgTPmaker(
		const MGLBRep* perim[4],///< Perimeter curves of surface to build.
		const std::vector<const MGFSurface*> faces[4]///<tangent faces a perimeter i in faces[i].
	);

	//! @brief   Makes an object of class %MGSBRepTP.
	//! @param   result    Output.
	void make_tp(MGSBRepTP& result);

	///@cond
	/*接続四辺面用に、細切れのTP曲線を計算する関数*/
	void make_subtp(MGSBRepVecTP& vectp);
	///@endcond

private:
	mgTPmaker(const mgTPmaker&)=delete;
	mgTPmaker& operator=(const mgTPmaker&)=delete;

	//! @brief Obtains the effective normal vector.
	//!
	//! Retrieves the effective normal vector at a vertex.
	//! It is used in constructing TP curve.
	const MGVector& normal(int vertex) const {	return m_cornerNormal[vertex];}

};
/** @}*/ //GEORelated

#endif // __mgTPmaker_H__
