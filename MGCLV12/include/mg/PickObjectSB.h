/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
/// MGPickObjectSB.h : MGPickObjectSB クラスの宣言およびインターフェイスの定義をします。
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef _MGPickObjectSB_HH_
#define _MGPickObjectSB_HH_

#include "mg/PickObject.h"
#include "mgGL/VBO.h"

class MGSurface;

/** @addtogroup MGObjectRelated
 *  @{
 */

///Is a MGPickObject that includes the perimeter information of a MGSurface.

///SB stands for surface boundary.
///MGPickObjectSB object is generated when users spedified 2-manifold and boundary
///selection, and the result is the boundary of a MGSurface.
/// MGPickObject is a class to locate where a picked object is in a group
/// hierarchy. Generally, A group includes other groups, and the included groups
/// include other groups. In that way the groups make a group hierachy.
/// MGPickObject represents this hierarcy, an MGObject or hierarchied MGGroup's.
/// When MGPickObject represents an MGObject, gel() returns MGObject
/// pointer and gel_is_object() returns true.
/// When MGPickObject represents an MGGroup, gel() returns MGGroup pointer,
/// and gel_is_object() returns false.
class MG_DLL_DECLR MGPickObjectSB:public MGPickObject{

public:

///Constructor.
MGPickObjectSB():MGPickObject(),m_perimeter(-1){;};
MGPickObjectSB(const MGPickObjectSB& psb);

///Conversion constructor from MGGelPosition and perimeter.
MGPickObjectSB(
	MGGelPosition& gelp,
	int perimeter
):MGPickObject(gelp),m_perimeter(perimeter){;};

///Conversion constructor from MGPickObject and start/end.
MGPickObjectSB(
	const MGPickObject& pobj,
	int perimeter
):MGPickObject(pobj),m_perimeter(perimeter){;};

///Assignment operator.
MGPickObjectSB& operator=(const MGPickObject& pobj);


///Generate a newed clone object.
MGPickObjectSB* clone()const;

///Highlightthe object using the display list of this object.
void hilight_using_display_list(
	int line_density	///<line density to draw a surface in wire mode.
)const;

///Return the edge pointer.
int perimeter()const{return m_perimeter;};

///Return the face of the edge.
MGSurface* surface();
const MGSurface* surface()const;

///Set the object pointer.
void set_perimeter(int perimeter){m_perimeter=perimeter;};

private:

	int m_perimeter;	///perimeter number of the MGSurface picked.
	mutable mgVBO m_vbo;//VBO to display the perimeter information.
};

/** @} */ // end of MGObjectRelated group

#endif
