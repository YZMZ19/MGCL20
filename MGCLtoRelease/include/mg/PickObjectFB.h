/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
/// MGPickObjectFB.h : MGPickObjectFB クラスの宣言およびインターフェイスの定義をします。
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef _MGPickObjectFB_HH_
#define _MGPickObjectFB_HH_

#include "mg/PickObject.h"
#include "mgGL/VBO.h"

class MGEdge;
class MGFace;

/** @addtogroup MGObjectRelated
 *  @{
 */

/// Is a class to locate where an object is in a group hierarchy.

/// Generally, A group includes other groups, and the included groups
/// include other groups. In this way, the groups make a group hierachy.
/// MGPickObject represents this hierarcy.
/// top_group() is the top MGGroup that includes
/// the object leaf_object() if m_Ghierarchy.size()==0. If m_Ghierarchy.size()>0,
/// top_group() includes m_Ghierarchy[0].
///	Let n=m_Ghierarchy.size(), then group m_Ghierarchy[i-1] includes
/// m_Ghierarchy[i] for n=0,...,n-2. m_Ghierarchy[n-1] includes leaf_object();
/// leaf_object() is the leaf MGObject pointer.
/// Although m_Ghierarchy[i] for i=0,...,n-2 are always MGGroup, m_Ghierarchy[n-1] may be
/// MGShell that includes MGFace. In this case, leaf_object() is the MGFace.
class MG_DLL_DECLR MGPickObjectFB:public MGPickObject{

public:

///Constructor.
MGPickObjectFB():MGPickObject(),m_edge(0){;};
MGPickObjectFB(const MGPickObjectFB& pfb);

///Conversion constructor from MGGelPosition and MGEdge.
MGPickObjectFB(
	MGGelPosition& gelp,
	const MGEdge* edge
):MGPickObject(gelp),m_edge(edge){;};

///Conversion constructor from MGPickObject and start/end.
MGPickObjectFB(
	MGPickObject& pobj,
	const MGEdge* edge
):MGPickObject(pobj),m_edge(edge){;};

///Assignment operator.
MGPickObjectFB& operator=(const MGPickObject& pobj);

///Generate a newed clone object.
MGPickObjectFB* clone()const;

///Return the edge pointer.
const MGEdge* edge()const{return m_edge;};

///Return the face of the edge.
MGFace* face();

///Highlightthe object using the display list of this object.
void hilight_using_display_list(
	int line_density	///<line density to draw a surface in wire mode.
)const;

///Set the object pointer.
void set_edge(const MGEdge* edge){m_edge=edge;};

private:

	const MGEdge* m_edge;	///MGEdge pointer of the face picked.
	mutable mgVBO m_vbo;//VBO to display the edge information.
};

/** @} */ // end of MGObjectRelated group
#endif
