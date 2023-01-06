/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"

#include "mg/PickObjectSB.h"
#include "topo/Face.h"
#include "topo/Edge.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//MGPickObjectSB is a MGPickObject that includes the perimeter information of
//a MGSurface.
//SB stands for surface boundary.
//MGPickObjectSB object is generated when users spedified 2-manifold and boundary
//selection, and the result is the boundary of a MGSurface.

/// MGPickObject is a class to locate where an object is in a group
/// hierarchy. Generally, A group includes other groups, and the included groups
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

MGPickObjectSB::MGPickObjectSB(const MGPickObjectSB& psb):MGPickObject(psb),
m_perimeter(psb.m_perimeter){
}

//Assignment operator.
MGPickObjectSB& MGPickObjectSB::operator=(const MGPickObject& pobj){
	MGPickObject::operator=(pobj);
	const MGPickObjectSB* psb=dynamic_cast<const MGPickObjectSB*>(&pobj);
	if(psb)
		m_perimeter=psb->m_perimeter;
	else
		m_perimeter=0;
	return *this;
}

////////////////オペレーション////////////

//Generate a newed clone object.
MGPickObjectSB* MGPickObjectSB::clone()const{
	return new MGPickObjectSB(*this);
}

//Return the face of the edge.
MGSurface* MGPickObjectSB::surface(){
	if(m_perimeter<0) return 0;
	return static_cast<MGSurface*>(leaf_object());
}

//Return the face of the edge.
const MGSurface* MGPickObjectSB::surface()const{
	if(m_perimeter<0) return 0;
	return static_cast<const MGSurface*>(leaf_object());
}
