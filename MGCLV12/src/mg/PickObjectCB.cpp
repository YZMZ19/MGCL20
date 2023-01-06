/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"

#include "mg/PickObjectCB.h"
#include "mg/Curve.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//MGPickObjectCB is a MGPickObject that includes the boundary information of
//a MGCurve.
//CB stands for curve boundary.
//MGPickObjectCB object is generated when users spedified 1-manifold and boundary
//selection.

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


MGPickObjectCB::MGPickObjectCB(const MGPickObjectCB& pcb):MGPickObject(pcb),
m_start_end(pcb.m_start_end){
}

//Assignment operator.
MGPickObjectCB& MGPickObjectCB::operator=(const MGPickObject& pobj){
	MGPickObject::operator=(pobj);
	const MGPickObjectCB* pcb=dynamic_cast<const MGPickObjectCB*>(&pobj);
	if(pcb)
		m_start_end=pcb->m_start_end;
	else
		m_start_end=0;
	return *this;
}

////////////////オペレーション////////////

//Generate a newed clone object.
MGPickObjectCB* MGPickObjectCB::clone()const{
	return new MGPickObjectCB(*this);
}

//Return the face of the edge.
const MGCurve* MGPickObjectCB::curve()const{
	if(m_start_end<0) return 0;
	return static_cast<const MGCurve*>(leaf_object());
}
