/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Default.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//  Implementation of class MGDefault
//
// Defines default values of each class.
//
extern MG_DLL_DECLR const MGPosition mgNULL_Pos(0);
extern MG_DLL_DECLR const MGPosition mgORIGIN(0.,0.,0.);
extern MG_DLL_DECLR const MGPosition mgORIGIN_2D(0.,0.);

extern MG_DLL_DECLR const MGVector mgX_UVEC(1.,0.,0.);
extern MG_DLL_DECLR const MGVector mgY_UVEC(0.,1.,0.);
extern MG_DLL_DECLR const MGVector mgZ_UVEC(0.,0.,1.);

extern MG_DLL_DECLR const MGVector mgX_UVEC_2D(1.,0.);
extern MG_DLL_DECLR const MGVector mgY_UVEC_2D(0.,1.);

extern MG_DLL_DECLR const MGInterval mgEMP_INTERV(MGINTERVAL_TYPE::MGINTERVAL_EMPTY);
extern MG_DLL_DECLR const MGInterval mgZERO_TO_DBLPAI(0.,mgDBLPAI);

extern MG_DLL_DECLR const MGBox mgNULL_BOX(0);
extern MG_DLL_DECLR const MGBox mgEMP_BOX = MGBox(
	MGInterval(),
	MGInterval(),
	MGInterval()
);
extern MG_DLL_DECLR const MGBox mgEMP_BOX_2D = MGBox(
	MGInterval(),
	MGInterval()
);

extern MG_DLL_DECLR const MGMatrix mgNULL_MATR(0);
extern MG_DLL_DECLR const MGMatrix mgUNIT_MATR(3,1.0);
extern MG_DLL_DECLR const MGMatrix mgUNIT_MATR_2D(2,1.0);

extern MG_DLL_DECLR const MGTransf mgNULL_TRANSF(0);
extern MG_DLL_DECLR const MGTransf mgID_TRANSF(MGMatrix(3,1.0), MGVector(0.,0.,0.));
extern MG_DLL_DECLR const MGTransf mgID_TRANSF_2D(MGMatrix(2,1.), MGVector(0.,0.));

extern MG_DLL_DECLR const MGVector mgNULL_VEC(0);
extern MG_DLL_DECLR const MGVector mgZERO_VEC(0.0, 0.0, 0.0);
extern MG_DLL_DECLR const MGVector mgZERO_VEC_2D(0.0, 0.0);

extern MG_DLL_DECLR const MGKnotVector mgNULL_KNOT_VECTOR(0);
extern MG_DLL_DECLR const MGKnotVector mgINFINITE_KNOT_VECTOR(2,2,-mgInfiniteVal, 2.*mgInfiniteVal);

//////////// Abstraction of Gell kind////////////
extern MG_DLL_DECLR const MGAbstractGel mgAll_Gell(MGALL_GELL,MGALL_TID);

extern MG_DLL_DECLR const MGAbstractGel mgAll_Group(MGTOP_KIND,MGGROUP_TID);
extern MG_DLL_DECLR const MGAbstractGel mgAll_Attrib(MGTOP_KIND,MGATTRIB_TID);

////////////////////////// MGObject ///////////////////
extern MG_DLL_DECLR const MGAbstractGel mgAll_Object(MGTOP_KIND,MGOBJECT_TID);

extern MG_DLL_DECLR const MGAbstractGel mgAll_0Manifold(MGMANIFOLD,MG0MANIFOLD);
extern MG_DLL_DECLR const MGAbstractGel mgAll_1Manifold(MGMANIFOLD,MG1MANIFOLD);
extern MG_DLL_DECLR const MGAbstractGel mgAll_2Manifold(MGMANIFOLD,MG2MANIFOLD);

extern MG_DLL_DECLR const MGAbstractGel mgAll_Geo(MGGEO_TOPO,MGGEOMETRY_TID);
extern MG_DLL_DECLR const MGAbstractGel mgAll_Topo(MGGEO_TOPO,MGTOPOLOGY_TID);
extern MG_DLL_DECLR const MGAbstractGel mgAll_STL(MGGEO_TOPO,MGSTL_TID);

extern MG_DLL_DECLR const MGAbstractGel mgAll_Point(MGGEO_KIND,MGPOINT_TID);
extern MG_DLL_DECLR const MGAbstractGel mgAll_Curve(MGGEO_KIND,MGCURVE_TID);
extern MG_DLL_DECLR const MGAbstractGel mgAll_Surface(MGGEO_KIND,MGSURFACE_TID);
extern MG_DLL_DECLR const MGAbstractGel mgAll_FSurface(MGFSURFACE_KIND,MGFSURFACE_TID);

extern MG_DLL_DECLR const MGAbstractGel mgAll_Straight(MGLEAF_KIND,MGSTRAIGHT_TID);
extern MG_DLL_DECLR const MGAbstractGel mgAll_Ellipse(MGLEAF_KIND,MGELLIPSE_TID);
extern MG_DLL_DECLR const MGAbstractGel mgAll_LBRep(MGLEAF_KIND,MGLBREP_TID);
extern MG_DLL_DECLR const MGAbstractGel mgAll_RLBRep(MGLEAF_KIND,MGRLBREP_TID);
extern MG_DLL_DECLR const MGAbstractGel mgAll_SurfCurve(MGLEAF_KIND,MGSRFCRV_TID);
extern MG_DLL_DECLR const MGAbstractGel mgAll_TrimmedCurve(MGLEAF_KIND,MGTRMCRV_TID);
extern MG_DLL_DECLR const MGAbstractGel mgAll_CompositeCurve(MGLEAF_KIND,MGCOMPCRV_TID);
extern MG_DLL_DECLR const MGAbstractGel mgAll_Plane(MGLEAF_KIND,MGPLANE_TID);
extern MG_DLL_DECLR const MGAbstractGel mgAll_SPhere(MGLEAF_KIND,MGSPHERE_TID);
extern MG_DLL_DECLR const MGAbstractGel mgAll_SBRep(MGLEAF_KIND,MGSBREP_TID);
extern MG_DLL_DECLR const MGAbstractGel mgAll_RSBRep(MGLEAF_KIND,MGRSBREP_TID);
extern MG_DLL_DECLR const MGAbstractGel mgAll_Cylinder(MGLEAF_KIND,MGCYLINDER_TID);
extern MG_DLL_DECLR const MGAbstractGel mgAll_PVertex(MGLEAF_KIND,MGPVERTEX_TID);
extern MG_DLL_DECLR const MGAbstractGel mgAll_BVertex(MGLEAF_KIND,MGBVERTEX_TID);
extern MG_DLL_DECLR const MGAbstractGel mgAll_Edge(MGLEAF_KIND,MGEDGE_TID);
extern MG_DLL_DECLR const MGAbstractGel mgAll_Face(MGLEAF_KIND,MGFACE_TID);
extern MG_DLL_DECLR const MGAbstractGel mgAll_Loop(MGLEAF_KIND,MGLOOP_TID);
extern MG_DLL_DECLR const MGAbstractGel mgAll_Shell(MGLEAF_KIND,MGSHELL_TID);
