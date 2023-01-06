/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGGEL_TID_HH_
#define _MGGEL_TID_HH_

#include <utility>

/** @addtogroup GelRelated
 *  @{
 */

///Type id of subclasses of MGGel.

///The mm in 0xmm000000L of the following TID's are subclass id of MGGel.
///Let 0x0010nnxxL is a type number, then 
///nn of is the manifold dimension, xx is the name id.
///
///When a new subclass is necessary to add, this mm will have a new number.
enum MGGEL_TID{
	MGALL_TID =			0x00000000L,	///<all of specified kind of MGGEL_KIND.

	MGATTRIBEDGEL_TID=  0x00000000L,	///<all of MGAttribedGel.
	MGATTRIB_TID =		0x02000000L,	///<all of MGAttrib.
	MGPCELL_TID =       0x03200000L,	///<all of MGPCell.
	MGBCELL_TID =       0x04200000L,	///<all of MGBCell.

	MG0MANIFOLD=		0x00000000L,	///<all of 0 manifold, a point.
	MG1MANIFOLD=		0x00000100L,	///<all of 1 manifold, a curve.
	MG2MANIFOLD=		0x00000200L,	///<all of 2 manifold, a surface.
	MG3MANIFOLD=		0x00000300L,	///<all of 3 manifold, a solid.

//        //*************Define MGOBJECT_TID***************
//	Following nn is the manifold dimension, xx is the name id.
//						0x0010nnxxL   :mm=00
    MGOBJECT_TID =      0x00000000L,	///<all of MGAttribedGel.
    MGGEOMETRY_TID =	0x00100000L,	///<all of MGGeometry.

	MGPOINT_TID =		0x00100001L,	///<MGPoint type id.

	MGCURVE_TID =		0x00100100L,	///<all of MGCurve.
	MGSTRAIGHT_TID=		0x00100101L,	///<MGStraight type id.
	MGELLIPSE_TID =		0x00100102L,	///<MGEllipse type id.
	MGLBREP_TID =		0x00100103L,	///<MGLBRep type id.
	MGRLBREP_TID =		0x00100104L,	///<MGRLBRep type id.
	MGSRFCRV_TID =		0x00100105L,	///<MGSurfCurve type id.
	MGTRMCRV_TID =		0x00100106L,	///<MGTrimmedCurve type id.
	MGCOMPCRV_TID =		0x00100107L,	///<MGCompositeCurve type id.
	MGBSUMCRV_TID =		0x00100108L,	///<MGBSumCurve type id.

	MGSURFACE_TID =		0x00100200L,	///<all of MGSurface.
	MGPLANE_TID =		0x00100201L,	///<MGPlane type id.
	MGSPHERE_TID =		0x00100203L,	///<MGSphere type id.
	MGSBREP_TID =		0x00100205L,	///<MGSBRep type id.
	MGRSBREP_TID =		0x00100206L,	///<MGRSBRep type id.
	MGCYLINDER_TID =	0x00100207L,	///<MGCylinder type id.
	MGBSUMSURF_TID =	0x00100208L,	///<MGBSumSurf type id.
	MGPLANEIMAGE_TID=	0x00100209L,	///<MGPlaneImage type id.

/////////////////////////////////////////////

    MGTOPOLOGY_TID =   0x00200000L,	///<all of MGTopology related,
	                    //related means includes MGPCell, MGBCell.

//	Following nn is the manifold dimension, xx is the name id,
//						0x0k20nnxxL
    MGPVERTEX_TID =     0x03200001L,	///<MGPVertex type id.
    MGBVERTEX_TID =     0x04200001L,	///<MGBVertex type id.

//	Following nn is the manifold dimension, xx is the name id,
//	and m is 0 for cells, and m is 1 for complexes.
//						0x002mnnxxL
	MGCELL_TID=			0x00200000L,///< is a cell.
	MGCOMPLEX_TID=		0x00210000L,///< is a complex.

	MGEDGE_TID =		0x00200101L,	///<MGEdge type id.
	MGFACE_TID =		0x00200201L,	///<MGFace type id.
	MGSOLID_TID =		0x00200301L,	///<MGSolid type id(not yet implemented).

	MGLOOP_TID =		0x00210101L,	///<MGLoop type id.
	MGSHELL_TID =		0x00210201L,	///<MGShell type id.

///////////////////////////////////
//id for MGSurface or MGFace.
	MGFSURFACE_TID =	0x00000200L,	///<MGFSurface type id.

///////////////////////////////////
//id for MGStl.
    MGSTL_TID =         0x00300200L,	///<MGStl type id.

///////////////////////////////////
//id for MGroup.
    MGGROUP_TID =       0x00400000L,	///<all of MGGroup.


/////////////////////////////////////////////
//  MGAttrib id.

	MGGLATTRIBUTE_TID=	0x02010000L,
	MGAPPEARANCE_TID =	0x02010001L,	///<MGAppearance(attributes).

	MGCONTEXT_TID =		0x02010010L,	///<MGContext type id.
	MGLIGHTS_TID =		0x02010020L,	///<MGLights type id.
	MGLIGHT_TID =		0x02010030L,	///<MGLight type id.
	MGDIRECTIONAL_LIGHT_TID =
						0x02010031L,	///<MGDirectionalLight type id.
	MGPOINT_LIGHT_TID = 0x02010032L,	///<MGPointLight type id.
	MGSPOT_LIGHT_TID =	0x02010033L,	///<MGSpotLight type id.

	MGFOG_TID =			0x02010040L,	///<MGFog type id.

	MGMATERIAL_TID =	0x02010050L,	///<MGMaterial type id.
	MGALPHA_FUNC_TID =	0x02010060L,	///<MGAlphaFunc type id.
	MGBLEND_FUNC_TID =	0x02010070L,	///<MGBlenFunc type id.
	MGCOLOR_TID =		0x02010080L,	///<MGColor type id.
	MGCOLOR_MASK_TID =	0x02010090L,	///<MGColorMask type id.
	MGDEPTH_FUNC_TID =	0x020100A0L,	///<MGDepthFunc type id.
	MGDEPTH_MASK_TID =	0x020100B0L,	///<MGDepthMask type id.
	MGLIGHT_ENABLE_TID=	0x020100C0L,	///<MGLightEnable type id.
	MGLINE_STIPPLE_TID=	0x020100D0L,	///<MGLineStipple type id.
	MGLINE_WIDTH_TID =	0x020100E0L,	///<MGLineWidth type id.
	MGPOLYGON_MODE_TID=	0x020100F0L,	///<MGPolygonMode type id.
	MGRENDER_ATTR_TID=	0x02010100L,	///<MGRenderAttr type id.
	MGSHADE_MODEL_TID=	0x02010110L,	///<MGShade type id.
	MGTRANSP_MODE_TID=	0x02010120L,	///<MGTransp type id.
	MGTEXTURE_TID=		0x02010200L,	///<MGTexture type id.
	MGNAME_TID=			0x02010300L,	///<MGName type id.
};

///MGGEL_KIND_TID is used to specify what kind of group is used to identify gels.
enum MGGEL_KIND{
	MGALL_GELL =	0x00000000L,///<all of the gels
	MGTOP_KIND=		0xff000000L,///<subkind is MGOBJECT_TID, MGGROUP_TID, MGATTRIB_TID
	MGMANIFOLD=		0xff00ff00L,///<subkind is MG0MANIFOLD, MG1MANIFOLD, MG2MANIFOLD, MG3MANIFOLD
	MGFSURFACE_KIND=0xff0fff00L,///<subkind is MGFSURFACE_TID(MGFace or MGSurface, MG2MANIFOLD that is not
		///<MGShell)
	MGGEO_TOPO=		0xfff00000L,///<MGGeometry, or MGTopology.
		///<subkind is MGGEOMETRY_TID, MGTOPOLOGY_TID
	MGGEO_KIND=		0xffffff00L,///<subkind is MGPOINT_TID0, MGCURVE_TID, MGSURFACE_TID
	MGLEAF_KIND=	0xffffffffL	///<subkind is all of the leaf class MGGEL_TID, that is:
		///<MGPOINT_TID,MGSTRAIGHT_TID,MGELLIPSE_TID,MGLBREP_TID,MGRLBREP_TID,
		///<MGSRFCRV_TID,MGTRMCRV_TID,MGCOMPCRV_TID,MGPLANE_TID,MGSPHERE_TID,
		///<MGSBREP_TID,MGRSBREP_TID,MGCYLINDER_TID,MGPVERTEX_TID,
		///<MGBVERTEX_TID,MGEDGE_TID,MGFACE_TID,MGLOOP_TID,MGSHELL_TID
};

///MGAbstractGel is a class to specify what kind of abstract gel group.

///Let MGAbstractGel agel(gel_kind, gel_tid), then
///gel_kind is either MGALL_GELL, MGTOP_KIND, MGMANIFOLD, MGGEO_TOPO,
///MGGEO_KIND, or MGLEAF_KIND. And gel_tid is the specific type of
///the gel_kind. For the value of gel_tid, see MGGEL_KIND above.
///Possible combinations are defined in MGDefault.h(as mgAll_xxxx).
///See the definition.
typedef std::pair<MGGEL_KIND,MGGEL_TID> MGAbstractGel;

/** @} */ // end of GelRelated group
#endif
