#ifndef _mgTL2Face_HH_
#define _mgTL2Face_HH_

#include <vector>
#include <map>
#include "Tl2/TL2parameter.h"
#include "Tl2/TL2Triangles.h"

/****************************************************************/
/*   Copyright (c) 2019 by System fugen G.K.                */
/*                       All rights reserved.                   */
/****************************************************************/

class MGCurve;
class MGFSurface;
class MGEdge;
class MGFace;
class MGLBRep;
class mgTL2Polyline;
class mgTL2LPline;
class mgTL2PlBridge;
class mgTLInputParam;

/** @file */
/** @addtogroup UseTessellation
 *  @{
 */

///UniqueTL2PlBridge definition.
using UniqueTL2PlBridge = std::unique_ptr<mgTL2PlBridge>;

///mgTL2Face is a proprietry class for Face tessellation.
class mgTL2Face{
private:
	mgTL2parameter m_param;///<Parameters of the tessellation.
	std::unique_ptr<MGFace> m_face;///<MGFace to tessellate that has only MGLoops without geometry.
							///The geometry is MGSurface of m_param.

public:

///Output to stream.
friend std::ostream& operator<< (std::ostream& out, const mgTL2Face& face);

//////////// constructor ///////////////
MG_DLL_DECLR mgTL2Face(
	const MGDrawParam& param,///<parameter for the tessellation.
	const MGSurface& face	///<テセレーションするフェイス
);
MG_DLL_DECLR mgTL2Face(
	const MGDrawParam& param,///<parameter for the tessellation.
	const MGFSurface& face,	///<テセレーションするフェイス
							///<Must be MGFace or MGSurface.
	const std::vector<SHLL_COM_EDGES>* polylines=0
		///< Input polygonized polylines of world coordinates for the face boundaries.
		///< polylines[i][j] is a j-th edge's polyline for face.loop(i),
		///< must be MGLBRep of order 2.
		///< polylines[i][j]=0 indicates loop i's edge j can be face's boundary and
		///< does not have any common edges.
		///< **polylines[i][j] must be the same direction as the faces's parameter edge.
);
MG_DLL_DECLR mgTL2Face(
	const MGFSurface& face,	///<テセレーションするフェイス
							///<Must be MGFace or MGSurface.
	double crvTol,			///<バウンダリのトレランス
	double surfTol,			///<平面とみなすトレランス
	double max_edge_len=-1.,///<when max_edge_len<=0, this means no limits on an edge length.
	const std::vector<SHLL_COM_EDGES>* polylines=0
		///< Input polygonized polylines for the face boundaries.
		///< polylines[i][j] is a j-th edge's polyline for face.loop(i),
		///< must be MGLBRep of order 2.
		///< polylines[i][j]=0 indicates loop i's edge j can be face's boundary and
		///< does not have any common edges.
		///< **polylines[i][j] must be the same direction as the faces's parameter edge.
);
MG_DLL_DECLR mgTL2Face(
	const mgTLInputParam& param,///<parameter for the tessellation.
	const MGFSurface& face	///<テセレーションするフェイス
);

mgTL2Face(const mgTL2Face& face);///<Copy constructor.


///Obtain the surface& of the tessellation target.
const MGSurface& surface()const{return m_param.get_surface();};

///Get mgTL2parameter.
mgTL2parameter& TL2param(){return m_param;};
const mgTL2parameter& TL2param()const{return m_param;};

///Get common boudary edges.
const std::vector<SHLL_COM_EDGES>* Bpoly()const{return m_param.Bpoly();};

///Do perform the tessellation.

///The result will be appended onto triangles.
///When triangles.is_uv()=false, all of the element of the triangle position data has normal data as
///(x,y,z,xn,yn,zn). Here (x,y,z) is the position data and (xn,yn,zn) is the normal vector
///at the position (x,y,z).
///When triangles.is_uv()=true, all of the element of the triange position data are (u,v).
void MG_DLL_DECLR tessellate(
	mgTL2Triangles& triangles	//Tessellated triangles will be output.
);

private:

	//////////// operator overload ////////////

mgTL2Face& operator= (const mgTL2Face& face)=delete;

///Polygonize all the boundaries of the target face,
///and make an MGFace that has the polygonized boundaries.
///The face made is m_face.
///The face does not have surface geometery, only has bounfaries.
void polygonizeBoundaries();

///Polygonize MGSurface boundaries,
///and make an MGFace that has the polygonized boundaries.
///The face made is m_face.
///The face does not have surface geometery, only has bounfaries.
///The target must be MGSurface.
void polygonizeSurfaceBoundaries();

};

///Express a splitting line for a face for tessellation.

///Private class for tessellation. Express a splitting line for a face,
///that starts from an MGLoop and that ends on an MGLoop.
///The both ends are expressed by m_le_start and m_le_end.
class mgTL2PlBridge{

friend std::ostream& operator<< (std::ostream& out, const mgTL2PlBridge& bridge);

public:
	mgTL2PlBridge(
		const MGLEPoint& le_start,//line's start point's MGLEPoint.
		const MGLEPoint& le_end//line's end point's MGLEPoint.
	);
	MGEdge* copy_edge();//Make a cop of m_polyLine and make MGEdge.
		//Function's return value is newed MGEdge pointer.
	MGEdge* release_edge();//Release m_polyLine object and make MGEdge.
		//Function's return value is newed MGEdge pointer.
		//After release_edge, m_polyLine is null.
	MGLEPoint& le_start(){return m_le_start;};
	MGLEPoint& le_end(){return m_le_end;};

private:
	std::unique_ptr<mgTL2Polyline> m_polyLine;//mgTL2Polyline whose start point is
		//m_le_start and whose end point is m_le_end;
	MGLEPoint m_le_start;//m_polyLine's start point's MGLEPoint.
	MGLEPoint m_le_end;//m_polyLine's end point's MGLEPoint.

};

///make_Edge() makes a polyline edge(parameter edge) of edge edgeuv which are wholly
///on a curve of m_param.Bpoly().
///All of the points of m_param.Bpoly()[id[0]][id[1]] will be converted to
///surface (u,v) parameter. These (u,v) representation makes the polyline edge.
MGEdge* make_Edge(
	const mgTL2parameter& tlpara,
	const MGEdge& edgeuv,
	short id[3],	//pointer to uniform LBRep of m_param.Bpoly(), must not be null.
					//On return, id of the starting point will be return.
	mgTL2Polyline*& poly		//generated mgTL2Polyline* for the edge will be returned.
);

///split face at the parameter t.

///Function's return value is true: if split was executed.
///false: if not.
bool splitTl(
	const MGFace& face,///<Target face to split.
	const mgTL2parameter& tlparam,
	std::stack<UniqueFace>& faces///<All of the splitted faces will be output.
);

//Split face obtaining bridges using getBridgeAtMiddlePara().
//Function's return value is true, if split was executed, false, if not.
bool splitByBridges(
	int dir, //=0: divide at u param, 1: at v param
	const MGFace& face,//Target face to split.
	const mgTL2parameter& tlparam,
	std::stack<UniqueFace>& faces//All of the splitted faces will be appended.
);

/** @} */ // end of UseTessellation group

#endif
