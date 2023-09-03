#ifndef _mgTL2parameter_HH_
#define _mgTL2parameter_HH_

#include <vector>
#include "mg/MGCL.h"
#include "mg/drawParam.h"
#include "topo/Face.h"

/****************************************************************/
/*   Copyright (c) 2019 by System fugen G.K.                */
/*                       All rights reserved.                   */
/****************************************************************/

class MGObject;
class MGFace;
class MGFSurface;
class MGSurface;
class MGDrawParam;

/** @file */
/** @addtogroup UseTessellation
 *  @{
 */

#define EDGE_LENGTH_DENOM 8.///<Default edge length parameter, which is the denominator of
		///the face object box length. Used in compute_max_edge_len().
#define NEAR_PARAM 0.002///< Parameter of the tessellation to avoid the same point division.
		///Data within parameter_span*NEAR_PARAM is regarded same point. Used in find_where().
#define MAX_DEVIATION_FROM_MIDDLE .26///<find_where() parameter to avoid extraneous edge split.
		///Used when split parameter value(u or v) is computed.
		///Split at the edge vertex value when the parameter value is 
		///within span_length*MAX_DEVIATION_FROM_MIDDLE from the middle.

#define SIN5DEGREE 0.087156	///= sin(5.0 degree)
#define LOOSE_ZERO_ANGLE  0.20		///= 11.45 degree

///Concavity data of mgTLPLines' m_concavity.
#define STRICT_CONCAVITY 0.0603	/// = 1-cos(20 degree)(160 degree open out of 360).
#define COS_LOOSE_ZERO_ANGLE_M1 -.0000061/// = cos(.2 degree)-1(180.2 degree open).
#define COS_LOOSE_ZERO_ANGLE2_M1 -.003805/// = cos(5 degree)-1(185. degree open).
#define COS_SHARPANGLE_M1 -1.866025	     /// = cos(150 degree)-1(320 degree open).

#define LARGE_RATIO 3.
#define SMALL_RATIO 2.
	/// Used in splitTl, defines ratio of u and v-param range.
	/// When the ratio is greater than LARGE_RATIO, the greater param range 
	/// is dvided.

#define NotObtainedConcav 3.

///Holds necessary parameter data for face tessellation.

///mgTL2parameter is a proprietry class for Face tessellation.
///In the constructor of mgTL2parameter(const MGFace&, double, double),
///all the parameters are initialized.
class mgTL2parameter{
public:

friend std::ostream& operator<< (std::ostream& out, const mgTL2parameter& para);

///Default Constructor.
mgTL2parameter()=default;

///Constructor that specifies each parameter.
mgTL2parameter(
	const MGFSurface& obj,	///<テセレーションするフェイス
							///<Must be MGFace or MGSurface.
	double crvTol,	///<Tessellation curve tolerance.
	double surfTol,///<Tessellation surface tolerance.
	const std::vector<SHLL_COM_EDGES>* polylines=0,
		///< Input polygonized polylines for the face boundaries.
		///< polylines[i][j] is a j-th edge's polyline for face.loop(i),
		///< must be MGLBRep of order 2.
		///< polylines[i][j]=0 indicates loop i's edge j can be face's bounday and
		///< has any common edges.
		///< **polylines[i][j] must be the same direction as the faces's parameter edge.
	double max_edge_len=-1.///<Maximum edge length of the final subdivided rectangles.
);

///Constructor of MGDrawParam.
mgTL2parameter(
	const MGFSurface& obj,	///<テセレーションするフェイス
							///<Must be MGFace or MGSurface.
	const MGDrawParam& param,///<parameter for the tessellation.
	const std::vector<SHLL_COM_EDGES>* polylines=0
		///< Input polygonized polylines for the face boundaries.
		///< polylines[i][j] is a j-th edge's polyline for face.loop(i),
		///< must be MGLBRep of order 2.
		///< polylines[i][j]=0 indicates loop i's edge j can be face's bounday and
		///< has any common edges.
		///< **polylines[i][j] must be the same direction as the faces's parameter edge.
);

mgTL2parameter(const mgTL2parameter& param2);

const std::vector<SHLL_COM_EDGES>* Bpoly()const{return m_Bpolylines;};
const MGFace& get_face()const{return *m_face;};
const MGSurface& get_surface()const{return *m_surface;};
double get_max_edge_len()const{return m_max_edge_len;};
double get_max_edge_len_sqr()const{return m_max_edge_len_sqr;};
double get_tess_crvError()const {return m_tess_crvError;};
double get_tess_srfError()const {return m_tess_srfError;};
double get_UError()const {return m_puerror;};
double get_VError()const {return m_pverror;};
double get_UVError()const{return m_uverror;};
bool target_is_face()const{ return m_face!=0;};

//Update the curve error. Returned is the old error.
double set_tess_crvError(double error);

//Update the surface error. Returned is the old error.
double set_tess_srfError(double error);

/////////Operator oveload/////////

///Assignment.
mgTL2parameter& operator=(const mgTL2parameter&);

private:

	const MGFace* m_face{ nullptr };	///<Original face to tessellate.
	const MGSurface* m_surface{ nullptr };///<Original surface to tessellate.
							///<If m_face!=0, m_surface=m_face->surface();
	const std::vector<SHLL_COM_EDGES>* m_Bpolylines{ nullptr };
		///< Input polygonized polylines for the face boundaries.
		///< polylines[i][j] is a j-th edge's polyline for face.loop(i),
		///< must be MGLBRep of order 2.
		///< polylines[i][j]=0 indicates loop i's edge j can be face's bounday and
		///< has any common edges.
		///< **polylines[i][j] must be the same direction as the faces's parameter edge.
	double m_puerror, m_pverror;///<Parameter error used for intersection computation.
	double m_uverror;///is sqrt(m_puerror*m_puerror+m_pverror*m_pverror);
	double m_tess_crvError;	///<Tessellation curve tolerance.
	double m_tess_srfError;	///<Tessellation surface tolerance.
	double m_max_edge_len;	///<Maximum edge length of the tessallated rectangles.
	double m_max_edge_len_sqr;	///<Square of m_max_edge_len.

void build_parameter(
	const MGFSurface& srf,
	double crvTol,
	double surfTol,
	double max_edge_len
);

};

///Compute maximum edge length for the tessellation from an object, twoManifold.
///twoManifold must be MGFSurface or MGShell.
double compute_max_edge_len(const MGObject& twoManifold);

/// <summary>
/// Get start(start=true, at the point of Vrtx+) or
/// end(start=false, at the point of Vrtx-) direction vector of the point Vrtx.
/// </summary>
const MGVector& getSEDirection(const MGLEPoint& Vrtx, bool start=true);

///Get the knotvector to guarantee the maximum line segment length square
///is less than maxElen2.
void getXYZline_ensuring_max_edge_length(
	double maxElen2,	///<square of maximum edge length.
	const MGLBRep& xyzpolyline,///<Input original LBRep of (x,y,z) of order 2.
	MGLBRep& xyzpolylineOut///<Output LBRep of order 2 whose knot vector is:
			///t(i)=i-1 for i=1,...,n and t(0)=0 and t(n+1)=n-1.
);

///Construct polyline MGLBRep whose maximum edge length is para.get_max_edge_len_sqr().
void getUVline_ensuring_max_edge_length(
	const mgTL2parameter& para,///<Input parameter.
	const MGCurve& uvline,///<Input original (u,v) line of the surface.
	MGLBRep& xyzpolylineOut///<Output LBRep of order 2 whose knot vector is:
			///<t(i)=i-1 for i=1,...,n and t(0)=0 and t(n+1)=n-1.
);

///test if concav is loosely concave or not. sideNum is the number of the sides.
inline
bool isLooselyConcave(double concav){return concav>=COS_LOOSE_ZERO_ANGLE2_M1;};

///test if angle is loosely concave or not.
inline
bool isConcave(double concav){return concav>=COS_LOOSE_ZERO_ANGLE_M1;};

///test if angle is sharp or not.
inline
bool isSharp(double concav){return concav<COS_SHARPANGLE_M1;};

//Find the 1st concave vertex of the outer loop lp from the edge number startEdge.
//Function's return value is the edge num found. When not found, -1 is returned.
int findConcaveVertex(
	const MGLoop& lp,//Target loop.
	int startEdge=0//Start edge number to find.
);

/** @} */ // end of UseTessellation group
#endif
