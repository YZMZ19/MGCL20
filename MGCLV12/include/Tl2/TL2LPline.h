#ifndef _mgTL2LPline_HH_
#define _mgTL2LPline_HH_

#include "Tl2/TL2parameter.h"
#include "Tl2/TL2Polyline.h"

/****************************************************************/
/*   Copyright (c) 2021 by System fugen G.K.                */
/*                       All rights reserved.                   */
/****************************************************************/

class MGPosition;
class mgTL2LPPoint;

/** @addtogroup UseTessellation
 *  @{
 */

///mgTL2LPline is limitted subinterval of mgTL2Polyline

/// , is a proprietry class for Face tessellation,
/// used to save copy of mgTL2Polyline, limitted subinterval of mgTL2Polyline,
/// holds the starting id of mgTL2Polyline and the number of vertices.
class mgTL2LPline{
	friend class mgTL2LPlines;
private:
	const mgTL2Polyline* m_polyline;//Original mgTL2Polyline pointer.
		//if(m_sharedLine) m_polyline=m_sharedLine.get(),
		//if(!m_sharedLine) m_polyline=mgTL2Face's Loop's edge's mgTL2Polyline pointer.
	std::shared_ptr<mgTL2Polyline> m_sharedLine;
	int m_idS;//starting id. Not relative, but absolute value of m_polyline.
	int m_nV;	//Number of vertices.
				//When m_nv>0, this mgTL2LPline's direction is the same as m_polyline's.
				//When m_nv<0, this mgTL2LPline's direction is opposite to m_polyline's.
		
	///Let m_polyline' vertices are (V0, V1, ... , Vn-1),
	///then (Vm_idS, Vm_idS+1, ... , Vm_ids+m_nv-1) are mgTL2LPline's polyline when m_nV>0,
	///     (Vm_idS, Vm_ids-1, ... , Vm_ids+m_nv+1) when m_nV<0.
	
	///when concave, 1 minus cangle, and when convex cangle minus 1.
	///cangle is cosine angle around Z-axis of (U,V) 2D coordinates plane
	///from pre-line's vector to aft-line's vector at the vertx i.
	///m_concavity's value is from -2 to 2.
	///2 is most concave(open), and -2 is most convex(closed).
	/// -2 means 180 degree convex, -1:90 degree convex, 0:flat
	/// 1:90 degree concave, 2:180 degree concave.
	mutable double m_concavity;

///string stream function
friend std::ostream& operator<< (std::ostream& ot, const mgTL2LPline& poly){
	return poly.toString(ot);
};

public:

///Default constructor.
mgTL2LPline():m_polyline(0),m_idS(0),m_nV(0),m_concavity(NotObtainedConcav){;};

///Copy constructor.
mgTL2LPline(
	const mgTL2LPline& lpline
):m_polyline(lpline.m_polyline), m_sharedLine(lpline.m_sharedLine)
,m_idS(lpline.m_idS),m_nV(lpline.m_nV)
, m_concavity(lpline.m_concavity){;};

//Constructor of whole mgTL2Polyline.
mgTL2LPline(
	const mgTL2Polyline* polyline//polyline must be the one of mgTL2Face's Loop's edge.
):m_polyline(polyline),m_idS(0),m_nV(polyline->bdim()), m_concavity(NotObtainedConcav){;};

mgTL2LPline(
	const mgTL2Polyline* polyline,//polyline must be the one of mgTL2Face's Loop's edge.
	int idS, //starting id. Not relative, but absolute value of m_polyline.
	int nV
):m_polyline(polyline),m_idS(idS),m_nV(nV), m_concavity(NotObtainedConcav) {
	assert(m_idS<m_polyline->bdim() && m_idS+m_nV+1>= 0);};

//Constructor of whole mgTL2Polyline.
mgTL2LPline(
	std::shared_ptr<mgTL2Polyline>& polyline//polyline must be the one of mgTL2Face's Loop's edge.
) :m_polyline(polyline.get()), m_sharedLine(polyline)
, m_idS(0), m_nV(polyline->bdim()), m_concavity(NotObtainedConcav) {
	;
};

mgTL2LPline(
	std::shared_ptr<mgTL2Polyline>& polyline,//polyline must be the one of mgTL2Face's Loop's edge.
	int idS, //starting id. Not relative, but absolute value of m_polyline.
	int nV
) :m_polyline(polyline.get()), m_sharedLine(polyline)
, m_idS(idS), m_nV(nV), m_concavity(NotObtainedConcav) {
	assert(m_idS < m_polyline->bdim() && m_idS + m_nV + 1 >= 0);
};

///Construct from subinterval of input lpline.
mgTL2LPline(
	const mgTL2LPline& lpline,
	int idS, //starting id of lpline
	int nV
);

///Construct from subinterval of input lpline.
mgTL2LPline(
	const mgTL2LPPoint& lpPoint,
	int nV
);

//////////// Member Function ///////////////

//Change this TL2LPline's point id to TL2Polyline's parameter t.
double changeIdToPolylineParameter(int id)const;

//Test if id is the id of the end point.
bool isEndPoint(int id)const{ return id==number_of_points()-1;};

//Test if id is the id of mid points.
bool isMidPoint(int id)const{ return id>0 && id<number_of_points()-1;};

///Construct new curve object by copying to newed area.
///User must delete this copied object by "delete".
mgTL2LPline* clone()const;

//Evaluation of this with the normalized parameter value t from 0. to 1.
//is provided. t=0. for the start point, t=1. for the end point.
MGVector eval(double t, int nderi=0)const;

//Get id of m_Bpolylines of m_param of m_polyline from vertex id i.
//Function's return value is
//true if the point of t is a boundary point, false if not.
bool get_id_from_VertexID(int i, short id[3])const;

///Get concavity of this edge.
///The concavity is obtained from the differece of two vectors,
///at the start and at the end point point tangent of this.
//Concavity's value is from -2 to 2. 2 is most concave, and
// -2 means 180 degree convex(most convex), -1:90 degree convex, 0:flat
// 1:90 degree concave, 2:180 degree concave.
double getConcavity()const;

///get the (u,v) parameter box.
MGBox getUvBox()const;

///Test if this is a null LPline.
bool is_null()const{return m_polyline==0;};
void setNull();

//Compute the intersections of sl and this mgTL2LPline.
//Function's return value is true if obtained.
// the (nearest) intersection vertex id of lp if >=0.
bool isectSlTl(
	const MGStraight& sl,
	double&  t//parameter value of the intersection.
)const;

//Update this by limiting the parameter range of the curve.
///Limitting is done at the knot parameter for both start and end.
void limit(
	int idS,	//start point id of this mgTL2LPline.
	int nV	//Number of vertices.
);

//Get the number of points of this closed polygon.
int number_of_points()const;

//Polygonize the (u,v) space straight line from this->uv(id1V) to pline2.uv(id2V).
//The direction of the output is from id1V to id2V.
//polygonizeSL does ensure the deviation from the surface within the surface
//tolerance.
void polygonizeSL(
	const mgTL2LPline& pline2,
	int id1V,	//id of this vertex.
	int id2V,	//id of pline2's vertex.
	mgTL2Polyline& uvpolyline
)const;

///Debug Function
std::ostream& toString(std::ostream& ostrm)const;

//Reverse the direction.
mgTL2LPline& negate();

//Subdivide at the id.
void subdivide(
	int id,	//Relative one that start from 0 even for opposite direction.
	mgTL2LPline& lp1,	//devided 1st mgTL2LPline.
	mgTL2LPline& lp2	//devided 2nd mgTL2LPline.
)const;

//Obtain the mid point of this line.
void mid(MGPosition& uvmid);

//Set the start point and the number of points to the whole line of m_polyline.
void setWholeLine();

std::shared_ptr<mgTL2Polyline>& sharedLine() { return m_sharedLine; };
const std::shared_ptr<mgTL2Polyline>& sharedLine() const { return m_sharedLine; };

const mgTL2Polyline* TL2Polyline()const{return m_polyline;};
const mgTL2parameter& TL2param()const{return m_polyline->TL2param();};


//Get i-th point surface parameter (u,v) of this polyline
MGPosition uv(int i)const;

//Get i-th point(x,y,z) of this polyline
MGPosition xyz(int i, bool need_normal=false)const;

void extractAsLBRep(MGLBRep& lb)const;

};

class mgTL2LPPoint {
public:
	const mgTL2LPline& m_line;
	const int m_id;
	mgTL2LPPoint(const mgTL2LPline& line, int id) :m_line(line), m_id(id) { ; }
	MGPosition uv()const { return m_line.uv(m_id); }
};

/// <summary>
/// Test if the value concavity is concave or not.
/// </summary>
inline bool isUndefinedConcavity(double concavity) {
	return concavity >= NotObtainedConcav;
};

/** @} */ // end of UseTessellation group
#endif
