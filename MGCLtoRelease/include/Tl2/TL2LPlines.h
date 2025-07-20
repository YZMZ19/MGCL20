#include "StdAfx.h"
#include "Tl2/TL2LPline.h"

class MGLoop;
class mgTL2Triangles;
class mgTL2Polyline;

/****************************************************************/
/*   Copyright (c) 2021 by System fugen G.K.                    */
/*                       All rights reserved.                   */
/****************************************************************/

//Utility class for mgTL2LPlines to indicate a point on a edge of a mgTL2LPlines.
class mgTLEdgePoint{
public:
	int m_edgeID; //Edge id 
	int m_pointID;//Point id in m_EdgeID

public:
mgTLEdgePoint():m_edgeID(-1), m_pointID(0){;};//-1 means undefined.
explicit mgTLEdgePoint(int edgeID, int pointID=0):m_edgeID(edgeID), m_pointID(pointID){;};
bool operator==(const mgTLEdgePoint& epid)const{
	return m_edgeID==epid.m_edgeID && m_pointID==epid.m_pointID;};
bool undefined()const{return m_edgeID==-1;};
bool defined()const { return m_edgeID != -1; };
};

//Utility class for mgTL2LPlines, to analyze concavity at vertices of mgTL2LPlines.
class mgTLConcavPid{
public:
	double concavity;//concavity at the point epID, see m_concavity of mgTL2LPlines.
	mgTLEdgePoint epID;

mgTLConcavPid(double cnocav=3., int eid=0, int pid=0)
	:concavity(cnocav),epID(eid,pid){;};
bool operator<(const mgTLConcavPid& cpid2)const{return concavity<cpid2.concavity;};

//Compare this with the input 3 mgTLConcavPid's, then if this is min or max, replace.
void getMinMax(
	mgTLConcavPid& maxConcav,//most concave point.
	mgTLConcavPid& max2Concav,//2ndly most concave point.
	mgTLConcavPid& maxConvex//most convex point.
);

};

/// class of a vector of mgTL2LPline.
/// To tessellate, construct this, then invoke tessellate.
class mgTL2LPlines{

public:

#define MAX_LINE_NUM 5
	enum SHARPorCONCAVE {
		CONCAVE = -1,
		SHARP = 1,
		OTHER = 0
	};

	std::stack<std::unique_ptr<mgTL2LPlines>>& m_stack;
	mgTL2Triangles& m_triangles;//mgTL2Triangles for the tessellated triangles to output.

private:
	mgTL2LPline m_plines[MAX_LINE_NUM];
	mutable short m_nlines;//Initially -1, then turned to side num.
	mutable short m_npoints;//Initially -1, then turned to num of point number.
	mutable short m_minConcavVid, m_maxConcavVid;//minimum conavity vertex id , and maximum id.
		//When =-1, not obtained yet.
	mutable double m_concavity[MAX_LINE_NUM];//concavity at the vertex i in [i].
		//vertex i is the start point of the edge m_plines[i].
		//Concavity is:
		//when concave, 1 minus cangle, and when convex cangle minus 1.
		//cangle is cosine angle around Z-axis of (U,V) 2D coordinates plane
		//from pre-line's vector to aft-line's vector at the vertx i.
		//m_concavity's value is from -2 to 2.
		//2 is most concave(open), and -2 is most convex(closed).
		// -2 means 180 degree convex, -1:90 degree convex, 0:flat
		// 1:90 degree concave, 2:180 degree concave.

public:
//////////// constructor ///////////////

//default constructor.
mgTL2LPlines()=delete;

mgTL2LPlines(
	std::stack<std::unique_ptr<mgTL2LPlines>>& stack,
	mgTL2Triangles& triangles//mgTL2Triangles for the tessellated triangles to output.
);

///Construct from an outer loop that has less than or equal to MAX_LINE_NUM edges.
mgTL2LPlines(
	std::stack<std::unique_ptr<mgTL2LPlines>>& stack,
	const MGLoop& oloop,
		//Outer boundary loop of the face of mgTL2Plygon
		//to tessellate whose number of edges is MAX_LINE_NUM at most.
	mgTL2Triangles& triangles	//mgTL2Triangles for the tessellated triangles to output.
);

const mgTL2LPline& operator[](int i)const{return m_plines[i];};
mgTL2LPline& operator[](int i){return m_plines[i];};

///Analyze the numbers of this, getting the number of 2 point edge, maximum and
///minimum point number edges.
void analyzePointNumber(
	int eids[5],
		//Following data are output:
		//eids[0]:num of 2 point edges,
		//eids[1], [2];//minimum number of vertices and the edge id
		//eids[3], [4];//Maximum vertex number and the edge id
	int np[5]//number of points of edge i is returned in np[i].
)const;

///Check the concavities of all the vertices, concavity or sharp angle.
///Function's return value is SHARPorCONCAVE:
//       CONCAVE if concavity found.
///      SHARP if concavity not found, and sharp angle found.
///      OTHER if both concavity ro sharp angle not found.
///concaveVID and convexVid are always returned regardless of functon's return value.
SHARPorCONCAVE analyzeConcavity(
	int& concaveVID,//most concave vertex id returned,
	int& convexVid	//most convex vertex id returned.
)const;

//Analyse concavity of this polygon,
//not only the vertices but the intermediate points also .
//Function's return value is the number of concave points.
int analyzeConcavityAllPoints(
	mgTLConcavPid& pidMostConcave,//most concave concavity and the point.
	mgTLConcavPid& pid2ndConcave,//2ndly most concave concavity and the point
	mgTLConcavPid& pidMostConvex//most convex concavity and the point
)const;

//Nullify this.
void setNull();

//Obtain how many m_plines[] are effective.
int getNumberOfLines()const{
	return m_nlines<0 ? getNumberOfLinesSub():m_nlines;
};

//Obtain how many points in this closed polygon..
int numberOfPoints()const{
	return m_npoints<0 ? getNumberOfPointsSub():m_npoints;
};

//Obtain how many points in m_plines[i].
int nPLine(int i)const{return m_plines[i].number_of_points();};

//Convert point id to epid(line id, id in the line).
//Here, point id is the point number from the start point of m_plines[0] to
//the point before the end point of the last of m_plines[]. 
void convertPointIDToLinePointID(int pid, mgTLEdgePoint& epid)const;

//Normalize the point id epid(eid, Pid), i.e. 
//if is the end point id, change it to the start point id of the next side.
void normalizePointID(mgTLEdgePoint& epid)const;

//Compare vid1's and vid2's openess(concavity),
//then if vid1's is more open, return true.
bool isMoreOpen(int vid1, int vid2)const;
bool isMoreClosed(int vid1, int vid2)const{return isMoreOpen(vid2,vid1);};

bool isUndefinedConcav(int vid)const{
	return ::isUndefinedConcavity(m_concavity[vid]);
};//>=NotObtainedConcav?

//Test if the edge eid is concave or not.
bool isConcaveEdge(int eid)const;

//Test if the vertex vid is loosely flat(true), or not.
bool isLooselyConcaveVertex(int vid)const;

/// <summary>
/// get 1st concave edge out of m_plines.
/// Function's return value is:
///  -1 if not found.
/// edge id if found.
/// Concavity is tested for all edge's start and end point tangent difference.
/// (Not concavity at a vertex)
/// </summary>
int get1stConcaveEdge()const;

//Get the most open vetex id.
int getMostOpenVid()const{return getMaximumConcavVid();};

//Return true if epid is the next point of epid2 around this closed polygon.
bool isNextPoint(const mgTLEdgePoint& epid, const mgTLEdgePoint& epid2)const;

//Return true if epid is the previous point of epid2 around this closed polygon.
bool isPrePoint(const mgTLEdgePoint& epid, const mgTLEdgePoint& epid2)const;

// Get the previous(decremented) point of epid.
// Function's return value is the decremented point around this mgTL2LPlines polygons.
mgTLEdgePoint decrement(
	const mgTLEdgePoint& epid,//input mgTLEdgePoint to decrement.
	int num=1
)const;

// Get the next(incremented) point of epid.
// Function's return value is the incremented point around this mgTL2LPlines polygons.
mgTLEdgePoint increment(
	const mgTLEdgePoint& epid,//input the pivot point.
	int num=1
)const;

///get the (u,v) parameter box.
MGBox getUvBox()const;

//Evaluate all the points into uv[], then
//evaluate all the point uv dirrerence data into dir[].
//Function's return value is the number of points evaluated.
int evalUVVec(
	MGPosition uv[],//(u.v) data
	MGVector dir[]//uv[i+1]-uv[i] vectors are output in [i].
)const;

///Get concavity of an edge eid, which is not at a vertex.
///That is, concavity is obtained from the differece of two vectors,
///at the start and at the end point point tangent of the same edge eid.
//Concavity's value is from -2 to 2. 2 is most concave, and
// -2 means 180 degree convex(most convex), -1:90 degree convex, 0:flat
// 1:90 degree concave, 2:180 degree concave.
double getEdgeConcavity(int eid)const{return m_plines[eid].getConcavity();};

///Get concavity at the vertex vid, which is not an edge's concavity.
///That is, the concavity is obtained from the difference of two vectors
///at the edge id vid's start point and at the previous edge's end point.
//Concavity's value is from -2 to 2. 2 is most concave, and
// -2 means 180 degree convex(most convex), -1:90 degree convex, 0:flat
// 1:90 degree concave, 2:180 degree concave.
double getVertexConcavity(int vid)const;

//Get the vertx id of minimum concavity(most convex, or most closed).
int getMinimumConcavVid()const;

//Get the vertx id of maximum concavity(most concave, or most open).
int getMaximumConcavVid()const;

//Obtain from point to subdivide.
//When nlines is 2 or 3, there are no 2 point edges.
//When nlinesis 4, 2point edge number is at most 1.
mgTLEdgePoint getFromPointToSubdivide(
	const int eids[5]
		//Following data are input(that are obtained by analyzePointNumber):
		//eids[0]:num of 2 point edges,
		//eids[1], [2];//minimum number of vertices and the edge id
		//eids[3], [4];//Maximum vertex number and the edge id
)const; 

//Get subdivide to-point, providing from point.
mgTLEdgePoint getToPointToSubdivide(
	const mgTLEdgePoint& from
	//indicates the edge and the point id where to subdivide,
)const;

//Print out this.
std::ostream& toString(std::ostream& ostrm)const;

///Let lpFrom=m_pline[from.eid], lpTo=m_pline[to.eid], then
///Subdivide this rectangle from the point from.pid of lpFrom
///to the point to.pid of lpTo.
///Then generate the two subdivided rectangles LPlines1 and LPlines2.
///When nlines =1, PidFrom(=from.m_pointID) or PidTo(=to.m_pointID) must be 0.
void subdivideFromTo(
	const mgTLEdgePoint& from,//id of edge and point that indicates which edge be subdivided.
	const mgTLEdgePoint& to,//id of edge and point that indicates which edge be subdivided.
	mgTL2LPlines& LPlines1,//1st subdivided rectangle, that is the after part
		//of lpFrom,  PidFrom to the end. The 1st edge is eidFrom(or a part of it).
	mgTL2LPlines& LPlines2,//2nd subdivided rectangle, that is the previous part of PidFrom.
		//The 1st edge is eidFrom(if Pidfrom>0), or subdividing bridge(if PidFrom==0).
	std::shared_ptr< mgTL2Polyline>& bridge //input or output bridge. When bridge is null,
		//subdivideFromTo makes it for output.
)const;
void subdivideFromTo(
	const mgTLEdgePoint& from,//id of edge and point that indicates which edge be subdivided.
	const mgTLEdgePoint& to,//id of edge and point that indicates which edge be subdivided.
	mgTL2LPlines& LPlines1,//1st subdivided rectangle, that is the after part
		//of lpFrom,  PidFrom to the end. The 1st edge is eidFrom(or a part of it).
	mgTL2LPlines& LPlines2//2nd subdivided rectangle, that is the previous part of PidFrom.
		//The 1st edge is eidFrom(if Pidfrom>0), or subdividing bridge(if PidFrom==0).
)const;

///Make strip data from pline0(=m_plines[eidMax]) and pline2(=m_plines[eidOpo]).
///The number of vertices of pline0 is greater or equal to the one of pline2,
///and the differecne must be at most 1.
///Neighbor edges of pline0 must be 2 points edges.
void makeStrip(
	int eidMax//edge of the maximum number of points.
)const;

///Make strip data from n*m points Suv and this boundary m_plines.
///The number of vertices of the edge id is greater or equal to the one of opposite,
///and the differecne must be at most 1. Furthermore, 
///Neighbor edges of pline0 must be the same number point edges.
void makeStrips(
	int edgeid,/// edge id of p
	const MGSPointSeq& Suv
)const;

///Make a fan of 1 triangle data from the pivot and the start point from.
///The number of points except pivot is numPoints.
void makeFan(
	const mgTLEdgePoint& pivot,//input pivot point
	const mgTLEdgePoint& from,//input the start point of the fan
	int numpoints//The number of points except pivot.
)const;

/// <summary>
/// Proprietry subprogram of tessellateShapr, makes a fan that includes the sharp vertex,
/// then tellellate the remainder.
/// The neighbor edge of the sharpID vertex must be concave.
/// </summary>
/// <param name="sharpID">the sharp vertex id is input</param>
void makeFanOfSharpTessellate(
	int sharpID,
	const mgTLEdgePoint& pivot
)const;

//Evaluate all the points of this polygon from the pivot(lineID,pivot)
//to build a fan, which is pushed back to m_triangles.
void makeFanAllPoints(
	const mgTLEdgePoint& pivot//input the pivot point.
)const;

///Make a fan of 1 triangle data. This point number is 3.
void makeFan3Points()const;

///Make fans from 4 point rectangle.
///This must have just 4 points.
///When nlines=4, each edge has (2, 2, 2, 2) points.
///When nlines=3, each edge has (2, 2, 3) points.
///When nlines=2, each edge has (3, 3), (2, 4) points.
void makeFan4Points(
	int eidMax // edge number whose point number is maximum).
)const;

//Make triangles. This must be 5 point polygon. That is
//When nlines=2, each edge has (3, 4), (2, 5) points.
//When nlines=3, each edge has (2, 3, 3), (2,2,4) points.
//When nlines=4, each edge has (2, 2, 2, 3) points.
void makeFan5Points(
	const int eids[5] //input the output of analyzePointNumber
				//(e.g. edge number whose point number is maximum).
)const;

///Make fan(s) from a rectangle that has 4 edges, whose number of vertices are (2,2,2,n).
///Here n>=3. The fan is pushed back to m_triangles.
void makeFan222n(
	int eidMax//edge id whose point number is more than 2.
)const;

///Do perform the tessellation for a rectangle whose points num is(2,m,n),
///that is a rectangel of 3 edges.
///Here m, n are greater than or equal to 3.
///Tessellated triangles are appended onto m_triangles.
void tessellate2mn(
	int eidTwo	//Edge id whose point number is 2.
)const;

///Tessellate a rectangle that has 2 continuous edges of only 2 points.
void tessellate22mn(
	int eidTwo
		//id of edge that has 2 vertices. Next edge of eidTwo also have only 2 vertices.
)const;

void tessellate2m2n(
	int eidMax	//Edge id whose points numuber is maximum.
)const;

///Do perform the tessellation for a concave rectangle.
void tessellateConcave(
	int concaveID	//vertiex id at which vertex concavity is detected.
)const;

///Do perform the tessellation for a sharp vertex rectangle.
void tessellateSharp(
	int sharpID	//vertiex id at which sharpness is detected.
)const;

///Do perform the tessellation for the mgTL2Plygon that has at most 4 edges outer loop,
///and that has no inner loops.
///The result will be appended onto triangles.
///When triangles.is_uv()=false, all of the element of the triangle position data
///has normal data as (x,y,z,xn,yn,zn). Here (x,y,z) is the position data and
///(xn,yn,zn) is the normal vector at the position (x,y,z).
///When triangles.is_uv()=true, all of the element of the triange position data are (u,v).
virtual void tessellate4()const;

/// <summary>
/// tessellation of a equal side 4-sided rectangle.
/// Let n=nPLine(idMax), nOpo=nPLine(idOpo), and m=nPLine(idPre), 
/// then nPline(idNxt)=m, and 0<=(n-nOpo)<=1.
/// </summary>
/// <param name="idMax">
/// edge id whose point number is equal to or greater thanthe oposite one,
///  and the difference must be less than or equal to 1.
/// </param>
void makeStripsE4(
	int idMax	/// edge id whose point number is equal to or greater than the oposite one,
				/// and the difference must be less than or equal to 1.
)const;

//Proprietry subprogram of makeStripsE4, handles when one of the edge's point number is 2.
void makeStripsE4Sub2(
	int id2	/// edge id whose point number is 2. The opposite edge's number must be3.
			//m or md must be >=3, and the difference of m and md must be at most 1.
)const;

//Proprietry subprogram of makeStripsE4, handles edges of any number of points.
void makeStripsE4Sub(
	int idMax	/// edge id whose point number is equal to or greater than the oposite one,
				/// and the difference must be less than or equal to 1.
)const;

/// <summary>
/// Tessellate if makeStripsE4 or tessellateEOpo is applicable.
/// Function's return value is true if tessellated using makeStripsE4 or tessellateEOpo,
/// false if not.
/// </summary>
/// <param name="id2">
/// edge id whose point number is 2.
/// </param>
bool t4makeStripsE4OrtessellateEOpo(
	const int eids[5],
	const int np[MAX_LINE_NUM]
)const;

/// <summary>
/// make strips from 4-sided rectangle if makeStripsE4 is applicable to this.
/// Function's retulrn value is true if makeStripsE4 was applied.
/// Let n=nPLine(id2), nOpo=nPLine(idOpo), m=nPLine(idPre), and md=nPLine(idNxt),
/// then n=2, nOpo>=2, and m, md>2.
/// </summary>
/// <param name="id2">
/// edge id whose point number is 2.
/// </param>
bool t4GeneralProcess()const;

/// <summary>
/// Tessellation of 4side rectangle that has 2 or 3 2 point edge.
/// </summary>
void t4Side22(
	int numEdgeTwoPoint,//number of edges whose points are only 2.
	int eidMin,//edge whose point number is minimum.
	int eidMax //edge whose point number is maximum.
)const;

/// <summary>
/// tessellation of 4-sided rectangle.
/// Let n=nPLine(idMin), nOpo=nPLine(idOpo), and m=nPLine(idPre), 
/// then n>=4, m>=3, nPline(idNxt)=m, and nOpo>=(n+2).
/// </summary>
void tessellateEOpo(
	int idMin	/// edge id whose point number is >=4.
				/// and nOpo(the one) of the idMax edge, nOpo>=(n+2);
)const;

/// <summary>
/// tessellation of 4-sided rectangle.
/// Let me=nPLine(id), m=nPLine(idOpo), nd=nPLine(idPre), n=nPLine(idNxt),
/// ndif=(nd-n), and mdif=(me-m),
/// then min(me,m,md,n)>=3, me>m, nd>n, and (mdif-ndif)>=2.
/// </summary>
void tessellate4General(int id)const;
//void tessellate4General2()const;

/// <summary>
/// Subdivide this at fromPoint by subdividing using from,
/// then tessellate both mgTL2LPlines.
/// </summary>
void tessellateBySubdivideAtPoint(
	const mgTLEdgePoint& from
)const;

/// <summary>
/// Subdivide this at mgTLEdgePoint(eid,pointId) by subdividing using from,
/// then tessellate both mgTL2LPlines.
/// </summary>
void tessellateBySubdivideAtPoint(
	int eid, int pointId
)const {
	tessellateBySubdivideAtPoint(mgTLEdgePoint(eid, pointId));
};

/// <summary>
/// Subdivide this from fromPoint to toPoint,
/// then tessellate both mgTL2LPlines.
/// </summary>
void tessellateBySubdivideFromTo(
	const mgTLEdgePoint& fromPoint,
	const mgTLEdgePoint& toPoint
)const;

private:

//Obtain minimum and maximum concavity vertex id.
void computeMinMaxConcavVid()const;

//get side length of m_plines[i] in sideLen[i].
//This is an approximation.
void getSideLength(double sideLen[4])const;

//This is a proprietry func of getToPointToSubdivide,
//gets to-point to subdivide by isectSlTl.
//Functin's return value is true if obtained.
bool getToPointByIsect(
	const mgTLEdgePoint& from,//id of edge and point that indicates which edge be subdivided.
	mgTLEdgePoint& to//to point is output.
)const;

//Change TL2Polyline's parameter t to mgTLEdgePoint.
//Midpoint value is changed to the best  nearest point id.
void changePolylineParameterToId(
	const mgTLEdgePoint& from,//id of edge and point that indicates which edge be subdivided.
	double t,//is a parameter value of  m_lines[to.edgeID] obtained by isectSlTl.
	mgTLEdgePoint& to//to.m_edgeid is input, and to.m_pointID will be updated
		//to indicate best point id from parameter t of m_lines[to.edgeID].
)const;

//Update (eidTo, PidTo) data to avoid
// (1) to generate more than 5 line mgTL2LPlines.
// (2) to point to be the same edge of from edge.
// (3) to be an end point.
void avoidAnomaly(
	const mgTLEdgePoint& from,//id of edge and point that indicates which edge be subdivided.
	mgTLEdgePoint& to//id of edge and point that indicates which edge be subdivided.
)const;

void avoidEndPoint(
	//const mgTLEdgePoint& from,//id of edge and point that indicates which edge be subdivided.
	mgTLEdgePoint& to//id of edge and point that indicates which edge be subdivided.
)const;

//Obtain how many lines(m_plines[]) are effective.
int getNumberOfLinesSub()const;

//Obtain how many points are in this.
int getNumberOfPointsSub()const;

///Check flatness of the edge eid and at the both ends of eid.
///When flat, make a fan putting the next or previous point as the pivot.
///Functions's return value is true if made.
bool makeFanIfFlatEdge(
	int eid//edge to test the flatness.
)const;

/// <summary>
/// make strips from 4-sided rectangle if makeStripsE4 is, tessellate2m2n, or
/// tessellate22mn,  applicable to this.
/// Function's retulrn value is true if tessellation was applied.
/// Let n=nPLine(id2), nOpo=nPLine(idOpo), m=nPLine(idPre), and md=nPLine(idNxt),
/// then n=2, nOpo>=2, and m, md>2.
/// </summary>
/// <param name="id2">
/// edge id whose point number is 2.
/// </param>
bool makeStrips2(
	int id2	// edge id whose point number is 2.
);

};

///Debug Function
inline
std::ostream& operator<< (std::ostream& ostrm, const mgTL2LPlines& lines){
	lines.toString(ostrm);
	return ostrm;
};

//Mix line1 and line2 by line1d*r1+line2d*(1-r1). Here,  line1d is transformed line1
//as line1's start point becomes P0, and the end point, P1.
//line2d is the same.
//Generated mgTL2Polyline's point number is nNew.If nNew is 0, line1's number is employed.
void getMixedLPline(
	const mgTL2LPline& line1, double r1,
	const mgTL2LPline& line2,
	const MGPosition& P0, const MGPosition& P1,
	mgTL2LPline& line,
	int nNew=0
);
void getMixedLPline(
	const mgTL2LPline& line1, double r1,
	const mgTL2LPline& line2,
	const MGPosition& Ps,//Start point.
	const mgTL2LPPoint& Pe,//End point.
	mgTL2LPline& line,
	int nNew = 0
);
void getMixedLPline(
	const mgTL2LPline& line1, double r1,
	const mgTL2LPline& line2,
	const mgTL2LPPoint& Ps,//Start point.
	const mgTL2LPPoint& Pe,//End point.
	mgTL2LPline& line,
	int nNew = 0
);
