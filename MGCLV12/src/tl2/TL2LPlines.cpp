#include "StdAfx.h"
#include "mg/Tolerance.h"
#include "mg/Straight.h"
#include "topo/Loop.h"
#include "Tl2/TL2parameter.h"
#include "Tl2/TL2Triangles.h"
#include "Tl2/TL2LPlines.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/****************************************************************/
/*   Copyright (c) 2021 by System fugen G.K.                    */
/*                       All rights reserved.                   */
/****************************************************************/

mgTL2LPlines::mgTL2LPlines(
	std::stack<std::unique_ptr<mgTL2LPlines>>& stack,
	mgTL2Triangles& triangles	//mgTL2Triangles for the tessellated triangles to output.
):m_stack(stack),m_triangles(triangles),m_nlines(-1),m_npoints(-1),
m_minConcavVid(-1), m_maxConcavVid(-1){
	for(int i=0; i<MAX_LINE_NUM; i++)
		m_concavity[i]=NotObtainedConcav;//the value that indicates not obtained yet.
}

///Construct from outer loop that has less than or equal to 5 edges.
///The outer boundary loop of the face must have edges whose number is up to 5.
mgTL2LPlines::mgTL2LPlines(
	std::stack<std::unique_ptr<mgTL2LPlines>>& stack,
	const MGLoop& oloop,
	mgTL2Triangles& triangles	//mgTL2Triangles for the tessellated triangles to output.
):m_stack(stack), m_triangles(triangles),m_nlines(-1), m_npoints(-1),
m_minConcavVid(-1), m_maxConcavVid(-1){
	MGComplex::const_iterator ei=oloop.pcell_begin();
	int nedges=oloop.number_of_edges();	assert(nedges<= MAX_LINE_NUM);

	for(int i=0; i<nedges; ei++, i++){
		const mgTL2Polyline* pli=TL2Polyline(ei);
		int nVi=pli->number_of_points();
		m_plines[i]=mgTL2LPline(pli,0,nVi);
		m_concavity[i]=NotObtainedConcav;//the value that indicates not obtained yet.
	}
}

/// <summary>
/// Subdivide this at mgTLEdgePoint(eid,pointId) by subdividing using from,
/// then tessellate both mgTL2LPlines.
/// </summary>
void mgTL2LPlines::tessellateBySubdivideAtPoint(
	const mgTLEdgePoint& from
)const {
	mgTLEdgePoint to = getToPointToSubdivide(from);
	tessellateBySubdivideFromTo(from, to);
}

/// <summary>
/// get 1st concave edge out of m_plines.
/// Function's return value is:
///  -1 if not found.
/// edge id if found.
///Concavity is tested for all edge's start and end point tangent difference.
///(Not concavity at a vertex)
/// </summary>
int mgTL2LPlines::get1stConcaveEdge()const {
	int nedges = getNumberOfLines();
	int eid = 0;
	for (; eid < nedges; eid++) {
		if (nPLine(eid) == 2)
			continue;
		if (isConcaveEdge(eid))
			break;
	}
	if (eid < nedges) //When a concave edge found.
		return eid;
	return -1;
}

/// <summary>
/// Subdivide this from fromPoint to toPoint,
/// then tessellate both mgTL2LPlines.
/// </summary>
void mgTL2LPlines::tessellateBySubdivideFromTo(
	const mgTLEdgePoint& fromPoint,
	const mgTLEdgePoint& toPoint
)const {
	mgTL2LPlines* LPlines1 = new mgTL2LPlines(m_stack, m_triangles);
	mgTL2LPlines* LPlines2 = new mgTL2LPlines(m_stack, m_triangles);
	subdivideFromTo(fromPoint, toPoint, *LPlines1, *LPlines2);
	m_stack.emplace(LPlines1);
	m_stack.emplace(LPlines2);
}

#define MIN_EDGDE 4
//Find an edge whose point num is too large compared to the opposite one.
//Function's return value is:
//-1 if not found,
//edge number whose point num is too large(0 - 3).
int edgeTooLargeDiff(int np[4]) {
	for (int i = 0; i < 4; i++) {
		int nm1 = np[i] - 1;
		int nOpom1 = np[(i + 2) % 4] - 1;
		if (nm1 == 1 || nOpom1 == 1)
			continue;
		if (nm1 > nOpom1 * LARGE_RATIO) {
			int npre = np[(i + 3) % 4], nnext = np[(i + 1) % 4];
			if (npre >= MIN_EDGDE && nnext >= MIN_EDGDE)
				return i;
		}
	}

	return -1;
}

///Do perform the tessellation for the mgTL2Plygon that has at most 4 edges outer loop,
///and that has no inner loops.
///The result will be appended onto triangles.
///When triangles.is_uv()=false, all of the element of the triangle position data
///has normal data as (x,y,z,xn,yn,zn). Here (x,y,z) is the position data and
///(xn,yn,zn) is the normal vector at the position (x,y,z).
///When triangles.is_uv()=true, all of the element of the triange position data are (u,v).
void mgTL2LPlines::tessellate4()const {
	int nlines = getNumberOfLines();
	if (nlines == 0)
		return;

	if (nlines == 1) {
		tessellateBySubdivideAtPoint(0, 0);
		return;
	}

	int eids[5], np[MAX_LINE_NUM];//np[i] is m_lines[i]'s number of points.
	analyzePointNumber(eids, np);
	int numEdgeTwoPoint = eids[0];//number of edges whose points are only 2.
	int eidMin = eids[2];//minimum number of vertices and the edge id will be stored.
	int eidMax = eids[4];//Maximum vertex number and the id of m_plines[i] will be stored.

	if (nlines == 2 && numEdgeTwoPoint) {
		m_triangles.makeFan(m_plines[eidMax], m_plines[eidMin]);
		return;
	}//this assures nMin>=3.

	int nPoints = numberOfPoints();//Total number of points of this rectangle.
	switch (nPoints) {
		case 3: makeFan3Points(); return;
		case 4:	makeFan4Points(eidMax); return;
		case 5:	makeFan5Points(eids); return;
		default:;
	}

	if (nlines == 3 && numEdgeTwoPoint) {
		if (numEdgeTwoPoint == 2)
			m_triangles.makeFan(m_plines[eidMax], m_plines[(eidMax + 2) % 3]);
		else
			tessellate2mn(eidMin);	
		return;
	}

	int concaveVID, convexVID;
	SHARPorCONCAVE concavity = analyzeConcavity(concaveVID, convexVID);
	if (concavity == CONCAVE) {
		if(nlines!=2)////If nlines==2, it is a illegal data, neglected.
			tessellateConcave(concaveVID);
		return;
	}
	if (nlines == 4 && numEdgeTwoPoint>=2) {
		t4Side22(numEdgeTwoPoint, eidMin, eidMax);
		return;
	}
	if (nlines<=4 && concavity == SHARP) {
		tessellateSharp(convexVID);
		return;
	}

	//Here, numEdgeTwoPoint<=1(at most 1).
	if (nlines == 4) {
		if (t4makeStripsE4OrtessellateEOpo(eids, np))
			return;

		int eDiff = edgeTooLargeDiff(np);
		if (eDiff >= 0) {//When too large difference edges found.
			int eNext = (eDiff + 1) % 4, ePre = (eDiff + 3) % 4;
			int eToDivide = np[eNext] >= np[ePre] ? eNext : ePre;
			tessellateBySubdivideAtPoint(eToDivide, np[eToDivide]/2);
			return;
		}

		//Check concavity of all edges.
		int eidConcave = get1stConcaveEdge();
		if (eidConcave >= 0) {//When a concave edge found.
			int eOpo = (eidConcave + 2) % nlines;
			if (np[eOpo] >= 3) {
				tessellateBySubdivideAtPoint(eOpo, np[eOpo]/2);
				return;
			}
		}

		//Here no opposite pairs are equal.
		int concavId = getMaximumConcavVid();
		if (!isLooselyConcave(getVertexConcavity(concavId))) {
			int eid= (eidMin + 2) % 4;//eid when numEdgeTwoPoint>0. Note that nPLine(eid)>=3.
			if (numEdgeTwoPoint==0) {
				int eidPre = (eidMax + 3) % 4, eidNxt = (eidMax + 1) % 4;
				int dif1 = np[eidMax] - np[(eidMax + 2) % 4], dif2 = np[eidPre] - np[eidNxt];
				int absdif2 = (dif2 < 0) ? -dif2 : dif2; assert(dif1 >= 1 && absdif2 >= 1);
				eid = eidNxt;
				if (dif1 > absdif2)
					eid = eidMax;
				else if (dif2 > 0)
					eid = eidPre;
			}
			if (np[(eid + 2) % 4] == 2 || !isConcaveEdge(eid)) {
				tessellate4General(eid);
				return;
			}
		}
	}

	tessellateBySubdivideAtPoint(getFromPointToSubdivide(eids));
}

///Do perform the tessellation for a concave rectangle.
void mgTL2LPlines::tessellateConcave(
	int concaveID	//vertiex id at which vertex concavity is detected.
)const{
	mgTLEdgePoint from(concaveID,0);
	int nlines=getNumberOfLines();
	if(numberOfPoints()==6){
		mgTLEdgePoint to=increment(from,3);
		int& PidTo=to.m_pointID;
		int& eidTo=to.m_edgeID;
		if(PidTo==0){
			if(eidTo==(concaveID+1)%nlines)
				PidTo++;
			else{
				int eidPre=(concaveID+nlines-1)%nlines;
				if(eidTo==eidPre){
					eidTo=(eidPre+nlines-1)%nlines;;
					PidTo=nPLine(eidTo)-2;
				}
			}
		}
		tessellateBySubdivideFromTo(from, to);
	}
	else
		tessellateBySubdivideAtPoint(from);
}

///Do perform the tessellation for a sharp vertex rectangle.
void mgTL2LPlines::tessellateSharp(
	int sharpID	//vertiex id at which sharpness is detected.
)const{
	int nlines=getNumberOfLines();
	int eidPre=(sharpID+nlines-1)%nlines, eidAft=(sharpID+1)%nlines;
	int n = nPLine(sharpID), nPre=nPLine(eidPre);

	mgTLEdgePoint from, to;//will be tested if undefined.
	if (nlines <= 3 && nPre >= 3 && n >= 3) {
		from = nPre > n ? mgTLEdgePoint(sharpID, 1) : mgTLEdgePoint(eidPre, nPre - 2);
	}else if(nlines>=4){
		if(isConcaveEdge(sharpID)){
			from = mgTLEdgePoint(eidPre, nPre - 2);
		}else if(isConcaveEdge(eidPre)){
			from = mgTLEdgePoint(sharpID, 1);
		}
	}

	if (from.defined()) {
		makeFanOfSharpTessellate(sharpID, from);
		return;
	}

	if(nlines>=4){
		if(isLooselyConcaveVertex(eidAft))
			from=mgTLEdgePoint(eidAft,0);
		else if(isLooselyConcaveVertex(eidPre))
			from=mgTLEdgePoint(eidPre,0);
	}
	if(from.undefined() && (nlines<=4 || nPre==2 || n==2)){
		//|| nPre==2 || n==2 is the condition that does
		//not increment the edge number of this when nlines=5.
		from= n==2 ? mgTLEdgePoint(eidAft,0) : mgTLEdgePoint(sharpID, 1);
		to=mgTLEdgePoint(eidPre,nPre-2);
	}
	if(from.undefined())
		from = isMoreOpen(eidAft, eidPre) ? mgTLEdgePoint(eidAft,0) : mgTLEdgePoint(eidPre,0);

	if(to.undefined())
		to= getToPointToSubdivide(from);
	tessellateBySubdivideFromTo(from, to);
}

///Tessellate a rectangle that has 2 continuous edges of only 2 points.
void mgTL2LPlines::tessellate22mn(
	int eidTwo
		//id of edge that has 2 vertices. Next edge of eidTwo also have only 2 vertices.
)const{
	int eidAft=(eidTwo+1)%4;
	int eidPre=(eidTwo+3)%4;
	if(isLooselyConcaveVertex(eidPre)){
		makeFanAllPoints(mgTLEdgePoint(eidAft,0));
	}else{
		tessellateBySubdivideAtPoint(getMaximumConcavVid(), 0);
	}
}

///Do perform the tessellation for a rectangle whose points num is(2,m,2,n).
///Here m, n are greater or equal to 3.
///Tessellated triangles are appended onto m_triangles
void mgTL2LPlines::tessellate2m2n(
	int eidMax//edge id whose point number is maximum.
)const {
	const mgTL2LPline& eMax = m_plines[eidMax];
	int n = eMax.number_of_points();

	int eidOpo = (eidMax + 2) % 4;;
	const mgTL2LPline& eOpo = m_plines[eidOpo];
	int nOpo = eOpo.number_of_points(); assert(nOpo <= n);
	if (n - nOpo <= 1) {
		makeStrip(eidMax);
		return;
	}
	if (nOpo == 2) {
		makeFan222n(eidMax);
		return;
	}

	mgTLEdgePoint from(eidOpo, nOpo / 2);
	mgTLEdgePoint to(eidMax, (n - 1) / 2);
	tessellateBySubdivideFromTo(from, to);
}

///Make a fan from a rectangle that has 4 edges, whose number of vertices are (2,2,2,n).
///Here n>=3. The fan is pushed back to m_triangles.
void mgTL2LPlines::makeFan222n(
	int eidMax//edge id whose point number is more than 2.
)const{
	if(!makeFanIfFlatEdge(eidMax)){
		const mgTL2LPline& eMax=m_plines[eidMax];
		int n= eMax.number_of_points();
		int eidOpo =(eidMax+2)%4;
		const mgTL2LPline& eOpo =m_plines[eidOpo];
		MGPosition Ps= eOpo.uv(0), Pe= eOpo.uv(1);

		MGPosition uviPre= eMax.uv(0);
		int i=1, nm1=n-1;
		for(; i<nm1; i++){
			const MGPosition& uvi= eMax.uv(i);
			MGVector Vi=uvi-uviPre;
			MGVector Vsi=Ps-uvi;
			MGVector Vei=Pe-uvi;
			double sangS=Vi.sangleSigned2D(Vsi);
			double sangE=Vi.sangleSigned2D(Vei);
			if(sangS>=sangE)//When pe is more away from uvi.
				break;
			uviPre=uvi;
		}
		int nHalf1=i;
		int nHalf2=n-nHalf1+1;
		mgTL2LPline pl1(eMax,0,nHalf1);
		mgTL2LPline pl2(eMax,nHalf1-1,nHalf2);
		if(nHalf1>=2)
			m_triangles.makeFan(pl1, eOpo,false,true);
		m_triangles.makeFan(pl2, eOpo,true);
	}
}

///Do perform the tessellation for a rectangle whose points num is(2,m,n),
///that is a rectangel of 3 edges.
///Here m, n are greater than or equal to 3.
///Tessellated triangles are appended onto m_triangles
void mgTL2LPlines::tessellate2mn(
	int eidTwo	//Edge id whose point number is 2.
)const{
	int eidAft=(eidTwo+1)%3, eidPre=(eidTwo+2)%3;
	mgTLEdgePoint from(eidAft, nPLine(eidAft) -2), to(eidPre,1);
	tessellateBySubdivideFromTo(from, to);
}


//get side length of m_plines[i] in sideLen[i].
//This is an approximation.
void mgTL2LPlines::getSideLength(double sideLen[4])const {
	assert(getNumberOfLines() == 4);

	const MGSurface& srf = *(m_triangles.surface());
	MGPosition corners[4];
	for (int i = 0; i < 4; i++) 
		corners[i] = srf.eval(m_plines[i].uv(0));

	for (int i = 0; i < 4; i++) {
		MGPosition midPoint = srf.eval(m_plines[i].eval(.5));

		sideLen[i] = (midPoint - corners[i]).len();
		sideLen[i] += (corners[(i + 1) % 4]- midPoint).len();
	}
}
