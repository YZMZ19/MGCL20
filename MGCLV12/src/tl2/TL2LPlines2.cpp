#include "StdAfx.h"
#include "mg/Position.h"
#include "Tl2/TL2Triangles.h"
#include "Tl2/TL2Fans.h"
#include "Tl2/TL2LPlines.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//Compare this and input 3 mgTLConcavPid's, then if this is min or max, replace.
void mgTLConcavPid::getMinMax(
	mgTLConcavPid& maxConcav,//most concave point.
	mgTLConcavPid& max2Concav,//2ndly most concave point.
	mgTLConcavPid& maxConvex//most convex point.
){
	mgTLConcavPid& concavNew=*this;
	if(maxConcav<concavNew)
		maxConcav=concavNew;
	else if(max2Concav<concavNew)
		max2Concav=concavNew;
	if(concavNew<maxConvex)
		maxConvex=concavNew;
}

//Evaluate all the points(uv data) into P[].
//Function's return value is the number of points evaluated.
int mgTL2LPlines::evalUVVec(
	MGPosition uv[],//(u.v) data
	MGVector dir[]//uv[i+1]-uv[i] vectors are output in [i].
)const{
	int nlines=getNumberOfLines();
	int nPoints=0;
	for(int i=0; i<nlines; i++){
		const mgTL2LPline& lpi=m_plines[i];
		int nPim1=lpi.number_of_points()-1;
		for(int j=0; j<nPim1; j++)
			uv[nPoints++]=lpi.uv(j);
	}
	for(int i=0, ip1=1; ip1<nPoints;ip1++){
		dir[i]=uv[ip1]-uv[i];
		i=ip1;
	}
	int nm1=nPoints-1;
	dir[nm1]=uv[0]-uv[nm1];
	return nPoints;
}

//analyse concavity of this polygon.
//Function's return value is the number of concave points.
int mgTL2LPlines::analyzeConcavityAllPoints(
	mgTLConcavPid& pidMostConcave,//most concave concavity and the point.
	mgTLConcavPid& pid2ndConcave,//2ndly most concave concavity and the point
	mgTLConcavPid& pidMostConvex//most convex concavity and the point
)const{
	//Find a concave vertex.
	pidMostConcave.concavity=-3., pid2ndConcave.concavity=-3.;
	pidMostConvex.concavity=3.;
	int concaveCount=0;

	int nlines=getNumberOfLines(), j=0;//j is edge index.
	for(; j<nlines; j++){
		const mgTL2LPline& lpj=m_plines[j];
		mgTLConcavPid concav0(getVertexConcavity(j),j,0);
		concav0.getMinMax(pidMostConcave,pid2ndConcave,pidMostConvex);
		if(isConcave(concav0.concavity)){//IF loose concave.
			concaveCount++;
		}

		int nlpjm1=lpj.number_of_points()-1;	assert(nlpjm1<=4);
		if(nlpjm1<=1)
			continue;

		MGPosition uv[5];
		for(int i=0; i<=nlpjm1; i++)
			uv[i]=lpj.uv(i);

		MGVector Vpre=uv[1]-uv[0];
		for(int i=1; i<nlpjm1; i++){//i is point id in edge j.
			MGVector Vaft=uv[i+1]-uv[i];
			mgTLConcavPid concavityi(MGCL::concavity(Vpre,Vaft),j,i);
			concavityi.getMinMax(pidMostConcave,pid2ndConcave,pidMostConvex);
			if(isConcave(concavityi.concavity)){//IF loose concave.
				concaveCount++;
			}
			Vpre=Vaft;
		}
	}
	return concaveCount;
}

//Convert point id to epid(line id, id in the line).
//Here, point id is the point number from the start point of m_plines[0] to
//the point before the end point of the last of m_plines[]. 
void mgTL2LPlines::convertPointIDToLinePointID(int pid, mgTLEdgePoint& epid)const{
	int& edge=epid.m_edgeID;
	int& point=epid.m_pointID;

	int nPoints=0;//Total number of points summed so far.
	for(edge=0; edge<MAX_LINE_NUM; edge++){
		const mgTL2LPline& lpi=m_plines[edge];
		if(lpi.is_null())
			break;
		int nlpim1=lpi.number_of_points()-1;
		nPoints+=nlpim1;
		if(nPoints>pid){
			point=pid-(nPoints-nlpim1);
			break;
		}
	}
}

//Evaluate all the points of this polygon from the pivot(lineID,pivot)
//to build a fan, which is pushed back to m_triangles.
void mgTL2LPlines::makeFanAllPoints(
	const mgTLEdgePoint& pivot//input the pivot point.
)const{
	int nPoints=numberOfPoints();
	mgTLEdgePoint from=increment(pivot);
	makeFan(pivot,from,nPoints-1);
}

///Make a fan of 1 triangle data from the pivot and the start point from.
///The number of points except pivot is numPoints.
void mgTL2LPlines::makeFan(
	const mgTLEdgePoint& pivot,//input pivot point
	const mgTLEdgePoint& from,//input the start point of the fan
	int numpoints//The number of points except pivot.
)const{
	mgTL2Triangle* fanP=new mgTL2Triangle(numpoints+1, mgTESTRIANG::mgTESTRIANG_FAN);
	mgTL2Triangle& fan=*fanP;
	bool isUV=m_triangles.is_uv(), normalRequired=m_triangles.need_normal();
	auto PonEdge = [isUV, normalRequired](const mgTL2LPline& edge, int i) {
		return isUV ? edge.uv(i) : edge.xyz(i, normalRequired);
	};

	//1. set pivot.
	int i=0, lid=pivot.m_edgeID, j=pivot.m_pointID;
			//i is fan's point id to store,
			//j is point id in the edge lid.
	fan[i++] = PonEdge(m_plines[lid],j);

	//2. set fan data.
	int nlines=getNumberOfLines();
	lid=from.m_edgeID, j=from.m_pointID;
	while(i<=numpoints){
		const mgTL2LPline& lp=m_plines[lid++];
		int nPm1=lp.number_of_points()-1;
		for(; j<nPm1 && i<=numpoints; j++)
			fan[i++] = PonEdge(lp,j);
		lid=lid%nlines;//Next line
		j=0;//Start from the start.
	}
	m_triangles.push_back(fanP);
}

///Check flatness of the edge eid and at the both ends of eid.
//When flat, make a fan putting the next or previous point as the pivot.
///Functions's return value is true if made.
bool mgTL2LPlines::makeFanIfFlatEdge(
	int eid//edge to test the flatness.
)const{
	int nlines=getNumberOfLines();
	mgTLEdgePoint pivot(-1,0);
	if(isLooselyConcaveVertex(eid)){//If at previous is flat
		int eidPre=(eid+nlines-1)%nlines;
			pivot=decrement(mgTLEdgePoint(eidPre,0));
	}else{
		int eidAft=(eid+1)%nlines;
		if(isLooselyConcaveVertex(eidAft)){
			int nAft=nPLine(eidAft);
			pivot=increment(mgTLEdgePoint(eidAft,nAft-1));
		}
	}
	if(pivot.m_edgeID>=0){//When two edges neighbor to eid are flat.
		makeFanAllPoints(pivot);
		return true;
	}
	return false;
}

///Make a strip data from pline0(=m_plines[eidMax]) and pline2(=m_plines[eidOpo]).
///The number of vertices of pline0 is greater or equal to the one of pline2,
///and the differecne must be at most 1.
///Neighbor edges of pline0 must be 2 points edges.
/// pline2 may be one point edge.
void mgTL2LPlines::makeStrip(
	int eidMax//edge of the maximum number of points.
)const{
	int nlines = getNumberOfLines(); assert(nlines == 4);
	int eidOpo = (eidMax + 2) % nlines, eidNxt = (eidMax + 1) % nlines,
		eidPre = (eidMax + 3) % nlines;
	const mgTL2LPline& pline0 = m_plines[eidMax];
	const mgTL2LPline& pline2 = m_plines[eidOpo];

	int n0 = pline0.number_of_points(), n2 = pline2.number_of_points();
	assert(n0 >= n2 && (n0 - n2) <= 1);
	bool needToMakeFan = false;
	int eidFrom = eidMax;
	if (n2>=2) {//if pline2 is not one point edge.
		int eidOpenest = getMostOpenVid();
		if (n0 > n2) {
			if (eidOpenest == eidMax)//eidMax's start point's angle is max.
				eidFrom = eidOpo;
			else if (eidOpenest == eidNxt) {
				needToMakeFan = true;
			}
		}
		else {//case that n0==n2.
			if (eidOpenest == eidPre) {
				eidFrom = eidOpo;
				needToMakeFan = true;
			}
			else if (eidOpenest == eidNxt) {
				needToMakeFan = true;
			}
		}
	}

	bool isUV = m_triangles.is_uv(), normalRequired = m_triangles.need_normal();
	auto PonEdge = [isUV, normalRequired](const mgTL2LPline& edge, int i) {
		return isUV ? edge.uv(i) : edge.xyz(i, normalRequired);
	};

	const mgTL2LPline* plBase = &m_plines[eidFrom];
	const mgTL2LPline* plTo = &m_plines[(eidFrom + 2) % nlines];

	int nBase = plBase->number_of_points(), nTo = plTo->number_of_points();
	int nBm1 = nBase - 1, nTm1 = nTo - 1;
	mgTL2LPline lineTo;
	if (needToMakeFan) {
		mgTL2Triangle& fan = *(new mgTL2Triangle(3, mgTESTRIANG::mgTESTRIANG_FAN));
		fan[0] = PonEdge(*plTo, 1);//plTo->uv(1);
		fan[1] = PonEdge(*plBase, nBm1);//plBase->uv(nBm1);
		fan[2] = PonEdge(*plTo, 0);//plTo->uv(0);
		m_triangles.push_back(&fan);

		lineTo = mgTL2LPline(*plTo, 1, nTm1);
		nTo = nTm1;
		plTo = &lineTo;
	}

	mgTL2Triangle& strip = *(new mgTL2Triangle(nBase + nTo, mgTESTRIANG::mgTESTRIANG_STRIP));
	int j = 0;
	int nMin = nTo < nBase ? nTo : nBase;
	for (int i = 0; i < nMin; i++) {
		strip[j++] = PonEdge(*plBase, nBm1 - i);//plBase->uv(nBm1 - i);
		strip[j++] = PonEdge(*plTo, i);//plTo->uv(i);
	}
	for (int i = nMin; i < nBase; i++)
		strip[j++] = PonEdge(*plBase, nBm1 - i);//plBase->uv(nBm1 - i);
	for (int i = nMin; i < nTo; i++)
		strip[j++] = PonEdge(*plTo, i);//plTo->uv(i);
	m_triangles.push_back(&strip);
}

///Make a fan of 1 triangle data from. This point number is 3.
void mgTL2LPlines::makeFan3Points()const{
	mgTL2Triangle* fanP=new mgTL2Triangle(3, mgTESTRIANG::mgTESTRIANG_FAN);
	mgTL2Triangle& fan=*fanP;
	bool isUV=m_triangles.is_uv(), normalRequired =m_triangles.need_normal();
	auto PonEdge = [isUV, normalRequired](const mgTL2LPline& edge, int i) {
		return isUV ? edge.uv(i) : edge.xyz(i, normalRequired);
	};
	if(m_plines[2].is_null()){//When 2 edges.
		int eid=(nPLine(0)==3) ? 0:1;
		for(int i=0; i<3; i++)
			fan[i]= PonEdge(m_plines[eid],i);
	}else{
		for(int i=0; i<3; i++)
			fan[i]= PonEdge(m_plines[i],0);
	}
	m_triangles.push_back(fanP);
}

///Make fans from 4 point rectangle.
///This must have just 4 points.
///When nlines=4, each edge has (2, 2, 2, 2) points.
///When nlines=3, each edge has (2, 2, 3) points.
///When nlines=2, each edge has (3, 3), (2, 4) points.
void mgTL2LPlines::makeFan4Points(
	int eidMax // edge number whose point number is maximum).
)const{
	int nlines=getNumberOfLines();
	if(nlines>=3){
		mgTLEdgePoint epid(getMaximumConcavVid(), 0);
		if(nlines==3){
			double concvityEidMax=getEdgeConcavity(eidMax);
			if(isLooselyConcave(concvityEidMax)){
				m_triangles.makeFan(m_plines[eidMax],m_plines[(eidMax+2)%3]);
				return;
			}
			if(concvityEidMax>getVertexConcavity(epid.m_edgeID)){
				epid.m_edgeID=eidMax;
				epid.m_pointID=1;
			}
		}
		makeFanAllPoints(epid);
		return;
	}

	mgTLConcavPid Pid, Pid2, Pid3;
	int concaveCount=analyzeConcavityAllPoints(Pid,Pid2,Pid3);
	makeFanAllPoints(Pid.epID);
}

//Make triangles. This must be 5 point polygon. That is
//When nlines=2, each edge has (3, 4), (2, 5) points.
//When nlines=3, each edge has (2, 3, 3), (2,2,4) points.
//When nlines=4, each edge has (2, 2, 2, 3) points.
void mgTL2LPlines::makeFan5Points(
	const int eids[5] //input the output of analyzePointNumber
				//(e.g. edge number whose point number is maximum).
)const{
	int eidMax=eids[4];//Maximum vertex number.
	int nlines=getNumberOfLines();
	if(nlines==3){
		int numEdgeTwoPoint = eids[0];//number of edges whose points are only 2.
		if(numEdgeTwoPoint==2){//Case of (2,2,4)
			m_triangles.makeFan(m_plines[eidMax],m_plines[(eidMax+2)%3]);
			return;
		}else{//case of (2,3,3), (3,2,3), (3,3,2)
			const int& eidMin=eids[2];//minimum number of vertices and the edge id will be stored.
			if(makeFanIfFlatEdge(eidMin))//Check 1st 3 point flatness.
				return;

			if(makeFanIfFlatEdge((eidMin+1)%3))//Check 2nd 3 point flatness.
				return;
		}
	}else if(nlines==4 && makeFanIfFlatEdge(eidMax)){//Check flatness.
		return;
	}

	mgTLConcavPid Pid, Pid2, Pid3;
	int concaveCount=analyzeConcavityAllPoints(Pid,Pid2,Pid3);

	//assert(concaveCount<=2);
	if(concaveCount>=2){
		mgTLEdgePoint Pid4=increment(Pid.epID,2);
		mgTLEdgePoint Pid5=increment(Pid2.epID,2);
		if(Pid2.epID==Pid4 || Pid.epID==Pid5){
			tessellateBySubdivideFromTo(Pid.epID, Pid2.epID);
			return;
		}
	}
	makeFanAllPoints(Pid.epID);
}

/// <summary>
/// Proprietry subprogram of tessellateShapr, makes a fan that includes the sharp vertex,
/// then tellellate the remainder.
/// The neighbor edge of the sharpID vertex must be concave.
/// </summary>
/// <param name="sharpID">the sharp vertex id is input</param>
void mgTL2LPlines::makeFanOfSharpTessellate(
	int sharpID,
	const mgTLEdgePoint& pivot
)const{
	//Find the number of points of the fan.
	mgTLEdgePoint from, to;//the start and end point to subdivide this.
	int nlines = getNumberOfLines();
	int eidPre = (sharpID + nlines - 1) % nlines, eidAft = (sharpID + 1) % nlines;
	int pivotEdgeId = pivot.m_edgeID;
	int concaveEdgeId = (sharpID == pivotEdgeId) ? eidPre : sharpID;
	const mgTL2LPline& concaveEdge = m_plines[concaveEdgeId];
	const mgTL2LPline& pivotEdge = m_plines[pivotEdgeId];
	MGPosition pivotUv = pivotEdge.uv(pivot.m_pointID);
	int npConcave = concaveEdge.number_of_points(); assert(npConcave >= 3);

	mgTLEdgePoint fromFan;//the start point of the fan to make.
	int nFanP;//number of points of the fan.
	if (sharpID == pivotEdgeId) {//When pre edge is concave.
		MGVector pivo2nm2 = concaveEdge.uv(npConcave - 2) - pivotUv;
		double lenPivo2nm2 = pivo2nm2 % pivo2nm2;
		int i = npConcave - 3;
		while (i >= 0) {
			MGVector pivo2i = concaveEdge.uv(i--) - pivotUv;
			double lenPivo2i = pivo2i % pivo2i;
			if (lenPivo2i >= lenPivo2nm2)
				break;
		}
		i++;
		fromFan = from = mgTLEdgePoint(concaveEdgeId, i), to=pivot;
		nFanP = npConcave - i;
	}
	else {//When sharpID edge is concave.
		MGVector pivo21 = concaveEdge.uv(1) - pivotUv;
		double lenPivo21 = pivo21 % pivo21;
		nFanP = 2;
		int nm1 = npConcave - 1;
		while (nFanP <=nm1) {
			MGVector pivo2i = concaveEdge.uv(nFanP++) - pivotUv;
			double lenPivo2i = pivo2i % pivo2i;
			if (lenPivo2i >= lenPivo21)
				break;
		}
		from =pivot,  to= mgTLEdgePoint(sharpID, nFanP-1);
		fromFan = mgTLEdgePoint(sharpID, 0);
	}
	makeFan(pivot, fromFan, nFanP);

	std::shared_ptr<mgTL2Polyline> bridge=
		std::make_shared< mgTL2Polyline>(*this, to, from);
	mgTL2LPlines dummy(m_stack, m_triangles);
	mgTL2LPlines* LPlines = new mgTL2LPlines(m_stack, m_triangles);
	subdivideFromTo(from, to, dummy, *LPlines, bridge);
	m_stack.emplace(LPlines);
}
