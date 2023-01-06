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
/*   Copyright (c) 2018 by System fugen G.K.                    */
/*                       All rights reserved.                   */
/****************************************************************/

//Nullify this.
void mgTL2LPlines::setNull(){
	for(int i=0; i<MAX_LINE_NUM; i++)
		m_plines[i].setNull();
	m_nlines=0;
}

//Test if the vertex vid is loosely flat(true), or not.
bool mgTL2LPlines::isLooselyConcaveVertex(int vid)const{
	return isLooselyConcave(getVertexConcavity(vid));
}

//Test if the edge eid is concave or not.
bool mgTL2LPlines::isConcaveEdge(int eid)const{
	return m_plines[eid].getConcavity() >=STRICT_CONCAVITY;
}

//Compare vid1's and vid2's openess(concavity),
//then if vid1's is more open, return true.
bool mgTL2LPlines::isMoreOpen(int vid1, int vid2)const{
	double concav1=getVertexConcavity(vid1), concave2=getVertexConcavity(vid2);
	return concav1>concave2;
}

///Get concavity at the vertex vid, which is not edge's concavity.
///That is, the concavity is obtained from the difference of two vectors
///at the edge id vid's start point and at the previous edge's end point.
//Concavity's value is from -2 to 2. 2 is most concave, and
// -2 means 180 degree convex(most convex), -1:90 degree convex, 0:flat
// 1:90 degree concave, 2:180 degree concave.
//Concavity's value is from -2 to 2. 2 is most concave, and
// -2 means 180 degree convex(most convex), -1:90 degree convex, 0:flat
// 1:90 degree concave, 2:180 degree concave.
double mgTL2LPlines::getVertexConcavity(int vid)const{
	double& concav=m_concavity[vid];
	if(isUndefinedConcav(vid)){//If not obtained yet.
		int nlines=getNumberOfLines();
		int preID=(vid+nlines-1)%nlines;
		const mgTL2LPline& preL=m_plines[preID];
		int nPre=preL.number_of_points();
		const mgTL2LPline& currentL=m_plines[vid];

		MGPosition Pvid=currentL.uv(0);
		MGVector Tpre=Pvid-preL.uv(nPre-2);
		MGVector Tnow=currentL.uv(1)-Pvid;
		concav=MGCL::concavity(Tpre,Tnow);
	}
	return concav;
}

//Obtain minimum and maximum concavity vertex id.
void mgTL2LPlines::computeMinMaxConcavVid()const{
	int nlines=getNumberOfLines();
	double minConcav=NotObtainedConcav, maxConcav=-NotObtainedConcav;
	for(int i=0; i<nlines; i++){
		double concav=getVertexConcavity(i);
		if(concav<minConcav){
			minConcav=concav;
			m_minConcavVid =i;
		}
		if(concav>maxConcav){
			maxConcav=concav;
			m_maxConcavVid =i;
		}
	}
}

//Get the vertx id of minimum concavity.
int mgTL2LPlines::getMinimumConcavVid()const{
	if(m_minConcavVid ==-1){
		computeMinMaxConcavVid();
	}
	return m_minConcavVid;
}

//Get the vertx id of maximum concavity.
int mgTL2LPlines::getMaximumConcavVid()const{
	if(m_maxConcavVid ==-1){
		computeMinMaxConcavVid();
	}
	return m_maxConcavVid;
}

///Check the concavities of all the vertices, concavity or sharp angle.
///Function's return value is SHARPorCONCAVE:
//       CONCAVE if concavity found.
///      SHARP if sharp angle found.
///      OTHER if both concavity ro sharp angle not found.
///concaveVID and convexVid are always returned regardless of functon's return value.
mgTL2LPlines::SHARPorCONCAVE mgTL2LPlines::analyzeConcavity(
	int& concaveVID,//most concave vertex id returned,
	int& convexVid	//most convex vertex id returned.
)const{
	SHARPorCONCAVE retval=OTHER;
	concaveVID= getMaximumConcavVid();
	convexVid= getMinimumConcavVid();
	if(isLooselyConcave(m_concavity[concaveVID])){
			retval=CONCAVE;
	}else if(isSharp(m_concavity[convexVid])){
		retval=SHARP;
	}
	return retval;
}

//Obtain from point to subdivide.
//When nlines is 2 or 3, there must be no 2 point edges.
//When nlinesis 4, 2point edge number is at most 1.
mgTLEdgePoint mgTL2LPlines::getFromPointToSubdivide(
	const int eids[5]
		//Following data are input(that are obtained by analyzePointNumber):
		//eids[0]:num of 2 point edges,
		//eids[1], [2];//minimum number of vertices and the edge id
		//eids[3], [4];//Maximum vertex number and the edge id
)const{
	mgTLEdgePoint from; //(eid=edge id, vid=point id) is returned, 
	int& eidFrom=from.m_edgeID;
	int nlines=getNumberOfLines();
	if(nlines==5){
		eidFrom=getMostOpenVid();
	}else{
		const int& eidMin=eids[2];
		const int& eidMax=eids[4];
		eidFrom=eidMax;
		if(nlines==4){
			assert(eids[0]<=1);//2point edge number must be at most 1.
			int eidOpo=(eidMax+2)%4;
			if(eids[0]==1){//When only one of the edges has only 2 vertices.
				if(eidMin==eidOpo){//When the oposite edge has only 2points.
					//We avoid 2points edge's opposite one to divide.
					eidFrom=(eidMax+1)%4, eidOpo=(eidMax+3)%4;
				}
			}else{
				MGPosition uv[4]={
					m_plines[0].uv(0), m_plines[1].uv(0),
					m_plines[2].uv(0), m_plines[3].uv(0) };//vertex i's (u,v) data.

				double spanLen[4];//each span length of edge i.
				for(int i=0; i<=3; i++){
					int ip1 = (i + 1) % 4;
					MGVector Vi=uv[ip1]-uv[i];
					int nspani=m_plines[i].number_of_points()-1;
					spanLen[i]=Vi%Vi/double(nspani*nspani);//square of the length.
				}

				double maxRatio2= LARGE_RATIO* LARGE_RATIO;
				int iTooLarge=0;
				for(; iTooLarge<4; iTooLarge++){
					int iOpo=(iTooLarge+2)%4;
					if(spanLen[iTooLarge]>spanLen[iOpo]* maxRatio2)
						break;
				}
				if(iTooLarge<4){
					eidFrom=(iTooLarge+1)%4, eidOpo=(iTooLarge+3)%4;
				}
			}
			if(nPLine(eidFrom)>nPLine(eidOpo))
				eidFrom=eidOpo;
		}
		from.m_pointID=nPLine(eidFrom)/2;
	}
	return from;
}

//Normalize the point id epid(eid, Pid), i.e. 
//if is the end point id, change it to the start point id of the next side.
void mgTL2LPlines::normalizePointID(
	mgTLEdgePoint& epid
)const{
	int& eid=epid.m_edgeID;
	int& Pid=epid.m_pointID;
	if(m_plines[eid].isEndPoint(Pid)){
		int nlines=getNumberOfLines();
		//normalize.
		eid=(eid+1)%nlines;
		Pid=0;
	}
}

//Return true if epid is the next point of epid2 around this closed polygon.
bool mgTL2LPlines::isNextPoint(
	const mgTLEdgePoint& epid, const mgTLEdgePoint& epid2
)const{
	mgTLEdgePoint epidTemp=increment(epid2);
	return (epid==epidTemp);
}

//Return true if epid is the previous point of epid2 around this closed polygon.
bool mgTL2LPlines::isPrePoint(
	const mgTLEdgePoint& epid, const mgTLEdgePoint& epid2
)const{
	mgTLEdgePoint epidTemp=increment(epid);
	return (epidTemp==epid2);
}

// Get the next(incremented) point of epid.
// Function's return value is the incremented point around this mgTL2LPlines polygons.
mgTLEdgePoint mgTL2LPlines::increment(
	const mgTLEdgePoint& epid,//input the target point.
	int num
)const{
	const int& Pid=epid.m_pointID;
	const int& eid=epid.m_edgeID;

	mgTLEdgePoint epid2;
	int& eid2=epid2.m_edgeID;
	int& Pid2=epid2.m_pointID;
	int nlines=getNumberOfLines();
	Pid2=Pid+num;
	eid2=eid;
	int npim1=nPLine(eid)-1;
	while(Pid2>=npim1){
		Pid2-=npim1;
		eid2=(eid2+1)%nlines;
		npim1=nPLine(eid2)-1;
	}
	return epid2;
}

// Get the previous(decremented) point of epid.
// Function's return value is the decremented point around this mgTL2LPlines polygons.
mgTLEdgePoint mgTL2LPlines::decrement(
	const mgTLEdgePoint& epid,//input mgTLEdgePoint to decrement.
	int num
)const{
	mgTLEdgePoint epid2;//decremented mgTLEdgePoint.
	int& eid2=epid2.m_edgeID;
	int& Pid2=epid2.m_pointID;

	const int& Pid=epid.m_pointID;
	const int& eid=epid.m_edgeID;
	int nlines=getNumberOfLines();
	Pid2=Pid-num;
	eid2=eid;
	while(Pid2<0){
		eid2=(eid2+nlines-1)%nlines;//previous edge.
		int npim1=nPLine(eid2)-1;
		Pid2+=npim1;
	}
	return epid2;
}

//Obtain how many lines(m_plines[]) are effective.
int mgTL2LPlines::getNumberOfLinesSub()const{
	for(m_nlines=0; m_nlines<MAX_LINE_NUM; m_nlines++){
		if(m_plines[m_nlines].is_null())
			break;
	}
	return m_nlines;
}

//Obtain how many points are in this.
int mgTL2LPlines::getNumberOfPointsSub()const{
	m_npoints=0;
	for(int i=0; i<MAX_LINE_NUM; i++){
		const mgTL2LPline& pli=m_plines[i];
		if(pli.is_null())
			break;
		m_npoints+=short(pli.number_of_points()-1);
	}
	return m_npoints;
}

///get the (u,v) parameter box.
MGBox mgTL2LPlines::getUvBox()const{
	MGBox uvbox;
	int n=getNumberOfLines();
	for(int i=0; i<n; i++){
		const mgTL2LPline& pli=m_plines[i];
		if(pli.is_null())
			break;
		uvbox|=pli.getUvBox();
	}
	return uvbox;
}

///Analyze the numbers of this, getting the number of 2 point edge, maximum and
///minimum point number edges.
void mgTL2LPlines::analyzePointNumber(
	int eids[5],
		//Following data are output:
		//eids[0]:num of 2 point edges,
		//eids[1], [2];//minimum number of vertices and the edge id
		//eids[3], [4];//Maximum vertex number and the edge id
	int np[5]//number of points of edge i is returned in np[i].
)const{
	int& num2=eids[0]=0;
	int &nMin=eids[1]=-1, &eidMin=eids[2]=0;//minimum number of vertices and the edge id will be stored.
	int &nMax=eids[3]=0, &eidMax=eids[4]=0;//Maximum vertex number and the id of m_plines[i] will be stored.
	m_npoints=0;
	int nLines=getNumberOfLines();
	for(int i=0; i<nLines; i++){
		int nVi=np[i]=nPLine(i);
		m_npoints+=short(nVi)-1;

		if(nVi==2)
			num2++;
		if(nMin<0 || nVi<nMin){
			nMin=nVi;
			eidMin=i;
		}
		if(nVi>nMax){
			nMax=nVi;
			eidMax=i;
		}
	}
	int nlm1=nLines-1;
	int nLast= np[nlm1];
	if(eidMin==0 && nLast==nMin)
		eidMin=nlm1;
	if(eidMax==0 && nLast==nMax)
		eidMax=nlm1;
}
void mgTL2LPlines::subdivideFromTo(
	const mgTLEdgePoint& from,//id of edge and point that indicates which edge be subdivided.
	const mgTLEdgePoint& to,//id of edge and point that indicates which edge be subdivided.
	mgTL2LPlines& LPlines1,//1st subdivided rectangle, that is the after part
		//of lpFrom,  PidFrom to the end. The 1st edge is eidFrom(or a part of it).
	mgTL2LPlines& LPlines2//2nd subdivided rectangle, that is the previous part of PidFrom.
		//The 1st edge is eidFrom(if Pidfrom>0), or subdividing bridge(if PidFrom==0).
)const {
	std::shared_ptr<mgTL2Polyline> bridge;
	subdivideFromTo(from, to, LPlines1, LPlines2, bridge);
}

///Let lpFrom=m_pline[from.eid], lpTo=m_pline[to.eid], then
///Subdivide this rectangle from the point from.pid of lpFrom
///to the point to.pid of lpTo.
///Then generate the two subdivided rectangles LPlines1 and LPlines2.
///When nlines =1, PidFrom(=from.m_pointID) or PidTo(=to.m_pointID) must be 0.
void mgTL2LPlines::subdivideFromTo(
	const mgTLEdgePoint& from,//id of edge and point that indicates which edge be subdivided.
	const mgTLEdgePoint& to,//id of edge and point that indicates which edge be subdivided.
	mgTL2LPlines& LPlines1,//1st subdivided rectangle, that is the after part
		//of lpFrom,  PidFrom to the end. The 1st edge is eidFrom(or a part of it).
	mgTL2LPlines& LPlines2,//2nd subdivided rectangle, that is the previous part of PidFrom.
		//The 1st edge is eidFrom(if Pidfrom>0), or subdividing bridge(if PidFrom==0).
	std::shared_ptr<mgTL2Polyline>& bridge //input or output bridge.
		//When bridge is null, subdivideFromTo makes it for output.
)const{
	const int& eidFrom=from.m_edgeID;
	const int& PidFrom=from.m_pointID;
	int eidTo=to.m_edgeID;
	int PidTo=to.m_pointID;
	assert(!m_plines[eidFrom].isEndPoint(PidFrom));

	int nlines=getNumberOfLines();
	if(m_plines[eidTo].isEndPoint(PidTo)){
		eidTo =(eidTo +1)%nlines;
		PidTo=0;
	}
	const mgTL2LPline& lpFrom=m_plines[eidFrom];
	const mgTL2LPline& lpTo=m_plines[eidTo];
	
	//Build bridge line.
	if(!bridge)
		bridge = std::shared_ptr<mgTL2Polyline>(new mgTL2Polyline);
	if( bridge->is_null())
		lpFrom.polygonizeSL(lpTo,PidFrom,PidTo, *bridge);

	mgTL2LPline saveToLp1;
	int id1=0, id2=0;//Used to count the LPlines1,2 area id.

	//build LPlines1,2
	if(PidFrom){
		lpFrom.subdivide(PidFrom,LPlines2[id2++],LPlines1[id1++]);
		LPlines2.m_concavity[0]=m_concavity[eidFrom];
	}else if(eidFrom!= eidTo){
		LPlines1[id1++]=lpFrom;
	}

	LPlines2[id2++]=mgTL2LPline(bridge);
	if(PidTo){
		lpTo.subdivide(PidTo,saveToLp1,LPlines2[id2++]);
	}else
		LPlines2[id2++]=lpTo;

	//Complete LPlines2.
	for(int idE=(eidTo +1)%nlines;idE!=eidFrom; idE=(idE+1)%nlines){
		LPlines2[id2]=m_plines[idE];
		LPlines2.m_concavity[id2++]=m_concavity[idE];
	}

	//build LPlines1.
	if(eidFrom!= eidTo){
		for(int idE=(eidFrom+1)%nlines;idE!= eidTo; idE=(idE+1)%nlines){
			LPlines1[id1]=m_plines[idE];
			LPlines1.m_concavity[id1++]=m_concavity[idE];
		}
	}
	if(PidTo){
		LPlines1[id1]=saveToLp1;
		LPlines1.m_concavity[id1++]=m_concavity[eidTo];
	}
	LPlines1[id1]=mgTL2LPline(bridge);
	LPlines1[id1++].negate();

	assert(id1<=MAX_LINE_NUM && id2<=MAX_LINE_NUM);
	if(LPlines1.getNumberOfLines()<=2 || LPlines2.getNumberOfLines()<=2){
		//This case occurs when illegal boudary data, which is neglected.
		LPlines1.setNull();
		LPlines2.setNull();
	}
}

//This is a proprietry func of getToPointToSubdivide, get to point to subdivide
//by isectSlTl.
//Functin's return value is true if obtained.
bool mgTL2LPlines::getToPointByIsect(
	const mgTLEdgePoint& from,//id of edge and point that indicates which edge be subdivided.
	mgTLEdgePoint& to//to point is output.
)const{
	const int& eidFrom=from.m_edgeID; const int& PidFrom=from.m_pointID;

	int nlines=getNumberOfLines();assert(2<nlines &&nlines<=5);
	int eidNext=(eidFrom+1)%nlines;
	int eidOpo= nlines>=4 ? (eidFrom+2)%nlines : eidNext;
	int eidPre=(eidFrom+nlines-1)%nlines;
	int eidPre2 = nlines==5 ? (eidFrom + 3)%5 : eidPre;
	const mgTL2LPline& lp=m_plines[eidFrom];
	const mgTL2LPline& lpPre=m_plines[eidPre];
	int nPre=lpPre.number_of_points();

	MGPosition P=lp.uv(PidFrom);
	MGPosition Ppre=(PidFrom==0 ? lpPre.uv(nPre-2):lp.uv(PidFrom-1));
	MGVector dir1=P-Ppre;
	MGVector dir2=lp.uv(PidFrom+1)-P;
	MGStraight slMid(dir1,dir2,P);//subdividing straight.

	int id[5]={ eidOpo, eidNext, eidPre2,  eidFrom, eidFrom};//ordered from most probable edge.
	if (nlines <= 3) {
		id[nlines - 1] = eidFrom;
		if (nlines == 3)
			id[1] = eidPre;
	}else if (nlines == 5)
		id[4] = eidPre;

	//get an intersection with all sides of this polygon.
	for(int i=0; i<nlines; i++){
		double t;
		if (m_plines[id[i]].isectSlTl(slMid, t)) {//If intersection found.
			//change the parameter value to indicate a vetex point of the polygon.
			to.m_edgeID = id[i];
			changePolylineParameterToId(from,t, to);
			return true;
		}
	}
	return false;
}

//Change TL2Polyline's parameter t to mgTLEdgePoint.
//Midpoint value is changed to the best  nearest point id.
void mgTL2LPlines::changePolylineParameterToId(
	const mgTLEdgePoint& from,//id of edge and point that indicates which edge be subdivided.
	double t,//is a parameter value of  m_lines[to.edgeID] obtained by isectSlTl.
	mgTLEdgePoint& to//to.m_edgeID is input, and to.m_edgeID and to.m_pointID will be updated
		//to indicate best point id of this mgTL2LPlines
		//from parameter t of m_lines[to.edgeID].
)const{
	const mgTL2LPline& lto = m_plines[to.m_edgeID];
	t -= lto.m_polyline->param_s();
	t = lto.m_nV > 0 ? t - double(lto.m_idS) : double(lto.m_idS) - t;
	double delta = t - double(int(t));

	int& toPointID = to.m_pointID;
	int id = int(t);
	int idp1 = id + 1;
	if (delta <= 0.333333)
		toPointID = id;
	else if (delta >= 0.666666)
		toPointID = idp1;
	else{
		const mgTL2LPline& lfrom = m_plines[from.m_edgeID];
		toPointID =
			double(from.m_pointID)<double(lfrom.number_of_points())*.5 ? idp1 : id;
	}
	avoidAnomaly(from, to);
}

void mgTL2LPlines::avoidEndPoint(
	//const mgTLEdgePoint& from,//id of edge and point that indicates which edge be subdivided.
	mgTLEdgePoint& to//id of edge and point that indicates which edge be subdivided.
)const {
	const mgTL2LPline& line = m_plines[to.m_edgeID];
	if (line.number_of_points() >= 3) {
		if (line.isEndPoint(to.m_pointID))
			to = decrement(to);
		else if (to.m_pointID == 0)
			to = increment(to);
	}
}

//Update (eidTo, PidTo) data to avoid
// (1) to generate more than 5 line mgTL2LPlines.
// (2) 'to' point to be the same edge of from edge.
// (3) to be an end point.
void mgTL2LPlines::avoidAnomaly(
	const mgTLEdgePoint& from,//id of edge and point that indicates which edge be subdivided.
	mgTLEdgePoint& to//id of edge and point that indicates which edge be subdivided.
)const{
	/*std::cout << std::endl <<
		"avoidMoreThan5Lines, to=("<< to.m_edgeID << "," << to.m_pointID << ")" << std::endl;*/
	int nlines=getNumberOfLines();
	const int& eidFrom=from.m_edgeID; const int& PidFrom=from.m_pointID;

	int eidNext = (eidFrom + 1) % nlines;
	int eidOpo = nlines >= 4 ? (eidFrom + 2) % nlines : eidNext;
	int eidPre = (eidFrom + nlines - 1) % nlines;

	mgTLEdgePoint to2(to);
	normalizePointID(to2);
	int& eidTo = to2.m_edgeID; int& PidTo = to2.m_pointID;
	switch (nlines) {
	case(2): to = mgTLEdgePoint(eidOpo, nPLine(eidOpo)/2); break;

	case(4):
		if (eidTo == eidFrom){
			if (PidTo == 0)
				to = mgTLEdgePoint(eidPre, nPLine(eidPre) - 2);
			else if (PidFrom) {
				if (PidTo < PidFrom)
					to = mgTLEdgePoint(eidFrom, 0);
				else
					to = mgTLEdgePoint(eidNext, 0);
			}
		}else if (eidTo == eidNext){
			if(PidFrom == 0 && PidTo == 0)
				to = mgTLEdgePoint(eidNext, 1);
		}else if (eidTo == eidOpo) {//If opo.
			if (PidTo == 0)
				avoidEndPoint(to);
		}else //(eidTo == eidPre)
			if (PidFrom == 0 && PidTo == 0)//If from and to are on the same edge.
				to = mgTLEdgePoint(eidOpo, nPLine(eidOpo) - 2);
		break;
	case(3):
	case(5):
		avoidEndPoint(to); break;
	}
	normalizePointID(to);
	/*std::cout <<
		"            ,updated to=(" << to.m_edgeID << "," << to.m_pointID << ")" << std::endl;*/
}

//Get subdivide to-point, providing from point.
mgTLEdgePoint mgTL2LPlines::getToPointToSubdivide(
	const mgTLEdgePoint& from
	//indicates the edge and the point id where to subdivide,
)const{
	const int& eidFrom=from.m_edgeID; const int& PidFrom=from.m_pointID;
	int nlines=getNumberOfLines();assert(nlines<=5);

	int eidOpo= nlines>=4 ? (eidFrom+2)%nlines : (eidFrom+1)%nlines;
	int nOpo=m_plines[eidOpo].number_of_points();

	assert(!(nlines==2 && PidFrom==0)); 
	assert(!(nlines==2 && nOpo<=2)); 
	assert(!(nlines==3 && PidFrom==0 && nOpo<=2)); 
	assert(!(nlines==5 && PidFrom)); 
	mgTLEdgePoint to(eidOpo, nOpo / 2);//default for nlines == 2.
	if (nlines > 2 && !getToPointByIsect(from,to)){
		//This is a very illegal case, that logically does not take place.
		int& PidTo = to.m_pointID;
		if(nlines<=4){
			int nOpom2=nOpo-2;
			PidTo=(nOpom2<=PidFrom) ? nOpom2:PidFrom;
		}else
			PidTo=0;
	}
	return to;
}

///Debug Function
std::ostream& mgTL2LPlines::toString(std::ostream& ostrm)const{
	ostrm<<"mgTL2LPlines::"<<this;
	ostrm<<", m_polylines"<<std::endl;
	int nLines = getNumberOfLines();
	for(int i=0; i< nLines; i++)
		ostrm<<"["<<i<<"]="<<m_plines[i];
	ostrm << " m_concavity(";
	for (int i = 0; i < nLines; i++) {
		isUndefinedConcav(i) ?
			ostrm << "Undefined" : ostrm << m_concavity[i];
		ostrm << ",";
	}
	ostrm << ")" << std::endl;
	return ostrm;
}
