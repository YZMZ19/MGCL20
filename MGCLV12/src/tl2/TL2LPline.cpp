#include "StdAfx.h"
#include "mg/Tolerance.h"
#include "mg/Straight.h"
#include "mg/Surface.h"
#include "topo/Face.h"
#include "Tl2/TL2parameter.h"
#include "TL2/TL2LPline.h"
#include "TL2/TL2Polyline.h"

/****************************************************************/
/*   Copyright (c) 2021 by System fugen G.K.                */
/*                       All rights reserved.                   */
/****************************************************************/


///mgTL2LPline is a proprietry class for Face tessellation.
///mgTL2LPline is limitted subinterval of mgTL2Polyline,
///holds the starting id of mgTL2Polyline, and the number of vertices.

//	enum polyline_type{
//		WHOLE_INNER=0,	//whole points are free inner points.
//		START_BOUNDARY,	//only start point is connected to a boundary line.
//		END_BOUNDARY,	//only end point is connected to a boundary line.
//		START_END_BOUNDARY,//Only start and end points are connected to a boundary.
//						//Inner points except start and end are not on a boundary.
//		WHOLE_BOUNDARY	//All of the points are on a boundary.
//	};


//////////// constructor ///////////////

///Construct from subinterval of input lpline.
mgTL2LPline::mgTL2LPline(
	const mgTL2LPline& lpline,
	int idS, //starting id of lpline
	int nV
):m_polyline(lpline.m_polyline), m_sharedLine(lpline.m_sharedLine)
,m_idS(lpline.m_idS),m_nV(lpline.m_nV),
m_concavity(NotObtainedConcav){
	limit(idS,nV);
	assert(m_idS<m_polyline->bdim());
}

///Construct from subinterval of input lpline.
mgTL2LPline::mgTL2LPline(
	const mgTL2LPPoint& lpPoint,
	int nV
) :mgTL2LPline(lpPoint.m_line, lpPoint.m_id, nV){
	;
}

void mgTL2LPline::setNull(){
	m_polyline=nullptr;
	m_sharedLine.reset();
	m_idS=0;
	m_nV=0;
}

//Set the start point and the number of points to the whole line of m_polyline.
void mgTL2LPline::setWholeLine() {
	m_idS = 0;
	m_nV = m_polyline->bdim();
	m_concavity = NotObtainedConcav;
}

///Construct new curve object by copying to newed area.
///User must delete this copied object by "delete".
mgTL2LPline* mgTL2LPline::clone()const{
	return new mgTL2LPline(*this);
}

//Evaluation of this with the normalized parameter value t from 0. to 1.
//is provided. t=0. for the start point, t=1. for the end point.
MGVector mgTL2LPline::eval(double t, int nderi)const{
	assert(0.<=t && t<=1.);
	const mgTL2Polyline& pl=*(TL2Polyline());
	const MGKnotVector& kv=pl.knot_vector();
	double ts=kv(m_idS+1), te;
	if(m_nV>0){
		te=kv(m_idS+m_nV);
	}else{
		te=kv(m_idS+2+m_nV);
	}
	MGVector result=pl.eval(ts+(te-ts)*t,nderi);
	if(m_nV<0 && nderi%2)
		result.negate();
	return result;
}

///Get concavity of this edge.
///The concavity is obtained from the differece of two vectors,
///at the start and at the end point point tangent of this.
//Concavity's value is from -2 to 2. 2 is most concave, and
// -2 means 180 degree convex(most convex), -1:90 degree convex, 0:flat
// 1:90 degree concave, 2:180 degree concave.
double mgTL2LPline::getConcavity(
)const{
	if (isUndefinedConcavity(m_concavity)) {
		int n = number_of_points();
		MGVector Vpre = uv(1) - uv(0);
		MGVector Vaft = uv(n - 1) - uv(n - 2);
		m_concavity = MGCL::concavity(Vpre, Vaft);
	}
	return m_concavity;
}

//Get id of m_Bpolylines of m_param of m_polyline from vertex id i.
//Function's return value is
//true if the point of t is a boundary point, false if not.
bool mgTL2LPline::get_id_from_VertexID(int i, short id[3])const{
	assert(i<number_of_points());
	if(m_nV>0)
		return m_polyline->get_id_from_VertexID(i+m_idS, id);
	return m_polyline->get_id_from_VertexID(m_idS-i, id);
}

//Change this TL2LPline's point id to TL2Polyline's parameter t.
double mgTL2LPline::changeIdToPolylineParameter(int id)const{
	double t=m_nV>0 ? m_idS+id:t=m_idS-id;
	return t+m_polyline->param_s();
}

//Get i-th point surface parameter (u,v) of this polyline.
//i is relative one that start from 0 even for opposite direction.
MGPosition mgTL2LPline::uv(int i)const{
	assert(i<number_of_points());
	if(m_nV>0)
		return m_polyline->uv(i+m_idS);
	return m_polyline->uv(m_idS-i);
}

//Get i-th point(x,y,z,xn,yn,zn) of this polyline.
//Here (x,y,z) is the position data, and (xn,yn,zn) is the unit normal at (x,y,z).
//i is relative one that start from 0 even for opposite direction.
MGPosition mgTL2LPline::xyz(int i, bool need_normal)const{
	assert(i<number_of_points());
	if(m_nV>0)
		return m_polyline->xyz(i+m_idS,need_normal);
	return m_polyline->xyz(m_idS-i,need_normal);
}

//Update this by limiting the parameter range of the curve.
///Limitting is done at the knot parameter for both start and end.
void mgTL2LPline::limit(
	int idS,	//start point id of this mgTL2LPline.
				//idS is relative one that starts from 0 even for opposite direction.
	int nV	//Number of vertices.
){
	assert(nV<=number_of_points());
	if(m_nV>0){
		assert(idS<int(m_nV) && idS+nV<=int(m_nV));
		m_idS+=idS;
		m_nV=nV;
	}else{
		assert(int(idS)<-m_nV && int(idS+nV)<=-m_nV);
		m_idS-=idS;
		m_nV=-int(nV);
	}
	assert(m_idS<m_polyline->bdim() && m_idS + m_nV +1>=0);
	m_concavity = NotObtainedConcav;
}

//Obtain the mid point of this line.
void mgTL2LPline::mid(MGPosition& uvmid){
	int nmid=number_of_points()/2;
	uvmid=uv(nmid);
}

//Get the number of points of this closed polygon.
int mgTL2LPline::number_of_points()const{
	if(m_nV>0)
		return m_nV;
	return -m_nV;
}

//Reverse the direction.
mgTL2LPline& mgTL2LPline::negate(){
	if(m_nV>0){
		m_idS+=m_nV-1;
		m_nV=-m_nV;
	}else{
		m_idS=int(m_idS+1+m_nV);
		m_nV=-m_nV;
	}
	assert(m_idS<m_polyline->bdim());
	if (!isUndefinedConcavity(m_concavity))
		m_concavity *= -1.;
	return *this;
}

//Subdivide at the id.
void mgTL2LPline::subdivide(
	int id, //Relative one that start from 0 even for opposite direction.
	mgTL2LPline& lp1,	//devided 1st mgTL2LPline.
	mgTL2LPline& lp2	//devided 2nd mgTL2LPline.
)const{
	assert(m_idS<m_polyline->bdim());
	if(id<=0){
		lp2 = *this;
		lp1.setNull();
	}else if(id>=number_of_points()-1){
		lp1=*this;
		lp2.setNull();
	}else{
		lp1.m_sharedLine = m_sharedLine;
		lp2.m_sharedLine = m_sharedLine;
		lp1.m_polyline=lp2.m_polyline=m_polyline;
		lp1.m_idS=m_idS;
		if(m_nV>0){
			assert(int(id)<m_nV);
			lp1.m_nV=int(id+1);
			lp2.m_idS=m_idS+id;
			lp2.m_nV=m_nV-int(id);
		}else{
			assert(int(id)<(-m_nV));
			lp1.m_nV=-int(id+1);
			lp2.m_idS=m_idS-id;
			lp2.m_nV=m_nV+int(id);
		}
		lp1.m_concavity = lp2.m_concavity = NotObtainedConcav;
	}
}

//Compute the intersections of sl and this mgTL2LPline.
//Function's return value is the (nearest) intersection vertex id of lp if >=0.
//If return value is minus, intersection not found.
bool mgTL2LPline::isectSlTl(
	const MGStraight& sl,
	double& t//parameter value of the intersection.
)const{
	const mgTL2parameter& tlparam=TL2param();
	const MGSurface& suface=tlparam.get_surface();
	double serror = suface.parameter_error();
	const MGBox& prange=suface.param_range();

	double uError = tlparam.get_UError(), vError = tlparam.get_VError();
	mgTolSetWCZero wczeroSet(serror);//Set&save the error.

	const MGPosition& Ps=sl.root_point();
	const mgTL2Polyline& ipoly=*TL2Polyline();
	double ts=ipoly.param_s()+double(m_idS);
	double te=ts+double(m_nV-1);
	if(m_nV<0){
		te=ts;
		ts=te+double(m_nV+1);
	}
	double mzero=MGTolerance::mach_zero();
	ts-=mzero; te+=mzero;
	MGCCisects tlist=ipoly.isect(sl);
	MGCCisects::iterator i=tlist.begin(), iend=tlist.end();
	for(; i!=iend; i++){
		MGCCisect& ti=isectCast<MGCCisect>(i);
		MGPosition& Pi=ti.point();
		if(fabs(Pi[0]-Ps[0])<=uError  && fabs(Pi[1]-Ps[1])<=vError)
			continue;

		t=ti.param1();
		if(ts<=t && t<=te)
			return true;
	}
	return false;
}


//Polygonize the (u,v) space straight line from this->uv(id1V) to pline2.uv(id2V).
//The direction of the output is from id1V to id2V.
//polygonizeSL does ensure the deviation from the surface within the surface
//tolerance.
void mgTL2LPline::polygonizeSL(
	const mgTL2LPline& pline2,
	int id1V,	//id of this vertex.
	int id2V,	//id of pline2's vertex.
	mgTL2Polyline& uvpolyline
)const{
	MGPosition uvS=uv(id1V);
	MGPosition uvE=pline2.uv(id2V);
	MGStraight sl(uvE,uvS);
	uvpolyline=mgTL2Polyline(TL2param(),sl);
	uvpolyline.setBoundaryID(mgTL2LPPoint(*this, id1V),false);//For start of uvpolyline.
	uvpolyline.setBoundaryID(mgTL2LPPoint(pline2, id2V));
}

///Debug Function
std::ostream& mgTL2LPline::toString(std::ostream& ostrm)const{
	ostrm<<"mgTL2LPline::"<<this;
	ostrm<<", m_polyline="<<m_polyline<<", m_sharedLine="<<m_sharedLine
		<<", m_idS="<<m_idS<<", m_nV="<<m_nV<<std::endl;
	if(m_polyline)
		ostrm<<(*m_polyline);
	ostrm << ", concavity=";
	isUndefinedConcavity(m_concavity) ? ostrm << "Undefined" : ostrm << m_concavity;
	ostrm<<std::endl;
	return ostrm;
}

///get the (u,v) parameter box.
MGBox mgTL2LPline::getUvBox()const{
	MGBox uvbox;
	int nP=number_of_points();
	for(int i=0; i<nP; i++)
		uvbox|=uv(i);
	return uvbox;
}

///Extract this as a polygon MGLBRep whose parameter range is [0,abs(m_nV)-1]
void mgTL2LPline::extractAsLBRep(MGLBRep& poly)const{
	int idS, idE;//Knot vector id of this LBRep of order 2.
	if (m_nV < 0) {
		idE = m_idS+1;
		idS = m_idS + m_nV + 2;
	} else {
		idS = m_idS+1;
		idE = m_idS + m_nV;
	}
	m_polyline->shrinkToKnots(idS, idE, poly);
	if (m_nV < 0)
		poly.negate();
	if (idS >= 2)
		poly.change_range(0., double(idE - idS));
}
