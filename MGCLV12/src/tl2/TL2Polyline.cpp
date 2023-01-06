#include "StdAfx.h"
#include "mg/Tolerance.h"
#include "mg/KnotArray.h"
#include "mg/Surface.h"
#include "mg/SurfCurve.h"
#include "topo/LEPoint.h"
#include "topo/Edge.h"
#include "topo/Loop.h"
#include "topo/Face.h"
#include "TL2/TL2Polyline.h"
#include "TL2/TL2Face.h"
#include "TL2/TL2LPlines.h"

/****************************************************************/
/*   Copyright (c) 2019 by System fugen G.K.                */
/*                       All rights reserved.                   */
/****************************************************************/


//mgTL2Polyline is a proprietry class for Face tessellation.
//mgTL2Polyline holds (u,v) MGLBRep of order 2(polyline).

//	enum polyline_type{
//		WHOLE_INNER=0,	//whole points are free inner points.
//		START_BOUNDARY,	//only start point is connected to a boundary line.
//		END_BOUNDARY,	//only end point is connected to a boundary line.
//		START_END_BOUNDARY,//Only start and end points are connected to a boundary.
//						//Inner points except start and end are not on a boundary.
//		WHOLE_BOUNDARY	//All of the points are on a boundary.
//	};


///Construct polyline MGLBRep whose maximum edge length is square_root(maxElen2).
void getXYZline_ensuring_max_edge_length(
	double maxElen2,	///<square of maximum edge length.
	const MGLBRep& xyzpolyline,///<Input original LBRep of (x,y,z) of order 2.
	MGLBRep& xyzpolylineOut///<Output (x,y,x) LBRep of order 2 whose knot vector is:
			///<t(i)=i-1 for i=1,...,n and t(0)=0 and t(n+1)=n-1.
){
	assert(xyzpolyline.order()==2);
	if(maxElen2<=0.)
		maxElen2=0.;

	const MGKnotVector& tOld=xyzpolyline.knot_vector();
	int nOld=xyzpolyline.bdim();
	const MGBPointSeq& pOld=xyzpolyline.line_bcoef();
	MGVector Ppre=pOld(0), Pnow;
	MGNDDArray tau(nOld*2);//nOld*2=The guessed length of the updated knot vetor.
	int j=0;//Index of tau.
	double tpre=tau(j++)=tOld[1];
	for(int i=1; i<nOld; i++){
		Pnow=pOld(i);
		MGVector vec=Pnow-Ppre;
		double vec2=vec%vec;
		double tnext=tOld[i+1];
		if(vec2>maxElen2){
			int nadd=int(sqrt(vec2/maxElen2));
			double diff=(tnext-tpre)/double(nadd+1);
			for(int k=0; k<nadd; k++){
				tpre+=diff;
				tau.store_with_capacityCheck(j++,tpre);
			}
		}
		tpre=tnext;
		tau.store_with_capacityCheck(j++,tpre);
		Ppre=Pnow;
	}
	tau.set_length(j);

	MGBPointSeq& bp=xyzpolylineOut.line_bcoef();
	MGKnotVector& t=xyzpolylineOut.knot_vector();
	int nnew=tau.length();
	bp.resize(nnew,3);t.size_change(2,nnew);
	for(j=0; j<nnew; j++){
		bp.store_at(j,xyzpolyline.eval(tau(j)));
		t(j+1)=double(j);//Set the uniform knot configuration of uvpolyline.
	}
	t(0)=t[1];
	t(nnew+1)=t[nnew];
}


///Construct polyline MGLBRep whose maximum edge length is para.get_max_edge_len_sqr().
void getUVline_ensuring_max_edge_length(
	const mgTL2parameter& para,///<Input parameter.
	const MGCurve& uvline,//Input original (u,v) line of the surface.
	MGLBRep& uvlineOut//Output LBRep of order 2 whose knot vector is:
			///t(i)=i-1 for i=1,...,n and t(0)=0 and t(n+1)=n-1.
){
	const MGSurface& surf=para.get_surface();//Target surface
	MGSurfCurve scrv(surf,uvline);

	double error=para.get_tess_crvError();
		//Error in world coordinates to approximate the line by a polyline.
	mgTolSetLineZero lineZeroSet(error);
	MGLBRep xyzpolyline;
	scrv.approximate_as_LBRep(xyzpolyline,2);//Here xyzpolyline is (x,y,z).
	lineZeroSet.restore();

	const MGKnotVector& tOld=xyzpolyline.knot_vector();
	int nOld=xyzpolyline.bdim();
	const MGBPointSeq& pOld=xyzpolyline.line_bcoef();

	double maxElen2=para.get_max_edge_len_sqr();
	double maxElen21=maxElen2*1.1;
	MGVector Ppre=pOld(0), Pnow;
	MGNDDArray tau(nOld*2);//nOld*2=The guessed length of the updated knot vetor.
	int j=0;//Index of tau.
	double tpre=tau(j++)=tOld[1];
	for(int i=1; i<nOld; i++){
		Pnow=pOld(i);
		MGVector vec=Pnow-Ppre;
		double vec2=vec%vec;
		double tnext=tOld[i+1];
		if(vec2>maxElen2){
			int nadd=int(sqrt(vec2/maxElen2));
			double diff=(tnext-tpre)/double(nadd+1);
			for(int k=0; k<nadd; k++){
				double tpresave=tpre;
				tpre+=diff;
				MGPosition Pnow2=scrv.eval(tpre);
				MGVector Vdiff2=Pnow2-Ppre;
				double vec21=Vdiff2%Vdiff2;
				if(vec21>=maxElen21){
					double diff2=diff*.5;
					tau.store_with_capacityCheck(j++,tpresave+diff2);
				}
				tau.store_with_capacityCheck(j++,tpre);
				Ppre=Pnow2;
			}
		}
		tpre=tnext;
		tau.store_with_capacityCheck(j++,tpre);
		Ppre=Pnow;
	}
	tau.set_length(j);

	MGBPointSeq& bp=uvlineOut.line_bcoef();
	MGKnotVector& t=uvlineOut.knot_vector();
	int nnew=tau.length();
	bp.resize(nnew,2);t.size_change(2,nnew);
	for(j=0; j<nnew; j++){
		bp.store_at(j,uvline.eval(tau(j)));
		t(j+1)=double(j);//Set the uniform knot configuration of uvpolyline.
	}
	t(0)=t[1];
	t(nnew+1)=t[nnew];
}

//////////// constructor ///////////////

//copy constructor.
mgTL2Polyline::mgTL2Polyline(const mgTL2Polyline& pline2)
:MGLBRep(pline2),m_tlparam(pline2.m_tlparam),m_type(pline2.m_type),
m_uv_start(pline2.m_uv_start), m_dire_start(pline2.m_dire_start),
m_dire_end(pline2.m_dire_end){
	for(int i=0; i<3; i++){
		m_start[i]=pline2.m_start[i];
		m_end[i]=pline2.m_end[i];
	}
}

//move constructor.
mgTL2Polyline::mgTL2Polyline(mgTL2Polyline&& pline2)noexcept
:MGLBRep(std::move(pline2)),m_tlparam(pline2.m_tlparam),m_type(pline2.m_type),
m_uv_start(std::move(pline2.m_uv_start)), m_dire_start(std::move(pline2.m_dire_start)),
m_dire_end(std::move(pline2.m_dire_end)){
	for(int i=0; i<3; i++){
		m_start[i]=pline2.m_start[i];
		m_end[i]=pline2.m_end[i];
	}
}
//move assignment.
mgTL2Polyline& mgTL2Polyline::operator=(mgTL2Polyline&& other)noexcept{
	m_tlparam=other.m_tlparam, m_type=other.m_type;
	for(int i=0; i<3; i++){
		m_start[i]=other.m_start[i];
		m_end[i]=other.m_end[i];
	}
	m_uv_start=std::move(other.m_uv_start);
	m_dire_start = std::move(other.m_dire_start);
	m_dire_end = std::move(other.m_dire_end);
	MGLBRep::operator=(std::move(other));
	return *this;
}


//Construct mgTL2Polyline from (u,v) curve representation of the surface.
//The type of the boundary is WHOLE_INNER.
mgTL2Polyline::mgTL2Polyline(
	const mgTL2parameter& para,
	const MGCurve& uvline
):m_tlparam(&para), m_type(WHOLE_INNER) {
	getUVline_ensuring_max_edge_length(para, uvline, *this);
}

//Function's return value is
//true if the point i t is a boundary point, false if not.
bool get_LEpointID(const MGLEPoint& le, short id[3]){
	const mgTL2Polyline* ipoly=TL2Polyline(le.iterator());
	return ipoly->get_id(le.param(),id);
}

//Construct dummy mgTL2Polyline of lengh n.
mgTL2Polyline::mgTL2Polyline(
	const mgTL2parameter& para,
	int n //nubmer of points to 
):MGLBRep(n, 2, 2) //of length=n, order=2, space dimension =2
, m_tlparam(&para), m_type(WHOLE_INNER) {;}

//Construct mgTL2Polyline of straight line of point number n from uvS to uvE.
//The type of the boundary is WHOLE_INNER, and so necessary to set_endID/set_startID,
//if necessary.
mgTL2Polyline::mgTL2Polyline(
	const mgTL2parameter& para,
	int n, //nubmer of points to 
	const MGPosition& uvE,
	const MGPosition& uvS
):mgTL2Polyline(para,n) //of length=n, order=2, space dimension =2
{
	double dnm1 = double(n - 1);
	MGStraight sl(uvE, uvS, dnm1);//end parameter is dnm1.
	MGBPointSeq& bp = line_bcoef();
	MGKnotVector& t = knot_vector();
	t(0)=0.;
	for (int i = 0; i < n; i++) {
		double tip1=t(i + 1) = double(i);
		bp.store_at(i, sl.eval(tip1));
	}
	t(n+1) = dnm1;
}

//Construct mgTL2Polyline from the straight from pvS to pvE.
//The type of the boundary is set according the boundary infromation of pvS and pvE.
mgTL2Polyline::mgTL2Polyline(
	const mgTL2LPlines& lines,
	const mgTLEdgePoint& pvE,
	const mgTLEdgePoint& pvS
):mgTL2Polyline(lines[0].TL2param(),2,
	lines[pvE.m_edgeID].uv(pvE.m_pointID), lines[pvS.m_edgeID].uv(pvS.m_pointID)){
	setBoundaryID(mgTL2LPPoint(lines[pvE.m_edgeID], pvE.m_pointID));
	setBoundaryID(mgTL2LPPoint(lines[pvS.m_edgeID], pvS.m_pointID),false);
}

//Construct mgTL2Polyline from the straight of pvS to pvE.
//The straight is guaranteed to keep the deviation from the surface
//within tessellation curve error and within the maximum edge length.
//The type of the boundary is set according the boundary infromation of pvS and pvE.
mgTL2Polyline::mgTL2Polyline(
	const MGLEPoint& pvE,
	const MGLEPoint& pvS
):m_type(WHOLE_INNER) {
	MGPosition uv0=pvS.eval(), uv1=pvE.eval();
	MGStraight sl(uv1,uv0);

	const mgTL2Polyline* ipoly=TL2Polyline(pvE.iterator());
	const mgTL2parameter& para=ipoly->TL2param();
	m_tlparam = &para;
	getUVline_ensuring_max_edge_length(para, sl, *this);

	short id[3];
	if(get_LEpointID(pvS,id))
		set_startID(id);
	if(get_LEpointID(pvE,id))
		set_endID(id);
}


//Compute the angle of this curves's world rep and u=const iso-parameter line world rep
//at the parameter t. The parameter t is normalized value from 0.(start point)
//to 1.(end point).
double mgTL2Polyline::angle2Uline(double t)const{
	const MGSurface& srf=TL2param().get_surface();
	MGSurfCurve crvOnSrf(srf,*this);
	MGPosition uvt;
	MGVector Vcrv;
	double para_t=param_s();
	if(t==0.){
		uvt=uv_start();
		Vcrv=direWorld_start();
	}else{
		para_t+=(param_e()-param_s())*t;
		uvt=eval(para_t);
		Vcrv=crvOnSrf.eval(para_t,1);
	}
	MGVector Vu=srf.eval(uvt,1);
	MGUnit_vector N=srf.unit_normal(uvt);
	return Vu.angle2pai(Vcrv,N);
}

//Compute the angle of this curves's world rep and u=const iso-parameter line world rep
//at the point id i.
double mgTL2Polyline::angle2Uline_at_i(int i)const{
	const MGSurface& srf=TL2param().get_surface();
	MGSurfCurve crvOnSrf(srf,*this);
	MGPosition uvt;
	MGVector Vcrv;
	double para_t=param_s();
	if(i==0){
		uvt=uv_start();
		Vcrv=direWorld_start();
	}else{
		uvt=uv(i);
		if(i==number_of_points()-1)
			Vcrv=direWorld_end();
		else{
			para_t=param_s()+double(i);
			Vcrv=crvOnSrf.eval(para_t,1);
		}
	}
	MGVector Vu=srf.eval(uvt,1);
	MGUnit_vector N=srf.unit_normal(uvt);
	return Vu.angle2pai(Vcrv,N);
}

///Construct new curve object by copying to newed area.
///User must delete this copied object by "delete".
mgTL2Polyline* mgTL2Polyline::clone()const{
	return new mgTL2Polyline(*this);
}

//Evaluate the direction at t(not normalized ordinary parameter value) in its
//world representation.
MGVector mgTL2Polyline::direWorld_with_non_normalized(double t)const{
	const mgTL2parameter& para=TL2param();
	const MGSurface& srf=para.get_surface();
	MGSurfCurve crvAftOnsrf(srf,*this);
	return crvAftOnsrf.eval(t,1);
}

//Get the direction of the line in the world coordinate at the normalized parameter t.
//The parameter t is normalized value from 0.(start point) to 1.(end point).
MGVector mgTL2Polyline::direWorld(double t)const{
	if(t == 0.)
		return direWorld_start();
	if(t == 1.)
		return direWorld_end();

	double s=param_s();
	s=s+(param_e()-s)*t;
	return direWorld_with_non_normalized(s);
}

//Get the start or end point direction of the line in the world coordinate.
const MGVector& mgTL2Polyline::direWorld_start()const{
	if(m_dire_start.is_null()){
		double ts=param_s();
		m_dire_start=direWorld_with_non_normalized(ts);
	}
	return m_dire_start;
}
const MGVector& mgTL2Polyline::direWorld_end()const{
	if(m_dire_end.is_null()){
		double te=param_e();
		m_dire_end=direWorld_with_non_normalized(te);
	}
	return m_dire_end;
}

//Get id of m_Bpolylines of mgTL2Face from point id of this polyline.
//Function's return value is
//true if the point of idVertex is a boundary point, false if not.
bool mgTL2Polyline::get_id_from_VertexID(int idVertex, short id[3])const{
	int npm1=number_of_points()-1;
	assert(0<=idVertex && idVertex<=npm1);

	switch(m_type){
	case WHOLE_INNER: return false;
	case WHOLE_BOUNDARY:
		for(int i=0; i<3; i++)
			id[i]=m_start[i];
		if(m_end[2]>=m_start[2])
			id[2]+=idVertex;
		else
			id[2]-=idVertex;
		return true;
	case START_BOUNDARY:
		if(idVertex)
			return false;
		for(int i=0; i<3; i++)
			id[i]=m_start[i];
		return true;
	case END_BOUNDARY:
		if(idVertex==npm1){
			for(int i=0; i<3; i++)
				id[i]=m_end[i];
			return true;
		}
		return false;;
	case START_END_BOUNDARY:
		if(idVertex==0){
			for(int i=0; i<3; i++)
				id[i]=m_start[i];
			return true;
		}else if(idVertex==npm1){
			for(int i=0; i<3; i++)
				id[i]=m_end[i];
			return true;
		}
	}
	return false;;
}

//Get id of m_Bpolylines of mgTL2Face from the parameter t.
//Function's return value is
//true if the point i t is a boundary point, false if not.
bool mgTL2Polyline::get_id(double t, short id[3])const{
	short idt=short(t-param_s()+.5);
	return get_id_from_VertexID(idt,id);
}

void mgTL2Polyline::set_startID(const short start[3]){
	for(int i=0; i<3; i++){
		m_start[i]=start[i];
	}
	if(m_type==WHOLE_INNER)
		m_type=START_BOUNDARY;
	else if(m_type==END_BOUNDARY)
		m_type=START_END_BOUNDARY;
}
void mgTL2Polyline::set_endID(const short end[3]){
	for(int i=0; i<3; i++){
		m_end[i]=end[i];
	}
	if(m_type==WHOLE_INNER)
		m_type=END_BOUNDARY;
	else if(m_type==START_BOUNDARY)
		m_type=START_END_BOUNDARY;
}

void mgTL2Polyline::get_startID(short start[3]){
	for(int i=0; i<3; i++){
		start[i]=m_start[i];
	}
}
void mgTL2Polyline::get_endID(short end[3]){
	if(m_type==WHOLE_BOUNDARY){
		end[0]=m_start[0];
		end[1]=m_start[1];
		end[2]=m_end[2];
	}else{
		for(int i=0; i<3; i++)
			end[i]=m_end[i];
	}
}

void mgTL2Polyline::setBoundaryID(
	const mgTL2LPPoint& P, //End point.
	bool isEnd
) {
	short ids[3];
	if (P.m_line.get_id_from_VertexID(P.m_id, ids))
		isEnd ? set_endID(ids) : get_startID(ids);
}

///Change direction of the line.
void mgTL2Polyline::negate(){
	MGLBRep::negate();
	m_uv_start.set_null();
	m_dire_start.set_null();
	m_dire_end.set_null();
	if(m_type==START_BOUNDARY)
		m_type=END_BOUNDARY;
	else if(m_type==END_BOUNDARY)
		m_type=START_BOUNDARY;
	else if(m_type==WHOLE_INNER)
		return;
	for(int i=0; i<3; i++){
		short save=m_start[i];
		m_start[i]=m_end[i];
		m_end[i]=save;
	}
}

//Get i-th point surface parameter (u,v) of this polyline
MGPosition mgTL2Polyline::uv(int i)const{
	assert(i>=0 && i<number_of_points());
	if(i==0)
		return uv_start();
	const MGBPointSeq& uvbp=line_bcoef();
	return uvbp(i);
}

//Get start point surface parameter (u,v) of this polyline
const MGPosition& mgTL2Polyline::uv_start()const{
	if(m_uv_start.is_null())
		m_uv_start=line_bcoef()(0);
	return m_uv_start;
}

//Get i-th point(x,y,z,xn,yn,zn) of this polyline.
//Here (x,y,z) is the position data, and (xn,yn,zn) is the unit normal at (x,y,z).
MGPosition mgTL2Polyline::xyz(int i, bool need_normal)const{
	assert(i>=0 && i<number_of_points());

	const mgTL2parameter& para=TL2param();
	const MGSurface& srf=para.get_surface();
	MGPosition uvi=uv(i);
	MGPosition N;
	int sdim=need_normal ? 6:3;
	MGPosition xyzPN(sdim);
	if(need_normal){
		N=srf.unit_normal(uvi);
		xyzPN.store_at(3,N);
	}
	short ids[3];
	bool on_boundary=get_id_from_VertexID(i,ids);
	if(on_boundary){
		const std::vector<SHLL_COM_EDGES>& bpolyls=*(para.Bpoly());
		SHLL_COM_EDGES edges=bpolyls[ids[0]];
		const MGLBRep* lb=edges[ids[1]];
		xyzPN.store_at(0,lb->line_bcoef()(ids[2]));
	}else{
		xyzPN.store_at(0,srf.eval(uvi));
	}
	return  xyzPN;
}

//Update this by limiting the parameter range of the curve.
///Limitting is done at the knot parameter for both start and end.
void mgTL2Polyline::limit(const MGInterval& i1){
	short nOldm1=bdim()-1;
	MGKnotVector& t=knot_vector();
	double ts=param_s(), te=param_e();
	if(i1.includes(ts) && i1.includes(te))
		return;

	m_dire_end.set_null();
	m_dire_start.set_null();
	m_uv_start.set_null();

	MGInterval i1Limit=param_range();i1Limit&=i1;
	short idS=short(i1Limit.low_point()-ts+.5);
	short idE=short(i1Limit.high_point()-ts+.5);
	int n=idE-idS+1;//New vertices number.
	if(n<=1)
		n=2;

	MGBPointSeq& bp=line_bcoef();
	if(idS){
		for(int k=0; k<2; k++)
			for(int j=0; j<n; j++)
				bp(j,k)=bp(j+idS,k);
	}
	bp.reshape(n);
	t.size_change(2,n);
	double tsnew=ts+double(idS);
	for(int i=0; i<n; i++)
		t(i+1)=tsnew+double(i);
	short nm1=n-1;
	t(0)=t[1];
	t(n+1)=t[n];

	short sid=m_start[2];
	switch(m_type){
		case WHOLE_INNER: break;
		case WHOLE_BOUNDARY:
			if(sid<=m_end[2]){
				m_start[2]=sid+idS;
				m_end[2]=sid+idE;
			}else{
				m_start[2]=sid-idS;
				m_end[2]=sid-idE;
			}
			break;
		case START_BOUNDARY:
			if(idS)
				set_type(WHOLE_INNER);
			break;
		case END_BOUNDARY:
			if(idE<nOldm1)
				set_type(WHOLE_INNER);
			break;
		case START_END_BOUNDARY:
			if(idS)
				set_type(END_BOUNDARY);
			if(idE<nOldm1)
				if(m_type==END_BOUNDARY)
					set_type(WHOLE_INNER);
				else
					set_type(START_BOUNDARY);
	}
}

///Debug Function
std::ostream& mgTL2Polyline::toString(std::ostream& ostrm)const{
	ostrm<<"mgTL2Polyline::"<<this;
	ostrm<<", m_type=";
	switch(m_type){
		case WHOLE_INNER: ostrm<<"WHOLE_INNER";break;
		case WHOLE_BOUNDARY:
			ostrm<<"WHOLE_BOUNDARY";
			ostrm<<std::endl<<"m_start=("<<m_start[0]<<","<<m_start[1]<<","<<m_start[2]<<")";
			break;
		case START_BOUNDARY:
			ostrm<<"START_BOUNDARY";
			ostrm<<std::endl<<"m_start=("<<m_start[0]<<","<<m_start[1]<<","<<m_start[2]<<")";
			break;
		case END_BOUNDARY:
			ostrm<<"END_BOUNDARY";
			ostrm<<std::endl<<"m_end=("<<m_end[0]<<","<<m_end[1]<<","<<m_end[2]<<")";
			break;
		case START_END_BOUNDARY:
			ostrm<<"START_END_BOUNDARY";
			ostrm<<std::endl<<"m_start=("<<m_start[0]<<","<<m_start[1]<<","<<m_start[2]<<")";
			ostrm<<"m_end=("<<m_end[0]<<","<<m_end[1]<<","<<m_end[2]<<")";
			break;
	}
	ostrm<<","<<std::endl;
	MGLBRep::toString(ostrm);
	return ostrm;
}

//Update the point number. The parameter range is also updated according to the newN.
mgTL2Polyline& mgTL2Polyline::updatePointNumber(int newN) {
	double dnewNm1= double(newN - 1);
	change_range(0., dnewNm1);

	MGLBRep newLb(newN, 2,2);
	MGBPointSeq& newBp=newLb.line_bcoef();
	MGKnotVector& newT=newLb.knot_vector();
	newT(0) = 0.;
	for (int i = 0; i < newN; i++) {
		double tip1 = newT(i + 1) = double(i);
		newBp.store_at(i, eval(tip1));
	}
	newT(newN + 1) = dnewNm1;
	MGLBRep& lbThis=*this;
	lbThis = std::move(newLb);
	return *this;
}

//Get curve representation(MGLBRep of uniform B-Spline of order 2).
const mgTL2Polyline* TL2Polyline(const MGEdge* edg){
	return static_cast<const mgTL2Polyline*>(edg->base_curve());
}

const mgTL2Polyline* TL2Polyline(MGComplex::const_iterator ei){
	const MGEdge* edg=edge_from_iterator(ei);
	return TL2Polyline(edg);
}
