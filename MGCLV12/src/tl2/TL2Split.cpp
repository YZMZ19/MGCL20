/********************************************************************/
/* Copyright (c) 2021 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/

#include "StdAfx.h"
#include <iomanip>
#include "mg/Box.h"
#include "mg/Curve.h"
#include "mg/Tolerance.h"
#include "mg/CParam_list.h"
#include "mg/SurfCurve.h"
#include "topo/LEPoint.h"
#include "topo/CPointsVec.h"
#include "topo/LCisect_vector.h"
#include "topo/Loop.h"
#include "topo/Edge.h"
#include "topo/Face.h"
#include "Tl2/TL2Polyline.h"
#include "Tl2/TL2Face.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

mgTL2PlBridge::mgTL2PlBridge(
	const MGLEPoint& le_start,//line's start point's MGLEPoint.
	const MGLEPoint& le_end//line's end point's MGLEPoint.
) :m_polyLine(new mgTL2Polyline(le_end, le_start)), m_le_start(le_start), m_le_end(le_end) {
	;
}

//Make a cop of m_polyLine and make MGEdge.
//Function's return value is newed MGEdge pointer.
MGEdge* mgTL2PlBridge::copy_edge() {
	mgTL2Polyline* line = new mgTL2Polyline(*m_polyLine);
	return new MGEdge(line);
}

//Release m_polyLine object and make MGEdge.
//Function's return value is newed MGEdge pointer.
//After release_edge, m_polyLine is null.
MGEdge* mgTL2PlBridge::release_edge() {
	mgTL2Polyline* line = m_polyLine.release();
	return new MGEdge(line);
}

std::ostream& operator<< (std::ostream& out, const mgTL2PlBridge& bridge) {
	out << "mgTL2PlBridge:" << *(bridge.m_polyLine);
	out << "start=" << bridge.m_le_start << ", end=" << bridge.m_le_end
		<< std::endl << std::endl;
	return out;
}

const MGVector& getSEDirection(const MGLEPoint& Vrtx, bool start) {
	const MGEdge* edge = Vrtx.edge();
	return start ?
	 TL2Polyline(edge)->direWorld_start() :
	 TL2Polyline(edge->pre_edge())->direWorld_end();
}

//Split the face giving bridges over loops. Splitting is done by finding the smallest
//closed areas out of bridges.
bool split_by_bridges(
	const MGFace& face,
	const std::vector<UniqueTL2PlBridge>& bridges,
	std::stack<UniqueFace>& faces//Result trimmed face(s) will be appended.
){
	int newFaceNum=0;
	double rzero=MGTolerance::rc_zero();
	const MGBox& bxuv=face.box_param();
	double minimum_loop_area=bxuv[0].length();
	minimum_loop_area+=bxuv[1].length().value();
	minimum_loop_area*=rzero;

	mgCPointsVec cpointsVec(face);
		//Found loops that are connected to
		//an inner loop of loop id loop_id will be stored
		//in this cpointsVec[loop_id].

	int nbridges=(int)bridges.size();
	for(int ii=nbridges-1; ii>=0; ii--)
		cpointsVec.build_CPointsVec(*(bridges[ii]));
	int nboundaries=face.number_of_boundaries();
	for(int i=0; i<nboundaries; i++){
		mgCPoints& lveci=cpointsVec[i];
		lveci.sort();
	}

	for(int i=0; i<nboundaries; i++){
		mgCPoints& lveci=cpointsVec[i];
		int nlv=lveci.size();
		assert((nlv%2)==0);//nlv must be even(pair of coming_in and going_out).
		if(!nlv)
			continue;

		for(int j=0; j<nlv; j++){
			std::unique_ptr<MGFace> nface2(new MGFace);
			std::unique_ptr<MGLoop> loop;
			std::deque<MGTrimLoop*> tloops;
			mgCONNECT_POINT& fpoint=lveci[j];
			int extract_inf=cpointsVec.extract_loop(fpoint,loop,tloops);
			if(!extract_inf)
				continue;
		
			std::vector<bool> used_loops(nboundaries,false);
			size_t ntloop=tloops.size();
			for(size_t k=0; k<ntloop; k++){
				tloops[k]->set_used_loop_flag(used_loops);
				tloops[k]->set_null();
			}

			double area=loop->area();
			if(area<=minimum_loop_area && area>= -minimum_loop_area)
				continue;

			nface2->add_boundary(loop.release());
			for(int j=nboundaries-1; j>=0; j--){
				if(!used_loops[j]){
					const UniqueLoop& lpj=face.loop(j);
					MGPosition uvS=lpj->start_point();
					if(nface2->in_range(uvS))
						nface2->add_boundary(face.loop(j)->clone());
				}
			}
			faces.emplace(nface2.release());
			newFaceNum++;
		}
	}
	if (newFaceNum <= 1 && faces.size())
		faces.pop();
	return newFaceNum>=2;
}

//Extract MGTrimLoop & mgCPoints info, and build this mgCPointsVec data.
void mgCPointsVec::build_CPointsVec(
	mgTL2PlBridge& bridge	//original loops to trim.
){
	for(int i=0; i<2; i++){
		MGLoop* lp=new MGLoop;
		MGLEPoint leS, leE;
		MGEdge* edg;
		if(i==0){
			leS=bridge.le_start(), leE=bridge.le_end();
			edg=bridge.copy_edge();
		}else{
			leE=bridge.le_start(), leS=bridge.le_end();
			edg=bridge.release_edge();
			edg->negate();
		}
		lp->append(edg);

		int lidS=leS.loop()->get_loop_id_in_face();
		int lidE=leE.loop()->get_loop_id_in_face();
		MGTrimLoop* trloop=new MGTrimLoop(lp,-lidS,leS,-lidE,leE);
		m_trim_loops.emplace_back(trloop);
		push_at(lidS,mgCONNECT_POINT(trloop,mgCONNECT_POINT::going_out));
		push_at(lidE,mgCONNECT_POINT(trloop,mgCONNECT_POINT::coming_in));
	}
}

//Compute all the intersections of face boundaries and sl, which are all parameter
//representation.
MGLCisect_vector isectBoundary(
	const MGFace& face,//face to split
	const MGStraight& sl
){
	MGLCisect_vector pvec=face.loop(0)->isectLoop(sl);
	int n_inner=face.number_of_inner_boundaries();
	for(int k=0; k<n_inner; k++){
		const UniqueLoop& innerLoopk=face.loop(k);
		pvec.append(innerLoopk->isectLoop(sl));
	}
	return pvec;
}

//Compute the intersections of coordinate f (u or v value of
//the surface parameter representation) and loop lp.
void isect1DTl(
	const mgTL2parameter& tlparam,
	int coordinate,	// f coordinate kind, =0: u, =1:v
	double f,
	const MGLoop& lp,
	std::vector<MGLEPoint>& pvec, //intersections be appended.
	std::pair<MGLEPoint,bool>* limit
		//when limit->second==true, isect->first be employed as the isect.
){
	double error = coordinate ? tlparam.get_VError() : tlparam.get_UError(); 
	mgTolSetWCZero wczeroSet(error);//Set&save the error.
	
	double flimit=0.;
	bool employ_greater=true;
	int coordinate2=(coordinate+1)%2;
	if(limit){
		MGLEPoint uvLE=limit->first;
		pvec.push_back(uvLE);
		employ_greater=limit->second;
		MGPosition uv=uvLE.eval();
		flimit=uv[coordinate2];
		if(employ_greater)
			flimit+=error*10.;
		else
			flimit-=error*10.;
	}
	MGComplex::const_iterator ei=lp.pcell_begin(), ee=lp.pcell_end();
	for(; ei!=ee; ei++){
		const mgTL2Polyline* ipoly=TL2Polyline(ei);
		if(limit){
			const MGInterval& frange=ipoly->box()[coordinate2];
			if(employ_greater){
				if(frange<flimit)
					continue;
			}else{
				if(frange>flimit)
					continue;
			}
		}

		int nbd=ipoly->bdim();
		short idLast=short(nbd-1);
		MGCParam_list tlist=ipoly->isect_1D(f, coordinate);
		MGCParam_list::iterator ti=tlist.begin(), te=tlist.end();
		MGComplex::const_iterator ei2=ei;
		double ts=ipoly->param_s();

		for(;ti!=te; ti++){
			double  t=*ti;
			if(limit){
				MGPosition uv=ipoly->eval(t);
				double ft=uv[coordinate2];
				if(employ_greater){
					if(ft<flimit)
						continue;
				}else{
					if(ft>flimit)
						continue;
				}
			}
			//change the parameter value to indicate a vetex point of the polygon.
			short idt=short(t-ts+.5);
			MGLEPoint le(ei,ts+double(idt));
			if(idt==idLast)
				le=le.start_of_next_edge();
			size_t npvec=pvec.size();
			size_t i=0;
			for(; i<npvec; i++)
				if(pvec[i]==le)
					break;
			if(i==npvec)
				pvec.push_back(le);
		}
	}
}

//Compute minimum sine angle of crv and the edge at pv[] in their world representation.
//When sin angle is small enough crv is parallel to loop(boundary of the face) at
//the start or end point.
//Function's return value is min of sin angle at both end.
double angle_edge_crv(
	const MGSurface& srf,
	const MGCurve& crv,//(u,v) parameter line representation
	std::vector<MGLEPoint>::const_iterator pv
		//a point on the loop of the face at start and end of crv.
){
	MGSurfCurve crvOnsrf(srf,crv);//Build world curve rep.
	double t[2]={crvOnsrf.param_s(), crvOnsrf.param_e()};
	double angmin=3., ang1;
	for(int i=0; i<2; i++){//test for crv's start and end point.
		MGVector crvDir=crvOnsrf.eval(t[i],1);//crv's direction
		if(i)
			crvDir*=-1.;//When end point of crvOnsrf(to direct to inner of the loop)

		const MGLEPoint& pvi=pv[i];
		const MGEdge* edgi=pvi.edge();
		const mgTL2Polyline* edgeiCrv=TL2Polyline(edgi);
		if(pvi.is_Edge_start_point()){
			const MGVector& v0 = getSEDirection(pvi);
			ang1= MGCL::parallelism(crvDir,v0);

			MGVector Vpre= getSEDirection(pvi, false);
			Vpre*=-1.;//To direct away from the vertex.
			double ang2= MGCL::parallelism(crvDir,Vpre);
			if(ang2<ang1)
				ang1=ang2;
		}else{
			MGVector v0=edgeiCrv->direWorld_with_non_normalized(pvi.param());
			ang1= MGCL::parallelism(crvDir,v0);
		}
		if(ang1<angmin)
			angmin=ang1;
	}
	return angmin;
}

typedef std::pair<int, double> SplitParam;//defines split parameter of a loop.
//SplitParam.first=0: second is u-value, =1: second is v-value, =-1:undefined.
//SplitParam.second=parameter value.

//Get_split_network.
//Function's return value is true, if network is obtained, false if not.
//False means the parameter value t_to_divide was not adequate for the split.
bool get_split_bridge(
	int tryNum,//input how many try number this is, from 1.
	const SplitParam para[2],
	const MGFace& face,//Target face to split.
	const mgTL2parameter& tlparam,
	std::vector<UniqueTL2PlBridge>& bridges,
	std::pair<MGLEPoint,bool>* limit
){
	const MGSurface& sufaceOrg=tlparam.get_surface();
	int nloop=face.number_of_loops();
	double zeroAngle=LOOSE_ZERO_ANGLE/double(tryNum);

	bool obtained = false;
	for(int numtry=0; numtry<2; numtry++){

		if(para[numtry].first==-1)
			break;

		int kcod0= para[numtry].first;
		double t_to_divide=para[numtry].second;

		std::vector<MGLEPoint> pvec;
		for(int k=0; k<nloop; k++)
			isect1DTl(tlparam,kcod0,t_to_divide,*(face.loop(k)),pvec,limit);
		int nisecm1=(int)pvec.size()-1;
		if (nisecm1 > 1) {
			int kcod1 = (kcod0 + 1) % 2;
			auto pvComp= [kcod1](const MGLEPoint& le1, const MGLEPoint& le2) {
				return le1.eval()[kcod1] < le2.eval()[kcod1];
			};
			std::sort(pvec.begin(), pvec.end(), pvComp);
		}

		auto pveci=pvec.begin();
		for(int i=0;i<nisecm1; i++, pveci++){
			MGLEPoint& pvi=pvec[i]; MGLEPoint& pvip1=pvec[i+1];
			MGPosition uv0=pvi.eval(), uv1=pvip1.eval();
			MGPosition uvmid=(uv0+uv1)*.5;
			int in=face.in_range_with_on(uvmid);
			if(in!=2)
				continue;

			MGStraight sl(uv1,uv0);
			if(angle_edge_crv(sufaceOrg,sl,pveci)<=zeroAngle)//When nearly parallel
				break;

			const mgTL2Polyline* poly=TL2Polyline(pvi.edge());
			if(int(pvi.param()-poly->param_s()+.5)==poly->bdim()-2)
				if(pvi.start_of_next_edge()==pvip1)//If pvi and pvipi are the same.
				continue;

			MGLCisect_vector tempI=isectBoundary(face,sl);
			if(tempI.entries()>2)
				continue;

			bridges.emplace_back(new mgTL2PlBridge(pvi, pvip1));
			obtained=true;
		}
		if(obtained)
			break;
		else
			bridges.clear();

	}
	return obtained;
}

//Check the subdivision parameter direction.
//Function's return value is:
//true: dir direction is good for Vrtx,false: bad for Vrtx.
bool goodDirection(
	int dir, const MGLEPoint& Vrtx, const MGSurface& srf,
	int& isoIsParallelToEdge
		// = 1; Vpre is parallel to u=const paramline,
		// = 2; Vpre is parallel to v=const paramline, 
		// = 3; Vaft is parallel to u=const paramline,
		// = 4; Vaft is parallel to v=const paramline.
		// = 0: otherwise.
) {
	const MGVector& Vaft = getSEDirection(Vrtx);
	const MGVector& Vpre = getSEDirection(Vrtx, false);

	MGPosition uv = Vrtx.eval();
	MGVector Vu = srf.eval(uv, 1), Vv = srf.eval(uv, 0, 1);//
	double angupre = Vu.sangle(Vpre), angvpre = Vv.sangle(Vpre);
	double anguaft = Vu.sangle(Vaft), angvaft = Vv.sangle(Vaft);

	//Define is_u. We assume Vu is not parallel to Vv(not singular).
	isoIsParallelToEdge = 0;
	int dir2 = -1;
	if (angupre <= LOOSE_ZERO_ANGLE) {
		isoIsParallelToEdge = 1;//Vpre is parallel to Vu.
		if (angvaft > LOOSE_ZERO_ANGLE) 
			dir2 = 0;// must dir=0(cannot be dir=1).
	}else {
		if (angvpre <= LOOSE_ZERO_ANGLE) {
			isoIsParallelToEdge = 2;//Vpre is parallel to Vv.
			if (anguaft > LOOSE_ZERO_ANGLE)
				dir2 = 1;// must dir=1(cannot be dir=0).
		}else if (anguaft <= LOOSE_ZERO_ANGLE) {
			isoIsParallelToEdge = 3;//Vaft is parallel to Vu.
			dir2 = 0;// must dir=0(cannot be dir=1).
		}else if (angvaft <= LOOSE_ZERO_ANGLE) {
			isoIsParallelToEdge = 4;//Vaft is parallel to Vv.
			dir2 = 1;// must dir=1(cannot be dir=0).
		}
	}
	return dir==-1 || dir2==-1 || dir==dir2;
}

//Find the 1st concave vertex of the outer loop lp from the edge number startEdge.
//Function's return value is the edge num found. When not found, -1 is returned.
int findConcaveVertex(
	const MGLoop& lp,//Target loop.
	int startEdge//Start edge number to find.
){
	const MGSurface& srf = TL2Polyline(lp.last_edge())->TL2param().get_surface();

	int eNum = lp.number_of_edges();
	for (int j=startEdge; j<eNum; j++) {
		MGLEPoint Vrtx = MGLEPoint(lp,j, lp.edge(j)->param_s());
		const MGVector& Vpre = getSEDirection(Vrtx, false);
		const MGVector& Vaft = getSEDirection(Vrtx, true);
		MGVector N = srf.normal(Vrtx.eval());
		double concavityi = MGCL::concavity(Vpre, Vaft, N);
		if(isLooselyConcave(concavityi)){
			return j;//found, return the edge number j.
		}
	}
	return -1;//Not found.
}

//Find concave vertex of the outer loop lp.
//Function's return value is true if concave vertex is found.
bool concaveAtVertex(
	int dir,//=0: find u param, =1: v param, -1: either will do.
	const MGLoop& lp,//Target loop.
	std::pair<MGLEPoint, bool>& limitData,
		//.first=MGLEPoint of the concave vertex be returned,
		//       which is always the starting vertex of an edge.
		//.second= true if greater part of the para.second ' parameter value be employed.
	SplitParam& para //parameter to split face:
		//para.first= 0 if u value, 1 if v.
		//    .second= u or v parameter value to split.
) {
	MGLEPoint& Vrtx = limitData.first;
	const MGSurface& srf = TL2Polyline(lp.last_edge())->TL2param().get_surface();

	int j = 0;//Start edge number.
	int eNum = lp.number_of_edges();
	int isoIsParallelToEdge = 0;
	while (j < eNum) {
		int jFound = findConcaveVertex(lp, j);
		if (jFound == -1)
			return false;

		j = jFound;
		Vrtx = MGLEPoint(lp,j, lp.edge(j)->param_s());
		if (goodDirection(dir, Vrtx, srf, isoIsParallelToEdge))
			break;
		j++;
		if (j >= eNum)
			return false;
	}

	MGPosition uv=Vrtx.eval();//(u,v) surface parameter of vertex.
	MGVector N= srf.normal(uv);//Surface normal at vertex be stored.
	const MGVector& Vpre = getSEDirection(Vrtx, false);
	const MGVector& Vaft = getSEDirection(Vrtx, true);

	//Now concave vertex is found.
	//Define SplitParam& para. We assume Vu is not parallel to Vv(not singular).
	bool& employ_greater = limitData.second;
	bool is_u = (dir == -1 && isoIsParallelToEdge) ? isoIsParallelToEdge % 2 : dir != 1;
	double sangVpreVaft = Vpre.sangle(Vaft);
	MGVector Viso = is_u ? srf.eval(uv, 0, 1) : srf.eval(uv, 1);
	if (-SIN5DEGREE <= sangVpreVaft && sangVpreVaft <= SIN5DEGREE) {
		//When angle of concave is small.
		MGVector Visoaft = Viso * Vaft;
		employ_greater = Visoaft % N < 0.;
	}else {
		if (isoIsParallelToEdge) {
			if (isoIsParallelToEdge <= 2)
				employ_greater = Viso % Vaft < 0.;
			else
				employ_greater = Viso % Vpre > 0.;
		}else {
			double ang = Vaft.angle2pai(Viso, N);
			employ_greater = ang <= mgPAI;
		}
	}

	para.first = is_u ? 0 : 1;
	para.second = uv[para.first];
	return true;
}

double normalized_dif_angle(double angle, double base_angle) {
	double dif1 = angle - base_angle;
	double dif2 = dif1;
	if (dif1 < 0.)
		dif2 += mgDBLPAI;
	else
		dif1 -= mgDBLPAI;
	if (fabs(dif1) > fabs(dif2))
		dif1 = dif2;
	return dif1;
}

//Find edges that constitute concavity of the input angle from the outer loop lp.
//Let numEdge be function's return value, then numEdge continuous edges
//constitute the concavity if numEdge>0.
int find_concave_over_edges(
	double concavityAngle,//angle of concavity in radian.
	const MGLoop& lp,//Target face to split.
	const MGEdge*& edgeStart,//Starting edge that constitues the concavity will be returned.
		//numEdge continuous edges constitute the concavity.
	double& difAngle//detected concavity angle is output.
){
	double wzero2=MGTolerance::wc_zero_sqr();
	int numEdge=0;
	double angleEnd =difAngle=0.;

	int nedges=lp.number_of_edges();
	MGComplex::const_iterator i=lp.pcell_begin(), ie=lp.pcell_end();
	for(int j=0; j<nedges; j++,i++){
		const MGEdge* ei=edge_from_iterator(i);
		const mgTL2Polyline* crvi=TL2Polyline(ei);
		const MGVector& V0=crvi->direWorld_start();//
		if(V0%V0<=wzero2){
			numEdge=0;
			difAngle=0.;
			continue;
		}

		double angleStart=crvi->angle2Uline(0.);
		if(numEdge){
			double angle2NewEdge=normalized_dif_angle(angleStart, angleEnd);
			if(angle2NewEdge>=LOOSE_ZERO_ANGLE){
				if(difAngle < concavityAngle){
					return numEdge;
				}else{
					numEdge=0;
					difAngle=0.;
				}
			}else
				difAngle+=angle2NewEdge;
		}
		double angle1third=crvi->angle2Uline(0.3333);
		difAngle+=normalized_dif_angle(angle1third, angleStart);
		double angle2third=crvi->angle2Uline(0.6666);
		difAngle+=normalized_dif_angle(angle2third, angle1third);
		angleEnd=crvi->angle2Uline(1.);
		difAngle+=normalized_dif_angle(angleEnd, angle2third);
		if(difAngle>=LOOSE_ZERO_ANGLE){
			numEdge=0;
			difAngle=0.;
		}else{
			if(!numEdge)
				edgeStart=ei;
			numEdge++;
		}
	}
	if(difAngle >= concavityAngle)
		numEdge=0;
	return numEdge;
}

//Get the split u or v value and the parameter limit at uv of a concave edge.
void compute_split_limit(
	const MGLEPoint& Vrtx,//MGLEPoint data of uv.
	const MGSurface& srf,
	SplitParam& para,
	bool& employ_greater
){
	const mgTL2Polyline* crv = TL2Polyline(Vrtx.iterator());
	int idV = int(Vrtx.param() - crv->param_s() + .5);
	MGPosition uv = crv->uv(idV);
	MGVector Ve = crv->direWorld_with_non_normalized(Vrtx.param());

	bool is_u = para.first == 0;
	MGVector Vu=srf.eval(uv,1), Vv=srf.eval(uv,0,1);

	double uangl=Ve.sangle(Vu), vangl=Ve.sangle(Vv);
	if(uangl>vangl*2.)
		is_u=false;
	if(vangl>uangl*2.)
		is_u=true;
	const MGVector& Viso=is_u ? Vv:Vu;
	MGVector N=Vu*Vv;
	MGVector NVe=N*Ve;
	double direction=NVe%Viso;
	employ_greater=direction>0.;

	para.first= is_u ? 0:1;
	para.second= uv[para.first];
}

//Find the 1st vertex of an edge that exeeds the half of concave angle difAngle.
//The search is started from the 1st vertex of concaveEdge and is repeated
//to the numE edges.
int find_halfConcave_point(
	int numE,	//Number of edges to search from the edge concaveEdge.
	const MGEdge* concaveEdge,//1st edge that constitutes the concavity. numE edges constitute
				//the total concavity that have the concavity angle difAngle.
	double difAngle,//Total concavity angle.
	MGLEPoint& Vrtx//MGLEPoint data of uv.
){
	MGComplex::const_iterator i=concaveEdge->edge_iterator();
	double halfA=difAngle*.5;
	double  angleSummed=0., angle0=0.;
	for(int j=0; j<numE; j++,i++){
		const mgTL2Polyline* crv=TL2Polyline(i);
		int nv=crv->number_of_points();
		int k;
		for(k=0; k<nv; k++){
			double anglek=crv->angle2Uline_at_i(k);
			if(j || k){
				angleSummed+=normalized_dif_angle(anglek,angle0);
				if(angleSummed<=halfA)//When the concavity angle exceeds the half.
					break;
			}
			angle0=anglek;
		}
		if(angleSummed<=halfA){//When the concavity angle exceeds the half.
			int nv3=nv/3;
			if(k<=nv3 && j){
				k=0;
			}else if((nv-nv3)<=k && j<(numE-1)){
				k=0;
				crv=TL2Polyline(++i);
			}
			double t=crv->param_s()+double(k);
			Vrtx=MGLEPoint(i,t);
			return k;
		}
	}
	return -1;//This should not happen.
}

//Find at what parameter the loop be subdivided for the tessellation.
//Candidate params are stored in para.
void findBestMiddlePara(
	const MGLoop& lp,//Input loop.
	const MGBox& uvrange,
	int dir,	///<input indicates if t is u or v of the surface parameter (u,v).
		///=0: must divide u-direction, =1:must divide v-direction, =otherwise both are good.
		//Split is done around the middle of u or v parameter.
	int tryNum,//input how many try number this is, from 1.
	SplitParam para[2],//[0] is 1st candidate, [1] is 2nd.
	double* t_to_avoid//parameter to avoid. if input, employ other value.
){
	para[0].first=para[1].first=-1;//set undefined.
	int id= dir==1 ? 1:0;
	const MGInterval& trange=uvrange[id];
	double tspan=trange.length();
	double terror=tspan*NEAR_PARAM;
	double DeviMaxFromMid=tspan*MAX_DEVIATION_FROM_MIDDLE;//maximum deviation.
	double mid=trange.mid_point();//mid point of trange.

	double lengthMinFromMid, lengthMin2FromMid;//length from mid for para[0] and para[1].
	double& t0=para[0].second;//target parameter value stored(1st candidate).
	double& t1=para[1].second;//target parameter value stored(2nd candidate).

	MGComplex::const_iterator i=lp.pcell_begin(), ie=lp.pcell_end();
	for(; i!=ie; i++){
		const mgTL2Polyline* pli=TL2Polyline(i);
		MGPosition uvi=pli->uv(0);
		double ti=uvi[id];
		double leni=(ti<=mid ? mid-ti:ti-mid);
		if(leni> DeviMaxFromMid)//If ti is far away from mid.
			continue;
		if(t_to_avoid && *t_to_avoid-terror<=ti && ti<=*t_to_avoid+terror)
			continue;

		if(para[0].first==-1){
			para[0] = SplitParam(id, ti); lengthMinFromMid=leni;
		}else if(para[1].first==-1){
			if(leni<lengthMinFromMid){
				para[1]=para[0];lengthMin2FromMid=lengthMinFromMid;
				para[0]=SplitParam(id, ti); lengthMinFromMid=leni;
			}else{
				para[1] = SplitParam(id, ti); lengthMin2FromMid=leni;
			}
		}else{
			if(leni<lengthMinFromMid){
				para[1]=para[0]; lengthMin2FromMid=lengthMinFromMid;
				para[0] = SplitParam(id, ti); lengthMinFromMid=leni;
			}else if(leni<lengthMin2FromMid){
				para[1] = SplitParam(id, ti); lengthMin2FromMid = lengthMinFromMid;
			}
		}
	}
	if(para[0].first==-1){
		if(tryNum==1){
			para[0].first=id, t0=mid;
			para[1].first=id, t1=mid+ DeviMaxFromMid *.5;
		}else{
			para[0].first=id, t0=mid- DeviMaxFromMid *.5;
			para[1].first=id, t1=mid+ DeviMaxFromMid *.8;
		}
	}
}

/// <summary>
/// Find a split param para if lp has edge(s) that constitute concavity.
/// Function's return value is:
///    0 if concavite not found.
///    1 if found and para and limitData are set.
///	   2 if found and only para is set(limitData is invalid).
/// </summary>
int concaveOverEdges(
	const MGSurface& srf, const MGLoop& lp,
	SplitParam para[2], //parameter to split face is output:
		//para.first  = 0 if u value, 1 if v.
		//    .second = u or v parameter value to split.
	std::pair<MGLEPoint, bool>& limitData
) {
	MGLEPoint& Vrtx = limitData.first;
	bool& employ_greater = limitData.second;

	const MGEdge* concaveEdge;
	double difAngle;
	double concaveAngle = -mg30DEGREE;
	int numE = find_concave_over_edges(concaveAngle,lp, concaveEdge, difAngle);
	if (numE) {//If concave found over edges.
		if (numE == 1) {
			const MGBox& edgbx = concaveEdge->box();
			const MGBox& lpbx = lp.box();
			double r = MAX_DEVIATION_FROM_MIDDLE;
			double rpH = r+.5, Hmr=.5-r;//rpH+Hmr=1.
			for (int i = 0; i < 2; i++) {//i=0: u, =1: v value test.
				double t0 = lpbx[i].low_point(), t1 = lpbx[i].high_point();
				const MGInterval& range = edgbx[i];
				if (range < (rpH * t0 + Hmr *t1) || (Hmr *t0 + rpH * t1) < range) {
					findBestMiddlePara(lp, lpbx, i, 1, para, nullptr);
					return 2;
				}
			}
		}
		int idV = find_halfConcave_point(numE, concaveEdge, difAngle, Vrtx); assert(idV >= 0);
		compute_split_limit(Vrtx, srf, para[0], employ_greater);
		return 1;
	}
	return 0;
}

//Obtain bridges at the middle of uv range of lp to split face.
//Function's return value is true if obtained.
bool getBridgeAtMiddlePara(
	const MGFace& face, const mgTL2parameter& tlparam,	const MGLoop& lp,
	int dir, 
	std::vector<UniqueTL2PlBridge>& bridges
){
	const MGBox& uvrange = lp.box();
	SplitParam para[2];
	findBestMiddlePara(lp, uvrange, dir, 1, para, nullptr);
	if (get_split_bridge(1, para, face, tlparam, bridges, nullptr))
		return true;

	double tOld = para[0].second;
	findBestMiddlePara(lp, uvrange, dir, 2, para, &tOld);
	return get_split_bridge(2, para, face, tlparam, bridges, nullptr);
}

//Given bridges, split face. If bridges not given, obtained using getBridgeAtMiddlePara().
//Function's return value is true, if split was executed, false, if not.
bool splitByBridges(
	int dir, //=0: divide at u param, 1: at v param
	const MGFace& face,//Target face to split.
	const mgTL2parameter& tlparam,
	std::stack<UniqueFace>& faces//All of the splitted faces will be appended.
) {
	std::vector<UniqueTL2PlBridge> bridges;
	if (getBridgeAtMiddlePara(face, tlparam, *(face.loop(0)), dir, bridges)) {
		return split_by_bridges(face, bridges, faces);
	}else
		return false;
}

int divideDirSub(
	const MGBox& uvbox, const MGSurface& srf, double ratio
) {
	MGPosition uvmid = uvbox.mid();
	double ulen = (srf.eval(uvmid, 1,0)).len() * uvbox[0].length();
	double vlen = (srf.eval(uvmid, 0,1)).len() * uvbox[1].length();
	int dir = -1;
	if (ulen > vlen * ratio)
		dir = 0;
	else if (vlen > ulen * ratio)
		dir = 1;
	return dir;
}

//Decide which direction be employed to subdivide.
//=0: must divide u-direction, =1:must divide v-direction, =-1: both are good.
int divideDir(
	const MGFace& face,//Target face to split.
	const mgTL2parameter& tlparam,
	int loopID
){
	const MGSurface& srf = tlparam.get_surface();
	const MGLoop& loop = *face.loop(loopID);
	int dir = divideDirSub(loop.box(),srf, loopID ? SMALL_RATIO : LARGE_RATIO);
	if (loopID && dir == -1)//When inner loop and either direction was good.
		dir= divideDirSub((*face.loop(0)).box(), srf, LARGE_RATIO);
	return dir;
}

//Split face.
//Function's return value is true, if split was executed.
//false, if not.
bool splitTl(
	const MGFace& face,//Target face to split.
	const mgTL2parameter& tlparam,
	std::stack<UniqueFace>& faces//All of the splitted faces will be appended.
){
	//Define at what parameter face be split.
	const MGSurface& srf = tlparam.get_surface();
	int loopID=(face.number_of_loops() >1) ? 1:0;//=0 when no inner loops.

	int dir = divideDir(face, tlparam, loopID);
		//=0: must divide u-direction, =1:must divide v-direction, =-1: both are good.

	bool obtained = false;
	std::vector<UniqueTL2PlBridge> bridges;
	const MGLoop& lp = *(face.loop(loopID));

	if(loopID==0){//When face has no inner loops.
	// 1. get split bridge(s) at concave outer loop.

		std::pair<MGLEPoint, bool> limitData;
		std::pair<MGLEPoint, bool>* limit = &limitData;

		SplitParam para[2];//[0] is 1st, [1] is 2nd candidate to split.
		para[0].first = para[1].first = -1;//indicates undefined.
		if (concaveAtVertex(dir, lp, limitData, para[0]))
			//If concave found at a vertex.
			obtained = get_split_bridge(1, para, face, tlparam, bridges, limit);

		if (!obtained){
			int found = concaveOverEdges(srf, lp, para, limitData);
			if (found) {
				if(found==2)
					limit = nullptr;
				obtained = get_split_bridge(1, para, face, tlparam, bridges, limit);
			}
		}
	}

	if(!obtained){//When face has innner loops, or
				  //has only an outer loop and bridges failed to obtain,
	// 2. get split bridge at a middle point of u or v parameter.

		//1st try with dir.
		if (!getBridgeAtMiddlePara(face, tlparam, lp, dir, bridges)) {
			dir = dir != 1 ? 1 : 0;
			//2nd try by changing the direction.
			getBridgeAtMiddlePara(face, tlparam, lp, dir, bridges);
		}
	}

	bool splitted=false;
	if(bridges.size())
		splitted=split_by_bridges(face,bridges,faces);
	return splitted;
}
