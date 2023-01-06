/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include <iomanip>
#include "mg/Interval.h"
#include "mg/Box.h"
#include "mg/Position.h"
#include "mg/Curve.h"
#include "mg/SurfCurve.h"
#include "mg/Surface.h"
#include "mg/Tolerance.h"
#include "mg/TrimmedCurve.h"
#include "topo/CPointsVec.h"
#include "topo/LEPoint.h"
#include "topo/Loop.h"
#include "topo/PVertex.h"
#include "topo/Edge.h"
#include "topo/Face.h"
#include "topo/LCisect_vector.h"
#include "topo/LLisect_vector.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGLoop Class.
//MGLoop is a boundary of a face, a boundary of 2D manifold cell.
//MGLoop always accepts parameter space curve and world space curve
//of a boundary curve, and constructs a boundary of a face from the
//two types of curves.

//Input curve direction indicates which part of the face will be target
//part after trimed by the boundary. In 2D space (u,v) of the parameter
//space, LEFT side of the parameter curve along the curve's direction
//is the target part of face.

struct LCComp{
	bool operator()(const MGLCisect& lc1, const MGLCisect& lc2){return lc1.t()<lc2.t();};
};

//Merge with param_curve as a network loop.
//Function's return value is:
//true: merge is done and param_curve is processed into this loop as a network.
//false: merge is not performed since no intersection with this loop were found.
//When false is returned, this loop does not have the input param_curve information
//as an edge information.
//This loop must not be empty loop, and the kind is always changed to NETWORK.
bool MGLoop::merge_network(const MGCurve& param_curve, double error){
	m_kind=NETWORK;

	MGLCisect_vector lcis=isect_with_endpoints(param_curve);
	//isect of this loop and param_curve.
	int n=lcis.entries();
	if(n==0)
		return false;

	if(n>1)
		std::sort(lcis.begin(),lcis.end(), LCComp());

	MGPVertex* Vnew=0;

	const MGLCisect& lc_start=lcis[0];
	double tnow=param_curve.param_s(), tnext=lc_start.t()
		, tspan=param_curve.param_span();
	const MGLEPoint& le0=lc_start.lp();
	makeVertexWithLcis(le0,lcis,Vnew);

	MGEdge* edgeStart=0;
	if(!MGREqual_base(tnext,tnow,tspan)){
		edgeStart=new MGEdge(param_curve, MGInterval(tnow,tnext));
		edgeStart->connect_at_end(*Vnew);
		tnow=tnext;
	}

	MGEdge* edgeNewi=edgeStart;
	int nm1=n-1;
	for(int i=0; i<nm1; i++){
		//Process one edge of param_curve(the start and end point connect).
		const MGLCisect& lc_aft=lcis[i+1];
		tnext=lc_aft.t();
		edgeNewi=new MGEdge(param_curve, MGInterval(tnow,tnext));
		edgeNewi->connect_at_start(*Vnew);

		const MGLEPoint& leip1=lc_aft.lp();
		makeVertexWithLcis(leip1,lcis,Vnew);
		edgeNewi->connect_at_end(*Vnew);

		if(!edgeStart)
			edgeStart=edgeNewi;
		tnow=tnext;
	}

	MGEdge* edgeEnd=edgeNewi;
	const MGLCisect& lc_end=lcis[nm1];
	tnext=param_curve.param_e();
	if(!MGREqual_base(tnext,tnow,tspan)){
		edgeEnd=new MGEdge(param_curve, MGInterval(tnow,tnext));
		edgeEnd->connect_at_start(*Vnew);
	}
	if(edgeStart && edgeEnd){
		mgTolSetWCZero wczeroSet(error);//Set&save the error.
		if(edgeStart->start_point()==edgeEnd->end_point()){
			edgeStart->connect_at_id(0, edgeEnd,1);
		}
	}
	return true;
}

//Remove pendent edge.
bool MGLoop::remove_pendent_edge(const MGFSurface& face){
	bool deleted=false;
	iterator i=pcell_begin(), iend=pcell_end();
	for(;i!=iend;){
		iterator iSave=i;
		MGEdge* edge= dynamic_cast<MGEdge*>(i->get());i++;
		MGEdge* pre=edge->pre_edge();
		MGEdge* aft=edge->aft_edge();
		if(pre&&aft)
			continue;
		if(!pre){
			MGPosition uv=edge->start_point();
			int in=face.in_range_with_on(uv);
			if(0<=in && in<=2){
				deleted=true;
				erase(iSave);
				continue;
			}
		}
		if(!aft){
			MGPosition uv=edge->end_point();
			int in=face.in_range_with_on(uv);
			if(0<=in && in<=2){
				erase(iSave);
				deleted=true;
				continue;
			}
		}
	}
	return deleted;
}

///Remove garbage edge.
///Error is the parameter error allowed of the edge length to judge removal.
void MGLoop::remove_garbage_edge(double error){
	iterator i=pcell_begin(), iend=pcell_end();
	for(;i!=iend;){
		iterator iSave=i;//Save
		MGEdge* edge= dynamic_cast<MGEdge*>(i->get());i++;
		MGTrimmedCurve crv=edge->trimmed_curve();
		double edgeLen=crv.length(crv.param_s(), crv.param_e());
		if(edgeLen>error)
			continue;

		MGEdge* pre=edge->pre_edge();
		MGEdge* aft=edge->aft_edge();
		if(pre&&aft){
			if(edge->vertex_end()->binder_vertex()!=
				edge->vertex_start()->binder_vertex())
				continue;
		}
		erase(iSave);
	}
}

//Build a boundary loop by tracking tloops backward.
//Function's return value is:
//0: no loops were extracted
//1: a boundary loop was extracted from the first non-null tloops
int build_boundary_loop(
	std::deque<MGTrimLoop*>& tloops,//used MGTrimLoop to extract loop is input.
	std::unique_ptr<MGLoop>& loop,//the loop will be output.
	std::vector<bool>& used_loops//vector of length face.number_of_boundaries().
		//used loops to build loop will be set true.
){
	if(tloops.empty())
		return 0;
	MGTrimLoop* tloop=nullptr;
	while(!tloops.empty()){
		tloop=tloops.back(); tloops.pop_back();
		if(!tloop->is_null())
			break;
	}
	if(!tloop || tloop->is_null())
		return 0;

	int nboundaries=(int)used_loops.size();
	for(int j=0; j<nboundaries; j++)
		used_loops[j]=false;
	tloop->set_used_loop_flag(used_loops);

	int first_point_loopid=tloop->end_loopid();
	int next_point_loopid=tloop->start_loopid();
	MGLEPoint first_point=tloop->end_lep();
	MGLEPoint ts=tloop->start_lep();
	loop=std::unique_ptr<MGLoop>(tloop->release_loop());
	loop->negate();
	while(next_point_loopid!=first_point_loopid && !tloops.empty()){
		//while next_point_loopid!=first_point_loopid, tloops must not be empty.
		tloop=tloops.back(); tloops.pop_back();
		MGLEPoint te=tloop->end_lep();
		assert(ts.loop()==te.loop());
		std::unique_ptr<MGLoop> loop2=trim_out_subloop(ts,te);
		loop->join(false,loop2);

		next_point_loopid=tloop->start_loopid();used_loops[next_point_loopid]=true;
		ts=tloop->start_lep();
		MGLoop* loop3=tloop->release_loop();
		loop3->negate();
		loop->join(false,loop3);
	}
	std::unique_ptr<MGLoop> loop4=trim_out_subloop(ts,first_point);
	loop->join(false,loop4);
	loop->make_close();

	return 1;
}

//Trim the face giving networks loops. Trimming is done by removing the smallest
//closed area out of networks that includes the parameter position uv(u,v).
void MGFace::trim(
	const std::vector<UniqueLoop>& networks,///<network to trim the face.
	const MGPosition& uv,
	std::vector<UniqueFace>& faces//Result trimmed face(s) will be appended.
)const{
	mgCPointsVec cpointsVec(*this,uv);
		//Found loops that are connected to
		//an inner loop of loop id loop_id will be stored
		//in this cpointsVec[loop_id].

	std::vector<UniqueLoop> iloops;
	int nnetworks=(int)networks.size();
	for(int i=nnetworks-1; i>=0; i--){
		std::vector<UniqueLoop> oloops;
		int closed_loop_obtained
			=cpointsVec.extract_trim_loop_info(*(networks[i]),oloops,iloops);
		if(closed_loop_obtained){
			oloops[0]->negate();
			MGFace* face2=new MGFace(*this);
			face2->add_boundary(oloops[0].release());
			faces.emplace_back(face2);
			return;
		}
	}

	std::unique_ptr<MGLoop> loop;
	std::deque<MGTrimLoop*> used_tloops;
	int eout=cpointsVec.extract_uv_boundary(loop,used_tloops);
	if(!eout){
		int niloops=(int)iloops.size();
		for(int i=0; i<niloops; i++){
			MGLoop* iloop=iloops[i].release();
			iloop->negate();
			MGFace* face2=new MGFace(*this);
			face2->add_boundary(iloop);
			faces.emplace_back(face2);
		}
		return;
	}

	int nboundaries=number_of_boundaries();
	std::vector<bool> used_loops(nboundaries);
	std::unique_ptr<MGLoop> new_loop;
	while(build_boundary_loop(used_tloops,new_loop,used_loops)){
		MGFace* face2=new MGFace(*this);
		for(int j=nboundaries-1; j>=0; j--){
			if(used_loops[j])
				face2->erase_boundary(j);
		}
		face2->add_boundary(new_loop.release());
		faces.emplace_back(face2);
	}
}

void set_used_loop(
	const std::deque<MGTrimLoop*>& used_tloops,
	std::vector<bool> used_loops
){
	int ntloop=(int)used_tloops.size();
	for(int i=0; i<ntloop; i++){
		used_tloops[i]->set_used_loop_flag(used_loops);
	}
}

//Extract sub face that is bounded by network loops.
//Extracted sub face is the smallest closed part of this face bounded by
//the networks that includes the parameter position uv(u,v).
void MGFace::extract_sub_face(
	const std::vector<UniqueLoop>& networks,///<(u,v) representation networks.
	const MGPosition& uv,
	std::unique_ptr<MGFace>& face//Result extracted face will be output.
)const{
	face.reset(0);
	mgCPointsVec cpointsVec(*this,uv);
		//Found loops that are connected to
		//an inner loop of loop id loop_id will be stored
		//in this cpointsVec[loop_id].

	std::vector<UniqueLoop> iloops;
	int nnetworks=(int)networks.size();
	for(int i=nnetworks-1; i>=0; i--){
		std::vector<UniqueLoop> oloops;
		int closed_loop_obtained
			=cpointsVec.extract_trim_loop_info(*(networks[i]),oloops,iloops);
		if(closed_loop_obtained){
			face=std::unique_ptr<MGFace>(new MGFace(*this));
			face->add_boundary(oloops[0].release());
			return;
		}
	}

	std::unique_ptr<MGLoop> loop;
	std::deque<MGTrimLoop*> used_tloops;
	int eout=cpointsVec.extract_uv_boundary(loop,used_tloops);
	if(eout){
		face=std::unique_ptr<MGFace>(new MGFace(*this));
		int nboundaries=number_of_boundaries();
		std::vector<bool> used_loops(nboundaries);
		set_used_loop(used_tloops,used_loops);
		for(int j=nboundaries-1; j>=0; j--){
			if(used_loops[j])
				face->erase_boundary(j);
		}
		face->add_boundary(loop.release());
	}else{
		int niloops=(int)iloops.size();
		if(!niloops)
			return;
		face=std::unique_ptr<MGFace>(new MGFace(*this));
		for(int i=0; i<niloops; i++){
			face->add_boundary(iloops[i].release());
		}
	}
}

//Split the face giving networks loops. Splitting is done by finding the smallest
//closed areas out of networks.
void MGFace::split(
	const std::vector<UniqueLoop>& networks,///<Network to split.
	std::vector<UniqueFace>& faces//Result trimmed face(s) will be appended.
)const{
	double rzero=MGTolerance::rc_zero();
	const MGBox& bxuv=box_param();
	double minimum_loop_area=bxuv[0].length();
	minimum_loop_area+=bxuv[1].length().value();
	minimum_loop_area*=rzero;

	bool handled=false;
	mgCPointsVec cpointsVec(*this);
		//Found loops that are connected to
		//an inner loop of loop id loop_id will be stored
		//in this cpointsVec[loop_id].

	std::vector<UniqueLoop> iloops,oloops;
	int nnetworks=(int)networks.size();
	for(int ii=nnetworks-1; ii>=0; ii--)
		cpointsVec.extract_trim_loop_info(*(networks[ii]),oloops,iloops);

	int nolp=(int)oloops.size();
	for(int i=0; i<nolp; i++){
		MGFace* nface=new MGFace(*this);
		nface->add_boundary(oloops[i].release());
		faces.emplace_back(nface);
		handled=true;
	}

	std::unique_ptr<MGFace> face2(new MGFace(*this));
	int nilp=(int)iloops.size();
	for(int i=0; i<nilp; i++){
		face2->add_boundary(iloops[i].release());
		handled=true;
	}

	int nboundaries=number_of_boundaries();
	for(int i=0; i<nboundaries; i++){
		mgCPoints& lveci=cpointsVec[i];
		int nlv=lveci.size();
		assert((nlv%2)==0);//nlv must be even(pair of coming_in and going_out).
		if(!nlv)
			continue;

		for(int j=0; j<nlv; j++){
			std::unique_ptr<MGLoop> loop;
			std::deque<MGTrimLoop*> tloops;
			mgCONNECT_POINT& fpoint=lveci.m_cpoints[j];
			int extract_inf=cpointsVec.extract_loop(fpoint,loop,tloops);
			if(!extract_inf)
				continue;
		
			std::vector<bool> used_loops(nboundaries,false);
			int ntloop=(int)tloops.size();
			for(int k=0; k<ntloop; k++){
				tloops[k]->set_used_loop_flag(used_loops);
				tloops[k]->set_null();
			}

			loop->get_kind();
			if(loop->m_area<=minimum_loop_area && loop->m_area >= -minimum_loop_area)
				continue;

			MGFace* nface2=face2->clone();
			for(int j=nboundaries-1; j>=0; j--){
				if(used_loops[j])
					nface2->erase_boundary(j);
			}
			nface2->add_boundary(loop.release());
			faces.emplace_back(nface2);
			handled=true;
		}
	}

	if(!handled)
		faces.emplace_back(face2.release());
}

//Get the clone of ts.loop()(==te.loop()), and trim the clone from ts to te.
//Function's return value is the trimmed loop.
std::unique_ptr<MGLoop> trim_out_subloop(const MGLEPoint& ts, const MGLEPoint& te){
	//assert(!ts.is_end_point() && !te.is_start_point());
	assert(ts.loop()==te.loop());
	const MGLoop* loop1=ts.loop();
	std::unique_ptr<MGLoop> loop2(new MGLoop(*loop1));
	MGLoop& loop2ref=*loop2;
	int ens=ts.edge_num();
	int ene=te.edge_num();
	MGLEPoint ts2(*loop1,ts,loop2ref);
	MGLEPoint te2(*loop1,te,loop2ref);
	loop2ref.trim(ts2,te2);
	return loop2;
}

const int ORDER_NEW=4;
void edgeChangeParam(
	double error,
	int parameter_normalization,
	MGEdge& edge
){
	MGCurve* bcrv=edge.base_curve();
	if(parameter_normalization
		|| edge.param_s()!=bcrv->param_s() || edge.param_e()!=bcrv->param_e()){
		mgTolSetWCLineZero wczeroSet(error,error);//Set&save the error.

		MGTrimmedCurve tci=edge.trimmed_curve();
		std::unique_ptr<MGLBRep> lbNew(new MGLBRep);
		tci.approximate_as_LBRep(*lbNew,ORDER_NEW,parameter_normalization);
		double ts=lbNew->param_s(), te=lbNew->param_e();
		edge.resetBinder();
		edge.set_extent(std::move(lbNew));//Replace the curve rep.
		edge.set_only_param_range(ts,te);//Set start and end parame of eip1.
	}
}

///Join adjacent two C1 continuous edges to one edge.
///When the number of edges is equal to or less than 4, join is not executed.
void MGLoop::join_C1_edges(
	int parameter_normalization
		///Specify how rebuilt edes's parameterization be done.
		//=0: no parameter normalization.
		//=1: normalize to range=(0., 1.);
		//=2: normalize to make the average length of the 1st derivative 
		//    is as equal to 1. as possible.
){
	MGFace* f=face();
	if(!f)
		return;

	int nedges=number_of_edges();
	if(nedges<=4)
		return;

	if(parameter_normalization && parameter_normalization!=1)
		parameter_normalization=2;

	double err=error();//Error for edge representation.
	MGSurface* srf=f->surface();

	bool EiJoined=false;
	MGEdge* ei=first_edge();
	MGTrimmedCurve tci=ei->trimmed_curve();
	MGSurfCurve worldCrvi(*srf,tci);
	MGVector taniE=worldCrvi.eval(worldCrvi.param_e(),1);
	int nem1=nedges-1;
	for(int j=0; j<nedges; j++){
		MGEdge* eip1=ei->aft_edge();
		if(!eip1)
			break;

		MGTrimmedCurve tcip1=eip1->trimmed_curve();
		MGSurfCurve worldCrvip1(*srf,tcip1);
		MGVector tanip1S=worldCrvip1.eval(worldCrvip1.param_s(),1);
		if(j<nem1 && taniE.parallel(tanip1S) && taniE%tanip1S>0.){
			//Join ei and eip1 to one edge.
			//ei will be deleted and will be included in eip1.

			std::unique_ptr<MGLBRep> lbip1New(new MGLBRep);
			if(!EiJoined)//If ei is not already processed to join to 1 edge.
				edgeChangeParam(err,parameter_normalization,*ei);
			lbip1New.reset(static_cast<MGLBRep*>(ei->free_extent()));

			MGLBRep lbip1;// will be (u,v) parameter rep of eip1;
		mgTolSetWCLineZero wczeroSet(err,err);//Set&save the error.
			tcip1.approximate_as_LBRep(lbip1,ORDER_NEW);
			lbip1New->connect(1,2,lbip1);
			lbip1New->approximate_as_LBRep(*lbip1New,ORDER_NEW,parameter_normalization);
		wczeroSet.restore();

			///delete redundunt knots by approximate_as_LBRep.
			eip1->resetBinder();
			double tsnew=lbip1New->param_s(), tenew=lbip1New->param_e();
			eip1->set_extent(std::move(lbip1New));//Replace the curve rep.
			eip1->set_only_param_range(tsnew,tenew);//Set start and end parame of eip1.
			MGEdge* eim1=ei->pre_edge();
			delete ei;
			if(eim1){
				if(eim1==eip1){
					make_close();
					break;
				}else
					eim1->join(false,eip1);
			}
			tcip1=eip1->trimmed_curve();
			worldCrvip1=MGSurfCurve(*srf,tcip1);
			EiJoined=true;
		}else{
			edgeChangeParam(err,parameter_normalization,*ei);
			EiJoined=false;
		}

		taniE=worldCrvip1.eval(worldCrvip1.param_e(),1);
		ei=eip1;
		tci=tcip1;
	}
}
