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
#include "mg/Surface.h"
#include "mg/Tolerance.h"
#include "topo/CPointsVec.h"
#include "topo/Loop.h"
#include "topo/Edge.h"
#include "topo/LCisect_vector.h"
#include "topo/LLisect_vector.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//Debug Function
std::ostream& operator<< (std::ostream& out, const mgCONNECT_POINT& cpoint){
	out<<"mgCONNECT_POINT::m_trim_loop="<<*(cpoint.m_trim_loop);
	if(cpoint.m_outin==mgCONNECT_POINT::coming_in){
		out<<",coming_in";
	}else{
		out<<",going_out";
	}
	return out;
}

std::ostream& operator<< (std::ostream& out, const mgCPointsVec& cpvec){
	out<<"mgCPointsVec::"<<std::endl;
	int n=(int)cpvec.m_trim_loops.size();
	for(int i=0; i<n; i++)
		out<<" m_trim_loops["<<i<<"]:"<<*(cpvec.m_trim_loops[i]);
	int n1=(int)cpvec.m_cpoints_vec.size();
	for(int i=0; i<n1; i++)
		out<<" m_cpoints_vec["<<i<<"]:"<<cpvec.m_cpoints_vec[i];
	return out;
}

bool mgCONNECT_POINT::operator< (const mgCONNECT_POINT& cpoint2)const{
	MGLEPoint le1=lep();
	MGLEPoint le2=cpoint2.lep();
	if(le1 == le2){
		MGVector deri1=eval_deriv();
		if(m_outin==coming_in)
			deri1.negate();
		MGVector deri2=cpoint2.eval_deriv();
		if(cpoint2.m_outin==coming_in)
			deri2.negate();
		if(deri1.parallel(deri2)){
			if(deri1%deri2>0.){
				if(m_outin==going_out)
					return true;
				return false;
			}
		}
		MGVector X=le1.eval(1);
		double angle1=X.angle2pai(deri1,mgZ_UVEC);
		double angle2=X.angle2pai(deri2,mgZ_UVEC);
		return angle1>angle2;
	}
	return le1<le2;
}
bool mgCONNECT_POINT::operator== (const mgCONNECT_POINT& cpoint2)const{
	if(m_trim_loop!=cpoint2.m_trim_loop)
		return false;
	if(m_outin!=cpoint2.m_outin)
		return false;
	return true;
}

MGVector mgCONNECT_POINT::eval_deriv()const{
	if(m_outin==going_out){
		return m_trim_loop->eval_start_deriv();
	}else{
		return m_trim_loop->eval_end_deriv();
	}
}

//obtain the next point of this connect point.
//The next point may be in the same loop, or in a different loop.
//function's return value is the mgCONNECT_POINT of the next point.
mgCONNECT_POINT mgCONNECT_POINT::next_point(
	mgCPointsVec& cpointsVec
){
	MGTrimLoop* tloop=trim_loop();
	int loopid=tloop->end_loopid();
	mgCPoints& cpoints=cpointsVec[loopid];
	int id=cpoints.find_trim_loop(tloop,coming_in);
	if(is_coming_in()){
		//The next point of the end point of tloop.
		return cpoints.next_point(id);
	}else{
		//The next point of the start point of trim_loop
		//(that is the end point of tloop).
		return cpoints.m_cpoints[id];
	}
}

//obtain the previous point of this connect point.
//The previous point may be in the same loop, or in a different loop.
//function's return value is the id of cpoints of the previous point.
mgCONNECT_POINT mgCONNECT_POINT::prev_point(
	mgCPointsVec& cpointsVec
){
	MGTrimLoop* tloop=trim_loop();
	int loopid=tloop->start_loopid();
	mgCPoints& cpoints=cpointsVec[loopid];
	int id=cpoints.find_trim_loop(tloop,going_out);
	if(is_coming_in()){
		//The previous point of the end point of trim_loop
		//(that is the start point of tloop).
		return cpoints.m_cpoints[id];
	}else{
		//The previous point of the start point of tloop.
		return cpoints.prev_point(id);
	}
}

//extract the loop whose start point is first_point.
//Function's return value is:
//0: no loop was extracted since first_point was null.
//1: a loop was extracted into loop.
int mgCPointsVec::extract_loop(
	mgCONNECT_POINT& first_point,
	std::unique_ptr<MGLoop>& loop, //the loop will be output.
	std::deque<MGTrimLoop*>& tloops//used MGTrimLoop to extract loop will be output.
){
	if(first_point.is_null())
		return 0;

	mgCONNECT_POINT cpoint, npoint;
	if(first_point.is_coming_in()){
		cpoint=first_point;
		loop=std::unique_ptr<MGLoop>(new MGLoop);
	}else{		
		loop=first_point.loop_clone();
		tloops.push_back(first_point.trim_loop());
		cpoint=first_point.next_point(*this);
	}

	//construct loop to the forward direction.
	int ntloops=(int)m_trim_loops.size();
	for(int counter=0; counter<ntloops; counter++){//This counter is to avoid an infinite loop.
		//assert(cpoint.is_coming_in());
		if(cpoint.is_going_out())
			return 0;
		npoint=cpoint.next_point(*this);
		//assert(npoint.is_going_out());
		if(npoint.is_coming_in())
			return 0;

		if(npoint.is_null())
			return 0;//This must not happen.

		MGLEPoint ts=cpoint.lep(), te=npoint.lep();
		assert(ts.loop()==te.loop());
		if(ts!=te){
			std::unique_ptr<MGLoop> loop2 = trim_out_subloop(ts, te);
			loop->join(false,loop2);
		}
		if(npoint==first_point){
			loop->make_close();
			return 1;
		}else{
			std::unique_ptr<MGLoop> npLoop = npoint.loop_clone();
			loop->join(false, npLoop);
			tloops.push_back(npoint.trim_loop());
			cpoint=npoint.next_point(*this);
		}
		if(cpoint==first_point){
			loop->make_close();
			return 1;
		}
	}
	return 0;//This must not happen.
}

//get the id of m_cpoints that includes tloop
int mgCPoints::find_trim_loop(const MGTrimLoop* tloop, mgCONNECT_POINT::OUTIN out_in)const{
	int n=(int)m_cpoints.size();
	for(int i=0; i<n; i++){
		const mgCONNECT_POINT& cpi=m_cpoints[i];
		if(tloop==cpi.trim_loop())
			if(cpi.outin()==out_in)
				return i;
	}
	return -1;
}

//Debug Function
std::ostream& operator<< (std::ostream& out, const mgCPoints& cpoints){
	int n=(int)cpoints.m_cpoints.size();
	out<<"mgCPoints, size="<<n<<",";
	for(int i=0; i<n; i++){
		out<<std::endl<<i<<"::"<<cpoints.m_cpoints[i]<<std::endl;
	}
	return out;
}

//Extract MGTrimLoop & mgCPoints info, and build this mgCPointsVec data.
//While processing, a closed outer loop or closed inner loop found will be stored in
//the oloops or iloops.
//When m_uv is not null, oloops contains the only closed outer loop that contain m_uv,
//and return function's return code as 1. If the closed loop was not found,
//generate m_trim_loops(the vector of MGTrimLoop).
//When m_uv is null, put all the outer loops into oloops, and inner loops into iloops.
//Function's return value is
//  1: when m_uv is not null, a closed loop that includes m_uv was extracted into oloops.
//  0: when m_uv is not null, no closed outer loop that includes m_uv was found
//     and trim loop info were output into m_trim_loops and
//     mgCONNECT_POINT's into m_cpoints_vec. Closed inner loops may be output to iloops.
//When m_uv is null, function's return value is always 0.
int mgCPointsVec::extract_trim_loop_info(
	const MGLoop& network,	//original loops to trim.
	std::vector<UniqueLoop>& oloops,//When m_uv is not nul, the closed outerboundary loop
		//that includes m_uv is output and the return code ==1.
		//When m_uv is null, all the detected closed outerboundary loops will be output,
		//and the function's return value is 0.
	std::vector<UniqueLoop>& iloops//closed inner boundary loop that
							//includes uv is output when return code =0 or 1.
){
	int nedge=network.number_of_edges();
	if(!nedge)
		return 0;//Any loop was not found.

	enum SFLAG{
		NOT_YET=0,
		SAME=1,
		OPPOSITE=2,
		BOTH=3
	};
	std::vector<int> searched(nedge,NOT_YET);
		//NOT_YET: not searched yet for the network.edge(i),
		//SAME: searched for the same direction of the original edge,
		//OPPOSITE: searched for the opposite direction of the original
		//BOTH(=SAME+OPPOSITE): searched for the both directions.

	//Serach for the same direction.
	bool handled=true;
	while(handled){
	handled=false;
	for(int i=0; i<nedge; i++){
		int ei_flag=searched[i];
		if(ei_flag>=BOTH)//If searched in the both directions.
			continue;

		std::unique_ptr<MGLoop> loop2(new MGLoop);
		handled=true;
		const MGEdge* first_edge=network.edge(i);//save the 1st edge to judge loop to the 1st.
		bool first_was_same_direction;
			//True means the target loop has the same direction as the previous edge.

		if(ei_flag==OPPOSITE || ei_flag==NOT_YET){//If search was not done at the same
			first_was_same_direction=true;
		}else{//If search was not done at the opposite
			first_was_same_direction=false;
		}

		const MGEdge* edgei=first_edge;
		bool is_same_direction=first_was_same_direction;
		int edge_id=i;
		const MGEdge* next_edge=0;
		do{//Search is performed toward the loop direction.
			MGEdge* ei=edgei->clone();
			//MGEdge* ei = new MGEdge(*edgei);

			bool atEnd = true;
			SFLAG flag = SAME;
			if(!is_same_direction){			
				ei->negate();
				atEnd = false;//aft edge at start point of the edge.
				flag = OPPOSITE;
			}
			searched[edge_id] += flag;
			loop2->append(ei);
			next_edge = edgei->aft_edge(atEnd);

			if(next_edge && next_edge!=first_edge){
				if(is_same_direction){
					if(edgei->is_connected_and_same_direction(false,*next_edge))
						is_same_direction=true;
					else
						is_same_direction=false;
				}else{
					if(edgei->is_connected_and_same_direction(true,*next_edge))
						is_same_direction=false;
					else
						is_same_direction=true;
				}
				edgei=next_edge;
				edge_id=network.edge_num(edgei);
			}
		}while(next_edge!=0 && next_edge!=first_edge);

		if(next_edge==first_edge){
			loop2->make_close();
			if(no_uv()){
				if(loop2->is_outer_boundary()){
					oloops.emplace_back(loop2.release());
				}else
					iloops.emplace_back(loop2.release());
			}else if(loop2->inside(m_uv)){
				if(loop2->is_outer_boundary()){
					oloops.emplace_back(loop2.release());
					return 1;
				}else
					iloops.emplace_back(loop2.release());
			}
		}else{
			assert(!next_edge);
			//When edge reached to the boundary of the original face,
			//opposite direction searching must be done.
			if(first_was_same_direction){
				next_edge=first_edge->pre_edge();
			}else{
				next_edge=first_edge->pre_edge(false);//pre edge at the end point of 1st edge.
			}
			edgei=first_edge;
			is_same_direction=first_was_same_direction;
			int iEdge=0;//Counter to prevent non-stop loop.
			while(next_edge && iEdge++<nedge){//Search is performed against the loop direction.
				if(is_same_direction){
					if(edgei->is_connected_and_same_direction(true,*next_edge))
						is_same_direction=true;
					else
						is_same_direction=false;
				}else{
					if(edgei->is_connected_and_same_direction(false,*next_edge))
						is_same_direction=false;
					else
						is_same_direction=true;
				}
				edgei=next_edge;
				edge_id=network.edge_num(edgei);
				MGEdge* ei = new MGEdge(*edgei);
				bool atEnd = true;
				SFLAG flag = SAME;
				if(!is_same_direction){
					flag=OPPOSITE;
					ei->negate();
					atEnd =false;//pre edge at the end point of the edge.
				}
				searched[edge_id] += flag;
				loop2->prepend(ei);
				next_edge = edgei->pre_edge(atEnd);
			}
			MGPosition uvS=loop2->start_point();
			MGPosition uvE=loop2->end_point();
			int inS=m_face.in_range_with_on(uvS),
				inE=m_face.in_range_with_on(uvE);
			if(inS>=0 && inS!=4)
				continue;
			if(inE>=0 && inE!=4)
				continue;
			int lidS=0; if(inS<0) lidS=-inS;
			int lidE=0; if(inE<0) lidE=-inE;
			double d;
			MGLEPoint leS=m_face.loop(int(lidS))->closest(uvS,d),
					  leE=m_face.loop(int(lidE))->closest(uvE,d);
			MGLoop* lp=loop2.release();
			MGTrimLoop* trloop=new MGTrimLoop(lp,inS,leS,inE,leE);
			m_trim_loops.emplace_back(trloop);
			push_at(lidS,mgCONNECT_POINT(trloop,mgCONNECT_POINT::going_out));
			push_at(lidE,mgCONNECT_POINT(trloop,mgCONNECT_POINT::coming_in));
		}

	}
	}

	int nboundaries=size();
	for(int i=0; i<nboundaries; i++){
		mgCPoints& lveci=(*this)[i];
		lveci.sort();
	}

	return 0;
}

//Extract a boundary out of networks that includes m_uv.
//m_uv is a parameter of m_face.
//Function's return value is
//  0: no loops were extracted,
//  2: Loops to trim face were output into used_tloops.
//     If trimming is performed by used_tloops[i] for i=0,...,n-1 one by one,
//     face will be trimmed by the part of the face that includes uv.
int mgCPointsVec::extract_uv_boundary(
	std::unique_ptr<MGLoop>& loop,	//the loop that includes m_uv will be output
		//when return code=2.
	std::deque<MGTrimLoop*>& used_tloops
){
	int nboundaries=size();
	for(int i=0; i<nboundaries; i++){
		mgCPoints& lveci=(*this)[i];
		int nlv=lveci.size();
		assert((nlv%2)==0);//nlv must be even(pair of coming_in and going_out).
		if(!nlv)
			continue;

		for(int j=0; j<nlv; j++){
			mgCONNECT_POINT& fpoint=lveci.m_cpoints[j];
			int extract_inf=extract_loop(fpoint,loop,used_tloops);
			if(!extract_inf)
				continue;
			int ntloop=(int)used_tloops.size(), k;
			if(!loop->inside(m_uv)){
				for(k=0; k<ntloop; k++)
					used_tloops[k]->set_null();
				continue;//if loop does not includes uv.
			}
			return 2;
		}
	}
	return 0;
}
