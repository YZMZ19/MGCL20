/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Curve.h"
#include "mg/Tolerance.h"
#include "topo/Edge.h"
#include "topo/Loop.h"
#include "topo/LLisect_vector.h"
#include "topo/LCisect_vector.h"

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

//Input curves direction indicate which part of the face will be target
//part after trimed by the boundary. In 2D space (u,v) of the parameter
//space, LEFT side of the parameter curve along the curve's direction
//is the target part of face.

//Compute common range of two loops, this and loop2.
//Function's return value is number of ranges obtained.
//In variable pranges1,2 and branges1,2, common ranges are output as:
//  Let n be output of the function, then pranges1.size()=pranges2.size()=2*n.
//  pranges1[2*i+0] and pranges1[2*i+1] are parameter values of this loop, and
//  pranges2[2*i+0] and pranges2[2*i+1] are parameter values of loop2.
//  Although pranges1[2*i+0] < pranges1[2*i+1] always holds, 
//  pranges2[2*i+0] < pranges2[2*i+1], or pranges2[2*i+0] > pranges2[2*i+1].
//  Let f1() be this loop, and f2() be loop2, then
//  f1(pranges1[j]) and f2(pranges2[j]) represent the same point in star
//  Face world for 0<=j<n*2 .
//  branges1,2 are parameter values of curves_world() of pranges1,2 each.
//  On return from common, every edge of pranges1,2 has curves_world.
//This and loop2 must have each star faces.
int MGLoop::common(
	const MGLoop& loop2,
	std::vector<MGLEPoint>& pranges1,
	std::vector<double>& branges1,
	std::vector<MGLEPoint>& pranges2,
	std::vector<double>& branges2
)const{

	assert(surface());		//Surface of the face must be attached.
	assert(loop2.surface());//Surface of the face must be attached.

	pranges1.clear(); branges1.clear();
	pranges2.clear(); branges2.clear();

	std::vector<UniqueCurve> curves1=curves_world();
	std::vector<UniqueCurve> curves2=loop2.curves_world();

	double lzero=MGTolerance::line_zero();//Save the tolerance.
	double lerror=lzero*.5;//The tolerance is line_zero().
	MGBox box1, box2;
	std::vector<MGBox> boxes1, boxes2;

	auto crv1=curves1.begin(); auto crv1end=curves1.end();
	for(; crv1!=crv1end; crv1++){
		MGBox boxi=(*crv1)->box(); boxi.expand(lerror);
		box1 |=boxi;
		boxes1.push_back(boxi);
	}
	auto crv2=curves2.begin(); auto crv2end=curves2.end();
	for(; crv2!=crv2end; crv2++){
		MGBox boxi=(*crv2)->box(); boxi.expand(lerror);
		box2 |=boxi;
		boxes2.push_back(boxi);
	}

	if((box1&box2).empty()) return 0;//If no overrideness, return.

	double lzero3=lzero*3.;//wider tolerance is necessary for common below.
	mgTolSetLineZero lineZeroSet(lzero3);
	int i,j, m=(int)boxes1.size(), n=(int)boxes2.size();
	MGComplex::const_iterator ip=pcell_begin(), jp;
	for(i=0; i<m; i++, ip++){
	
	jp=loop2.pcell_begin();
	for(j=0; j<n; j++, jp++){
		if(!(boxes1[i]&boxes2[j]).empty()){
			std::vector<double> rangesij;
			int ncom=curves1[i]->common(*(curves2[j]),rangesij);
			if(ncom>=1){
				const MGEdge* pei=edge_from_iterator(ip);	//Parameter edge of i.
				const MGEdge* pej=edge_from_iterator(jp);	//Parameter edge of j.
				int inc, now;	//Used to sort pranges in assending order.
				if(pei->equal_direction_to_binder()){
					inc=1; now=0;
				}else{
					inc=-1; now=4*ncom-3;
				}
				for(int k=0; k<ncom; k++){
					double bt10=rangesij[now], bt11=rangesij[now+inc];
					pranges1.push_back(MGLEPoint(ip,pei->param_pcell(bt10)));
					pranges1.push_back(MGLEPoint(ip,pei->param_pcell(bt11)));
					branges1.push_back(bt10); branges1.push_back(bt11);
					double bt20=rangesij[now+2], bt21=rangesij[now+inc+2];
					pranges2.push_back(MGLEPoint(jp,pej->param_pcell(bt20)));
					pranges2.push_back(MGLEPoint(jp,pej->param_pcell(bt21)));
					branges2.push_back(bt20); branges2.push_back(bt21);
					now+=(4*inc);
				}
			}
		}
	}

	}
	int num=(int)(pranges1.size()/2);
	return num;
}

//Make this loop as closed.
//This loop's 1st edge's start point must be the same as the last edge's end point.
//However, this is not tested in make_close.
void MGLoop::make_close(){
	first_edge()->connect_at_id(0, last_edge(),1);
	assert(closed());
}

//Make a vertex at lp and subdivide the edge into two edges.
//Returned is true if subdivision is done and false if no subdivision
//is done since lp was one of existed vertex.
//When function's return value is true, pre is always the same edge as
//lp's edge and the iterator is unchanged.
bool MGLoop::make_vertex(
	const MGLEPoint& lp,//point to subdivide of this loop.
	MGEdge*& pre,	//pre and aft-edge of the lp will be output
	MGEdge*& aft,	//after make_vertex's execution, may be null.
	MGPVertex** pvertex ///MGPVertex of the made vertex which is the boundary of lp.
){
	bool made=false;
	double err=error()*4.;
	MGEdge* e=lp.edge_to_update();
	double ps=e->param_s(), pe=e->param_e(), t=lp.param();
	int vetexID=1;
	if(t-ps<=err){
		aft=e; pre=e->pre_edge();
		vetexID=0;
	}else if(pe-t<=err){
		pre=e; aft=e->aft_edge();
	}else{
		int aft2ID;
		MGEdge* aft2=e->aft_edge(true,&aft2ID);	//Save the current after edge.
		pre=e; aft=new MGEdge(*pre);
		pre->trim_end(t);
		//Add the new cell into the loop complex in the order of the edges.
		iterator pcellp=pcellIterator(lp.edge_num()); pcellp++;
		add_cell(aft, pcellp);
		aft->trim_start(t);//Trim is done here after add_cell
						//since binder of the edge can be also trimmed.
		pre->join(false,aft);
		if(aft2)
			aft->connect_at_id(1,aft2,aft2ID);
		made=true;
	}
	if(pvertex)
		*pvertex=e->vertex(vetexID).get();
	return made;
}

///Make a vertex at lp and subdivide the edge into two edges.
///Returned is true if subdivision is done and false if no subdivision
///is done since lp was one of existed vertex.
///When function's return value is true, pre is always the same edge as
///lp's edge and the iterator is unchanged.
bool MGLoop::makeVertexWithLcis(
	const MGLEPoint& lp,///<point to subdivide of this loop.
	MGLCisect_vector& lcis,///<MGLCisect_vector before makeVertex,
				///will be so updated according to the vertex make.
	MGPVertex*& pvertex ///MGPVertex of the made vertex.
){
	MGEdge* pre; MGEdge* aft;
	bool made=make_vertex(lp,pre,aft,&pvertex);
	if(made)
		lcis.update_lepoint(lp);
	return made;
}


//Subdivide this loop so that one parameter range in rangesin becomes one edge.
//In rangesin, parameter range of this loop is stored as:
//Let n=rangesin.size()/2, then i-th span from ranges[2*i] to ranges[2*i+1] is one
//parameter span that is supposed to be one edge for 0<=i<n.
//Returned are new edge pointers that correspond to the i-th span
//after subdivided.
//****Currently this does not conform to non_manifold model.
//That is, this loop must not have partner edges already.
std::vector<MGEdge*> MGLoop::subdivide(
	const std::vector<MGLEPoint>& rangesin
){
	int nle=(int)rangesin.size();
	int n=nle/2;
	if(!n) return std::vector<MGEdge*>();

	//subdivide the edges.
	//Subdivision is done once for all the MGLEPoints of one same edge.
	MGEdge *epre, *eaft;
	std::vector<MGEdge*> edges(n);
	int i=nle-1;
	int ei;
	while(i>0){
		ei=(i+1)/2;	//ei=(id of edges to store)+1.
		std::vector<MGLEPoint> ranges;
		ranges.push_back(rangesin[i]); ranges.push_back(rangesin[i-1]);
		int incr=-1; if(rangesin[i]>rangesin[i-1]) incr=1;
			//When LEPoints in renges[] are stored in decending order, incr=1,
			//and when in ascending order, incr=-1.
		i-=2;
		while(i>0 && rangesin[i+1].equal_edge(rangesin[i])){
			 ranges.push_back(rangesin[i]); ranges.push_back(rangesin[i-1]);
			 i-=2;
		}

		int mle=(int)ranges.size();//num of LMEPoints in the same edge.
		int m=mle/2;		  //num of new edges in the same edge.
		int id;		//id of ranges, will be processed in descending order.
		if(incr==1){//In this case, ranges are stored in descending order.
			id=0; ei--;
		}else{      //In this case, ranges are stored in ascending order.
			id=mle-1; ei-=m;
		}
		for(int j=0; j<m; j++){
			make_vertex(ranges[id], epre, eaft); id+=incr;
			make_vertex(ranges[id], epre, eaft); id+=incr;
			edges[ei]=eaft; ei-=incr;
		}
	}
	return edges;
}

//le1, and 2 must be of the same edge.
//new edge between le1 and le2 will be output.
MGEdge* MGLoop::subdivide(
	MGLEPoint& le1, MGLEPoint& le2
		//parameter ranges of binder edges, not loop's parameter edge.
){
	assert(le1.equal_edge(le2));//Both must be of the same edge.

	//subdivide the edges.
	MGEdge *epre, *eaft;
	if(le1>le2){
		make_vertex(le1, epre, eaft);
		make_vertex(le2, epre, eaft);
	}else{
		make_vertex(le2, epre, eaft);
		make_vertex(le1, epre, eaft);
	}
	return eaft;
}
