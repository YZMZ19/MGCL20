/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD508.
//!	@author System fugen

#include "StdAfx.h"
#include "mg/CompositeCurve.h"
#include "topo/BVertex.h"
#include "topo/Edge.h"
#include "topo/Loop.h"
#include "mgiges/IgesGsec.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesPD508.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace MGCL;
using namespace MGIges;
//!	@brief MGIgesPD508 is the class for Iges parameter data type 508(LOOP Entity).

// Constructors.

//! Constructs an object of class MGIgesPD508.
MGIgesPD508::MGIgesPD508(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(LOOP,DEpointer), m_edges(1){
}

//Read in parameter data from string stream data.
void MGIgesPD508::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	int num_edges;
	get_integer(pDelimeter,pdstream,num_edges);
	m_edges.resize(num_edges+1);
	for(int i=1; i<=num_edges; i++){
		MGIges508Edge* edge=new MGIges508Edge;
		edge->read_in(pDelimeter,pdstream);
		m_edges[i].reset(edge);
	}
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD508::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	std::vector<std::string>& plines ///<output plines.
)const{
	int num_edges=(int)(m_edges.size()-1);
	put_integer(num_edges,gsec,plines);
	for(int i=1; i<=num_edges; i++){
		const MGIges508Edge& edgei=*(m_edges[i]);
		edgei.write_out_into_string(gsec,plines);
	}
}

//Read in parameter data from string stream data.
void MGIges508Edge::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	get_integer(pDelimeter,pdstream,m_type);
	get_DEpointer(pDelimeter,pdstream,m_edgeListDE);
	get_integer(pDelimeter,pdstream,m_edge);
	get_integer(pDelimeter,pdstream,m_orientation);
	int num_parameter_curve;
	get_integer(pDelimeter,pdstream,num_parameter_curve);
	m_isoparameterics.resize(num_parameter_curve+1);
	m_pcurves.resize(num_parameter_curve+1);
	for(int i=1; i<=num_parameter_curve; i++){
		int iso;
		get_integer(pDelimeter,pdstream,iso);
		m_isoparameterics[i]=iso ? 1:0;
		int pcurveDE;//pointer to the DE of the underlying parameter curve.
		get_DEpointer(pDelimeter,pdstream,pcurveDE);
		m_pcurves[i]=pcurveDE;
	}
}

///Write out this MGIges508Edge as MGIgesParamLine's(into plines).
///Except for string data, one integer or double data is output
///into one MGIgesParamLine, not striding over more than one line.
///Only when string data is output(to Holleris string), the data
///may stride over more than one lines.
///plines[i] for 0<=i<plines.size() are valid.
void MGIges508Edge::write_out_into_string(
	const MGIgesGSec& gsec,	///<Input gsec to input delimeter_param and delimeter_record;
	std::vector<std::string>& plines ///<output plines.
)const{
	put_integer(m_type,gsec,plines);
	put_DEpointer(m_edgeListDE,gsec,plines);
	put_integer(m_edge,gsec,plines);
	put_integer(m_orientation,gsec,plines);
	int num_parameter_curve=(int)(m_pcurves.size()-1);
	put_integer(num_parameter_curve,gsec,plines);
	for(int i=1; i<=num_parameter_curve; i++){
		int iso=m_isoparameterics[i] ? 1:0;
		put_integer(iso,gsec,plines);
		put_DEpointer(m_pcurves[i],gsec,plines);
	}
}

///Get the edge pointer(newed object).
MGBVertex* MGIges508Edge::get_BVertex(const MGIgesIfstream& ifs)const{
	return ifs.m_vertexListMap.get_BVertex(m_edgeListDE,m_edge);
}

void buildPEdgeFromWorldCurve(
	const MGSurface& srf,//The surface on which binder edge bedge lies.
	MGEdge* bedge,       //Binder edge(The curve is of world space).
	std::unique_ptr<MGLoop>& loop,//loop to append the built parameter edge into.
	std::vector<double>& pspan,//When bedge lies on a perimeter of srf, the parameter span is input.
	int peri_num,        //perimeter number when bedge lies on a perimeter.
	bool oppositeDirection//true if bedge is opposite to the parameter edge to build.
){
	const MGCurve* wcurve = bedge->base_curve();
	double tLast = wcurve->param_s(), terror = wcurve->param_error();
	MGEdge* eAppended=loop->append_edge_from_crvWorld(
		srf, *wcurve, tLast, terror, pspan, peri_num,oppositeDirection
	);
	eAppended->shareBinderIfHadPartner(bedge);
}

//Convert de(type=508: LOOP entry) to MGLoop.
//Returned is a newed MGLoop object.
MGLoop* MGIgesIfstream::convert_loop(
	const MGIgesDirectoryEntry& de,	//directory entry of type 508.
	const MGSurface& srf //Base surface whose boundary this loop will make.
)const{
	const std::unique_ptr<MGIgesPD>& pd=de.paramData();
	const MGIgesPD508* pd508=static_cast<const MGIgesPD508*>(pd.get());
	int nedges=pd508->number_of_edges();
	if(!nedges)
		return 0;

	std::vector<double> pspan[2]; int peri_num[2];
	int startID=1;
	const MGIges508Edge& edge1=pd508->edge(1);
	MGEdge* bedge1=nullptr;
	MGEdge* pedge1 = nullptr;
	const MGCurve* wcurve1 = nullptr;

	std::unique_ptr<MGLoop> loop = std::unique_ptr<MGLoop>(new MGLoop);
	if(!edge1.is_vertex()){
		m_edgeListMap.get_edgePair(edge1, bedge1, pedge1);//This is a binder edge.
		if(!bedge1)
			return 0;

		if(pedge1){
			loop->append(pedge1);
			pedge1->shareBinderIfHadPartner(bedge1);
		} else{
			wcurve1 = bedge1->base_curve();
			int nCom1st = srf.getPerimeterCommon(*wcurve1, pspan, peri_num);
			if(nCom1st>=2)
				startID = 2;
			else{
				buildPEdgeFromWorldCurve(
					srf, bedge1, loop, pspan[0], peri_num[0], edge1.direction_is_opposite()
				);
			}
		}
	}

	std::vector<double> pspan2[2]; int peri_num2[2];
	std::vector<double>& pspan0=pspan2[0];
	std::vector<double>& pspan1=pspan2[1];
	for(int i=2; i<=nedges; i++){
		const MGIges508Edge& edgei=pd508->edge(i);
		if(edgei.is_vertex())
			continue;

		MGEdge* bedgei;
		MGEdge* pedgei;
		m_edgeListMap.get_edgePair(edgei,bedgei, pedgei);//This is a binder edge.
		if(!bedgei)
			return 0;

		if(pedgei){
			loop->append(pedgei);
			pedgei->shareBinderIfHadPartner(bedgei);
		} else{
			const MGCurve* wcurve = bedgei->base_curve();
			int nComi = srf.getPerimeterCommon(*wcurve, pspan2, peri_num2);
			if(nComi>=2){
				MGPosition uv0 = srf.perimeter_uv(peri_num2[0], pspan0[2]);
				MGPosition uv1 = srf.perimeter_uv(peri_num2[1], pspan1[2]);
				MGPosition uvE = loop->end_point();
				if((uvE-uv1).len()<(uvE-uv0).len()){
					pspan0 = pspan1; peri_num2[0] = peri_num2[1];
				}
			}
			buildPEdgeFromWorldCurve(
				srf, bedgei, loop, pspan0, peri_num2[0], edgei.direction_is_opposite()
			);
		}
	}

	if(startID==2){
		MGPosition uv0=srf.perimeter_uv(peri_num[0],pspan[0][2]);
		MGPosition uv1=srf.perimeter_uv(peri_num[1],pspan[1][2]);
		MGPosition uvE=loop->end_point();
		if((uvE-uv1).len()<(uvE-uv0).len()){
			pspan[0]=pspan[1]; peri_num[0]=peri_num[1];
		}
		buildPEdgeFromWorldCurve(
			srf, bedge1, loop, pspan[0], peri_num[0], edge1.direction_is_opposite()
		);
	}

	loop->make_close();
	return loop.release();
}

///String stream function
std::ostream& operator<< (std::ostream& ostrm, const MGIges508Edge& edge){
	int n=edge.number_of_pcurves();
	ostrm<<"MGIges508Edge::"<<&edge<<",m_type="<<edge.m_type;
	ostrm<<",m_edgeListDE="<<edge.m_edgeListDE<<",m_edge="<<edge.m_edge;
	ostrm<<",m_orientation="<<edge.m_orientation
		<<",num pcurves="<<n<<std::endl;
	for(int i=1; i<=n; i++){
		ostrm<<i<<"("<<edge.m_isoparameterics[i]<<","<<edge.m_pcurves[i]<<")";
		if(i<n)
			ostrm<<",";
		else
			ostrm<<std::endl;
	}
	return ostrm;
}

