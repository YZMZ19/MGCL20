/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

//! @file
//!	@brief  Declaration for class MGIges504EdgeListMap.
//!	@author System fugen

#include "StdAfx.h"
#include "mg/CompositeCurve.h"
#include "mgIges/IgesVertexListMap.h"
#include "mgIges/iges504EdgeListMap.h"
#include "mgIges/igespd508.h"
#include "mgIges/Igesifstream.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace std;

//!	@brief MGIges504EdgeListMap is the class to store MGEdge*(newed objects) that
//are generated for MGIges504 EDGE list.

//Obtain a pair of (binder MGEdge* Bedge, parameter MGEdge* Pedge) of edge508.
//Bedge is from the the DE of the EDGE List Entry(MGIgesPD504) and the indgex edge.
//Pedge is from m_pcurves, which may be null when not specified.
void MGIges504EdgeListMap::get_edgePair(
	const MGIges508Edge& edge508,
	MGEdge*& Bedge,
	MGEdge*& Pedge
){
	//1. Get Bedge.

	int DEid = edge508.m_edgeListDE;	//DE id of type 504
	int edge = edge508.m_edge;	//List index of the EDGE List Entryt DE edge_list
		
	const MGIgesFstream::UniqueDE& de=m_ifstream->directoryEntry(DEid);
	const std::unique_ptr<MGIgesPD>& pd=de->paramData();
	const MGIgesPD504& pd504=*(static_cast<const MGIgesPD504*>(pd.get()));//reference to MGIgesPD504.
	int vID=(int)m_EdgesVector.size();
	pair<PD504toEdgesMap::iterator, bool> insertR=
		m_PD504PtoVecID.insert(make_pair(&pd504,vID));
	if(insertR.second){
		//If the pair(&pd504,vID) is inserted, make std::vector<MGEdge*> for this pd504.
		int nedges=(int)pd504.m_edges.size();
		MGEdge* eddummy=0;
		m_EdgesVector.push_back(std::vector<MGEdge*>(nedges,eddummy));//insert dummy array.
	}else//If not.
		vID=insertR.first->second;//Already inserted &pd504's array index.

	assert(vID<int(m_EdgesVector.size()));
	std::vector<MGEdge*>& edges=m_EdgesVector[vID];
	assert(edge<int(edges.size()));
	Bedge =edges[edge];
	if(!Bedge){
		//This is the 1st refernece to pd504's edge. Generate MGEdge*.
		const MGIges504Edge& edge504=pd504[edge];
		MGCurve* crv= dynamic_cast<MGCurve*>(m_ifstream->convert_to_gel(edge504.m_curve_DE));
		if(crv){
			Bedge =new MGEdge(crv);
			edges[edge]= Bedge;
		}
	}

	//2. Get Pedge.
	Pedge = nullptr;
	int npcurves = edge508.number_of_pcurves();
	if(!npcurves)
		return;
	const std::vector<int>& pedges = edge508.m_pcurves;
	MGCurve* pcrv = dynamic_cast<MGCurve*>(m_ifstream->convert_to_gel(pedges[1]));
	pcrv->change_dimension(2);
	if(npcurves==1){
		Pedge = new MGEdge(pcrv);
		return;
	}
	MGCompositeCurve* pcompo = new MGCompositeCurve(pcrv);
	for(size_t i = 2; i<=npcurves; i++){
		MGCurve* pcrvi = dynamic_cast<MGCurve*>(m_ifstream->convert_to_gel(pedges[i]));
		pcompo->connect_to_end(pcrvi);
	}
	Pedge = new MGEdge(pcompo);
}
