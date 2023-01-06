#include "StdAfx.h"
#include "mg/Tolerance.h"
#include "mg/drawParam.h"
#include "mg/BPointSeq.h"
#include "mg/SurfCurve.h"
#include "topo/Edge.h"
#include "topo/LEPoint.h"
#include "topo/Loop.h"
#include "topo/Face.h"
#include "Tl2/TLInputParam.h"
#include "Tl2/TL2parameter.h"
#include "Tl2/TL2Fans.h"
#include "Tl2/TL2Polyline.h"
#include "Tl2/TL2LPlines.h"
#include "Tl2/TL2Face.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

/****************************************************************/
/*   Copyright (c) 2019 by System fugen G.K.                */
/*                       All rights reserved.                   */
/****************************************************************/


///////////Constructor//////////

//copy constructor.
mgTL2Face::mgTL2Face(const mgTL2Face& face)
:m_param(face.m_param){
}

mgTL2Face::mgTL2Face(
	const MGDrawParam& param,//parameter for the tessellation.
	const MGSurface& face	//テセレーションするフェイス
							//Must be MGFace or MGSurface.
):m_param(face,param){
	polygonizeBoundaries();
}

mgTL2Face::mgTL2Face(
	const MGDrawParam& param,///<parameter for the tessellation.
	const MGFSurface& face,
	const std::vector<SHLL_COM_EDGES>* polylines
		///< Input polygonized polylines for the face boundaries.
		///< polylines[i][j] is a j-th edge's polyline for face.loop(i),
		///< must be MGLBRep of order 2.
		///< polylines[i][j]=0 indicates loop i's edge j can be face's bounday and
		///< has any common edges.
		///< **polylines[i][j] must be the same direction as the faces's parameter edge.
):m_param(face,param,polylines){
	polygonizeBoundaries();
}

mgTL2Face::mgTL2Face(
	const MGFSurface& face,///<テセレーションするフェイス
						///<Must be MGFace or MGSurface.
	double crvTol,		///<バウンダリのトレランス
	double surfTol,		///<平面とみなすトレランス
	double max_edge_len,///<when max_edge_len<=0, this means no limits on an edge length.
	const std::vector<SHLL_COM_EDGES>* polylines
		///< Input polygonized polylines for the face boundaries.
		///< polylines[i][j] is a j-th edge's polyline for face.loop(i),
		///< must be MGLBRep of order 2.
		///< polylines[i][j]=0 indicates loop i's edge j can be face's bounday and
		///< has any common edges.
		///< **polylines[i][j] must be the same direction as the faces's parameter edge.
):m_param(face,crvTol,surfTol,polylines,max_edge_len){
	polygonizeBoundaries();
}

mgTL2Face::mgTL2Face(
	const mgTLInputParam& param,///<parameter for the tessellation.
	const MGFSurface& face	///<テセレーションするフェイス
) : mgTL2Face(face, param.crvTol(), param.surfTol(), param.max_edge_len()) {
	;
}

////////Member function//////////

///make_Edge() makes a polyline edge(parameter edge) of edge edgeuv which are wholly
///on a curve of Bpoly().
///All of the points of m_param.Bpoly()[id[0]][id[1]] will be converted to
///surface (u,v) parameter. These (u,v) representation makes the polyline edge.
MGEdge* make_Edge(
	const mgTL2parameter& tlpara,
	const MGEdge& edgeuv,
	short id[3],	//pointer to uniform LBRep of m_param.Bpoly(), must not be null
					//(id[0] and id[1]).
					//On return, id of the starting point will be return.
	mgTL2Polyline*& poly		//generated mgTL2Polyline* for the edge will be returned.
){
	const std::vector<SHLL_COM_EDGES>& Bpolylines=*(tlpara.Bpoly());
	const MGLBRep* edgexyz=Bpolylines[id[0]][id[1]];
	const MGSurface& srf=tlpara.get_surface();
	MGTrimmedCurve crvuv=edgeuv.trimmed_curve();
	MGSurfCurve plinexyz(srf,crvuv);
	const MGBPointSeq& bpxyz=edgexyz->line_bcoef();

	int nbd=edgexyz->bdim();
	bool equalDirection=edgeuv.equal_direction_to_binder();
	double s0=crvuv.param_s(), s1=crvuv.param_e();
	double sguess=s0;

	poly=new mgTL2Polyline(tlpara);
	poly->set_type(mgTL2Polyline::WHOLE_BOUNDARY);
	MGBPointSeq& uvbp=poly->line_bcoef();
	MGKnotVector& uvKnotv=poly->knot_vector();
	uvbp.resize(nbd,2);uvKnotv.size_change(2,nbd);
	for(int k=0; k<nbd; k++){
		int kp1=k+1;
		MGPosition xyz=equalDirection ? bpxyz(k):bpxyz(nbd-kp1);
		double tp;
		int pobtained=plinexyz.perp_guess(s0,s1,xyz,sguess,tp);
		if(!pobtained)
			plinexyz.on(xyz,tp);
		sguess=tp;
		MGVector uv=crvuv.eval(tp);
		uvbp.store_at(k,uv);
		uvKnotv(kp1)=double(k);
	}
	uvKnotv(0)=uvKnotv[1];
	uvKnotv(nbd+1)=uvKnotv[nbd];

	id[2]=equalDirection ? short(nbd-1):0;
	poly->set_endID(id);
	id[2]=equalDirection ? 0:short(nbd-1);
	poly->set_startID(id);
	return new MGEdge(poly);
}

///Polygonize all the boundaries of the target face,
///and make an MGFace that has the polygonized boundaries.
///The face made will be m_face.
///The face does not have surface geometery, only has bounfaries.
void mgTL2Face::polygonizeBoundaries(){
	size_t nloop=0;
	const MGFace& fOrigin=m_param.get_face();
	if(m_param.target_is_face())
		nloop=(size_t)fOrigin.number_of_loops();

	if(!nloop){
		polygonizeSurfaceBoundaries();
		return;
	}

	m_face=std::unique_ptr<MGFace>(new MGFace);
	size_t nBpolylines=0;
	const std::vector<SHLL_COM_EDGES>& Bpolylines=*(Bpoly());
	if(&Bpolylines)
		nBpolylines=Bpolylines.size();

	for(size_t i=0; i<nloop; i++){
		const UniqueLoop& lpi=fOrigin.loop((int)i);
		int nEdges_Bpolylinesi=0;
		const SHLL_COM_EDGES* Bpolylinesi=0;
		if(i<nBpolylines){
			Bpolylinesi=&(Bpolylines[i]);
			nEdges_Bpolylinesi=(int)Bpolylinesi->size();
		}

		short id[3]; id[0]=(short)i;//id[0] is loop number.
		short idtemp[3];
		int nedge=lpi->number_of_edges();
		MGLoop* lp=new MGLoop;
		const MGLBRep* edgepoly=0;
		mgTL2Polyline* polyPre=0;
		mgTL2Polyline* polyStart=0;
		mgTL2Polyline* poly=0;
		for(int j=0; j<nedge; j++){
			id[1]=j;//id[1] is endge number in the loop.
			MGEdge* eij=0;//edge j in loop i(polyline reresentation).
			const MGEdge& edgeuv=*(lpi->edge(j));//The original edge j in loop i.
			const MGLBRep* eijxyz=0;
			if(j<nEdges_Bpolylinesi){
				eijxyz=(*Bpolylinesi)[j];
				if(eijxyz){
					eij=make_Edge(m_param,edgeuv,id,poly);//id[2] was set in make_Edge().
					if(polyPre){
						if(polyPre->boundaryType()<mgTL2Polyline::START_END_BOUNDARY)
							polyPre->set_endID(id);
					}
				}
			}
			if(!eijxyz){
				poly=new mgTL2Polyline(m_param,edgeuv.trimmed_curve());
				eij= new MGEdge(poly);
				if(polyPre){
					if(polyPre->boundaryType()==mgTL2Polyline::WHOLE_BOUNDARY){
						polyPre->get_endID(idtemp);
						poly->set_startID(idtemp);
					}
				}
			}
			lp->append(eij);
			if(j==0)
				polyStart=poly;
			polyPre=poly;
		}
		if(polyStart->boundaryType()<mgTL2Polyline::START_END_BOUNDARY
			&& poly->boundaryType()==mgTL2Polyline::WHOLE_BOUNDARY){
			poly->get_endID(idtemp);
			polyStart->set_startID(idtemp);
		}else if(polyStart->boundaryType()==mgTL2Polyline::WHOLE_BOUNDARY
			&& poly->boundaryType()<mgTL2Polyline::START_END_BOUNDARY){
			polyStart->get_startID(idtemp);
			poly->set_endID(idtemp);
		}
		lp->make_close();
		m_face->append_boundary(lp);
	}
}

///Polygonize MGSurface boundaries,
///and make an MGFace that has the polygonized boundaries.
///The face made will be m_face.
///The face does not have surface geometery, only has bounfaries.
///The target must be MGSurface.
void mgTL2Face::polygonizeSurfaceBoundaries(){
	m_face=std::unique_ptr<MGFace>(new MGFace);
	const MGSurface& srf=m_param.get_surface();
	MGLoop* lp=new MGLoop;
	std::vector<UniqueCurve> bcurves=srf.outer_boundary_param();
	size_t n=bcurves.size();
	for(size_t i=0; i<n; i++){
		mgTL2Polyline* poly=new mgTL2Polyline(m_param,*(bcurves[i]));
		MGEdge* ei=new MGEdge(poly);
		lp->append(ei);
	}
	lp->make_close();
	m_face->append_boundary(lp);
}

ostream& operator<< (ostream& out, const mgTL2Face& face){
	out<<"mgTL2Face="<<(&face);
	out<<face.m_param<<endl;
	out<<"mgTL2Face::m_face="<<*(face.m_face)<<endl;
	return out;
}

///Do perform the tessellation.
///The result be appended onto triangles.
///When triangles.is_uv()=false, all of the element of the triangle position data has normal data as
///(x,y,z,xn,yn,zn). Here (x,y,z) is the position data and (xn,yn,zn) is the normal vector
///at the position (x,y,z).
///When triangles.is_uv()=true, all of the element of the triange position data are (u,v).
void mgTL2Face::tessellate(
	mgTL2Triangles& triangles	//Tessellated triangles will be output.
){
	const MGSurface& srf=m_param.get_surface();
	triangles.set_surface(&srf);
	double srftol=m_param.get_tess_srfError();
	double mlen2=m_param.get_max_edge_len_sqr();

	std::stack<UniqueFace> faceStack;
	faceStack.emplace(m_face.release());	///<MGFace sequence ordered.
	std::vector<UniqueFace> flatPolygons;
	while(!faceStack.empty()){
		UniqueFace faceP(faceStack.top().release()); faceStack.pop();
		MGFace& face = *faceP;//std::cout << face << std::endl;
		if(face.number_of_loops()==1){
			const MGLoop& oloop = *(face.loop(0));
			int nedges=oloop.number_of_edges();
			if(nedges<=4 && findConcaveVertex(oloop)==-1){
				std::stack<std::unique_ptr<mgTL2LPlines>> LPlinesStack;
				LPlinesStack.push(std::make_unique< mgTL2LPlines>(LPlinesStack,oloop, triangles));
				while (!LPlinesStack.empty()) {
					std::unique_ptr<mgTL2LPlines> LPlines(LPlinesStack.top().release());
					LPlinesStack.pop();
					//std::cout << *LPlines << std::endl;
					LPlines->tessellate4();
				}
				continue;
			}
			if(srf.is_flat_and_small(oloop.box(),srftol,mlen2)){
				flatPolygons.emplace_back(faceP.release());
				continue;
			}
		}

		if(!splitTl(face,m_param, faceStack))
			//When splitTl failed to split.
			flatPolygons.emplace_back(faceP.release());
	}

	//triangulate each polygon.
	for(auto& i: flatPolygons)
		triangulate(*(i->loop(0)),triangles);
}
