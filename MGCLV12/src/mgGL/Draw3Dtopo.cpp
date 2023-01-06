/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include <map>
#include "mg/Point.h"
#include "mg/Plane.h"
#include "mg/PickObjectFB.h"
#include "topo/PVertex.h"
#include "topo/BVertex.h"
#include "topo/Edge.h"
#include "topo/Loop.h"
#include "topo/Face.h"
#include "topo/Shell.h"
#include "Tl2/TL2Face.h"
#include "Tl2/TL2Triangles.h"
#include "Tl2/TLInputParam.h"
#include "mgGL/VBO.h"
#include "mgGL/Appearance.h"

using namespace std;

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//Implements the drawWire functions of all the classes.

///Draw the object in wire mode, in the world coordinates.
///The object is converted to curve(s) and is drawn.
void MGComplex::drawWire(
	mgVBO& vbo,
	int line_density	//line density to draw a surface in wire mode.
)const{
	mgLightModeSwitcher switcher(vbo, mgGLSL::NoShading);
	for(auto& pcelli:m_pcells){
		pcelli->drawWire(vbo, line_density);
	}
}

//Draw curve in the star face world coordinates.
///This is a boundary loop(MGPCell) of a face, and the curves are extracted
///from the binder edges of this and drawn.
void MGLoop::drawWire_in_star(
	mgVBO& vbo,
	int line_density	//line density to draw a surface in wire mode.
)const{
	if(!star()) return;

	mgLightModeSwitcher switcher(vbo, mgGLSL::NoShading);
	for(auto& pcelli:pcells()){
		const MGEdge& ei = *(dynamic_cast<const MGEdge*>(pcelli.get()));
		ei.drawWire_in_star(vbo, line_density);
	}
}

//////////////////////////////////////////////

//Draw 3D curve in the star face world coordinates.
///This is a boundary edge(MGPCell) of a face, and the curves are extracted
///from the binder edge and drawn.
void MGEdge::drawWire_in_star(
	mgVBO& vbo,
	int line_density	//line density to draw a surface in wire mode.
)const{
	if(!star()) return;

	mgLightModeSwitcher switcher(vbo, mgGLSL::NoShading);
	SharedBCell& bndr=make_binder_with_extent();
	const MGEdge* be = dynamic_cast<const MGEdge*>(bndr.get());
	be->drawWire(vbo,line_density);
}

void MGBVertex::drawVertex(
	mgVBO& vbo
)const{
	make_extent();
	const MGPoint* pnt=point();
	if(pnt){
		mgLightModeSwitcher switcher(vbo, mgGLSL::NoShading);
		const MGPosition& P=pnt->position();
		vbo.drawPoint(P[0],P[1],P[2]);
		return;
	}
}

void MGEdge::drawWire(
	mgVBO& vbo,
	int line_density	//line density to draw a surface in wire mode.
)const{
	const MGCurve* crv=base_curve();
	crv->drawSE(vbo, param_s(), param_e());
}

void MGEdge::drawVertex(mgVBO& vbo)const{
	mgLightModeSwitcher switcher(vbo, mgGLSL::NoShading);
	if(active_start()){
		MGPosition Ps=start_point();
		vbo.drawPoint(Ps[0],Ps[1],Ps[2]);
	}
	if(active_end()){
		MGPosition Pe=end_point();
		vbo.drawPoint(Pe[0],Pe[1],Pe[2]);
	}
}

///Draw vertices of this loop. The vertices coordinate kind is the same as
///the edges. That is, the vertices are MGBVertex.
void MGLoop::drawVertex(
	mgVBO& vbo
)const{
	std::vector<MGBCell*> bvec;
	get_all_boundary_binders(bvec);
	for(auto& bcel:bvec){
		auto bv = dynamic_cast<const MGBVertex*>(bcel);
		bv->drawVertex(vbo);
	}
}

void MGFace::drawVertex(
	mgVBO& vbo
)const{
	const MGBox& pbx=box_param();
	if(pbx.is_null()) return;

	mgLightModeSwitcher switcher(vbo, mgGLSL::NoShading);
	//Draw boundaries' vertexes.
	for(auto& loop : boundaries())
		loop->drawVertex(vbo);
}

static const float edgeColor[4]={1.,.5,.5,1.};//Edge color to hilight.

///Highlight the object using the display list of this object.
void MGPickObjectFB::hilight_using_display_list(
	int line_density	///<line density to draw a surface in wire mode.
)const{
	MGPickObject::hilight_using_display_list(line_density);
	if(!m_vbo.is_made()){
		m_vbo.setStaticAttribColor(edgeColor);//glColor4fv(edgeColor);
		const MGEdge& e=*(edge());
		e.drawWire(m_vbo);
	}
	m_vbo.draw(MGCL::HIGHLIGHT);
}

//Make a display list of this gel.
void MGShell::make_display_list(
	MGCL::VIEWMODE vmode
)const{
	mgVBO* vbo=dlist_name();
	vbo->initializeVBO(vmode);
	
	if(vmode==MGCL::WIRE || vmode==MGCL::WIRE_AND_SHADING){
		mgLightModeSwitcher switcher(*vbo, mgGLSL::NoShading);
		int nfaces=number_of_faces();
		//Make display list of the gel that has an object.
		for(int i=0; i<nfaces; i++){
			const MGFace& facei=*(face(i));
			if(facei.no_display())
				continue;
			vbo->drawGel(facei);//glCallList(name2);//call object.
		}
		vbo->setDirty(false,MGCL::WIRE);
	}
	if(vmode==MGCL::SHADING || vmode==MGCL::WIRE_AND_SHADING){
		mgLightModeSwitcher switcher(*vbo, mgGLSL::Shading);
		const MGDrawParam& para=mgVBOElement::getDrawParam();
		shade(*vbo,para);
		vbo->setDirty(false, MGCL::SHADING);
	}
}

//Map of MGEdge(binder edge) to the polygonized MGLBRep.
typedef std::map<const MGEdge*,MGLBRep*> ELMAP;

//Set up common edges world coordinate line data in polylines.
//polylines[i] is vector<SHLL_COM_EDGES> for face(i)
//for i=0,...,number_of_faces()-1.
void set_up_shell_shade(
	const MGShell& shell,//Target shell
   	MGDrawParam& para,///< input tessellation parameter.
				/// When para.maximum_edge_length_tess()<=0. is input, maximum_edge_length is
				///computed and set in para.
	std::vector<UniqueLBRep>& boudaries,///<Container of the boundary data.
		//Polylined boundaries will be held in boundaries.
	std::vector< std::vector<SHLL_COM_EDGES> >& polylines//polylines[i] holds boundary
		//data of shell.face(i)'s loops polylined boundaries data.
){
	//shell.ensure_BVertices_of_ModelEdges();
	double error=para.curve_tolerance_tess();
	double melen=para.maximum_edge_length_tess();
	if(melen<=0.){
		melen=compute_max_edge_len(shell);
		para.set_maximum_edge_length_tess(melen);
	}
	double melen2=melen*melen;

	//Map to find which poligonized MGLBRep original MGEdge corresponds to.
	ELMAP boundary_map;//The key is MGEdge

	MGComplex::const_iterator i=shell.pcell_begin();
	int nfaces=shell.number_of_faces();
	polylines.resize(nfaces);
	for(int j=0; j<nfaces; j++,i++){//iterate over all the faces.
		const MGFace& fj=*(shell.face(i));

		vector<SHLL_COM_EDGES>& fjpolylines=polylines[j];
		int nloop=fj.number_of_loops();
		fjpolylines.resize(nloop);
		for(int k=0; k<nloop; k++){//iterate over all the loops in face fj.
			const MGLoop& lpk=*(fj.loop(k));

			SHLL_COM_EDGES& lpkcomedges=fjpolylines[k];
			int nedges=lpk.number_of_edges();
			lpkcomedges.resize(nedges);
			for(int l=0;l<nedges; l++){//iterate over all the edges of loop lpk.
				const MGEdge& el=*(lpk.edge(l));
				const MGEdge* bel=el.make_binder_with_curve();
				ELMAP::iterator mapi=boundary_map.find(bel);
				MGLBRep* poly;
				if(mapi==boundary_map.end()){
					//If not found, the 1st reference to the edge.
					unique_ptr<MGCurve> belCrv=unique_ptr<MGCurve>(bel->curve_limitted());
					poly=new MGLBRep();
					MGLBRep polyTemp;
					belCrv->polygonize(error,polyTemp);
					getXYZline_ensuring_max_edge_length(melen2,polyTemp, *poly);
					boudaries.emplace_back(poly);
					boundary_map.insert(make_pair(bel,poly));
				}else{
					poly=mapi->second;
				}
				lpkcomedges[l]=poly;
			}
		}
	}
}

//Shade the object in world coordinates.
void MGFace::shade(
	mgVBO& vbo,
	const MGDrawParam& para,
	MGCL::DRAW_TARGET target
)const{
	mgTL2Face face(para,*this);
	mgTL2Triangles tris(MGCL::XYZNormal,surface());
	face.tessellate(tris);	
	vbo.drawShade(tris,target);
}

//Shade the object in world coordinates.
void MGFace::shade(
	mgVBO& vbo,
	const MGDrawParam& para,
	std::vector<SHLL_COM_EDGES>* polylines,
	MGCL::DRAW_TARGET target
)const{
	mgTL2Face face(para,*this,polylines);
	mgTL2Triangles tris(MGCL::XYZNormal,surface());
	face.tessellate(tris);	
	vbo.drawShade(tris,target);
}

///Triangulate this object(MGShell, MGFace, or MGSurface is the target).
void MGFace::triangulate(
	const MGDrawParam& para,
	MGCL::TL_DATA_KIND dkind,
	std::vector<mgTL2Triangles>& trisVec
)const{
	mgTL2Face face(para,*this);
	mgTL2Triangles tris(dkind,surface());
	face.tessellate(tris);	
	trisVec.push_back(std::move(tris));
}

//Shade the object in world coordinates.
void MGShell::shade(
	mgVBO& vbo,
	const MGDrawParam& para,
	MGCL::DRAW_TARGET target
)const{
	vector<UniqueLBRep> boundaries;//Container of the boundary data.
	vector< vector<SHLL_COM_EDGES> > polylines;
	MGDrawParam para2(para);
	set_up_shell_shade(*this,para2,boundaries,polylines);//shade

	MGComplex::const_iterator i=pcell_begin();
	int nfaces=number_of_faces();
	for(int j=0; j<nfaces; j++,i++){
		const MGFace& fj = *(face(i));
		mgTL2Face face(para2,fj,&(polylines[j]));
		mgTL2Triangles tris(MGCL::XYZNormal,fj.surface());
		face.tessellate(tris);	
		mgVBO* fjvbo=fj.dlist_name();
		fjvbo->drawShade(tris,target);
	}
}

///Triangulate this object(MGShell, MGFace, or MGSurface is the target).
void MGShell::triangulate(
	const MGDrawParam& para,
	MGCL::TL_DATA_KIND dkind,
	std::vector<mgTL2Triangles>& trisVec
)const{
	vector<UniqueLBRep> boundaries;//Container of the boundary data.
	std::vector< std::vector<SHLL_COM_EDGES> > polylines;
	MGDrawParam para2(para);
	set_up_shell_shade(*this,para2,boundaries,polylines);

	MGComplex::const_iterator i=pcell_begin();
	int nfaces=number_of_faces();
	for(int j=0; j<nfaces; j++,i++){
		const MGFace& fj=*(face(i));
		mgTL2Face face(para2,fj,&(polylines[j]));
		mgTL2Triangles tris(dkind,fj.surface());
		face.tessellate(tris);
		trisVec.push_back(std::move(tris));
	}
}
