/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"

#include "mg/Tolerance.h"
#include "mg/DnameControl.h"
#include "mg/Object.h"
#include "mg/Box.h"
#include "mg/LBRep.h"
#include "mg/SPointSeq.h"
#include "mg/FSurface.h"
#include "mg/Surface.h"
#include "mg/Plane.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/MGStl.h"
#include "mg/PickObjects.h"
#include "mg/PickObjectCB.h"
#include "mg/PickObjectSB.h"
#include "Tl2/TL2Triangles.h"
#include "Tl2/TL2Face.h"
#include "mgGL/Appearance.h"
#include "mgGL/Context.h"

#ifndef _CONSOLE

#include "mgGL/OpenGLView.h"
#include "mgGL/Lights.h"
#include "mgGL/Light.h"
#include "mgGL/PlaneImage.h"
#include "mgGL/DirectionalLight.h"

#endif //_CONSOLE

#include "mgGL/VBO.h"

//でばっぐよう

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

MGDrawParam::MGDrawParam(
	const MGContext& contx,
	double span_length_wire
) :m_span_length_wire(span_length_wire),
m_line_desity_wire_face(contx.line_density()),
m_maximum_edge_length_tess(contx.tess_maximum_edge_length()) {
	build_crv_srf_tolerance(contx.tess_curve_tolerance(), contx.tess_surface_tolerance());
}

///Triangulate this object(MGShell, MGFace, or MGSurface is the target).
void MGSurface::triangulate(
	const MGDrawParam& para,
	MGCL::TL_DATA_KIND dkind,
	std::vector<mgTL2Triangles>& trisVec
)const {
	mgTL2Face face(para, *this);
	mgTL2Triangles tris(dkind, this);
	face.tessellate(tris);
	trisVec.push_back(std::move(tris));
}

#ifndef _CONSOLE

//Make 2 types of display list of this gel(wire and shading).
void MGObject::make_display_list(
	MGCL::VIEWMODE vmode
)const{
	mgVBO* vbo=dlist_name();
	vbo->initializeVBO(vmode);

	const MGDrawParam& para=mgVBOElement::getDrawParam();
	int line_density=para.line_desity_wire_face();
	if(vmode!=MGCL::SHADING){
		mgLightModeSwitcher switcher(*vbo, mgGLSL::NoShading);
		drawWire(*vbo,line_density);
		vbo->setDirty(false, MGCL::WIRE);
	}
	if(vmode==MGCL::SHADING || vmode==MGCL::WIRE_AND_SHADING){
		if(manifold_dimension()>=2){
			mgLightModeSwitcher switcher(*vbo, mgGLSL::Shading);
			shade(*vbo,para, MGCL::SHADING);
			vbo->setDirty(false,MGCL::SHADING);
		}
	}
}

///Make a display list of this gel.
void MGPlane::make_display_list(
	MGCL::VIEWMODE vmode
)const{
	mgVBO* vbo=dlist_name();
	if(vmode==MGCL::SHADING){
		vbo->clearElements(MGCL::SHADING);
	}else{
		if(vmode!=MGCL::WIRE)
			vmode=MGCL::WIRE_AND_SHADINGVIEW;
		vbo->clearElements(MGCL::WIRE_AND_SHADING);
	}
	vbo->clearStaticAttributes();
	drawAttrib(*vbo,false);

	//Here vmode=WIRE, SHADING, or WIRE_AND_SHADING.
	const MGDrawParam& para=mgVBOElement::getDrawParam();
	int line_density=para.line_desity_wire_face();
	if(vmode!=MGCL::SHADING){
		drawWire(*vbo,line_density);
		vbo->setDirty(false,MGCL::WIRE);
	}
	if(vmode!=MGCL::WIRE){
		shade(*vbo,para,MGCL::SHADING);
		vbo->setDirty(false,MGCL::SHADING);
	}
}

///Process of draw or render attributes.
void MGGroup::drawAttrib(
	mgVBO& vbo,///<The target graphic object.
	bool no_color///<if true, color attribute will be neglected.
)const{
	const MGAppearance* app=appearance();
	if(app){
		app->drawAttrib(vbo,no_color);
	} else {
		const MGContext* pCtx = this->context();
		if(pCtx)
			pCtx->drawAttrib(vbo, no_color);
	}
}

//Make a display list of this gel.
void MGGroup::make_display_list(
	MGCL::VIEWMODE vmode
)const{
	mgVBO* vbo=dlist_name();
	vbo->initializeVBO(MGCL::WIRE_AND_SHADINGVIEW);

	//Make display list of the gel that has an object.
	std::vector<const MGGroup*> nested_groups;
	MGGroup::const_iterator i,is=begin(), ie=end();	
	for(i=is; i!=ie; i++){
		MGGel* geli=i->get();
		const MGAttribedGel* agel = dynamic_cast<const MGAttribedGel*>(geli);
		if(!agel)
			continue;
		if(agel->no_display())
			continue;
		mgVBO* vboi=agel->dlist_name();
		if(!vboi)//This is to judge MGAppearance(MGAppearance must be excluded).
			continue;
		vbo->drawGel(*agel);
		agel->setDirty(true);
		const MGGroup* grp=dynamic_cast<const MGGroup*>(agel);
		if(!grp)
			continue;
		nested_groups.push_back(grp);
	}
	
	///make display list of the nested MGGroup.
	size_t ngroups(nested_groups.size());
	for(size_t j=0; j<ngroups; j++){
		const MGGroup& grpj=*(nested_groups[j]);
		grpj.make_display_list(vmode);
	}
	vbo->setDirty(false);
}

void MGPlane::display_arrows(mgSysGL& sgl)const{
	MGVector U,V;
	get_uv_display_vector(U,V);
	const MGPosition& cen = center();
	MGBox box(cen-U,cen+U);
	box.expand(cen-V);
	box.expand(cen+V);
	MGPosition pos[10], uv=center_param();
	arrow(box,uv[0],uv[1],pos);

	const MGColor& ucolor=MGColor::get_instance(MGColor::Red);
	const MGColor& vcolor=MGColor::get_instance(MGColor::Green);
	const MGColor& white=MGColor::get_instance(MGColor::White);
	ucolor.exec(sgl);
	sgl.drawArrow(&pos[0]);

	pos[3] = pos[0];
	vcolor.exec(sgl);
	sgl.drawArrow(&pos[3]);

	pos[6] = pos[0];			
	white.exec(sgl);
	sgl.drawArrow(&pos[6]);
}

//Shade the object in world coordinates.
void MGStl::shade(
	mgVBO& vbo,
	const MGDrawParam& para,
	MGCL::DRAW_TARGET target
)const{
	mgLightModeSwitcher switcher(vbo, mgGLSL::Shading);
	vbo.drawSTL(*this,target,GL_FILL);
}

//Draw 3D curve in world coordinates.
//The object is converted to curve(s) and is drawn.
void MGStl::drawWire(
	mgVBO& vbo,
	int line_density	//line density to draw a surface in wire mode.
)const{
	mgLightModeSwitcher switcher(vbo, mgGLSL::NoShading);
	vbo.drawSTL(*this, MGCL::WIRE,GL_LINE);
}

// 三角形ごとの法線ベクトルを表示する
void MGStl::display_arrows(mgSysGL& sgl)const{
	// 三角形の数を取得
	size_t nTriangle(m_vecNormlTriang.size());
	// 矢印描画のための座標値の配列
	MGPosition pos[4];
	// 色を示すインスタンスを生成
	const MGColor& white=MGColor::get_instance(MGColor::White);

	for(int i = 0; size_t(i) < nTriangle; i++){
	// それぞれの三角形の法線方向に矢印を描画
		// 三角形の中点を取得
		int i3 = i*3;
		pos[0] = (m_vecPos[m_indices[i3]]+m_vecPos[m_indices[i3+1]]+m_vecPos[m_indices[i3+2]])/3;

		// 各辺の長さを取得
		double dist[3];
		dist[0] = m_vecPos[m_indices[i3]].distance(m_vecPos[m_indices[i3+1]]);
		dist[1] = m_vecPos[m_indices[i3+1]].distance(m_vecPos[m_indices[i3+2]]);
		dist[2] = m_vecPos[m_indices[i3]].distance(m_vecPos[m_indices[i3+2]]);
		
		// 三角形の中から最長の辺の長さを取得し
		// その1/2の値を矢印の軸の長さに用いる
		double max = dist[0];
		for(int j = 0; j < 2; j++){
			if(dist[j] < dist[j+1]){
				max = dist[j+1];
			}
		}
		double len = max/2;

		// 矢印の先端の座標を計算
		const MGVector& vecX = m_vecNormlTriang[i] * len;
		pos[1] = pos[0] + vecX;

		// 矢印の両端の座標を求める処理
		// 矢印の先端座標に加えるベクトルを計算
		const MGVector& head_rootx = vecX * .3;
		// 三角形の任意の辺のベクトルを取得し、面の法線べクトルとの積算を行う
		MGUnit_vector arrowVec = (m_vecPos[m_indices[i3+1]]- m_vecPos[m_indices[i3]]).normalize();
		arrowVec *= m_vecNormlTriang[i];
		// 矢印の先端座標に加えるもう１つのベクトルを計算
		const MGVector& head_rooty = arrowVec*(.5*.3*vecX.len());
		// 矢印の両端の座標を計算
		pos[2]=pos[1]-head_rootx+head_rooty;
		pos[3]=pos[1]-head_rootx-head_rooty;

		// 矢印の描画を行う
		white.exec(sgl);
		sgl.drawArrow(pos);
	}
}

//Draw 3D curve in world coordinates.
//The object is converted to curve(s) and is drawn.
void MGFSurface::drawWireFS(
	mgVBO& vbo,
	int line_density	//line density to draw a surface in wire mode.
)const{
	mgLightModeSwitcher switcher(vbo, mgGLSL::NoShading);

	vbo.LineWidth(1.);//glLineWidth(1.);
	std::vector<UniqueCurve> ilines=inner_skeleton(line_density);
	for (const auto& ilinesi : ilines)
		ilinesi->drawWire(vbo, line_density);

	vbo.LineWidth(2.);//glLineWidth(2.);
	std::vector<UniqueCurve> bndries=get_all_boundaries();
	for(const auto& bndriesi: bndries)
		bndriesi->drawWire(vbo, line_density);
}

//Draw 3D curve in world coordinates.
//The object is converted to curve(s) and is drawn.
void MGFSurface::drawWireFS_to_highlight(
	mgVBO& vbo,
	int line_density	//line density to draw a surface in wire mode.
)const{
	mgLightModeSwitcher switcher(vbo, mgGLSL::NoShading);
	std::vector<UniqueCurve> ilines=skeleton(line_density);
	for(const auto& line:ilines)
		line->drawWire(vbo, line_density);
}

///Display direction arrows on the surface.
void MGFSurface::display_arrowsFS(mgSysGL& sgl,int udiv, int vdiv)const{

	MGPosition uv(2), pos[10];
	const MGBox& box = box_param2();
	double us = box[0].low_point(), ue = box[0].high_point(),
		vs = box[1].low_point(), ve = box[1].high_point();
	const MGColor& ucolor=MGColor::get_instance(MGColor::Red);
	const MGColor& vcolor=MGColor::get_instance(MGColor::Green);
	const MGColor& white=MGColor::get_instance(MGColor::White);
	for(int i = 0; i <= udiv; i++){
		uv(0) = (us*(udiv-i)+ue*i)/udiv;
		for(int j = 0; j <= vdiv; j++){
			uv(1) = (vs*(vdiv-j)+ve*j)/vdiv;
			arrow(uv, pos);

			ucolor.exec(sgl);
			sgl.drawArrow(&pos[0]);
			pos[3] = pos[0];
			
			vcolor.exec(sgl);
			sgl.drawArrow(&pos[3]);
			pos[6] = pos[0];
			
			white.exec(sgl);
			sgl.drawArrow(&pos[6]);
		}
	}
	//::glColor3f(0.0f, 0.0f, 0.0f);
}

///Shade the object in world coordinates.
void MGFSurface::shadeFS(
	mgVBO& vbo,
	const MGDrawParam& para,
	MGCL::DRAW_TARGET target
)const{
	const MGObject* obj=object_pointer();
	obj->shade(vbo,para,target);
}

///Delete the mgVBO of the i-th element.
void MGGroup::delete_displayList(
	const_iterator x
)const{
	mgVBO* vbo=m_VBO.get();
	if(vbo){
		const MGGel* geli=x->get();
		const MGAttribedGel* ageli=dynamic_cast<const MGAttribedGel*>(geli);
		if(ageli)
			vbo->deleteGel(*ageli);
	}
}

///Delete display list of the sequence [first, last).
void MGGroup::delete_displayList(
	const_iterator first, const_iterator last
)const{
	mgVBO* vbo=m_VBO.get();
	if(vbo){
		const_iterator i=first;
		for(; i!=last; i++){
			const MGGel* geli=i->get();
			const MGAttribedGel* ageli=dynamic_cast<const MGAttribedGel*>(geli);
			if(ageli)
				vbo->deleteGel(*ageli);
		}
	}
}

///Delete the mgVBO of gels_to_delete
void MGGroup::delete_displayList(
	const std::vector<const MGGel*>& gels_to_delete
)const{
	mgVBO* vbo=m_VBO.get();
	if(vbo){
		std::vector<const MGGel*>::const_iterator i,ie;
		i=gels_to_delete.begin();
		ie=gels_to_delete.end();

		for(; i!=ie; i++){
			const MGGel* geli=*i;
			const MGAttribedGel* ageli=dynamic_cast<const MGAttribedGel*>(geli);
			if(!ageli)
				continue;
			vbo->deleteGel(*ageli);
		}
	}
}

static const float edgeColor[4]={1.,.5,.5,0.};//Edge color to hilight.
static const float white[4]={1.,1.,1.,0.};//Highlight back color.
static const float endPointColor[4]={.5,1.,.5,0.};//Start/End point color to hilight.

///Highlightthe object using the display list of this object.
void MGPickObject::hilight_using_display_list(
	int line_density	///<line density to draw a surface in wire mode.
)const{
	mgVBO* nm=top_object()->dlist_name();
	nm->highlight();
}

///Highlight the object using the display list of this object.
void MGPickObjectSB::hilight_using_display_list(
	int line_density	///<line density to draw a surface in wire mode.
)const{
	MGPickObject::hilight_using_display_list(line_density);
	if(!m_vbo.is_made()){
		m_vbo.setStaticAttribColor(edgeColor);//glColor4fv(edgeColor);
		std::unique_ptr<MGCurve> e(surface()->perimeter_curve(perimeter()));
		e->drawWire(m_vbo);
	}
	m_vbo.highlight();
}

///Highlight the object using the display list of this object.
void MGPickObjectCB::hilight_using_display_list(
	int line_density	///<line density to draw a surface in wire mode.
)const{
	MGPickObject::hilight_using_display_list(line_density);
	const MGCurve* c=curve();
	MGPosition P=m_start_end ? c->end_point() : c->start_point();
	if(!m_vbo.is_made()){
		m_vbo.setStaticAttribColor(edgeColor);//glColor4fv(edgeColor);
		m_vbo.drawPoint(P);
	}
	m_vbo.highlight();
}

//////display member function.
#define NDIV 2
void MGCurve::display_arrows(mgSysGL& sgl)const{
	double ts=param_s(), te=param_e();
	MGPosition pos[4];
	for(int i = 0; i <= NDIV; i++){
		double param = (ts*(NDIV-i)+te*i)/NDIV;
		arrow(param, pos);
		sgl.drawArrow(pos);
	}
}
void MGCurve::display_break_points(mgSysGL& sgl)const{
	const MGKnotVector& t=knot_vector();
	double ts=param_s(), te=param_e();
	int k=t.order(), n=t.bdim();
	for(int i=k-1; i<=n; i++){
		double tau=t[i];
		if(i>=k && tau==t[i-1]) continue;
		if(tau<ts || tau>te) continue;
		MGVector P=eval(tau);
		sgl.drawPoint(P[0],P[1],P[2]);
	}
}
void MGLBRep::display_control_polygon(mgSysGL& sgl)const{
	const MGBPointSeq& bp=line_bcoef();
	sgl.drawPointSeq(bp);
}
void MGRLBRep::display_control_polygon(mgSysGL& sgl)const{
	MGBPointSeq bp = non_homogeneous_bcoef();
	sgl.drawPointSeq(bp);
}

void MGCurve::display_curvatures(
	mgSysGL& sgl,
	int		density,//densitiy of the graph.
	bool	use_radius,//true:radius display, false:curvature display.
	double	scaleRelative	//scaling of the graph. =1. is defalut length.
)const{
	sgl.drawCurvaGraph(*this,density,use_radius, scaleRelative);
}

//////display member function.
void MGSurface::display_arrows(mgSysGL& sgl)const{
	display_arrowsFS(sgl);
}

//Display control polygons using mgVBO::MGDrawPointSeq(sp)
void MGSBRep::display_control_polygon(mgSysGL& sgl)const{
	const MGSPointSeq& sp=surface_bcoef();
	sgl.drawPointSeq(sp);
}

void MGRSBRep::display_control_polygon(mgSysGL& sgl)const{
	MGSPointSeq sp=non_homogeneous_bcoef();
	sgl.drawPointSeq(sp);
}

//////display member function.
void MGFace::display_arrows(mgSysGL& sgl)const{
	display_arrowsFS(sgl);
}
void MGFace::display_control_polygon(mgSysGL& sgl)const{
	sgl.setLineStipple(2,0x5555);
	surface()->display_control_polygon(sgl);
}

//////display member function.
void MGShell::display_arrows(mgSysGL& sgl)const{
	int n=number_of_faces();
	for(int i=0; i<n; i++)
		face(i)->display_arrows(sgl);
}
void MGShell::display_control_polygon(mgSysGL& sgl)const{
	sgl.setLineStipple(2,0x5555);
	int n=number_of_faces();
	for(int i=0; i<n; i++)
		face(i)->display_control_polygon(sgl);
}

///Judge if the display list for vmode is made or not.
bool MGAttribedGel::displayList_is_made(MGCL::VIEWMODE vmode)const{
	mgVBO* vbo=m_VBO.get();
	if(vbo)
		return vbo->is_made(vmode);
	return false;
}

///Get the number of shading elements of m_VBO.
int MGAttribedGel::getVBOElementsNumber()const{
	mgVBO* vbo=dlist_name();
	if(!vbo)
		return 0;
	return (int)vbo->m_elements.size();
}

///Get the number of shading elements of m_VBO.
int MGAttribedGel::getVBOShaderElementsNumber()const{
	mgVBO* vbo=dlist_name();
	if(!vbo)
		return 0;
	return (int)vbo->m_elementsShade.size();
}

//Shade the object in world coordinates.
void MGSurface::shade(
	mgVBO& vbo,
	const MGDrawParam& para,
	MGCL::DRAW_TARGET target
)const{
	mgTL2Face face(para,*this);
	mgTL2Triangles tris(MGCL::XYZNormal,this);
	face.tessellate(tris);	
	vbo.drawShade(tris,target);
}

///Set no display for this vector of MGPickObject.
void MGPickObjects::setNoDisplay()const{
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		const MGPickObject& pobji=**i;
		const MGObject* obji=pobji.top_object();
		mgVBO* vbo=obji->dlist_name();
		if(vbo)
			vbo->set_no_display();
	}
}

///Set no display for this vector of MGPickObject.
void MGPickObjects::setDisplay()const{
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		const MGPickObject& pobji=**i;
		const MGObject* obji=pobji.top_object();
		mgVBO* vbo=obji->dlist_name();
		if(vbo)
			vbo->set_display();
	}
}

#else //_CONSOLE

using namespace std;
void MGObject::make_display_list(MGCL::VIEWMODE vmode)const {
	cout<< "MGObject::make_display_list, vmode="<<vmode<<endl;}
void MGPlane::make_display_list(MGCL::VIEWMODE vmode)const {
	cout << "MGPlane::make_display_list, vmode=" << vmode << endl;
}
void MGGroup::drawAttrib(mgVBO& vbo,bool no_color)const {
	cout << "MGGroup::drawAttrib, vbo="
		<< &vbo << ", co_color=" << no_color << endl;
}
void MGGroup::make_display_list(MGCL::VIEWMODE vmode)const {
	cout << "MGGroup::make_display_list, vmode=" << vmode << endl;
}
void MGPlane::display_arrows(mgSysGL& sgl)const {
	cout << "MGPlane::display_arrows, sgl=" << &sgl << endl;
}
void MGStl::shade(mgVBO& vbo, const MGDrawParam& para, MGCL::DRAW_TARGET target)const {
	cout << "MGStl::shade, vbo=" << &vbo<<", para="<<&para<<", target="<<target << endl;
}
void MGStl::drawWire(mgVBO& vbo, int line_density)const {
	cout << "MGStl::drawWire, vbo=" << &vbo <<", line_density="<< line_density << endl;
}
void MGStl::display_arrows(mgSysGL& sgl)const {
	cout << "MGStl::display_arrows, sgl=" << &sgl << endl;
}
void MGFSurface::drawWireFS(mgVBO& vbo, int line_density)const {
	cout << "MGFSurface::drawWireFS, vbo=" << &vbo
		<< ", line_density=" << line_density << endl;
}
void MGFSurface::drawWireFS_to_highlight(mgVBO& vbo, int line_density)const {
	cout << "MGFSurface::drawWireFS_to_highlight, vbo=" << &vbo
		<<", line_density=" << line_density << endl;
}
void MGFSurface::display_arrowsFS(mgSysGL& sgl, int udiv, int vdiv)const {
	cout << "MGFSurface::display_arrowsFS, sgl=" << &sgl
		<< ", udiv=" << udiv<<", vdiv="<<vdiv << endl;
}
void MGFSurface::shadeFS(mgVBO& vbo,const MGDrawParam& para,MGCL::DRAW_TARGET target
)const {
	cout << "MGFSurface::shadeFS, vbo=" << &vbo
		<< ", para=" << &para << ", target=" << target << endl;
}
void MGGroup::delete_displayList(const_iterator x)const {
	cout << "MGGroup::delete_displayList, const_iterator"<< endl;
}
void MGGroup::delete_displayList(const_iterator first, const_iterator last)const {
	cout << "MGGroup::delete_displayList, const_iterator(first, last)" << endl;
}
void MGGroup::delete_displayList(
	const std::vector<const MGGel*>& gels_to_delete
)const {
	cout << "MGGroup::delete_displayList, gels_to_delete="<<&gels_to_delete << endl;
}
void MGPickObject::hilight_using_display_list(int line_density)const {
	cout << "MGPickObject::hilight_using_display_list, line_density=" << line_density << endl;
}
void MGPickObjectSB::hilight_using_display_list(int line_density)const {
	cout << "MGPickObjectSB::hilight_using_display_list, line_density=" << line_density << endl;
}
void MGPickObjectCB::hilight_using_display_list(int line_density)const {
	cout << "MGPickObjectCB::hilight_using_display_list, line_density=" << line_density << endl;
}
void MGCurve::display_arrows(mgSysGL& sgl)const {
	cout << "MGCurve::display_arrows, sgl=" << &sgl << endl;
}
void MGCurve::display_break_points(mgSysGL& sgl)const {
	cout << "MGCurve::display_break_points, sgl=" << &sgl << endl;
}
void MGLBRep::display_control_polygon(mgSysGL& sgl)const {
	cout << "MGLBRep::display_control_polygon, sgl=" << &sgl << endl;
}
void MGRLBRep::display_control_polygon(mgSysGL& sgl)const {
	cout << "MGRLBRep::display_control_polygon, sgl=" << &sgl << endl;
}
void MGCurve::display_curvatures(
	mgSysGL& sgl,int density,bool use_radius,double scaleRelative)const {
	cout << "MGCurve::display_curvatures, sgl=" << &sgl
		<<", dinsity="<<density<<", use_radius="<< use_radius 
		<<", scaleRelative="<< scaleRelative << endl;
}
void MGSurface::display_arrows(mgSysGL& sgl)const {
	cout << "MGSurface::display_arrows, sgl=" << &sgl << endl;
}
void MGSBRep::display_control_polygon(mgSysGL& sgl)const {
	cout << "MGSBRep::display_control_polygon, sgl=" << &sgl << endl;
}
void MGRSBRep::display_control_polygon(mgSysGL& sgl)const {
	cout << "MGRSBRep::display_control_polygon, sgl=" << &sgl << endl;
}
void MGFace::display_arrows(mgSysGL& sgl)const {
	cout << "MGFace::display_arrows, sgl=" << &sgl << endl;
}
void MGFace::display_control_polygon(mgSysGL& sgl)const {
	cout << "MGFace::display_control_polygon, sgl=" << &sgl << endl;
}
void MGShell::display_arrows(mgSysGL& sgl)const {
	cout << "MGShell::display_arrows, sgl=" << &sgl << endl;
}
void MGShell::display_control_polygon(mgSysGL& sgl)const {
	cout << "MGShell::display_control_polygon, sgl=" << &sgl << endl;
}
void MGSurface::shade(mgVBO& vbo,const MGDrawParam& para,MGCL::DRAW_TARGET target)const {
	cout << "MGSurface::shade, vbo=" << &vbo<<", para="<<para<<", target="<<target << endl;
}
void MGPickObjects::setNoDisplay()const {
	cout << "MGPickObjects::setNoDisplay" << endl;
}
void MGPickObjects::setDisplay()const {
	cout << "MGPickObjects::setDisplay" << endl;
}

#endif //_CONSOLE
