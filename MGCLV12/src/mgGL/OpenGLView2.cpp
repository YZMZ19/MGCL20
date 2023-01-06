/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/drawParam.h"
#include "mgGL/Context.h"
#include "mgGL/OpenGLView.h"
#include "mgGL/SysGLList.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//Attach to Command Drawer.
//Function's return value is the number of drawer's attached after attached.
int MGOpenGLView::attach_drawer(
	mgVBO* drawer
	,bool common///<if true, drawer is attachd as common picture for all the sibling window.
		///<if flase, drawer is attached to this specific pictures.
){
	std::list<mgVBO*>* vbos;
	if(common)
		vbos=&m_command_drawersCommon;
	else
		vbos=&m_command_drawersSpecific;

	vbos->push_back(drawer);
	return (int)vbos->size();
}

//Detach from the Command Drawer..
void MGOpenGLView::detach_drawer(
	mgVBO* drawer
	,bool common///<if true, drawer is attachd as common picture for all the sibling window.
		///<if flase, drawer is attached to this specific pictures.
){
	std::list<mgVBO*>* vbos;
	if(common)
		vbos=&m_command_drawersCommon;
	else
		vbos=&m_command_drawersSpecific;
	vbos->remove(drawer);
}

///Set line density for a surface to draw in wire mode.
void MGOpenGLView::set_line_density(int line_density){
	draw_param().set_line_density(line_density);
}

//Set line density and tessellation parameter of draw_param.
void MGOpenGLView::setDrawParam(const MGContext& ctx){
	draw_param()=MGDrawParam(ctx,span_length());
}

///Set m_eyeP and m_up_vector(the eye position and view-up-vector of the OpenGL).
void MGOpenGLView::setEyePositionUpVector(
	const MGPosition& eyeP,
	const MGVector& upVector
){
	m_viewAttrib.setEyePositionUpVector(eyeP,upVector);
}

///Import draw attrib data from MGContext.
///Draw attribs are colors(Bcolor, Gcolor, Hcolor),
///line approximation smoothness, pick aperture, and MGDrawParam.
void MGOpenGLView::importDrawAttribFromContext(
	const MGContext& ctx///<Target context.
){
	setBcolor(ctx.Bcolor());
	setGcolor(ctx.Gcolor());
	mgVBOElement::setHilightColor(ctx.Hcolor());
	set_smooth(ctx.smooth());
	set_pick_aperture(ctx.pick_aperture());
	setDrawParam(ctx);
}

///Import MGContext's mgGLViewAttrib data, which are
///view mode, Construction Plane, and viewing transfromation data.
void MGOpenGLView::importViewAttribFromContext(
	const MGContext& ctx///<Target context.
){
	const MGglViewAttrib& vw=ctx.theView();
	copy(vw);	
}

///Import construction plane's grid data from a MGContext, which are
///  1)(u,v) grid numbers, 2)(u,v) grid spans. Colors are not imported.
void MGOpenGLView::importGridAttrib(const MGContext& ctx){
	MGConstructionPlane& cpl=cplane();
	cpl.importGridAttrib(ctx);
}

///Import MGOpenGLView data from MGContext.
void MGOpenGLView::import_context(
	const MGContext& ctx
){
	importDrawAttribFromContext(ctx);
	if(!get_parent_OpenGLView())
		importViewAttribFromContext(ctx);
}

///Set view mode
void MGOpenGLView::setViewMode(MGCL::VIEWMODE vmode){
	if(vmode==MGCL::DONTCARE)
		vmode=MGCL::WIRE_AND_SHADINGVIEW;
	else if(vmode==MGCL::HIGHLIGHT)
		vmode=MGCL::WIREVIEW;
	m_viewAttrib.m_viewMode=vmode;
}

///set the smooth factor of this view.
void MGOpenGLView::set_smooth(float smooth){
	MGOpenGLView* pglv=get_parent_OpenGLView();
	if(!pglv)
		pglv=this;
	pglv->m_smooth=smooth;
}

///Get the smooth factor of this view.
float MGOpenGLView::smooth()const{
	const MGOpenGLView* pglv=get_parent_OpenGLView();
	if(!pglv)
		pglv=this;
	return pglv->m_smooth;
}

///Get the pick aperture.
float MGOpenGLView::pick_aperture()const{
	const MGOpenGLView* pglv=get_parent_OpenGLView();
	if(!pglv)
		pglv=this;
	return pglv->m_pick_aperture;
}

void MGOpenGLView::set_pick_aperture(double pick_aperture){
	set_pick_aperture((float)pick_aperture);
}
void MGOpenGLView::set_pick_aperture(float pick_aperture){
	MGOpenGLView* pglv=get_parent_OpenGLView();
	if(!pglv)
		pglv=this;
	pglv->m_pick_aperture=pick_aperture;
}

//Set background color;
void MGOpenGLView::setBcolor(const MGColor& color){
	MGOpenGLView* pglv=get_parent_OpenGLView();
	if(!pglv)
		pglv=this;

	pglv->m_Bcolor=color;
}

//Set default object color;
void MGOpenGLView::setGcolor(const MGColor& color){
	MGOpenGLView* pglv=get_parent_OpenGLView();
	if(!pglv)
		pglv=this;

	pglv->m_Gcolor=color;
}

const MGColor& MGOpenGLView::Bcolor()const{
	const MGOpenGLView* pglv=get_parent_OpenGLView();
	if(!pglv)
		pglv=this;
	return pglv->m_Bcolor;
}

const MGColor& MGOpenGLView::Gcolor()const{
	const MGOpenGLView* pglv=get_parent_OpenGLView();
	if(!pglv)
		pglv=this;
	return pglv->m_Gcolor;
}

///Enable when bEnabled=true, else disable grid snap.
///When bEnabled=true, cpnlane is enbaled.
void MGOpenGLView::enable_grid_snap(bool bEnabled) {
	MGConstructionPlane& cpl = cplane();
	if (!cpl.valid())
		return;

	if (bEnabled) {
		cpl.set_bind_to_grid_enable();
	}
	else {
		cpl.set_bind_to_grid_disable();
	}
}

// grid snap is enabled
bool MGOpenGLView::is_enabled_grid_snap()const {
	if (cplane().disabled()) return false;
	return cplane().is_bind_to_grid();
}

const MGDrawParam& MGOpenGLView::draw_param()const {
	return mgVBOElement::getDrawParam();
}
MGDrawParam& MGOpenGLView::draw_param() {
	return mgVBOElement::getDrawParam();
}

//get viewport of OpenGL.
void MGOpenGLView::get_viewport(
	int vp[4]
	///<(vp[0],vp[1]) is (left,bottom) coordinates.
	///<(vp[2],vp[3]) is (width, height) of the viewport.
)const {
	for (int i = 0; i < 4; i++) {
		vp[i] = m_viewPort[i];
	}
}

///get window(m_width, m_height);
void MGOpenGLView::get_window(int& width, int& height)const {
	width = m_width;
	height = m_height;
}
