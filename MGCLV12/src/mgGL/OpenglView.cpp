/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
// MGOpenGLView.cpp : インプリメンテーション ファイル
//
#include "StdAfx.h"
#include "mg/Tolerance.h"
#include "mg/Position.h"
#include "mg/AttribedGel.h"
#include "mg/DnameControl.h"
#include "mg/CParam_list.h"
#include "mg/CSisect.h"
#include "mg/Group.h"
#include "mgGL/Context.h"
#include "mgGL/glslprogram.h"
#include "mgGL/OpenGLView.h"
#include "mgGL/GLAttrib.h"
#include "mgGL/SysGLList.h"
#include "mgGL/glViewAttrib.h"
#include "mgGL/VBO.h"

using namespace std;

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

namespace{
	int DISPNAME_BASE_ID=1;
}
unsigned OpenGLStartDisplayName(){
	return DISPNAME_BASE_ID;
}


/////////////////////////////////////////////////////////////////////////////
// MGOpenGLView

MGOpenGLView::MGOpenGLView(
	bool perspective	//indicates if the view is pespective or not.
):m_width(m_viewPort[2]),m_height(m_viewPort[3]),
m_viewAttrib(perspective){
	draw_param().set_span_length(span_length());
}

//Construct from MGglViewAttrib.
MGOpenGLView::MGOpenGLView(
	const MGglViewAttrib& glatr
):m_width(m_viewPort[2]),m_height(m_viewPort[3]),
m_viewAttrib(glatr){
	draw_param().set_span_length(span_length());
}

//Copy the informations of glview2 into this.
//Data that is not copied from glview2 are:
//m_sysgllist must be made by invoking openGL's
//display list generation functions.
void MGOpenGLView::copy(const MGOpenGLView& glview2){
	m_Bcolor=glview2.m_Bcolor;
	m_Gcolor=glview2.m_Gcolor;
	m_smooth=glview2.m_smooth;
	m_pick_aperture=glview2.m_pick_aperture;
	m_viewAttrib=glview2.m_viewAttrib;
	setLookAtMat();
}
void MGOpenGLView::copy(const MGglViewAttrib& glatr){
	m_viewAttrib=glatr;
	setLookAtMat();

	const MGVector& upV=glatr.view_up_vector();
	MGUnit_vector XAxis=upV*eye_position();
	for(int i=0; i<3; i++){
		m_XAxis_current[i]=(float)XAxis[i];
		m_center_current[i]=(float)glatr.m_center[i];
		m_up_vector_current[i]=(float)upV[i];
	}
	draw_param().set_span_length(span_length());
}

//Return display list name.
mgVBO* MGOpenGLView::display_list(){
	if(!has_parent_OpenGLView())
		return m_display_list;
	return get_parent_OpenGLView()->display_list();
}
//Draw the scene defined in this view including the current objects as hilighted.
void MGOpenGLView::drawScene(const MGPickObjects* pobjs){
	if(m_width <= 0 || !has_display_list())
		return;

	setupDrawEnv(Bcolor());
	int& x = m_viewPort[0]; int& y = m_viewPort[1];//(left, bottom)
	glViewport(x, y, m_width, m_height);

	mgGLSLProgram* glsl=mgGLSLProgram::getCurrentGLSLProgram();
	glsl->setFuncType(mgGLSL::standard);

//1. Construction plane drawing
	execDefaultStaticAttrib();
	MGConstructionPlane& cpl=m_viewAttrib.cplane();
	if(cpl.enabled())
		cpl.draw();
	
//2. Target objects drawing
	execDefaultStaticAttrib();
	MGOpenGLView* gltarget= m_parent_glView ? m_parent_glView:this;
	gltarget->m_display_list->draw(viewMode());

//3. the system generated display drawing.
	execDefaultStaticAttrib();
	gltarget->m_sysgllist.draw_list(viewMode());

//4.draw command specific pictures of parent OpenGLView.
	drawCommandDrawer(gltarget->m_command_drawersCommon);

//5.draw command specific pictures of this OpenGLView.
	drawCommandDrawer(m_command_drawersSpecific);

//6. current objects highlighting.
	highlight(pobjs);
	glFinish();
}

///Draw std::list<mgVBO*>.
void MGOpenGLView::drawCommandDrawer(
	std::list<mgVBO*>& drawers
){
	execDefaultStaticAttrib();
	std::list<mgVBO*>::iterator i=drawers.begin(), ie=drawers.end();
	for(; i!=ie; i++){
		mgVBO* drawer=*i;
		if(drawer)
			drawer->draw();
	}
}

//Highlight pobjs.
void MGOpenGLView::highlight(const MGPickObjects* pobjsP){
	if(!pobjsP)
		return;

	int ld=line_density();
	GLboolean depthEnabled=glIsEnabled(GL_DEPTH_TEST);
	GLboolean blendEnabled=glIsEnabled(GL_BLEND);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);

	const MGPickObjects& pobjs=*pobjsP;
	size_t nHL=pobjs.size();
	for(size_t i=0; i<nHL; i++){
		pobjs[i].hilight_using_display_list(ld);
	}

	if(depthEnabled)
		glEnable(GL_DEPTH_TEST);
	if(blendEnabled)
		glEnable(GL_BLEND);
}

//Get the surface parameter value uv(u,v) where screen coordinate (sx,sy) is projected on.
//If no projection points are found, the nearest point to a perimeter of surf will
//be returned.
bool MGOpenGLView::get_surface_parameter_glv(
	const MGFSurface& surf,
	int sx, int sy,	//Screen coordinates. (left, bottom) is (0,0).
	MGPosition& uv	//surface parameter (u,v) where (sx,sy) is projected on will be returned.
)const{
	MGStraight sl;
	unproject_to_sl_glv(sx,sy,sl);
	MGCSisects csis=surf.isectFS(sl);
	if(csis.size()){
		auto is=static_cast<const MGCSisect*>(csis.front().get());
		uv=is->param_surface();
		return true;
	}
	uv=surf.closest_on_boundary(sl);
	return false;
}

//Test if this has a display list to draw.
bool MGOpenGLView::has_display_list()const{
	if(m_display_list)
		return true;
	if(!has_parent_OpenGLView())
		return false;
	return get_parent_OpenGLView()->has_display_list();
}

//locate the screen coordinate (x,y) in the  3D world coordinate.
//(x, y)'s origin is (left, bottom) of the screen.
MGPosition MGOpenGLView::locate_glv(int x, int y, MGPosition* uv)const{
	MGStraight sl;
	unproject_to_sl_glv(x,y,sl);
	const MGConstructionPlane& pl=cplane();
	MGPosition xyz,uv2;
	xyz=pl.locate(sl,uv2);
	if(uv)
		*uv=uv2;
	return xyz;
}

//construct the construction plane along with its display list.
void MGOpenGLView::make_construction_plane(
	const MGPosition& mid,	//center of the construction plane.
	const MGVector& uderi,	//u-axis vector of the construction plane.
	const MGVector& vderi,	//v-axis vector of the construction plane.
	double uspan,			//span length between the lines along u-axis.
	double vspan,			//span length between the lines along v-axis.
	int ulnum,			//number of lines to draw along u-axis.
	int vlnum			//number of lines to draw along v-axis.
){
	m_viewAttrib.m_cplane.set_grid_data(MGPlane(uderi,vderi,mid), uspan,vspan,ulnum,vlnum);
}

//Make openGL display list in this glview.
void MGOpenGLView::make_display_list(const MGContext& ctx,const MGGroup& grp){
//Generate OpenGL display list.
	if(has_parent_OpenGLView())
		return;
	setDrawParam(ctx);
	m_display_list=grp.dlist_name();
	grp.make_display_list(viewMode());
}


//Rotate the current view by the angle along the vector(x,y,z),
//performs a counterclockwise rotation of angle angle about
//the vector from the origin through the point (x, y, z).
///The rotation matrix made is stored in m_PreCenterMat.
void MGOpenGLView::rotate(float angle, float x, float y, float z){
	const MGPosition& c=center();
	glm::vec3 cntr((float)c[0],(float)c[1],(float)c[2]);
	glm::mat4& preM=m_viewAttrib.m_PreCenterMat;
	preM=glm::translate(preM,cntr);
	preM=glm::rotate(preM,angle,glm::vec3(x,y,z));//Rotate around center.
	preM=glm::translate(preM,-cntr);
}

///Rotate the current view by the angle along vector around center.
///Performs a counterclockwise rotation.
///The rotation matrix made is stored in m_PreCenterMat.
void MGOpenGLView::rotate(float angle, float vector[3], const MGPosition& center){
	glm::vec3 vec(vector[0],vector[1],vector[2]);
	glm::vec3 cntr((float)center[0],(float)center[1],(float)center[2]);
	glm::mat4& preM=m_viewAttrib.m_PreCenterMat;
	preM=glm::translate(preM,cntr);
	preM=glm::rotate(preM,angle,vec);//Rotate around center.
	preM=glm::translate(preM,-cntr);
}

///Rotate the current view
///by the angle[0] around m_up_vector_current(Y-direction of the screen),
///and by the angle[1] around m_XAxis_current(view-up-vector*eye-vector=X-direction of the screen).
///Rotation is performed around m_center_current.
///The rotation matrix made is stored in m_PreCenterMat.
void MGOpenGLView::rotate(const float angle[2]){
	glm::mat4& preM=m_viewAttrib.m_PreCenterMat;
	double scl = m_viewAttrib.m_scale*.2;//.2 is adequate?
	preM=glm::translate(preM,m_center_current);
	preM=glm::rotate(preM,float(angle[1]/scl),m_XAxis_current);//Rotate around center.
	preM=glm::rotate(preM,float(angle[0]/scl),m_up_vector_current);//Rotate around center.
	preM=glm::translate(preM,-m_center_current);
}

void MGOpenGLView::set_center_current(int x, int y){
	MGStraight sl1, sl2;
	unproject_to_sl_glv(x,y,sl1);

	const MGPosition& C=center();
	MGPosition centrOnCursor=sl1.eval(sl1.param(C));
	for(int i=0; i<3; i++)
		m_center_current[i]=(float)centrOnCursor[i];

	unproject_to_sl_glv(x,y+200,sl2);
	MGPosition centerUp=sl2.eval(sl2.param(C));
	MGVector upV(3);
	for(int i=0; i<3; i++){
		upV(i)=centerUp[i]-m_center_current[i];
	}
	upV.set_unit();
	for(int i=0; i<3; i++){
		m_up_vector_current[i]=(float)upV[i];
	}

	const MGVector& eyeP=sl1.direction();
	MGVector X=upV*eyeP.normalize();//Screnn's X-direction
	X.set_unit();
	for(int i=0; i<3; i++)
		m_XAxis_current[i]=(float)X[i];
}

///Scale the current view by the factor and the view center.
void MGOpenGLView::scale(
	double factor,	///< saling factor.
	int* center///< scaling center of the screen coordinate (sx, sy) on the view plane,
		///where sx=center[0], sy=center[1] if center!=0.
){
	if(center){	
		double wx,wy;
		screen_to_world(m_viewPort +2,double(center[0]),double(center[1]),wx,wy);
		m_viewAttrib.m_cx=wx-(wx-m_viewAttrib.m_cx)/factor;
		m_viewAttrib.m_cy=wy-(wy-m_viewAttrib.m_cy)/factor;
	}
	m_viewAttrib.m_scale*=factor;
}

void MGOpenGLView::set_window(int width, int height){
	m_width=width;m_height=height;
}

//Set the parent MGOpenGLView.
void MGOpenGLView::set_parent_OpenGLView(MGOpenGLView* parent){
	m_parent_glView=parent;
}

//Convert the screen coordinate (sx, sy) to world coordinate (wx, wy) on the 
//view plane.
void MGOpenGLView::screen_to_world(
	int wh[2],	//width(wh[0]) and height(wh[1]) of the screen.
	double sx,double sy, double& wx, double& wy
)const{
	double sheight=double(wh[1]), swidth=double(wh[0]);
	double wsratio=view_volume_height()/sheight;
	wx=m_viewAttrib.m_cx+wsratio*(sx-swidth*.5);
	wy=m_viewAttrib.m_cy+wsratio*(sy-sheight*.5);
}

//Translate the current view by (dx, dy).
void MGOpenGLView::translate(double dx, double dy){
	m_viewAttrib.m_cx-=dx/m_viewAttrib.m_scale;
	m_viewAttrib.m_cy-=dy/m_viewAttrib.m_scale;
}

//Translate the current view by (dx, dy) without current scale.
void MGOpenGLView::translate_without_scale(double dx, double dy){
	m_viewAttrib.m_cx-=dx; m_viewAttrib.m_cy-=dy;
}

//Translate and scale the current view.
//(x0, y0) to (x1,y1) is the rectangle of screen coordinate whose origin is
//(left,bottom).
void MGOpenGLView::pan_zoom(int x0, int y0, int x1, int y1){
	int dx=x1-x0, dy=y1-y0;
	if(dx==0 && dy==0) return;

	if(dx<0) dx*=-1;if(dy<0) dy*=-1;
	double sdy=double(dy);

	double sheight=double(m_height), swidth=double(m_width);
	if(dx){
		double sdx=double(dx);
		double aspect=sheight/swidth;
		if(sdy/sdx < aspect) sdy=sdx*aspect;
	}

	double sxm=(double(x0)+double(x1))*.5;
	double sym=(double(y0)+double(y1))*.5;
	screen_to_world(m_viewPort+2,sxm,sym,m_viewAttrib.m_cx,m_viewAttrib.m_cy);
	double wsratio=view_volume_height()/sheight;
	m_viewAttrib.m_scale=diameter()/(wsratio*sdy);
}

//Translate and scale the current view.
//box is world coordinate's box cube.
void MGOpenGLView::pan_zoom(const MGBox& box){
	glm::mat4 mmat; get_model_matrix(mmat);

	MGPosition wc=box.mid();
	MGPosition sc; project(wc,sc,&mmat,0);

	//Set the viewing center position.
	screen_to_world(m_viewPort+2,sc[0],sc[1],m_viewAttrib.m_cx,m_viewAttrib.m_cy);

	//Set the scaling factor.
	float* mp=&mmat[0][0];
	double p=mp[2]*wc[0]+mp[6]*wc[1]+mp[10]*wc[2]+mp[14];
	p*=-.9;//was -1.0. Maybe scaled to too large object.
	double len=view_volume_near()*box.len()/p;
	m_viewAttrib.m_scale=diameter()/len;
}

//Convert the windows screen coordinate (x,y) to MGCL's straight line.
//and get the intersection of the straight line and the construction plane.
//The origin of the screen coordinate is left, bottom. Not left, top.
void MGOpenGLView::unproject(
	int x, int y,	//screen coordinate whose origin is (left, bottom).
	MGStraight& sl,	//The straight line of (x,y) will be returnred.
	MGCSisect& is	//the intersectio of the sl and the construction plane
					//will be returned.
)const{
	unproject_to_sl_glv(x,y,sl);
	if(m_viewAttrib.m_cplane.valid())
		sl.relation(m_viewAttrib.m_cplane.plane(), is);
	else
		is=MGCSisect(sl.root_point(), 0.,MGPosition(0.,0.));
}

//Convert the windows screen coordinate (x,y) to MGCL's straight line.
//The origin of the screen coordinate is left, bottom. Not left, top.
void MGOpenGLView::unproject_to_sl_glv(int x, int y, MGStraight& sl)const{
	glm::mat4 modelMatrix, projMatrix;
	get_projection_matrix(m_viewPort,projMatrix);
	get_model_matrix(modelMatrix);

	glm::ivec4 viewport(m_viewPort[0],m_viewPort[1],m_viewPort[2],m_viewPort[3]);
	glm::vec3 winPoint(x,y,0.);
	glm::vec3 P=glm::unProject(winPoint,modelMatrix,projMatrix,viewport);
	MGPosition origin(P[0],P[1],P[2]);
	winPoint[2]=100.;
	glm::vec3 Q=glm::unProject(winPoint,modelMatrix,projMatrix,viewport);
	MGPosition point(Q[0],Q[1],Q[2]);

	sl=MGStraight(MGSTRAIGHT_TYPE::MGSTRAIGHT_UNLIMIT, point-origin,origin);
	const MGConstructionPlane& pl=cplane();
	MGPosition xyz,uv2;
	xyz=pl.locate(sl,uv2);
	sl.update_root(xyz);
}

void MGOpenGLView::get_near_position(
	const MGCurve* crv, 
	const float center[2],///<screen coordinates whose origin is (left, bottom).
	double& t	//parameter value of the curve crv near to (sx,sy) will be returned.
){
	MGStraight sl;

	int sx=(int)center[0], sy=(int)center[1];
	unproject_to_sl_glv(sx,sy,sl);
	std::unique_ptr<MGCurve> crv2(crv->clone());
	const MGPosition rp=sl.root_point();
	*crv2-=rp;
	MGMatrix mat; mat.set_axis(sl.direction(),2);
	*crv2*=mat;
	t=crv2->closest2D(MGDefault::origin_2D());
}

//Initialize the viewing environment.
//When the eye position is necessary to update, setEyePositionUpVector() must be invoked
//before initializeViewingEnvironmentByBox since the current eye position direction is used.
//Initialization is done by the parameter box.
void MGOpenGLView::initializeViewingEnvironmentByBox(
	const MGBox& box
) {
	m_viewAttrib.compute_viewing_environment(box);
	setLookAtMat();
	const MGPosition& cntr = center();
	for (int i = 0; i < 3; i++)
		m_center_current[i] = (float)cntr[i];
	draw_param().set_span_length(span_length());
}

//Initialize the viewing environment to view objects projected onto a plane.
///The pespectiveness is unchanged.
//center of the object is the origin of the plane.
void MGOpenGLView::update_viewing_environment(
	const MGPlane& plane,
	double diameter//diameter of the view. This is set to m_diameter.
					///<diameter of the sphere that sorround the model.
) {
	const MGPosition& cntr = plane.root_point();
	for (int i = 0; i < 3; i++)
		m_center_current[i] = (float)cntr[i];

	if (diameter <= 0.)
		diameter = m_viewAttrib.diameter();
	MGPosition eyeP(0., 0., 1.);
	MGVector upVector(0., 1., 0.);
	setEyePositionUpVector(eyeP, upVector);
	m_viewAttrib.compute_viewing_environment(cntr, diameter);
	setLookAtMat();

	MGMatrix mat;
	mat.set_xy_axis(plane.u_deriv(), plane.v_deriv());

	glm::mat4 glmat;//double glmat[16];
	mat.convert_to_glMatrix(glmat);
	m_viewAttrib.m_modelViewMat *= glmat;//glMultMatrixd(glmat);
	m_viewAttrib.m_cplane.set_plane(plane);
}

//Set the viewing context of MGOpenGLView view to the context.
void MGOpenGLView::set_view_context(
	MGContext& ctxt
) {
	ctxt.set_line_density(line_density());
	ctxt.set_smooth(smooth());
	ctxt.set_pick_aperture(pick_aperture());
	ctxt.set_Bcolor(Bcolor());
	ctxt.set_Gcolor(Gcolor());
	ctxt.set_Hcolor(Hcolor());

	ctxt.set_view(viewAttrib());
}
