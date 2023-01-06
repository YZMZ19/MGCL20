/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#pragma once

class MGBox;
class MGGroup;
class MGObject;
class MGGLAttrib;
class MGAttribedGel;
class MGGelPositions;
class MGStraight;
class MGEdge;
class mgVBO;

#include <list>
#include <set>
#include <glm/glm.hpp>
#include "mg/Types.h"
#include "mg/PickObjects.h"
#include "mg/Position.h"
#include "mg/Straight.h"
#include "mg/Group.h"
#include "mgGL/Appearance.h"
#include "mgGL/Color.h"
#include "mgGL/sysGLList.h"
#include "mgGL/glViewAttrib.h"
#include "mgGL/ConstructionPlane.h"
#include "mgGL/glslprogram.h"

/** @addtogroup DisplayHandling
 *  @{
 */

///Defines OpenGL display class for MGCL objects.

///MGOpenGLView provides various functions to draw or render pictures
///of MGCL objects in OpenGL' windows.
class MG_DLL_DECLR MGOpenGLView{

friend class MGglViewAttrib;

public:

//////////////////////////////PUBLIC MEMBER DATA/////////////

///System display list manager.
mgSysGLList m_sysgllist;

//////////////////////////////////METHOD/////////////////////

MGOpenGLView(
	bool perspective=true	///<indicates if the view is pespective or not.
);

///Construct from MGglViewAttrib.
MGOpenGLView(
	const MGglViewAttrib& glatr
);

virtual ~MGOpenGLView()=default;

/////////////////////////////////////////////////////////////////////////////////

///Attach to VBO reference to Command Drawer list.
///Function's return value is the number of drawer's attached after attached.
int attach_drawer(
	mgVBO* drawer///<The target graphic object.
	,bool common=true///<if true, drawer is attachd as common picture for all the sibling window.
		///<if flase, drawer is attached to this specific pictures.
);

///Detach the atttached mgVBO.
void detach_drawer(
	mgVBO* drawer///<The target graphic object.
	,bool common=true///<if true, drawer is attachd as common picture for all the sibling window.
		///<if flase, drawer is attached to this specific pictures.
);

///Copy the informations of glview2 into this.
///Data that is not copied from glview2 are:
///m_sysgllist.
///m_sysgllist must be made by invoking openGL's display list generation functions.
void copy(const MGOpenGLView& glview2);

///copy the attributes in glatr into this.
void copy(const MGglViewAttrib& glatr);

///Get color values, Bcolor:background color, Gcolor:Object color,
///Hcolor:hilighted object color of this OpenGLView.
const MGColor& Bcolor()const;
const MGColor& Gcolor()const;
static const MGColor& Hcolor(){return mgVBOElement::getHilightColor();};

///return the center of the box of this view.
const MGPosition& center()const{return m_viewAttrib.center();};

const MGConstructionPlane& cplane()const{return m_viewAttrib.cplane();};
MGConstructionPlane& cplane(){return m_viewAttrib.cplane();};

///delete system display list by the function code.
bool DeleteDisplayList_by_function(int fc){
	return m_sysgllist.delete_lists_by_function_code(fc);
}

bool DeleteDisplayList_by_function_object_code(int fc,const MGGel* gel){
	return m_sysgllist.delete_lists_by_function_object_code(fc,gel);
}

///Obtainthe the diameter of the sphere that surround the whole model.
double diameter()const{return m_viewAttrib.diameter();};

const MGDrawParam& draw_param()const;
MGDrawParam& draw_param();

///Return line density for a surface to draw in wire mode.
int line_density()const{return draw_param().line_desity_wire_face();};

void execDefaultStaticAttrib();

///Draw the scene defined in this view including the current objects
///as hilighted on the curent view.
void drawScene(const MGPickObjects* pobjs=0);

///Enable when bEnabled=true, else disable grid snap.
///When bEnabled=true, cpnlane is enbaled.
void enable_grid_snap(bool bEnabled);

///Get the eye position.
const MGPosition& eye_position()const{return m_viewAttrib.eye_position();};

///get ModelView matrix of OpenGL.
void get_model_matrix(
	glm::mat4& modelMat	//double modelMat[16] ///<OpenGL's model matrix
)const;

///get projection matrix, given the viewport data 
void get_projection_matrix(
	const int vp[4],///<viewport data ={left, bottom, widht, height}
	glm::mat4& projMat	//double projMat[16]	///<OpenGL's projection matrix
)const;

/// <summary>
/// get model matrix modelM, projection matrix projM, and viewpot vp.
/// </summary>
void getModelViewProjectionMatrices(
	glm::mat4& modelM,//double modelM[16],
	glm::mat4& projM,//double projM[16],
	int vp[4]		//={left, bottom, width, height}.
)const;

///Get the surface parameter value uv(u,v) where screen coordinate (sx,sy) is projected on.
///If no projection points are found, the nearest point to a perimeter of surf will
///be returned.
bool get_surface_parameter_glv(
	const MGFSurface& surf,///<Target surface.
	int sx,	///<Screen coordinates x. (left, bottom) is (0,0).
	int sy,	///< y.
	MGPosition& uv	///<surface parameter (u,v) where (sx,sy) is projected on will be returned.
)const;

//SysGLリストのを取得する
mgSysGL* getSysGLByFunctionCode(int fc){return m_sysgllist.getSysGLByFunctionCode(fc);};

///get viewport of OpenGL.
void get_viewport(
	int vp[4]
		///<(vp[0],vp[1]) is (left,bottom) coordinates.
		///<(vp[2],vp[3]) is (width, height) of the viewport.
)const;

///get window(m_width, m_height);
void get_window(int& width, int& height)const;

///Set the parent MGOpenGLView.
MGOpenGLView* get_parent_OpenGLView(){return m_parent_glView;};
const MGOpenGLView* get_parent_OpenGLView()const{return m_parent_glView;};
bool has_parent_OpenGLView()const{return m_parent_glView!=0;};

///Test if this has a display list to draw.
bool has_display_list()const;

///Transform to the home position.
void setHomeMatrix(){m_viewAttrib.setHomeMatrix();};

///Import MGOpenGLView data from MGContext.
void import_context(
	const MGContext& ctx///<Target context.
);

///Import MGContext's mgGLViewAttrib data, which are
///view mode, Construction Plane, and viewing transfromation data.
void importViewAttribFromContext(
	const MGContext& ctx///<Target context.
);

///Import draw attrib data from MGContext.
///Draw attribs are colors(Bcolor, Gcolor, Hcolor),
///line approximation smoothness, pick aperture, and MGDrawParam.
void importDrawAttribFromContext(
	const MGContext& ctx///<Target context.
);

///Import construction plane's grid data from a MGContext, which are
///  1)(u,v) grid numbers, 2)(u,v) grid spans. Colors are not imported.
void importGridAttrib(const MGContext& ctx);

//Initialize the viewing environment.
//When the eye position is necessary to update, setEyePositionUpVector() must be invoked
//before initializeViewingEnvironmentByBox since the current eye position direction is used.
//Initialization is done by the parameter box.
void initializeViewingEnvironmentByBox(
	const MGBox& box
);

///Update the viewing environment to view objects projected onto a plane.
///The pespectiveness and the cplane are unchanged.
///center of the object is the origin of the plane.
void update_viewing_environment(
	const MGPlane& plane,///<Target plane.
	double diameter=-1.///<diameter of the view. This is set to m_diameter.
		///<diameter of the sphere that sorround the model.
		///<If diameter<=0. the current diameter is not updated.
);

//Update the center and the scale of the view.
///The pespectiveness and the cplane are unchanged.
void updateCenterScalle(
	const MGPosition& center,
	double diameter=-1.///<diameter of the view. This is set to m_diameter.
		///<diameter of the sphere that sorround the model.
		///<If diameter<=0. the current diameter is not updated.
);

bool is_enabled_grid_snap()const;

///Return if this is a perspective view or not.
bool is_perspective() const{return m_viewAttrib.is_perspective();};

///locate the screen coordinate (x,y) in the  3D world coordinate.
///(x, y)'s origin is (left, bottom) of the screen.
///In uv, the construction plane's parameter(u,v) coordinate will be returned.
MGPosition locate_glv(int x, int y, MGPosition* uv=0)const;

void make_construction_plane(
	const MGPosition& mid,	///<center of the construction plane.
	const MGVector& uderi,	///<u-axis vector of the construction plane.
	const MGVector& vderi,	///<v-axis vector of the construction plane.
	double uspan,			///<span length between the lines along u-axis.
	double vspan,			///<span length between the lines along v-axis.
	int ulnum,			///<number of lines to draw along u-axis.
	int lnum			///<number of lines to draw along v-axis.
);

///Make openGL display list of a standard mgl file in the input glview.
///A standard mgl file means that all the objects are read into a MGGroup object,
///and that the objects are the target to make display list.
///If this is not the case, make display list as follows:
///1. Name the display list as e.g. NAME, and generate the display list by glNewName.
///2. Store all the data to display in the display list NAME.
///(Sometimes glCallList may be useful.)
///At this time, use glPushName() whose list name is int(MGObject*).
///MGOpenGLView's pick() will return the MGObject* when pick() is invoked.
///3. Let MGOpenGLView know the list NAME by invoking set_display_list().
void make_display_list(const MGContext& ctx,const MGGroup& grp);

///Translate and scale the current view.
///(x0, y0) to (x1,y1) is the rectangle of screen coordinate whose origin is
///(left,bottom).
void pan_zoom(int x0, int y0, int x1, int y1);

///Translate and scale the current view.
///box is world coordinate's box cube.
void pan_zoom(const MGBox& box);

///Pick objects in the display list generated by make_display_list.
///Function's return value is MGPickObject vector.
///All the objects which were inside the pick aperture will be output.
///This data can be accessed using current_object() or current_PickObject().
///pick_glv needs a current view.
MGPickObjects pick_glv(
	const float centerAperture[4],
	///<specifies center([0], [1]) and pick aperture of x and y([2], [3]).
	const MGAbstractGels& objtype=mgAll_Object
			///<Target object kind. See MGGEL_KIND in "mg/types.h" or "mg/default.h"
);

///Determine if screen coordinate (sx,sy) is closer to the start point or to the end
///of the curve curve.
///Functin's return value is 0: if start point, 1: if end point.
int pick_start_end_glv(
	const MGCurve& curve,///<Target curve.
	int sx,	///<X of Screen coordinate. (left, bottom) is (0,0).
	int sy	///< Y.
);

///Pick a perimeter of the surface surf. That is, obtain the perimeter number
///that passes input (sx,sy) when drawn in the current view matrix.
///Function's return value is perimeter number picked.
///When no perimeters are picked, -1 will be returned.
int pick_perimeter_glv(
	const MGSurface& surf,///<Target surface.
	int sx,	///<X of Screen coordinate. (left, bottom) is (0,0).
	int sy,///< Y.
	MGPosition* uv=0,	///<surface parameter (u,v) nearest to (sx,sy) will be returned.
	float aperturex=-1.,///<specifies pick aperture of x and y.
	float aperturey=-1.///<When <=0. value is specified, default value(the value
			///<obtained by pick_aperture() will be used.
);

///Pick an edge of the face f. That is, obtain the edge number
///that passes input (sx,sy) when drawn in the current view matrix.
///Function's return value is the edge pointer picked.
///When no edges are picked, null will be returned.
const MGEdge* pick_edge_glv(
	const MGFace& f,///<Target face.
	int sx,	///<X of Screen coordinate. (left, bottom) is (0,0).
	int sy,///< Y.
	MGPosition* uv=0,	///<surface parameter (u,v) nearest to (sx,sy) will be returned.
	float aperturex=-1.,///<specifies pick aperture of x and y.
	float aperturey=-1.///<When <=0. value is specified, default value(the value
			///<obtained by pick_aperture() will be used.
);

///Get the pick aperture.
float pick_aperture()const;

///Function's return value is the number of hit objects.
int pick_to_select_buf(
	const float centrApertr[4],///<Screen coordinates. center(left, bottom)
	    ///<and (aperturex, aperturey).
	mgVBO* display_list,///<display list that includes pick objects.
	std::set<unsigned>& selected///<Selected data will be returned. This data consist of
			///<the data set by setSelectName.
);


///project world coordinates to OpenGL's screen coordinates.
///Generally, users of project are recommended to get modelMat, projlMat, or
///vp, and input them to project.
///Before use of project, SetContext() must be invoked.
void project(
	const MGPosition& world,///<Target world coordinate.
	MGPosition& screen,///<Screen coordinate will be output.
	const glm::mat4* modelMat=0,	///<OpenGL's model matrix
	const glm::mat4* projlMat=0	///<OpenGL's projection matrix
)const;

///Push back (function_code, object id) to the end or the beginning of
///system display list.
///Function's return value is the mgSysGL pointer pushed.
///sysgl can be handled(can be erased after the generation) by
/// (1)function code fc (2)object id oi(in other words, object name).
///If the system display list is to handle by the id(object name), oi must be
///meaningful. Usually oi is recommended to be the object pointer.
///If display list handling is not intended by the object id, set oi as null.
///int push_back_to_sysgl(int fc, MGGel* oi=0){return m_sysgllist.push_back(fc,oi);};
///int push_front_to_sysgl(int fc, MGGel* oi=0){return m_sysgllist.push_front(fc,oi);};
mgSysGL* push_back_to_sysgl(int fc,const MGGel* oi=0){return m_sysgllist.push_back(fc,oi);};
mgSysGL* push_front_to_sysgl(int fc,const MGGel* oi=0){return  m_sysgllist.push_front(fc,oi);};

///Push back (function_code, object id) to the end or the beginning of
///system display list. sysgl must be a newed object, and the ownership will be
///transfered to this.
void push_back_to_sysgl(mgSysGL* sysgl){m_sysgllist.push_back(sysgl);};

///Rotate the current view by the angle along the vector(x,y,z),
///performs a counterclockwise rotation of angle about
///the vector from the m_center.
///The rotation matrix made is stored in m_PreCenterMat.
void rotate(float angle, float x, float y, float z);

///Rotate the current view by the angle along vector around center.
///Performs a counterclockwise rotation.
///The rotation matrix made is stored in m_PreCenterMat.
void rotate(float angle, float vector[3], const MGPosition& center);

///Rotate the current view
///by the angle[0] around m_up_vector_current(Y-direction of the screen),
///and by the angle[1] around m_XAxis_current(view-up-vector*eye-vector=X-direction of the screen).
///Rotation is performed around m_center_current.
///(m_center_current,m_up_vector_current,m_XAxis_current) are et by set_center_current.
///The rotation matrix made is stored in m_PreCenterMat.
void rotate(const float angle[2]);

///Scale the current view by the factor and the view center.
void scale(
	double factor,	///< saling factor.
	int* center=0///< scaling center of the screen coordinate (sx, sy) on the view plane,
		///<where sx=center[0], sy=center[1] if center!=0.
);

void set_center(const MGPosition& pos){ m_viewAttrib.set_center(pos);}
void set_center_current(int x, int y);

///Set line density for a surface to draw in wire mode.
void set_line_density(int line_density=1);

//Set line density and tessellation parameter of draw_param.
void setDrawParam(const MGContext& ctx);

void set_output_dpi(float dpi){m_dpi = dpi;};

float get_output_dpi(){return m_dpi;};


///Set m_eyeP and m_up_vector(the eye position and view-up-vector of the OpenGL).
void setEyePositionUpVector(
	const MGPosition& eyeP,
	const MGVector& upVector
);

///Set the parent MGOpenGLView.
void set_parent_OpenGLView(MGOpenGLView* parent=0);

///Set background color;
void setBcolor(const MGColor& color);

void set_fovy(double fovy){m_viewAttrib.set_fovy(fovy);};

//void set_defalut_colors();

///Set default object color;
void setGcolor(const MGColor& color);

///Set hilight color;
static void setHcolor(const MGColor& color){mgVBOElement::setHilightColor(color);};

void set_pick_aperture(double pick_aperture);
void set_pick_aperture(float pick_aperture);

///set the smooth factor of this view.
void set_smooth(float smooth);

///Set the window as(0,0)-(width,height).
///This does not invoke glViewport().
void set_window(int width, int height);

///Get the smooth factor of this view.
float smooth()const;

///Get the draw span length(approximate line segment length to draw curves).
double span_length()const{return m_viewAttrib.diameter()*double(smooth())*.1;};

///Translate the current view by (dx, dy).
void translate(double dx, double dy);

///Translate the current view by (dx, dy) without current scale.
void translate_without_scale(double dx, double dy);

///Convert the windows screen coordinate (x,y) to MGCL's straight line.
///and get the intersection of the straight line and the construction plane.
///The origin of the screen coordinate is left, bottom. Not left, top.
///The direction of the sl is from the screen to the viewer.
void unproject(
	int x,///< X screen coordinate whose origin is (left, bottom).
	int y,///< Y.
	MGStraight& sl,	///<The straight line of (x,y) will be returnred.
	MGCSisect& is	///<the intersectio of the sl and the construction plane
					///<will be returned.
)const;

///Convert the windows screen coordinate (x,y) to MGCL's straight line.
///The origin of the screen coordinate is left, bottom. Not left, top.
///The direction of the sl is from the screen to the viewer.
void unproject_to_sl_glv(int x, int y, MGStraight& sl)const;

///Get the view up vector.
const MGVector& view_up_vector()const{return m_viewAttrib.view_up_vector();};

///compute the view volume far.
double view_volume_far()const{return m_viewAttrib.view_volume_far();}

///compute the view volume height.
double view_volume_height()const{return m_viewAttrib.view_volume_height();}

///compute the view volume near.
double view_volume_near()const{return m_viewAttrib.view_volume_near();}

double get_scale() {return m_viewAttrib.m_scale;};

///Compute the nearest point parameter value of crv to the point center.
///center's data is given by the screen coordinate of this current screen.
void get_near_position(
	const MGCurve* crv, ///< Target curve.
	const float center[2],///<screen coordinates whose origin is (left, bottom).
	double& t	///<parameter value of the curve crv near to (sx,sy) will be returned.
);

///Return display list name that shows all of the objects in this view.
mgVBO* display_list();

///Update display list name that shows all of the objects in this view.
void set_display_list(
	mgVBO* dlist///<The target graphic object.
){m_display_list=dlist;};

//Set the viewing context of MGOpenGLView view to the context.
void set_view_context(
	MGContext& ctxt
);

///Set view mode
void setViewMode(MGCL::VIEWMODE vmode=MGCL::WIREVIEW);

///Set up drawing environment. That is,
///(1) Invoke basic OpenGL functions(Color, Depth, etc)
///(2) Set necessary uniform variables to shader program(transform matrices, ...)
void setupDrawEnv(
	const MGColor& backColor,
	const float* centrApertu = nullptr//centrApertu = nullptr means standard draw,
	//and centrApertu != null means selection mode.
);

///Convert the screen coordinate (sx, sy) to world coordinate (wx, wy) on the 
///view plane.
void screen_to_world(
	int wh[2],	///<width(wh[0]) and height(wh[1]) of the screen.
	double sx,	///<Input screen x.
	double sy,	///<             y
	double& wx,///<World x will be output.
	double& wy///< y will be output.
)const;

const MGglViewAttrib& viewAttrib()const{return m_viewAttrib;};
MGglViewAttrib& viewAttrib(){return m_viewAttrib;};

MGCL::VIEWMODE viewMode()const{return m_viewAttrib.viewMode();};

protected:
/// アトリビュート
	int m_viewPort[4]{ 0,0,0,0 };//Currently m_viewPort[0] and [1] are fixed as 0.
	int& m_width; int& m_height;//Window's width and height. Alias of m_viewPort[2], [3].

///0.
	MGOpenGLView* m_parent_glView{ nullptr };///<If this OpenGLView has the parent, the pointer
		///<will be set. The parent means all of the display list are shared
		///<with the parent's, and on this drawScene's invocation, the parent's
		///<display list will be drawn.
	
	///Display list name of the document which is displayed on this glview.
	///m_display_list=0 indicates no display list is generated for this glview.
	mgVBO* m_display_list{ nullptr };

///1. Static atributes(usually not changed after this is initaialized)

	MGColor m_Bcolor;	///<Background color.
	MGColor m_Gcolor;	///<Object lines color.
	float m_smooth{0.001f};	///<Smoothness of the curves to draw.
	///< 1/smooth is the division number of a curve whose length is the window width.
	///< When smooth becomes small, smoothness increases.	
	float m_pick_aperture{5.f};///<Pick aperture. Number of pixels to allow picking.

	///Attributes that are saved to .mgl file.
	MGglViewAttrib m_viewAttrib;
	
	glm::mat4 m_lookAtMat;///<glm::lookAt matrix constructied from eye_position()
		///< and view_up_vector().

///2. Dynamic atributes(depend on what is current viewing environment).
	glm::vec3 m_XAxis_current;///<When set_center_current is invoked, 
						///<the current unprojected straight line of the input cursor is set.
	glm::vec3 m_center_current{ 0.f, 0.f, 0.f };///<Current center after transformations(a point on m_unprojectSL_current).
	glm::vec3 m_up_vector_current;

	float m_dpi{ 96.0f }; // output device independ value;

//Highlight pobjs.
void highlight(const MGPickObjects* pobjsP);

///Set if this view is a perspective view(true), or orthnormal view(false).
void set_perspective(bool pers,double fovy=45.){
	m_viewAttrib.set_perspective(pers,fovy);
};

private:

///Current command's picture drawer.
///MGOpenGLView will invoke the drawer last in all the darwings.
std::list<mgVBO*> m_command_drawersCommon;

///Current command's picture drawer.
///MGOpenGLView will invoke the drawer after m_command_drawersCommon.
std::list<mgVBO*> m_command_drawersSpecific;

///Extract selectionName data from the frame buffer drawn by selectionDraw();
void extractSelected(
	const int viewport[4],///Viewport of the selection target window.
	std::set<unsigned>& selected///Selected name data will be returned.
		/// This data consist of the data set by selectionDraw.
);

///Draw std::list<mgVBO*>.
void drawCommandDrawer(
	std::list<mgVBO*>& drawers
);

///Constructglm::lookAtMatrix from eye_position() and view_up_vector().
void setLookAtMat();

};

unsigned OpenGLStartDisplayName();


///////////////////////////////////////////////////////////////////////////////////////

/** @} */ // end of DisplayHandling group
