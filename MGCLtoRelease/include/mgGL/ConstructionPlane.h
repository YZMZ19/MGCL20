/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
// MGConstructionPlane.h : MGConstructionPlane クラスの宣言およびインターフェイスの定義をします。
#ifndef _MGConstructionPlane_HH_
#define _MGConstructionPlane_HH_

class MGBox;
class MGPosition;
class MGContext;

#include "mg/Plane.h"
#include "mgGL/Color.h"
#include "mgGL/VBO.h"

/** @file */
/** @addtogroup DisplayHandling
 *  @{
 */

///MGConstructionPlane defines a construction plane ton input 3D data.

///MGConstructionPlane provides a local working 2D coordinate system and
///a local 3D coordinate system.
///MGConstructionPlane has the right hand coordinate system (U,V,N), where
///U=uspan*m_plane.uderiv(), V=uspan*m_plane.vderiv(), and N=m_nspan*m_plane.normal().
///m_plane.uderiv(), m_plane.vderiv(), and m_plane.normal() are set to length one
///on the construction.
///MGConstructionPlane has the cplane coordinate system such that the coordinate (x, y, z) is
///the normal world coordinate (R+x*U, R+y*V, R+z*N). The cplane coordinate conversion
///utilities are provided as convert_to_world() or convert_from_world().
class MG_DLL_DECLR MGConstructionPlane: public mgVBO{

public:

//////////Constructor/////////

///Default dtor.
MGConstructionPlane();

///Ctor of each necessary data.
MGConstructionPlane(
	double origin[3],	///<Origin's coordinate value.
	double uaxis[3],	///<A vector value of the horizontal direction.
	double vaxis[3],	///<A vector value of the vertical direction.
	double uspan,		///<span length along u axis.
	double vspan,		///<span length along v axis.
	int uline_num,		///<number of lines along u axis.
	int vline_num,		///<number of lines along v axis.
	double nspan=1.		///<span length along normal axis.
);

///Ctor of MGPlane.
MGConstructionPlane(
	const MGPlane& plane,///<construction plane.
	double uspan,		///<span length along u axis.
	double vspan,		///<span length along v axis.
	int uline_num,		///<number of lines along u axis.
	int vline_num,		///<number of lines along v axis.
	double nspan=1.		///<span length along normal axis.
);

///Bind the input point uv(this construction plane's parameter value)
///to the nearest grid point of this construction plane.
void bind_to_grid(
	const MGPosition& uv,
	MGPosition& uvout
)const;

///Change origin point.
void change_origin(
	const MGPosition& new_origin
);

///Convert cplane coordinates to the normal world coordinates.
///Function's return is the world coordinates.
MGPosition convert_to_world(const MGPosition& cplane_coord)const;

///Convert to cplane coordinates from the normal world coordinates.
///Function's return is the cplane coordinates.
MGPosition convert_from_world(const MGPosition& world_coord)const;

///Draw this plane using OpenGL.
virtual void make_display_list(MGCL::VIEWMODE vmode=MGCL::DONTCARE);

///Return if locate point on this plane should be bind to grid point or not.
///true if should be bound to grid point.
bool is_bind_to_grid()const{return m_bind_to_grid;};

///Obtain the position data of the parameter (u,v).
MGVector eval(const MGPosition& uv)const{return m_plane.eval(uv);};
MGVector eval(double u, double v)const{return m_plane.eval(u,v);};

///Get line and axis colors
void get_colors(
	MGColor& lineColor,		///<Grid line color
	MGColor& uaxisColor,	///<u axis
	MGColor& vaxisColor		///<v axis
)const;

///Obtain the grid data of this.
void get_grid_data(
	double& uspan,		///<span length along u axis.
	double& vspan,		///<span length along v axis.
	int& uline_num,		///<number of lines along u axis.
	int& vline_num,		///<number of lines along v axis.
	double& nspan		///<span length along normal axis.
)const;

///locate a point on this plane, given straight line.
///the located point will be the intersection(or nearest bound) point
///of sl and the plane.
///Function's return value is the world coordinate located.
MGPosition locate(
	const MGStraight& sl,///<input the ray straight line of the cursor.
	MGPosition& uv		///<the plane's parameter value (u,v) will be output.
)const;

///Disable the construction plane.
bool disabled()const{return m_disabled;};

///Enable the construction plane.
bool enabled()const{return !m_disabled;};

///set bind_to_grid enable or disable.
void set_bind_to_grid_enable(){m_bind_to_grid=true;};
void set_bind_to_grid_disable(){m_bind_to_grid=false;};

///Import grid data from a MGContext, which are
///  1)(u,v) grid numbers, 2)(u,v) grid spans. Colors are not imported.
void importGridAttrib(const MGContext& ctx);

///Construct grid and the plane data from box and view.
///When view_num=2(x,y):(u,v)=(X-axis, Y-axis) where (u,v) is the direction of this plane,
///     view_num=3(y,z):(u,v)=(Y-axis, Z-axis),
///     view_num=4(z,x):(u,v)=(Z-axis, X-axis).
///The colors of the axes are set according to the kind of the axis.
///X's color=gridColors[1], Y's color=gridColors[2], Z's color=gridColors[3].
///Other line color is gridColor[0].
void setGridDataByBox(
	const MGBox& box,///<The target box.
	int view_num=1,	///<Standard view number:
	///<1: 3D perspective view, 2:(x,y) 2D view, 3:(y,z) 2D view, 4:(z,x) 2D view
	///<0: non standard view.
	const MGColor* gridColors=0,///<The colors for the grid, axes. When null, colors not set.
	double spanIn=-1.///< span length of the grid. If span<=0., default value from the box is used.
);

///Compute grid data and the plane from the plane and the grid span data.
void set_grid_data(
	const MGPlane& plane,///<construction plane.
	double uspan,		///<span length along u axis.
	double vspan,		///<span length along v axis.
	int uline_num,		///<number of lines along u axis.
	int vline_num,		///<number of lines along v axis.
	double nspan=1		///<span length along normal axis.
);

///Test if this is valid or not.
bool valid()const{return m_plane.sdim()>0;};

///Get the MGPlane info.
const MGPlane& plane()const{return m_plane;};

///Get the u-direction span length.
double uspan()const{return m_uspan;};

///Get the v-direction span length.
double vspan()const{return m_vspan;};

///Get the mesh number of v-direction.
int vnum()const{return m_vnum;};

///Get the mesh number of u-direction.
int unum()const{return m_unum;};

///set line, u-axis, and v-axis colors.
void set_colors(const MGColor colors[3]);///< [0]=line, [1]=u-axis, [2]=vaxis

///set colors by VID. colors[0] is the one of grid lines.
///colors[i] are the ones for axis lines for 1<=i<=3.
///As a standard one, colors[1]=x-axis, [2]=y-axis, [3]=z-axis.
///When 2<=vid<=4, (u,v) colors are:
///vid=2:(u,v)=(x,y), vid=3:(u,v)=(y,z), vid=4:(u,v)=(z,x).
///Other vid is treated as vid=2.
void set_colorsByViewID(
	int vid,
	const MGColor colors[4]
);

///set grid line color.
void set_line_color(const MGColor& color);

///set uaxis color.
void set_uaxis_color(const MGColor& color);

///set vaxis color.
void set_vaxis_color(const MGColor& color);

///Set disable.
void set_disable(){m_disabled=true;};

///Set enable.
void set_enable(){m_disabled=false;};

///Set span length of u and v direction as span(the same).
void set_span(double span);

///Set u-direction span length.
void set_uspan(double span);

///Set v-direction span length.
void set_vspan(double span);

///Set mesh number of u and v(the same).
void set_num(int line_num);

///Set mesh number of u-direction.
void set_unum(int unum);

///Set mesh number of v-direction.
void set_vnum(int vnum);

///Set MGPlane data.
void set_plane(const MGPlane& plane);

private:
	bool m_disabled;///<true if this construction plane should not be displayed.
					///<false if this construction plane should be displayed.
	bool m_bind_to_grid;///<true if locate point should be bind to grid point of this plane.
					///<false if not.
	MGPlane m_plane;///<construction plane expression.
	double m_vspan;///<the span length of the direction m_plane.v_deriv().
	double m_uspan;///<the span length of the direction m_plane.u_deriv().
	double m_nspan;///<the span length of the direction m_plane.normal();
	int m_vnum;	///<number of lines along v axis.
	int m_unum;	///<number of lines along u axis
	MGColor m_lineColor,m_uaxisColor,m_vaxisColor;///<color of grid lines, u-axis and v-axis, each.

	////construction plane's color definitions.
	//const static float* m_lColorV;//Grid's line color except (u,v) plus axis.

///Compute cplane parameter from 3D box.
MG_DLL_DECLR friend void MGcplane_parameter(
	const MGBox& box,///<The target box.
	double& span,	///<span length will be output.
	int& lnum,	///<number of lines along vertical and horizontal will be output.
	int& sdid,	///<maxmum area coordinate pair will be output.
					///<0:(x,y), 1:(y,z), 2:(z,x)
	MGPosition& mid	///<rounded mid point will be output.
);

///Debug Function.
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream& out, const MGConstructionPlane& pln);

/// Serialization.
MG_DLL_DECLR friend MGOfstream& operator<< (MGOfstream& buf, const MGConstructionPlane& cpl);
MG_DLL_DECLR friend MGIfstream& operator>> (MGIfstream& buf, MGConstructionPlane& cpl);

};

///Compute cplane parameter from 3D box.
MG_DLL_DECLR void MGcplane_parameter(
	const MGBox& box,///<The target box.
	double& span,	///<span length will be output.
	int& lnum,	///<number of lines along vertical and horizontal will be output.
	int& sdid,	///<maxmum area coordinate pair will be output.
					///<0:(x,y), 1:(y,z), 2:(z,x)
	MGPosition& mid	///<rounded mid point will be output.
);

/** @} */ // end of DisplayHandling group
#endif
