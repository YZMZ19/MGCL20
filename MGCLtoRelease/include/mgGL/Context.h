/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#ifndef _MGContext_HH_
#define _MGContext_HH_

#include "mg/Attrib.h"
#include "mgGL/Color.h"
#include "mgGL/Appearance.h"
#include "mgGL/SnapAttrib.h"
#include "mgGL/glViewAttrib.h"
#include "Tl2/TLInputParam.h"

///////////////////////////////////////////////////

class MGOpenGLView;
class MGOfstream;
class MGIfstream;

/** @addtogroup DisplayHandling
 *  @{
 */

/// MGContext defines the attributes of a document.

///MGContext includes snap attributes, colors(background, object, highlight),
///smoothness definition data to approximate spline data by a polyline,
///pick-aperture, viewing matrix data of the main view, construction plane's
///attributes, tollerances, and tessellation parameters, and other default appearance data.
class MG_DLL_DECLR MGContext: public MGAttrib{

public:

MGContext();
MGContext(const MGBox& bx);

MGContext(
	const MGSnapAttrib& snap_attrib,///<Snap data
	
	int line_density,		///<line density for a surface to draw in wire mode.
	const MGColor& Bcolor,	///<Background color.
	const MGColor& Gcolor,	///<Object lines color.
	const MGColor& Hcolor,	///<Object highlight color.
	float smooth,	///<Smoothness of the curves to draw.
		///< 1/smooth is the division number of a curve whose length is the window width.
		///< When smooth becomes small, smoothness increases.	
	float pick_aperture,///<Pick aperture. Number of pixels to allow picking.
	const MGglViewAttrib& view,///<MGglViewAttrib of the mainView of the document.
	const MGColor gridColor[4],///Construction plane's grid colors.

	const int gridNum[2], ///< Construction plane's grid number.
	const double gridSpan[2], ///< Construction plane's grid  span.
	
	const double torelance[6],///<MGCL Tolerance.
		///<[0]=wc_zero;
		///<[1]=rc_zero;
		///<[2]=mach_zero;
		///<[3]=line_zero;
		///<[4]=angle_zero;
		///<[5]=max_knot_ratio;
	const mgTLInputParam& tessellate_param,///<tessellation parameter.
	MGAppearance* appearance=0///<must be a newed object of MGAppearance.
		///<the ownership will be transfered to this MGContext.
);

///copy constructor
MGContext(const MGContext& context);

///////// Destructor.//////

~MGContext();

///Assignment operator
MGContext& operator=(const MGGel& gel2);
MGContext& operator=(const MGContext& gel2);

///comparison
bool operator<(const MGContext& gel2)const;
bool operator<(const MGGel& gel2)const;

///Generate copied gel of this gel.
///Returned is a newed object. User must delete the object.
MGContext* clone()const{return new MGContext(*this);};

///draw attribute data.
virtual void drawAttrib(
	mgVBO& vbo,///<The target graphic object.
	bool no_color=false	///<if true, color attribute will be neglected.
)const{	exec_draw_attributes(vbo);};

/// set/get snap attrib data.
const MGSnapAttrib& snap_attrib()const{return m_snap_attrib;}; //Snap data
MGSnapAttrib& snap_attrib(){return m_snap_attrib;}; //Snap data
void set_snap_attrib(const MGSnapAttrib& snap_attrib){m_snap_attrib=snap_attrib;};

///get line_density attrib data.
int line_density()const{return m_line_density;};

///set line_density attrib data.
void set_line_density(int line_density){m_line_density=line_density;};

/// set/get color data.

///Background color.
const MGColor& Bcolor()const{return m_Bcolor;};

///Background color.
MGColor& Bcolor(){return m_Bcolor;};

///Background color setter.
void set_Bcolor(const MGColor& Bcolor);

///Object color.
const MGColor& Gcolor()const{return m_Gcolor;};

///Object color.
MGColor& Gcolor(){return m_Gcolor;};

///Object color setter.
void set_Gcolor(const MGColor& Gcolor);

///Object highlight color.
const MGColor& Hcolor()const{return m_Hcolor;};

///Object highlight color.
MGColor& Hcolor(){return m_Hcolor;};

///Object highlight color setter.
void set_Hcolor(const MGColor& Hcolor);


///CPlane's grid colors.
///[0]=line, [1]=x-axis, [2]=y-axis, [3]=z-axis.
const MGColor* gridColors()const{return m_gridColor;};
const int* gridNum()const{return m_gridNum;};
const double* gridSpan()const{return m_gridSpan;};
void set_gridColors(const MGColor gridColor[4]);
void set_gridNum(const int gridNum[2]);
void set_gridSpan(const double gridSpan[2]);

// Grid color operations.
/// Grid color.
const MGColor& gridColor()const{return m_gridColor[0];};

/// Grid color.
MGColor& gridColor(){return m_gridColor[0];};

/// Grid color setter.
void set_gridColor(const MGColor& color){m_gridColor[0] = color;};

// X-axis(u-axis) color operations.
/// X-axis(u-axis) color.
const MGColor& axisXcolor()const{return m_gridColor[1];};

/// X-axis(u-axis) color.
MGColor& axisXcolor(){return m_gridColor[1];};

/// X-axis(u-axis) color setter.
void set_axisXcolor(const MGColor& color){m_gridColor[1] = color;};

// Y-axis(v-axis) color operations.
/// Y-axis(v-axis) color.
const MGColor& axisYcolor()const{return m_gridColor[2];};

/// Y-axis(v-axis) color.
MGColor& axisYcolor(){return m_gridColor[2];};

/// Y-axis(v-axis) color setter.
void set_axisYcolor(const MGColor& color){m_gridColor[2] = color;};

// Z-axis color operations.
/// Z-axis color.
const MGColor& axisZcolor()const{return m_gridColor[3];};

/// Z-axis color.
MGColor& axisZcolor(){return m_gridColor[3];};

/// Z-axis color setter.
void set_axisZcolor(const MGColor& color){m_gridColor[3] = color;};

/// set/get smooth data.
///The smooht data is used for MGOpenGLView. See the explanation.
///Smoothness of the curves to draw.
float smooth()const{return m_smooth;};	

void set_smooth(float smooth){m_smooth=smooth;};

///pick_aperture.
///The pick_aperture data is used for MGOpenGLView. See the explanation.
float pick_aperture()const{return m_pick_aperture;};
void set_pick_aperture(double pick_aperture){m_pick_aperture=float(pick_aperture);};
void set_pick_aperture(float pick_aperture){m_pick_aperture=pick_aperture;};

///Get the view attribute.
const MGglViewAttrib& theView()const{return m_view;};
MGglViewAttrib& theView(){return m_view;};

///Set the view data of the view.
void set_view(
	const MGglViewAttrib& view
);

/// Get the tolerance data.
double* tolerance(){return m_tolerance;};
const double* tolerance()const{return m_tolerance;};

///Set the tolerance data.
void set_tolerance(
	double wc_zero,
	double rc_zero,
	double mach_zero,
	double line_zero,
	double angle_zero,
	double max_knot_ratio
);

/// Get the tessellation parameter.
const mgTLInputParam& tessellate_param()const{return m_tessellate_param;};
mgTLInputParam& tessellate_param(){return m_tessellate_param;};

/// Set the tessellation parameter.
void set_tessellate_param(const mgTLInputParam& tessellate_param){
	m_tessellate_param=tessellate_param;};

///Get the tessellation paprameter maximum_edge_length.
double tess_maximum_edge_length()const{return m_tessellate_param.max_edge_len();};

///Get the tessellation paprameter surface_tolerance.
double tess_surface_tolerance()const{return m_tessellate_param.surfTol();};

///Get the tessellation paprameter curve_tolerance.
double tess_curve_tolerance()const{return m_tessellate_param.crvTol();};

/// Appearance data.
const MGAppearance* appearance()const{return m_appearance;};
MGAppearance* appearance(){return m_appearance;};

///set appearance.
void set_appearance(
	MGAppearance* appearance	///<appearance must be a newed object. The ownership will
								///<be transfered to this MGContext.
);

///Remove appearance data which is returned as the function return value.
std::unique_ptr<MGAppearance> remove_appearance();

/// Return This object's typeID
long identify_type() const{return MGCONTEXT_TID;};

/////////Attributes execution functions.///////

///Execution of MGCL tolerance. Set the MGCL tolerance of this context.
void exec_tolerance()const;

///Execution of drawing OpenGL attributes.
///Valid OpenGL rendering context must be made current.
void exec_draw_attributes(mgVBO& vbo)const;

///Execution of rendering OpenGL attributes.
///Valid OpenGL rendering context must be made current.
void exec_render_attributes(mgVBO& vbo)const;

///stream text output.
std::ostream& toString(std::ostream& ostrm) const;

///Get the name of the class.
std::string whoami()const{return "Context";};

///Read all member data.
void ReadMembers(MGIfstream& buf);

///Write all member data
void WriteMembers(MGOfstream& buf)const;

private:

	MGSnapAttrib m_snap_attrib;///<Snap data
	
	///////Folowing are MGOpenGLView attributes./////////
		int m_line_density;///<line density for a surface to draw in wire mode.
		MGColor m_Bcolor;	///<Background color.
		MGColor m_Gcolor;	///<Object lines color.
		MGColor m_Hcolor;	///<Object highlight color.
		float m_smooth;	///<Smoothness of the curves to draw.
			///< 1/smooth is the division number of a curve whose length is the window width.
			///< When smooth becomes small, smoothness increases.	
		float m_pick_aperture;///<Pick aperture. Number of pixels to allow picking.
		MGglViewAttrib m_view;///<MGglViewAttrib of the mainView of the document.
		MGColor m_gridColor[4];
			//[0]=grid color, [1]=X-axis color, [2]=Y-axis color, [3]=Z-axis color.
		int m_gridNum[2]; /// グリッドの本数
		double m_gridSpan[2];/// グリッドの間隔
	///////Above are MGOpenGLView attributes./////////

	double m_tolerance[6];///<MGCL Tolerance.
			///<[0]=wc_zero;
			///<[1]=rc_zero;
			///<[2]=mach_zero;
			///<[3]=line_zero;
			///<[4]=angle_zero;
			///<[5]=max_knot_ratio;
	mgTLInputParam m_tessellate_param;///<tessellation parameter.
	MGAppearance* m_appearance;///<Newed objects of MGAppearance.

	//void setDefaultTolerance();
};

/** @} */ // end of DisplayHandling group
#endif //#ifndef _MGContext_HH_
