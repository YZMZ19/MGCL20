/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#pragma once

#include "mg/Position.h"
#include "mgGL/ConstructionPlane.h"
#include <iosfwd>

class MGOfstream;
class MGIfstream;
class MGBox;

/** @file */
/** @addtogroup DisplayHandling
 *  @{
 */

///MGglViewAttrib is a class to serialize MGOpenGLView.

///MGglViewAttrib is used by MGOpenGLView and MGContext to define
///view-mode, construction plane's geometry data, and viewing transformation
///matrix data.
class MG_DLL_DECLR MGglViewAttrib{

public:
	
///Debug Function.
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream& out, const MGglViewAttrib& atr);

/// Serialization fucntion.
MG_DLL_DECLR friend MGOfstream& operator<< (MGOfstream& buf, const MGglViewAttrib& atr);
MG_DLL_DECLR friend MGIfstream& operator>> (MGIfstream& buf, MGglViewAttrib& atr);

///Ctor.
MGglViewAttrib(bool is_perspective=true);
MGglViewAttrib(const MGBox& box, const MGColor* gridColors=0);

///return the center of the box of this view.
const MGPosition& center()const{return m_center;};
MGPosition& center(){return m_center;};

///Compute the viewing environment the parameter box.
///compute_viewing_environment() uses (eye_position,view_up_vector) as input.
///They must be set before compute_viewing_environment.
///eye_position is used only to get the direction from the origin.
void compute_viewing_environment(
	const MGPosition& center,
	double diameter
);

///Compute the viewing environment the parameter box.
///compute_viewing_environment() uses (eye_position,view_up_vector) as input.
///They must be set before compute_viewing_environment.
void compute_viewing_environment(
	const MGBox& box///<Input the box of the target scene.
);

///Get the construction plane.
const MGConstructionPlane& cplane()const{return m_cplane;};
MGConstructionPlane& cplane(){return m_cplane;};

///Obtainthe the diameter of the sphere that surround the whole model.
double diameter()const{return m_diameter;};
double& diameter(){return m_diameter;};

///Get the eye position.
const MGPosition& eye_position()const{return m_eyeP;}
MGPosition& eye_position(){return m_eyeP;}

///Return if this is a perspective view or not.
bool is_perspective() const{return m_perspective;};

///Set the center of the objects.
void set_center(const MGPosition& pos){ m_center=pos;}

///Set m_eyeP and m_up_vector(the eye position and view-up-vector of the OpenGL).
void setEyePositionUpVector(
	const MGPosition& eyeP,
	const MGVector& upVector
);

///SDet the angle of top and bottom viewing pyramid in degrees.
void set_fovy(double fovy){m_fovy=fovy;};

///Set home matrix.
void setHomeMatrix();

///Set if this view is a perspective view(true), or orthnormal view(falsle).
void set_perspective(bool pers, double fovy=45.);

///Set the view-up vector.
void set_view_up_vector(const MGVector& up){m_up_vector=up;};

///Get view mode.
MGCL::VIEWMODE viewMode()const{return m_viewMode;};

///Get the view up vector.
const MGVector& view_up_vector()const{return m_up_vector;};
MGVector& view_up_vector(){return m_up_vector;};

///compute the view volume far.
double view_volume_far()const{return m_far;}

///compute the view volume height.
double view_volume_height()const{return m_diameter/m_scale;}

///compute the view volume near.
double view_volume_near()const{return m_near;}

private:
/// アトリビュート
	MGCL::VIEWMODE m_viewMode;
	MGConstructionPlane m_cplane;///<construction plane;

	bool m_perspective;	///<Indicate if this is perspective or orthographic view.
						///<true if perspective.
	double m_fovy;		///<angle of top and bottom viewing pyramid in degrees.
						///<See gluPerspective's fovy.

	double m_near, m_far;///<Viewing frustum's near, far, and
	MGPosition m_eyeP;	///<eye data of gluLookAt.
	MGVector m_up_vector;///<up vector data of gluLookAt.

	///m_center and m_diameter define a sphere whose center and diameter are m_center and
	///m_diameter. The sphere includes the whole model.
	MGPosition m_center;///<World coordinate of the center of the document. 
	double m_diameter;	///<diameter of the sphere that sorround the model.

	double m_scale;		///<Current scaling factor.
	double m_cx, m_cy;	///<center of the screen in world coordinate
			///<(when the center of the screen is supposed to be (0.,0.).
			///<or (-m_cx, -m_cy) is the current panning translation distance.


	glm::mat4 m_PreCenterMat;///Model matrix to apply before translate to m_center, used by get_model_matrix().
	glm::mat4 m_modelViewMat;

friend class MGOpenGLView;
};

///Initial scale for a construction plane that controls how many times bigger
///the plane area to be than the box area.
#define INITIAL_SCALE 6.

/** @} */ // end of DisplayHandling group
