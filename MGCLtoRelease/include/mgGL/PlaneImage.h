/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGPlaneImage_HH_
#define _MGPlaneImage_HH_

#include <memory>
#include "mg/Plane.h"
#include "mgGL/Image.h"
#include "mgGL/texture.h"

//
//Define MGPlaneImage Class.

class MGOfstream;
class MGIfstream;

/** @addtogroup GLAttrib
 *  @{
 */

///MGPlaneImage defines square image plane.

///This is displayed using texture mapping.
///The plane's root point is the (left, bottom) point of the image(m_image).
///m_uderiv and m_vderiv of the plane are always defined as unit_vectors.
///The pixel sizes are also defined as m_pixelSizeWidth and m_pixelSizeHeight.
///This implies width and hight can be different size.

class MG_DLL_DECLR MGPlaneImage:public MGPlane{

public:
	
MGImage m_image;///<Image to display on the plane.
	///The plane's root point corresponds to the (left, bottom) of m_image,
	/// and width direction of the image is uderiv() of the plane.

//////////////////Constructor//////////////////
MGPlaneImage():m_pixelSizeWidth(0.),m_pixelSizeHeight(0.),
m_rightTexCoord(0.),m_topTexCoord(0.){;};

MGPlaneImage(
	const MGVector& uderiv,///<image's width direction.
	const MGVector& vderiv,///<image's height direction.
		///vderiv is normal to uderiv. If this is not the case,
		///vderiv will be so transformed as normal to uderiv.
		///uderiv's direction is not changed.
	const MGPosition &origin,
		///(left,bottom) point of the image in the world coordinaes.
	MGImage& image,//Image data will be transfered to this.
	double pixelSizeWidth,///<one pixel size for the width direction
	double pixelSizeHeight=-1.///<one pixel size for the height direction.
		///When pixelSizeHeight<=0.,
		///pixelSizeHeight=pixelSizeWidth is assumed.
);

///Assignment
MGPlaneImage& operator=(const MGPlaneImage& gel2);

////////Member Function////////

///Obtain ceter coordinate of the geometry.
MGPosition center() const;

///Generate a newed clone object.
MGPlaneImage* clone()const;

///Compute the box from the scratch.
void compute_box(MGBox& bx) const;

///Draw the image in world coordinates using texture.
void drawWire(
	mgVBO& vbo,///<The target graphic object.
	int line_density=1	///<line density, not used.
)const{MGSurface::drawWire(vbo,line_density);};

///Shade the object in world coordinates.
void shade(
	mgVBO& vbo,///<The target graphic object.
	const MGDrawParam& para,///<Parameter to draw.
	MGCL::DRAW_TARGET target= MGCL::SHADING///<The target vbo element.
)const;

///Get image.
const MGImage& get_image()const{return m_image;};
MGImage& get_image(){return m_image;};

/// Return This object's typeID
long identify_type() const{return MGPLANEIMAGE_TID;};

//Return (width, height) of the image.
int image_width()const{return m_image.width();};
int image_height()const{return m_image.height();};

///Evaluate left_bottom point.
MGPosition left_bottom()const;

///Evaluate right_bottom point.
MGPosition right_bottom()const;

///Evaluate right_top point.
MGPosition right_top()const;

///Evaluate left_top point.
MGPosition left_top()const;

/// Return ending parameter value.
double param_e_u() const{return totalWidth();};
double param_e_v() const{return totalHeight();};


/// パラメータ範囲を返す。
///Return parameter range of the plane(Infinite box).
MGBox param_range() const;

/// Return starting parameter value.
double param_s_u() const{return 0.;};
double param_s_v() const{return 0.;};

/// Compute parameter curve.
///Returned is newed area pointer, and must be freed by delete.
MGCurve* parameter_curve(
	int is_u		///<Indicates x is u-value if is_u is true.
	, double x		///<Parameter value.
					///<The value is u or v according to is_u.
)const;

/// i must be < perimeter_num().
MGCurve* perimeter_curve(int i)const;

///Return how many perimeters this surface has.
int perimeter_num() const{return 4;};

/// Construct perimeter (u,v) parameter position.
/// i is perimeter number:
/// =0: v=min line, =1: u=max line, =2: v=max line, =3: u=min line
/// t is perimeter parameter line's parameter value of u or v.
MGPosition perimeter_uv(int i,double t) const;

///Obtain boundary and main parameter lines of the FSurface.
///skeleton includes boundary() and inner parameter lines.
///density indicates how many inner parameter lines are necessary
///for both u and v directions.
std::vector<UniqueCurve> skeleton(int density = 1)const;

///Get the name of the class.
std::string whoami()const{return "PlaneImage";};

///Read all member data.
void ReadMembers(MGIfstream& buf);
///Write all member data
void WriteMembers(MGOfstream& buf)const;

/// Output function.
std::ostream& toString(std::ostream&) const;

///Get the one pixel size in world coordinate for width.
double pixelSizeWidth()const{return m_pixelSizeWidth;};

///Get the one pixel size in world coordinate for height.
double pixelSeizeHeight()const{return m_pixelSizeHeight;};

///Get total image size in world coordnate for width.
double totalWidth()const;

///Get total image size in world coordnate for height.
double totalHeight()const;

private:

	mutable mgTexture m_texture;//Texture to draw this image.
	double m_pixelSizeWidth, m_pixelSizeHeight;
		///<size of one pixel of m_image in the world coordinates of the plane.

	double m_rightTexCoord, m_topTexCoord;///< Texture coordinate of the (right,top);

};

/** @} */ // end of GLAttrib group
#endif
