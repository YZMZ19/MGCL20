/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#include "StdAfx.h"
#include "mg/Straight.h"
#include "mgGL/PlaneImage.h"
#include "mgGL/VBO.h"
#include "mgGL/Texture.h"
#include "mgGL/GLSLProgram.h"

//
//Define MGPlaneImage Class.

///MGPlaneImage defines square image plane.
///This is displayed using texture mapping.
///The plane's root point is the (left, bottom) point of the image(m_image).
///m_uderiv and m_vderiv of the plane are always defined as unit_vectors.
///The pixel sizes are also defined as m_pixelSizeWidth and m_pixelSizeHeight.
///This implies width and hight can be different size.

//////////////////Constructor//////////////////

MGPlaneImage::MGPlaneImage(
	const MGVector& uderiv,///<image's width direction.
	const MGVector& vderiv,///<image's height direction.
		///vderiv is normal to uderiv. If this is not the case,
		///vderiv will be so transformed as normal to uderiv.
		///uderiv's direction is not changed.
	const MGPosition &origin,///(left,bottom) point of the image.
	MGImage& image,//Image data will be transfered to this MGPlaneImage.
	double pixelSizeWidth,//one pixel size in world coordinate for width.
	double pixelSizeHeight//one pixel size in world coordinate for height.
):MGPlane(uderiv,vderiv,origin), m_image(image),
m_pixelSizeWidth(pixelSizeWidth), m_pixelSizeHeight(pixelSizeHeight){
	MGPlane::normalize();///Change for the plane's (uderiv(), vderiv(), normal)
		///to construct a orthonormal system.

	assert(m_pixelSizeWidth>0.);
	if(m_pixelSizeHeight<=0.)
		m_pixelSizeHeight=m_pixelSizeWidth;
	m_rightTexCoord=1.;//m_rightTexCoord=double(image_width())/m_pixelSizeWidth;
	m_topTexCoord=1.;//m_topTexCoord=double(image_height())/m_pixelSizeHeight;
}

///Assignment
MGPlaneImage& MGPlaneImage::operator=(const MGPlaneImage& gel2){
	MGPlane::operator=(gel2);
	m_image=gel2.m_image;

	m_pixelSizeWidth=gel2.m_pixelSizeWidth;
	m_pixelSizeHeight=gel2.m_pixelSizeHeight;

	m_rightTexCoord=gel2.m_rightTexCoord;
	m_topTexCoord=gel2.m_topTexCoord;
	return *this;
}

////////Member Function////////

///Obtain ceter coordinate of the geometry.
MGPosition MGPlaneImage::center()const{
	double uc=param_e_u()*.5;
	double vc=param_e_v()*.5;
	return eval(uc,vc);
}

///Generate a newed clone object.
MGPlaneImage* MGPlaneImage::clone()const{
	MGPlaneImage* pimg=new MGPlaneImage;
	pimg->MGPlane::operator=(*this);
	std::unique_ptr<MGImage> image2(m_image.clone());
	pimg->m_image=*image2;

	pimg->m_pixelSizeWidth=m_pixelSizeWidth;
	pimg->m_pixelSizeHeight=m_pixelSizeHeight;
	pimg->m_rightTexCoord=m_rightTexCoord;
	pimg->m_topTexCoord=m_topTexCoord;
	return pimg;
}

///Draw 3D curve in world coordinates.
void MGPlaneImage::shade(
	mgVBO& vbo,
	const MGDrawParam& para,
	MGCL::DRAW_TARGET target
)const{
	if(!m_texture.getTextureID()){
		m_texture.set_image(m_image
			,true//mutable texture
			,false//Not PointSprite
			,GL_CLAMP,GL_NEAREST );
	}

	// MGPlaneImage‚ÍLight‚Ì‰e‹¿‚ðŽó‚¯‚È‚¢
	mgLightModeSwitcher switcher(vbo, mgGLSL::NoShading);

	vbo.Begin(GL_TRIANGLE_STRIP,target);
		MGPosition LT=left_top();//left top point
		vbo.TexCoord2d(0.,m_topTexCoord);vbo.Vertex3dv(LT.data());

		const MGPosition& LB=root_point();//Left bottom point.
		vbo.TexCoord2d(0.,0.); vbo.Vertex3dv(LB.data());

		MGPosition RT=right_top();//right top point
		vbo.TexCoord2d(m_rightTexCoord,m_topTexCoord);vbo.Vertex3dv(RT.data());

		MGPosition RB=right_bottom();//right bottom point.
		vbo.TexCoord2d(m_rightTexCoord,0.);vbo.Vertex3dv(RB.data());
		vbo.setTexture(&m_texture);
		vbo.setDrawType(mgGLSL::Texture);
	vbo.End();
}

///Evaluate left_bottom point.
MGPosition MGPlaneImage::left_bottom()const{
	return root_point();
}

///Evaluate right_bottom point.
MGPosition MGPlaneImage::right_bottom()const{
	return eval(param_e_u(),0.);
}

///Evaluate right_top point.
MGPosition MGPlaneImage::right_top()const{
	return eval(param_e_u(),param_e_v());
}

///Evaluate left_top point.
MGPosition MGPlaneImage::left_top()const{
	return eval(0.,param_e_v());
}

///Return parameter range of the plane(Infinite box).
MGBox MGPlaneImage::param_range() const{
	return MGSurface::param_range();
}

///Compute parameter curve.
///Returned is newed area pointer, and must be freed by delete.
MGCurve* MGPlaneImage::parameter_curve(
	int is_u		///<Indicates x is u-value if is_u is true.
	, double x		///<Parameter value.
					///<The value is u or v according to is_u.
)const{

	const MGVector* dir; const MGVector* dir2;
	double paramE;
	if(is_u){
		paramE=param_e_v();
		dir2=&u_deriv();
		dir=&v_deriv();
	}else{
		paramE=param_e_u();
		dir2=&v_deriv();
		dir=&u_deriv();
	}
	MGPosition origin(root_point());
	origin+=x*(*dir2);
	return new MGStraight(paramE,0.,*dir,origin);
}

/// i must be < perimeter_num().
///When perimeter_num()==0, this function is undefined.
MGCurve* MGPlaneImage::perimeter_curve(int i)const{
	return MGSurface::perimeter_curve(i);
}

/// Construct perimeter (u,v) parameter position.
/// i is perimeter number:
/// =0: v=min line, =1: u=max line, =2: v=max line, =3: u=min line
/// t is perimeter parameter line's parameter value of u or v.
MGPosition MGPlaneImage::perimeter_uv(int i,double t)const{
	return MGSurface::perimeter_uv(i,t);
}

///Obtain boundary and main parameter lines of the FSurface.
///skeleton includes boundary() and inner parameter lines.
///density indicates how many inner parameter lines are necessary
///for both u and v directions.
std::vector<UniqueCurve> MGPlaneImage::skeleton(int density)const{
	return MGFSurface::skeleton(density);
}

//Get total image size in world coordnate for width.
double MGPlaneImage::totalWidth()const{
	double w=m_pixelSizeWidth*double(image_width());
	return w;
}

//Get total image size in world coordnate for height.
double MGPlaneImage::totalHeight()const{
	double h=m_pixelSizeHeight*double(image_height());
	return h;
}

///Read all member data.
void MGPlaneImage::ReadMembers(MGIfstream& buf){
}
///Write all member data
void MGPlaneImage::WriteMembers(MGOfstream& buf)const{
}

/// Output function.
std::ostream& MGPlaneImage::toString(std::ostream& strm) const{
	strm<<"MGPlaneImage::"<<this<<",";
	strm<<"image's size=("<<image_width()<<","<<image_height()<<")"<<std::endl;
	MGPlane::toString(strm);
	return strm;
}

///Return minimum box that includes whole of the surface.
///Returned is a newed object pointer.
void MGPlaneImage::compute_box(MGBox& bx) const{
	double ue=param_e_u();
	double ve=param_e_v();
	bx=MGBox(eval(0.,0.), eval(ue,ve));
	bx.expand(eval(0.,ve));
	bx.expand(eval(ue,0.));
}
