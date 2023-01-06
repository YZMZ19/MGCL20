/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#if !defined( __VBOBYSCREEN_H__)
#define __VBOBYSCREEN_H__

#include "mgGL/VBO.h"

class MGPosition;

/** @file */
/** @addtogroup DisplayHandling
 *  @{
 */

///mgVBOByScreen is a VBO to draw pictures in screen coordinates of a constant z value. 

///All of the input coordinates(except z) are treated as screen coordinates
///whose screen data is set by set_viewport. The origin is (x,y) and the size is
///(width, height). Here (x,y,width,height) are the values
///set by set_viewport(). 
///All of the coordinates are converted to the normalized device coordinate(NDC)
///whose range is(-1,-1) to (1,1). Only z value given by the constructor or se_ZValue
///is the NDC.
class MG_DLL_DECLR mgVBOByScreen:public mgVBO{

public:

///construct given z(in NDC).
mgVBOByScreen(float z = -1.f):mgVBO(mgGLSL::NdcScreen), m_z(z){ ; };
	
///Convert screen coordinates to normalized world coordinates.
void convert2NormalizedWorld(
	const int xyS[2],	///<inout opject target to convert.
	float& x,			///<converted normalized world coordinates, x.
	float& y			///<converted normalized world coordinates, y.
);
	
///Convert screen coordinates to normalized world coordinates.
void convert2NormalizedWorld(
	float x0,	///<inout opject target to convert, x.
	float y0,	///<inout opject target to convert, y.
	float& x,	///<converted normalized world coordinates, x.
	float& y	///<converted normalized world coordinates, y.
);

///Get height.
float getHeight()const{return m_heightHalf*2.f;};

///Get width.
float getWidth()const{return m_widthHalf*2.f;};

///screen data をセットする
void set_viewport(const int vport[4]);

//Set the ormalized device coordinates z value (in NDC) to control depth test.

///All of the points' z value of Vertex() after set_ZValue is invoked
///will be this value.
void set_ZValue(float z){m_z=z;};

///頂点の座標値を指定する
void Vertex(const int v[2]);
void Vertex(int x, int y);
void Vertex(const MGPosition& v);
void Vertex(float x, float y, float z=0.0f);
void Vertex3d(double x, double y, double z=0.0);
void Vertex2fv(const float v[2]);
void Vertex3fv(const float v[3]);
void Vertex2dv(const double v[2]);
void Vertex3dv(const double v[3]);

protected:

	float m_centerX, m_centerY;///< Viewport data of the screen.
	float m_widthHalf, m_heightHalf;///<half of width and height.
	float m_z;///<Normalized device coordinates z value to control depth test.

};

/** @} */ // end of DisplayHandling group
#endif //__VBOBYSCREEN_H__
