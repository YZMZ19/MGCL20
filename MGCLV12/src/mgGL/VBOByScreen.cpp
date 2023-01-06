/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mgGL/VBOByScreen.h"
#include "mgGL/VBOLeaf.h"

///mgVBOByScreen is a VBO to draw pictures in screen coordinates of a constant
///z value. All of the input coordinates(except z) are treated as screen coordinates
///whose screen data is set by set_viewport. The origin is (x,y) and the size is
///(width, height). Here (x,y,width,height) are the values
///set by set_viewport(). 
///All of the coordinates are converted to the normalized device coordinate(NDC)
///whose range is(-1,-1) to (1,1). Only z value given by the constructor or se_ZValue
///is the NDC.
	
///screen data をセットする
void mgVBOByScreen::set_viewport(const int vport[4]){
	m_widthHalf=float(vport[2]*.5);
	m_heightHalf=float(vport[3]*.5);
	m_centerX=float(vport[0])+m_widthHalf;
	m_centerY=float(vport[1])+m_heightHalf;
}	

///頂点のscreen座標値を指定する
void mgVBOByScreen::Vertex(const int v[2]){
	Vertex(float(v[0]),float(v[1]));
}
void mgVBOByScreen::Vertex(int x, int y){
	Vertex(float(x),float(y));
}

void mgVBOByScreen::Vertex(const MGPosition& v){
	Vertex(float(v[0]),float(v[1]));
}

void mgVBOByScreen::Vertex(float x, float y, float z){
	if(!is_InBegin())
		return;

	vboFPoint fp;
	convert2NormalizedWorld(x,y,fp.m_x, fp.m_y);
	fp.m_z=m_z;
	m_builder->push_backVertex(fp);
}

void mgVBOByScreen::Vertex3d(double x, double y, double z){
	Vertex(float(x),float(y));
}

void mgVBOByScreen::Vertex2fv(const float v[2]){
	Vertex(v[0],v[1]);
}

void mgVBOByScreen::Vertex3fv(const float v[3]){
	Vertex(v[0],v[1]);
}

void mgVBOByScreen::Vertex2dv(const double v[2]){
	Vertex(float(v[0]),float(v[1]));
}

void mgVBOByScreen::Vertex3dv(const double v[3]){
	Vertex(float(v[0]),float(v[1]));
}

void mgVBOByScreen::convert2NormalizedWorld(
	const int xyS[2],	///<inout opject target to convert.
	float& x,			///<converted normalized world coordinates.
	float& y
){
	convert2NormalizedWorld(float(xyS[0]), float(xyS[1]),x,y);
}
	
///Convert screen coordinates to normalized world coordinates.
void mgVBOByScreen::convert2NormalizedWorld(
	float x0, float y0,	///<inout opject target to convert.
	float& x,			///<converted normalized world coordinates.
	float& y
){
	x=(x0-m_centerX)/m_widthHalf;
	y=(y0-m_centerY)/m_heightHalf;
}
