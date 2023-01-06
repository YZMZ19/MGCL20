/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#include "StdAfx.h"
#include "mgGL/StaticGLAttrib.h"
#include "mgGL/Color.h"
#include "mgGL/openglView.h"

//
//Define mgStaticGLAttrib Class.
///mgStaticGLAttrib defines MGColor and line width data of OpenGL.


mgStaticGLAttrib::mgStaticGLAttrib():m_lineWidth(1.f)
,m_stippleFactor(0),m_LineStipplePattern(0),m_lightMode(0){
	MGColor dcolor;
	const float* fc=dcolor.color();
	for(int i=0; i<4; i++)
		m_color[i]=fc[i];
}

mgStaticGLAttrib::mgStaticGLAttrib(const MGColor& color, float lineWidth)
:m_lineWidth(lineWidth)
,m_stippleFactor(0),m_LineStipplePattern(0),m_lightMode(0){
	const float* fc=color.color();
	for(int i=0; i<4; i++)
		m_color[i]=fc[i];
}

void mgStaticGLAttrib::setColor(const float color[4]){
	for(int i=0; i<4; i++)
		m_color[i]=color[i];
}

void mgStaticGLAttrib::setColor(const MGColor& color){
	const float* fc=color.color();
	for(int i=0; i<4; i++)
		m_color[i]=fc[i];
}

void mgStaticGLAttrib::getColor(MGColor& color){
	color=MGColor(m_color[0], m_color[1], m_color[2], m_color[3]);
}

///Line stipple属性をセットする。
///When factor=0 is input, line pattern is disabled. This means the line is solid.
///When factor<0, the stipple attribute is undefined. This means the attribute
///is defined by the environment.
///When factor<=0, pattern is unnecessary.
void mgStaticGLAttrib::setLineStipple(short int factor, GLushort pattern){
	m_stippleFactor=factor;
	m_LineStipplePattern=pattern;
}

void mgStaticGLAttrib::getLineStipple(short int& factor, GLushort& pattern)const{
	factor=m_stippleFactor;
	pattern=m_LineStipplePattern;
}
