/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno             */
/* All rights reserved.                                             */
/********************************************************************/

#if !defined(_MGVBOLEAFBUILDER__INCLUDED_)
#define _MGVBOLEAFBUILDER__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "StdAfx.h"
#include "mg/Position.h"
#include "mgGL/Color.h"

/////////////////////////////////////////////////////////////////////////////
// Vertex Buffer Object Class.

///@cond

class vboFPoint{
public:
	float m_x,m_y,m_z;
	vboFPoint():m_x(0.f),m_y(0.f),m_z(0.f){;};
	vboFPoint(float x, float y, float z):m_x(x),m_y(y),m_z(z){;};
	vboFPoint(double x, double y, double z):m_x((float)x),m_y((float)y),m_z((float)z){;};
	vboFPoint(const MGPosition& P):m_x((float)P[0]),m_y((float)P[1]),m_z((float)P[2]){;};
	vboFPoint(const MGVector& P):m_x((float)P[0]),m_y((float)P[1]),m_z((float)P[2]){;};
};
class vboFP2D{
public:
	float m_s,m_t;
	vboFP2D():m_s(0.f),m_t(0.f){;};
	vboFP2D(float s, float t):m_s(s),m_t(t){;};
	vboFP2D(double s, double t):m_s((float)s),m_t((float)t){;};
	vboFP2D(const MGPosition& P):m_s((float)P[0]),m_t((float)P[1]){;};
	vboFP2D(const MGVector& P):m_s((float)P[0]),m_t((float)P[1]){;};
};
class vboColor{
public:
	float m_color[4];
	vboColor(const float colr[4]){for(int i=0; i<4; i++) m_color[i]=colr[i];};
	vboColor(const MGColor& clr){clr.get_color(m_color);};
};



///Assistant class for mgVBO.
///mgVBOLeafBuilder holds a temporary data for mgVBO's Begin() and End()
///to send the data to OpenGL.
class mgVBOLeafBuilder{

friend class mgVBOLeaf;
private:

	mgGLSL::DrawType m_drawType;

	mutable mgTexture* m_texture;

	unsigned m_typeBegin;///<type specified by Begin(). After End() m_typeBegin is set to 0.
	MGColor m_colorStatic;///Color specified by setStaticAttribColor().
		///The color of mgVBOLeaf generated is set to this color.
	GLfloat m_sizeStatic=0.;///The size specified by setStaticAttribSize.
		///The size of mgVBOLeaf generated is set to this size.

	short int m_stippleFactor;///Line stipple factor.
		///If m_stippleFactor=0, line Stipple is disabled.
		///   m_stippleFactor<0, line stipple is undefined.
	GLushort m_LineStipplePattern;///m_LineStipplePatternindicates the pattern.

	///light modeはm_elementsShadeに対してのみ有効。m_elementsに対しては常にlightはオフ
	int m_lightMode;/// <0: undefined, =0:Light is disabled, >0:Light is enabled.

	std::vector<vboFPoint> m_VertexData;
	std::vector<vboColor> m_ColorData;
	std::vector<vboFPoint> m_NormalData;
	std::vector<vboFP2D> m_TextureData;

// オペレーション
public:

mgVBOLeafBuilder(GLenum type):m_typeBegin(type),m_drawType(mgGLSL::Primitive),m_texture(0)
,m_stippleFactor(-1),m_LineStipplePattern(0),m_lightMode(-1){;};

mgVBOLeafBuilder(GLenum type, const MGColor& colorStatic, GLfloat sizeStatic)
:m_typeBegin(type),m_colorStatic(colorStatic),m_sizeStatic(sizeStatic),
m_drawType(mgGLSL::Primitive),m_texture(0)
,m_stippleFactor(-1),m_LineStipplePattern(0),m_lightMode(-1){;};

bool is_null()const{return m_VertexData.size()==0;};

void setTypeBegin(GLenum typeBegin){m_typeBegin=typeBegin;};
unsigned typeBegin()const{return m_typeBegin;};

///Static attributeを設定する。
///begin() - end()の間であればそのmgVBOLeafに対しても適用され
///その後に作成されるbegin() - end()のmgVBOLeafすべてに適用される。
void setStaticAttribColor(const MGColor& color){m_colorStatic=color;};
void setStaticAttribColor(const float color[4]){
	m_colorStatic=MGColor(color[0],color[1],color[2],color[3]);};
void setStaticAttribColor(float r, float g, float b){
	m_colorStatic=MGColor(r,g,b);};
const MGColor& colorStatic()const{return m_colorStatic;};

void setStaticAttribSize(GLfloat size){m_sizeStatic=size;};
GLfloat sizeStatic()const{return m_sizeStatic;};

///Line stipple属性をセットする。
///When factor=0 is input, line pattern is disabled. This means the line is solid.
///When factor<0, the stipple attribute is undefined. This means the attribute
///is defined by the environment.
///When factor<=0, pattern is unnecessary.
void setLineStipple(short int factor,GLushort pattern){
	m_stippleFactor=factor;
	m_LineStipplePattern=pattern;
};

///Set light mode. mode=-1:undefined, =0:disabled, =1:enabled.
void setLightMode(int mode){m_lightMode=mode;};

void setDrawType(mgGLSL::DrawType drawType){m_drawType = drawType;};
mgGLSL::DrawType getDrawType()const{return m_drawType;};

void setTexture(mgTexture* texture){m_texture=texture;};
mgTexture* getTexture()const{return m_texture;};

void push_backVertex(const vboFPoint& v){m_VertexData.push_back(v);};
void push_backColor(const vboColor& c){m_ColorData.push_back(c);};
void push_backNormal(const vboFPoint& n){m_NormalData.push_back(n);};
void push_backTexture(const vboFP2D& t){m_TextureData.push_back(t);};

};

///@endcond

#endif // !defined(_MGVBOLEAFBUILDER__INCLUDED_)
