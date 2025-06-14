/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#include "StdAfx.h"
#include "mg/BPointSeq.h"
#include "mgGL/Texture.h"
#include "mgGL/OpenGLView.h"
#include "mgGL/glslprogram.h"
#include "mgGL/VBOLeafBuilder.h"
#include "mgGL/VBOLeaf.h"

namespace{static const float white[4]={1.,1.,1.,1.};}//Highlight back color.
namespace{	GLenum glErr;};

////////////////////////  mgVBOLeaf  //////////////////////////////

mgVBOLeaf::mgVBOLeaf()
:m_vertexArrayID(0),m_size(1.), // TODO Check!!
m_bufferID(0), m_count(0),m_primitiveMode(0),m_texture(0),m_drawType(mgGLSL::Primitive),
m_ColorSpecified(false),m_NormalSpecified(false),m_TextureSpecified(false)
,m_stippleFactor(-1),m_LineStipplePattern(0),m_lightMode(-1){
;
}

mgVBOLeaf::mgVBOLeaf(GLfloat size,MGColor& color)
:m_vertexArrayID(0),m_size(size),m_color(color),
m_bufferID(0), m_count(0),m_primitiveMode(0),m_texture(0),m_drawType(mgGLSL::Primitive),
m_ColorSpecified(false),m_NormalSpecified(false),m_TextureSpecified(false)
,m_stippleFactor(-1),m_LineStipplePattern(0),m_lightMode(-1){
;
}

mgVBOLeaf::mgVBOLeaf(const mgVBOLeaf& vbol)
:m_vertexArrayID(0),m_size(vbol.m_size),m_color(vbol.m_color),
m_bufferID(0),m_count(0),m_primitiveMode(0),m_texture(0),m_drawType(mgGLSL::Primitive),
m_ColorSpecified(false),m_NormalSpecified(false),m_TextureSpecified(false)
,m_stippleFactor(vbol.m_stippleFactor),m_LineStipplePattern(vbol.m_LineStipplePattern)
,m_lightMode(vbol.m_lightMode){
;
}

mgVBOLeaf::mgVBOLeaf(const mgVBOLeafBuilder& builder)
:m_vertexArrayID(0),m_size(builder.sizeStatic()),m_color(builder.colorStatic()),
m_bufferID(0),m_primitiveMode(builder.typeBegin()),
m_texture(builder.getTexture()),m_drawType(builder.getDrawType()),
m_ColorSpecified(false),m_NormalSpecified(false),m_TextureSpecified(false)
,m_stippleFactor(builder.m_stippleFactor)
,m_LineStipplePattern(builder.m_LineStipplePattern),m_lightMode(builder.m_lightMode){
	m_count=(unsigned)builder.m_VertexData.size();

	if(m_count==0)
		return;
	if(m_primitiveMode==GL_QUAD_STRIP)
		m_primitiveMode=GL_TRIANGLE_STRIP;
	glErr=glGetError();
	glGenBuffers(1,&m_bufferID);
	if((glErr=glGetError()) != GL_NO_ERROR){
		mgGLSL::printOpenGLError(glErr);
	}

	float* area=new float[m_count*4];//Maximum arear is reserved.
	unsigned baseSize=sizeof(float)*m_count;
	unsigned sizeV(baseSize*3), sizeC(0), sizeN(0), sizeT(0);
	unsigned sizeTotal=sizeV;
	if(builder.m_ColorData.size()>=m_count){
		m_ColorSpecified=true;
		sizeC=baseSize*4;
		sizeTotal+=sizeC;
	}
	if(builder.m_NormalData.size()>=m_count){
		m_NormalSpecified=true;
		sizeN=baseSize*3;
		sizeTotal+=sizeN;
	}
	if(builder.m_TextureData.size()>=m_count){
		m_TextureSpecified=true;
		sizeT=baseSize*2;
		sizeTotal+=sizeT;
	}
	glBindBuffer(GL_ARRAY_BUFFER,m_bufferID);
	glBufferData(GL_ARRAY_BUFFER,sizeTotal,0,GL_STATIC_DRAW);

	//vPosition
	for(unsigned i=0; i<m_count; i++){
		const vboFPoint& Pi=builder.m_VertexData[i];
		unsigned i3=i*3;
		area[i3++]=Pi.m_x; area[i3++]=Pi.m_y; area[i3]=Pi.m_z;
	}
	glBufferSubData(GL_ARRAY_BUFFER,0,sizeV,area);
	unsigned offset=sizeV;

	//vColor
	if(m_ColorSpecified){
		for(unsigned i=0; i<m_count; i++){
			const vboColor& Ci=builder.m_ColorData[i];
			unsigned i4=i*4;
			for(unsigned j=0;j<4; j++) area[i4+j]=Ci.m_color[j];
		}
		glBufferSubData(GL_ARRAY_BUFFER,offset,sizeC,area);
		offset+=sizeC;
	}

	//vNormal
	if(m_NormalSpecified){
		for(unsigned i=0; i<m_count; i++){
			const vboFPoint& N=builder.m_NormalData[i];
			unsigned i3=i*3;
			area[i3++]=N.m_x;area[i3++]=N.m_y;area[i3]=N.m_z;
		}
		glBufferSubData(GL_ARRAY_BUFFER,offset,sizeN,area);
		offset+=sizeN;
	}

	//vTexture
	if(m_TextureSpecified){
		for(unsigned i=0; i<m_count; i++){
			const vboFP2D& T=builder.m_TextureData[i];
			unsigned i2=i*2;
			area[i2++]=T.m_s;area[i2]=T.m_t;
		}
		glBufferSubData(GL_ARRAY_BUFFER,offset,sizeT,area);
	}
	glBindBuffer(GL_ARRAY_BUFFER,0);
	delete[] area;
}

mgVBOLeaf::~mgVBOLeaf(){
	if(m_bufferID)
		glDeleteBuffers(1,&m_bufferID);
	if(m_vertexArrayID)
		glDeleteVertexArrays(1,&m_vertexArrayID);
}

#define HILIGHT_LINE_LARGE_SIZE 2.f
#define HILIGHT_LINE_SMALL_SIZE 1.f
#define HILIGHT_POINT_LARGE_SIZE 7.f
#define HILIGHT_POINT_SMALL_SIZE 5.f
void mgVBOLeaf::execStaticAttrib(
	MGCL::VIEWMODE viewMode,
	bool selection
){
	const float* colr;
	GLfloat size=m_size;
	int vColorLoc=mgGLSLProgram::vColor;
	mgGLSLProgram* glslP=mgGLSLProgram::getCurrentGLSLProgram();
	if(!selection){
		glDisableVertexAttribArray(vColorLoc);
	}
	if(viewMode==MGCL::HIGHLIGHT){//When selection=true, viewMode!=HIGHLIGHT.
		mgGLSL::execLightMode(0);//Light Off.
		mgGLSL::execStaticColorAttrib(getHilightColor());
		if(m_primitiveMode==GL_POINTS){
			glDisable(GL_PROGRAM_POINT_SIZE);
			glPointSize(HILIGHT_POINT_LARGE_SIZE);
			size=HILIGHT_POINT_SMALL_SIZE;
		}else{
			mgGLSL::execStaticLineWidth(HILIGHT_LINE_LARGE_SIZE);
			size=HILIGHT_LINE_SMALL_SIZE;
		}
		glDrawArrays(m_primitiveMode,0,(GLsizei)m_count);
		colr=white;
	}else{
		colr=m_color.color();
		mgGLSL::execLightMode(m_lightMode);
	}
	
	if(!selection){
		if(viewMode==MGCL::HIGHLIGHT || m_color.defined())
			mgGLSL::execStaticColorAttrib(colr);
	}

	//Process of m_size.
	if(m_primitiveMode==GL_POINTS && size>0.){
	//Point size
		if(viewMode==MGCL::HIGHLIGHT){
			glDisable(GL_PROGRAM_POINT_SIZE);
			glPointSize(size);
		}else{
			glEnable(GL_PROGRAM_POINT_SIZE);
			glslP->setUniform(mgGLSLProgram::pointSize,size);
		}
	}else if(m_primitiveMode==GL_LINES ||
		m_primitiveMode==GL_LINE_STRIP || m_primitiveMode==GL_LINE_LOOP){
	//Line width
		mgGLSL::execStaticLineWidth(size);
		mgGLSL::execStaticLineStipple(m_stippleFactor,m_LineStipplePattern);
	}else{
	//Polygon mode(m_size defines what kinds of polygon mode be applied).
		if(m_size<0.)
			glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
		else if(m_size<=m_PolygonModePointSizeBase){
			glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
			glLineWidth(size);
		}else{
			glPolygonMode(GL_FRONT_AND_BACK,GL_POINT);
			glPointSize(m_size-m_PolygonModePointSizeBase);
		}
	}
}

///描画関数draw()は、!is_made()であれば、作成し、表示処理をする。
///is_made()(描画データ作成済み）であれば、すでに作成されたmgVBOElementの描画を行う。
///mgVBOLeafでは常に作成すみであり表示処理をおこなう。
void mgVBOLeaf::draw(MGCL::VIEWMODE viewMode)
{
	if(is_no_display())
		return;
	if(m_bufferID==0)
		return;

	if(m_vertexArrayID==0)
		glGenVertexArrays(1,&m_vertexArrayID);
	glBindVertexArray(m_vertexArrayID);

	mgGLSLProgram* glsl=mgGLSLProgram::getCurrentGLSLProgram();
	glsl->setUniform(mgGLSLProgram::DrawType, m_drawType);
	if(texture_is_bound())
		m_texture->use();

	glBindBuffer(GL_ARRAY_BUFFER, m_bufferID);
	assert(glIsBuffer(m_bufferID));
	int vPositionLoc=mgGLSLProgram::vPosition;
	glVertexAttribPointer(vPositionLoc,3,GL_FLOAT,GL_FALSE,0,(GLvoid*)0);
	glEnableVertexAttribArray(vPositionLoc);

mgGLSL::pushStaticGLAttrib();
	execStaticAttrib(viewMode,false);

	size_t offset=0;
	unsigned baseSize=sizeof(float)*m_count;
	offset+=baseSize*3;
	if(m_ColorSpecified){
		int vColorLoc=mgGLSLProgram::vColor;
		glVertexAttribPointer(vColorLoc,4,GL_FLOAT,GL_FALSE,0,(GLvoid*)offset);
		glEnableVertexAttribArray(vColorLoc);

		offset+=baseSize*4;
	}
	if(m_NormalSpecified){
		int vNormalLoc=glsl->getvNormalLocation();
		if(vNormalLoc>=0){
			glVertexAttribPointer(vNormalLoc,3,GL_FLOAT,GL_FALSE,0,(GLvoid*)offset);
			glEnableVertexAttribArray(vNormalLoc);
		}
		offset+=baseSize*3;
	}
	if(m_TextureSpecified){
		int vTextureLoc=glsl->getvTextureCoordLocation();
		if(vTextureLoc>=0){
			glVertexAttribPointer(vTextureLoc,2,GL_FLOAT,GL_FALSE,0,(GLvoid*)offset);
			glEnableVertexAttribArray(vTextureLoc);
		}
	}

	glDrawArrays(m_primitiveMode,0,(GLsizei)m_count);

mgGLSL::popStaticGLAttrib();
	glBindBuffer(GL_ARRAY_BUFFER,0);
};
	
///描画関数selectionDraw()は、Object選択のための表示処理をする。
///通常のdrawとの相違：///Colorとしてm_bufferIDを用い、size処理以外の
///attributesの処理（normal, texture, color)をしない。
void mgVBOLeaf::selectionDraw(MGCL::VIEWMODE viewMode){
	if(is_no_display())
		return;
	if(m_bufferID==0)
		return;

	if(m_vertexArrayID==0)
		glGenVertexArrays(1,&m_vertexArrayID);
	glBindVertexArray(m_vertexArrayID);

	glBindBuffer(GL_ARRAY_BUFFER,m_bufferID);
	assert(glIsBuffer(m_bufferID));
	int vPositionLoc=mgGLSLProgram::vPosition;
	glVertexAttribPointer(vPositionLoc,3,GL_FLOAT,GL_FALSE,0,(GLvoid*)0);
	glEnableVertexAttribArray(vPositionLoc);

	execStaticAttrib(viewMode,true);

	glDrawArrays(m_primitiveMode,0,(GLsizei)m_count);
	glBindBuffer(GL_ARRAY_BUFFER,0);
};

///m_primitiveMode=GL_QUAD_STRIP, GL_TRIANGLES, GL_TRIANGLE_STRIP, GL_TRIANGLE_FAN
///のときにのみ有効でそのPolygonMode（GL_POINT, GL_LINE, GL_FILL)を指定する.
///setPolygonModeしていないmgVBOLeafはGL_FILLとされる。
void mgVBOLeaf::setPolygonMode(GLenum mode){
	if(m_primitiveMode==GL_TRIANGLES ||
		m_primitiveMode==GL_TRIANGLE_STRIP || m_primitiveMode==GL_TRIANGLE_FAN){
		GLfloat size=m_size;
		if(size>m_PolygonModePointSizeBase)
			size-=m_PolygonModePointSizeBase;
		else if(size<0.)
			size*=-1.;

		if(mode==GL_FILL){
			m_size=-size;
		}else if(mode==GL_POINT){
			m_size=size+m_PolygonModePointSizeBase;
		}else
			m_size=size;
	}
}

void mgVBOLeaf::setLightMode(int mode){
	m_lightMode=mode;
}
///Line stipple属性をセットする。
///When factor=0 is input, line pattern is disabled. This means the line is solid.
///When factor<0, the stipple attribute is undefined. This means the attribute
///is defined by the environment.
///When factor<=0, pattern is unnecessary.
void mgVBOLeaf::setLineStipple(short int factor, GLushort pattern){
	m_stippleFactor=factor;
	m_LineStipplePattern=pattern;
}

void mgVBOLeaf::setStaticAttribLineWidth(
	GLfloat size///size<=0. はundefinedを示す
){
	m_size=size;
}
void mgVBOLeaf::setStaticAttribPointSize(
	GLfloat size///size<=0. はundefinedを示す
){
	m_size=size;
}

///作成済のmgVBOLeafのIDの位置のVertex dataをPのデータに置き換える。
///0<= ID <getVerticesNumber();
void mgVBOLeaf::updateVertex(unsigned ID, const MGPosition& P){
	assert(ID<m_count);
	//vPosition update.
	float Pi[3]={(float)P[0],(float)P[1],(float)P[2]};
	updateVertices(ID,1,Pi);
}

///作成済のmgVBOLeafのstartIDの位置からnum個のVertex dataをPsのデータに
///置き換える。 0<= startID , startID+num<=getVerticesNumber().
///Here num=Ps.size();
void mgVBOLeaf::updateVertices(unsigned startID, const std::vector<MGPosition>& Ps){
	int num=(int)Ps.size();assert(startID+num<=m_count);
	float* area=new float[num*3];
	//vPosition update.
	for(int i=0, i3=0; i<num; i++){
		const MGPosition& Pi=Ps[i];
		area[i3++]=(float)Pi[0];
		area[i3++]=(float)Pi[1];
		area[i3++]=(float)Pi[2];
	}
	updateVertices(startID,num,area);
	delete[] area;
}

///作成済のmgVBOLeafのstartIDの位置からnum個のVertex dataをPsのデータに
///置き換える。 0<= startID , startID+num<=getVerticesNumber().
///Here num=Ps.length();
void mgVBOLeaf::updateVertices(unsigned startID, const MGBPointSeq& Ps){
	int num=Ps.length();assert(startID+num<=m_count);
	float* area=new float[num*3];
	//vPosition update.
	for(int i=0, i3=0; i<num; i++){
		area[i3++]=(float)Ps(i,0);
		area[i3++]=(float)Ps(i,1);
		area[i3++]=(float)Ps(i,2);
	}
	updateVertices(startID,num,area);
	delete[] area;
}

///作成済のmgVBOLeafのstartIDの位置からnumVerticesのVertex dataをareaのデータに
///置き換える。 0<= startID , startID+numVertices<=getVerticesNumber().
///areaの配列の長さ＝numVertices*3となる
void mgVBOLeaf::updateVertices(unsigned startID, unsigned numVertices, const float* area){
	assert(startID+numVertices<=m_count);
	//vPosition update.
	glBindBuffer(GL_ARRAY_BUFFER,m_bufferID);
	glBufferSubData(GL_ARRAY_BUFFER,startID*12,numVertices*12,area);
	glBindBuffer(GL_ARRAY_BUFFER,0);
}

///Textureをsetする。
void mgVBOLeaf::setTexture(mgTexture* texture){
	m_texture=texture;
}

///Textureをgetする。
const mgTexture* mgVBOLeaf::getTexture()const{
	return m_texture;
}

///Textureをgetする。
mgTexture* mgVBOLeaf::getTexture(){
	return m_texture;
}

///Test if a texture is binded.
bool mgVBOLeaf::texture_is_bound()const{
	if(m_texture){
		if(m_texture->getTextureID()!=0)
			return true;
	}
	return false;
}
