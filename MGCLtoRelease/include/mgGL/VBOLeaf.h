/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#if !defined(_MGVBOLEAF__INCLUDED_)
#define _MGVBOLEAF__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "StdAfx.h"
#include "mgGL/glslprogram.h"
#include "mgGL/Color.h"
#include "mgGL/VBOElement.h"

class MGPosition;
class MGBPointSeq;
class MGCurve;
class mgVBO;
class mgVBOLeafBuilder;
class mgTexture;

/** @file */
/** @addtogroup DisplayHandling
 *  @{
 */


/////////////////////////////////////////////////////////////////////////////
// mgVBOLeaf

///mgVBO�ɑ΂��ĕ`��f�[�^�쐬��̏���ێ����邽�߂̃N���X.

///mgVBO��mgVBOLeaf��mgVBO���g��element�Ƃ��ĕێ�����B
///mgVBO��element�Ƃ���̂�MGGroup�܂���MGShell�̃����o�[�f�[�^�̕\������
///�ێ����邽��
class MG_DLL_DECLR mgVBOLeaf:public mgVBOElement{

///Maximum line width.

///Used to judge the polygon mode of the polygon if this mgVBOLeaf is a triangle.
///When polygon mode = GL_LINE, m_size=the line width.
///When polygon mode = GL_POINT, the point size=m_size-m_PolygonModePointSizeBase.
///That is, if m_size>m_PolygonModePointSizeBase, its polygon mode is GL_POINT and
///the point size=m_size-m_PolygonModePointSizeBase.
static const GLfloat m_PolygonModePointSizeBase;

private:

	mgGLSL::DrawType m_drawType;/// �`�惂�[�h mgGLSL::DrawType�̒l������

	GLuint m_bufferID;///<glGenBuffers�łō쐬���ꂽBuffer id�B���쐬�̂Ƃ���=0.
				///m_bufferID��m_vertexArrayID��Bind�����B

	///���_���Ƃ�(Static�łȂ��j�l���w�肳��Ă��邩�ۂ�������
	bool m_ColorSpecified:1, m_NormalSpecified:1, m_TextureSpecified:1;

	mutable GLuint m_vertexArrayID;///<glGenVertexArrays�ō쐬���ꂽVertexArray.���쐬�̂Ƃ���=0.

	///(m_primitiveMode,m_count)�͊ȈՕ\��(VIEWMODE==WIRE)�܂��͋����\��(VIEWMODE==HIGHLIGHT)�ɗ��p�����B
	///(m_primitiveMode,m_count)��glDrawArrays��m_bufferID��first=0�Ŏg�p�����B
	GLenum m_primitiveMode;///<glDrawArrays�ɓ��͂Ƃ���primitive���
	unsigned m_count;///<glDrawArrays�ɓ��͂Ƃ���Vertex�̐�

	mutable mgTexture* m_texture;///< Texture to apply. When m_texture.get()!=null, m_texture.use()
		///<will be invoked when draw().

	MGColor m_color;///<Static Vertex-Attributes for this VBO to apply.
		///<Color data.

	GLfloat m_size;///<Static Vertex-Attributes for this VBO to apply.
		///<Line-width, or point-size.
		///<When m_primitiveMode=GL_POINTS, m_size means point size.
		///<When m_primitiveMode=GL_LINEXXX, m_size means line width.
		///<When m_primitiveMode=GLTRIANGLEXXXX, m_size means following:
		///<  m_size<0. GL_FILL
		///<  1.<=m_size<=m_PolygonModePointSizeBase. GL_LINE and the line width=m_size.
		///<  m_PolygonModePointSizeBase<m_size GL_POINT and
		///<  the point size=(m_size-m_PolygonModePointSizeBase).

	short int m_stippleFactor;///<Line stipple factor.
		///<If m_stippleFactor=0, line Stipple is disabled.
		///<   m_stippleFactor<0, line stipple is undefined.
	GLushort m_LineStipplePattern;///<m_LineStipplePatternindicates the pattern.

	///light mode��m_elementsShade�ɑ΂��Ă̂ݗL���Bm_elements�ɑ΂��Ă͏��light�̓I�t
	int m_lightMode;///<  <0: undefined, =0:Light is disabled, >0:Light is enabled.

friend class mgVBO;

public:

mgVBOLeaf();
mgVBOLeaf(GLfloat size,MGColor& color);
mgVBOLeaf(const mgVBOLeaf& vbol);
mgVBOLeaf(const mgVBOLeafBuilder& builder);
~mgVBOLeaf();

///����mgVBOElement��null(���܂�draw/make_display_list()��������Ă��Ȃ��j���𔻒�

///mgVBO��m_gel�ɑ΂���make_display_list()�������Ȃ���Ă��Ȃ��Ƃ�false���Ԃ����
///mgVBOLeaf��make_display_list()���ꂽ�Ƃ��̂ݑ��݂��邽�ߏ��true���Ԃ�B
bool is_made(MGCL::VIEWMODE viewMode=MGCL::DONTCARE){return true;};

///Remake the display list.

///���łɍ쐬�ς݂ł����Ă������I�ɍč쐬���s���B
///mgVBOLeaf�͍쐬�ς݂̏ꍇ�̂ݑ��݂��邽�߉������Ȃ��B
void make_display_list(MGCL::VIEWMODE vmode=MGCL::DONTCARE){;};

///Draw.

///�`��֐�draw()�́Ais_made()�ł���΁A�쐬���A�\������������B
///!is_made()(�`��f�[�^�쐬�ς݁j�ł���΁A���łɍ쐬���ꂽmgVBOElement�̕`����s���B
///mgVBOLeaf�ł͏�ɍ쐬���݂ł���\�������������Ȃ��B
void draw(MGCL::VIEWMODE viewMode=MGCL::DONTCARE);

///Redraw.

///draw()��mgVBOLeaf���쐬�ς݁inot null)�ł���΍쐬�������s��Ȃ����A
///redraw()�͋����I�ɍč쐬���s���`�揈���������Ȃ��B
///mgVBOLeaf�ł͏�ɍ쐬���݂ł���\�������������Ȃ��B
void redraw(MGCL::VIEWMODE viewMode=MGCL::DONTCARE){draw(viewMode);};

///Draw for the Selection.

///�`��֐�selectionDraw()�́AObject�I���̂��߂̕\������������B
///�ʏ��draw�Ƃ̑���F///Color�Ƃ���m_bufferID��p���Asize�����ȊO��
///attributes�̏����inormal, texture, color)�����Ȃ��B
virtual void selectionDraw(MGCL::VIEWMODE viewMode=MGCL::DONTCARE);

///DrawType���Z�b�g����
void setDrawType(mgGLSL::DrawType drawType){m_drawType=drawType;};
mgGLSL::DrawType getDrawType()const{return m_drawType;};

///Texture��set����B
void setTexture(mgTexture* texture);

///Texture��get����B
const mgTexture* getTexture()const;
mgTexture* getTexture();

///Test if a texture is binded.
bool texture_is_bound()const;

///�쐬�ς�mgVBOLeaf��ID�̈ʒu��Vertex data��P�̃f�[�^�ɒu��������B
///0<= ID <getVerticesNumber();
void updateVertex(unsigned ID, const MGPosition& P);

///�쐬�ς�mgVBOLeaf��startID�̈ʒu����num��Vertex data��Ps�̃f�[�^��
///�u��������B 0<= startID , startID+num<=getVerticesNumber().
///Here num=Ps.size();
void updateVertices(unsigned startID, const std::vector<MGPosition>& Ps);

///�쐬�ς�mgVBOLeaf��startID�̈ʒu����num��Vertex data��Ps�̃f�[�^��
///�u��������B 0<= startID , startID+num<=getVerticesNumber().
///Here num=Ps.length();
void updateVertices(unsigned startID, const MGBPointSeq& Ps);

///�쐬�ς�mgVBOLeaf��startID�̈ʒu����numVertices�̒��_data��area�̃f�[�^��
///�u��������B 0<= startID , startID+numVertices<=getVerticesNumber().
///area�̔z��̒�����numVertices*3�ƂȂ�
void updateVertices(unsigned startID, unsigned numVertices, const float* area);

///Vertex�̐������߂�
unsigned getVerticesNumber()const{return m_count;};

///Obtain display list name.

///0(null) �͂���mgVBOElement��mgVBOLeaf�ł���A���O�������Ȃ����Ƃ������B
///���O��mgVBO����������
unsigned getDName()const{return 0;};

///Selection�ɐݒ肷�閼�O�����߂�B

///=0�̂Ƃ��A���O�̐ݒ菈�������Ȃ��B
///0(null) �͂���mgVBOElement��mgVBOLeaf�ł���A���O�������Ȃ����Ƃ������B
///���O��mgVBO����������
virtual GLuint getSelectionName()const{return 0;};

///m_primitiveMode=GL_QUAD_STRIP, GL_TRIANGLES, GL_TRIANGLE_STRIP, GL_TRIANGLE_FAN
///�̂Ƃ��ɂ̂ݗL���ł���PolygonMode�iGL_POINT, GL_LINE, GL_FILL)���w�肷��.
///setPolygonMode���Ă��Ȃ�mgVBOLeaf��GL_FILL�Ƃ����B
void setPolygonMode(GLenum mode);

///Set color.
void setStaticAttribColor(const MGColor& color){
	m_color=color;
};

///Set primitive mode.
void setPrimitiveMode(GLenum primitiveMode){m_primitiveMode=primitiveMode;};

///Set light mode. mode=-1:undefined, =0:disabled, =1:enabled.
void setLightMode(int mode);

///Line stipple�������Z�b�g����B

///When factor=0 is input, line pattern is disabled. This means the line is solid.
///When factor<0, the stipple attribute is undefined. This means the attribute
///is defined by the environment.
///When factor<=0, pattern is unnecessary.
void setLineStipple(short int factor, GLushort pattern);

///Set line width.
void setStaticAttribLineWidth(GLfloat size);///size<=0. ��undefined������

///Set Point size.
void setStaticAttribPointSize(GLfloat size);///size<=0. ��undefined������

protected:
	

private:
	void execStaticAttrib(MGCL::VIEWMODE viewMode, bool selection=false);

};

/** @} */ // end of DisplayHandling group
#endif // !defined(_MGVBOLEAF__INCLUDED_)
