/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#if !defined(_MGVBOELEMENT__INCLUDED_)
#define _MGVBOELEMENT__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "StdAfx.h"
#include "mg/drawParam.h"
#include "mgGL/Color.h"
#include "mgGL/GLAttrib.h"

class MGPosition;
class MGBox;
class MGBPointSeq;
class MGSPointSeq;
class MGCurve;
class MGSPointSeq;
class mgTL2Triangles;
class MGStl;
class MGColor;
class MGAttribedGel;
class MGComplex;
class MGCell;
class mgVBOLeafBuilder;
class mgVBOLeaf;
class mgTexture;
class mgVBOElement;

/** @file */
/** @addtogroup DisplayHandling
 *  @{
 */

///UniqueVBOElement definition.
using UniqueVBOElement = std::unique_ptr<mgVBOElement>;

/////////////////////////////////////////////////////////////////////////////
// mgVBOElement

///Interface class to include an element in mgVBO class's.
class MG_DLL_DECLR mgVBOElement{
friend class mgVBO;

public:

mgVBOElement():m_no_display(false){;};

virtual ~mgVBOElement(){;};

///����mgVBOElement��null(���܂�draw/make_display_list()��������Ă��Ȃ��j���𔻒�
///mgVBOPointer��m_vbo�ɑ΂���make_display_list()�������Ȃ���Ă��Ȃ��Ƃ�false���Ԃ����
virtual bool is_made(MGCL::VIEWMODE viewMode=MGCL::DONTCARE)=0;
	
///���łɍ쐬�ς݂ł����Ă������I�ɍč쐬���s���B
virtual void make_display_list(MGCL::VIEWMODE vmode=MGCL::DONTCARE)=0;

///�`��֐�draw()�́A!is_made()�ł���΁A�쐬���A�\������������B
///is_made()(�`��f�[�^�쐬�ς݁j�ł���΁A���łɍ쐬���ꂽmgVBOElement�̕`����s���B
virtual void draw(MGCL::VIEWMODE viewMode=MGCL::DONTCARE)=0;

///draw()��is_made(�쐬�ς�)�ł���΍쐬�������s��Ȃ����A
///redraw()�͋����I�ɍč쐬���s���`�揈���������Ȃ��B
virtual void redraw(MGCL::VIEWMODE viewMode)=0;

///�`��֐�selectionDraw()�́AObject�I���̂��߂̕\������������B
///�ʏ��draw�Ƃ̑���F///Color�Ƃ���m_bufferID��p���Asize�����ȊO��
///attributes�̏����inormal, texture, color)�����Ȃ��B
virtual void selectionDraw(MGCL::VIEWMODE viewMode=MGCL::DONTCARE)=0;

///Obtain display list name.
///0(null) �͂���mgVBOElement��mgVBOLeaf�ł���A���O�������Ȃ����Ƃ������B
///���O��mgVBO����������
virtual unsigned getDName()const=0;

///Selection�ɐݒ肷�閼�O�����߂�B=0�̂Ƃ��A���O�̐ݒ菈�������Ȃ��B
///0(null) �͂���mgVBOElement��mgVBOLeaf�ł���A���O�������Ȃ����Ƃ������B
///���O��mgVBO����������
virtual GLuint getSelectionName()const=0;

///When this is a mgVBOPointer, return the vbo pointer referenced.
virtual mgVBO* vboPointer(){return 0;};

///When this is a mgVBOLeaf, return the mgVBOLeaf pointer.
virtual const mgVBOLeaf* leafPointer(){return 0;};

///set_display/set_no_display controls if this mgVBO be displayed or not.
virtual void set_display(){	m_no_display=false;};
virtual void set_no_display(){	m_no_display=true;};

///Test if this is no display mode or not.

///True if no display, false if display mode.
virtual bool getNoDisplayMode()const{return m_no_display;};
bool is_no_display()const{return getNoDisplayMode();};
bool is_display()const{return !is_no_display();};

///Set the draw param. This is applied to all the make_display_list ofmgVBO
///after setDrawParam().
static void setDrawParam(const MGDrawParam& dpara){m_drawPara=dpara;};
static MGDrawParam& getDrawParam(){return m_drawPara;};

///Set/get hilight color.
static void setHilightColor(const MGColor& hcolor){m_hilightColor=hcolor;};
static const MGColor& getHilightColor(){return m_hilightColor;};

///Set/get default point size.
static void setDefaultPointSize(GLfloat psize){m_pointSize=psize;};
static GLfloat getDefaultPointSize(){return m_pointSize;};

protected:
	
	bool m_no_display:1;///<Controls if this mgVBO be displayed or not.
		///< =true: not display, false: display.

private:
	static MGDrawParam m_drawPara;///<draw parameter for MGAttribedGel' make_display_list.
	static MGColor m_hilightColor;///Color to hilight.
	static GLfloat m_pointSize;///<Outer Point size to draw.

};

/** @} */ // end of DisplayHandling group
#endif // !defined(_MGVBOELEMENT__INCLUDED_)
