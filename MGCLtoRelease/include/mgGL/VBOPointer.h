/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#if !defined(_MGVBOPOINTER__INCLUDED_)
#define _MGVBOPOINTER__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "StdAfx.h"
#include "mg/drawParam.h"
#include "mgGL/Color.h"
#include "mgGL/GLAttrib.h"
#include "mgGL/VBO.h"
#include "mgGL/VBOLeaf.h"

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

/** @file */
/** @addtogroup DisplayHandling
 *  @{
 */

/////////////////////////////////////////////////////////////////////////////
// mgVBOPointer

///MGGroup�Ȃǂ�member��MGGel�̕`���mgVBO��element�Ƃ��邽�߂̃N���X.

///MGAttribedGel�ɑ΂���`���MGAttribedGel��unique_ptr<mgVBO>�Ƃ��ĕێ������B
///����pointer��ΏۂƂ���mgVBOPointer���쐬��mgVBO�̃����o�[�Ƃ���B
///mgVBOPointer��vbo���Q�Ƃ���̂݁B
class MG_DLL_DECLR mgVBOPointer:public mgVBOElement{
public:

///MGAttribedGel�p��constructor.
///mgVBOPointer��vbo���Q�Ƃ���̂�
mgVBOPointer(mgVBO& vbo):m_vbo(&vbo){;};

~mgVBOPointer(){;};

///����mgVBOElement��null(���܂�draw/make_display_list()��������Ă��Ȃ��j���𔻒�
///mgVBOPointer��m_vbo�ɑ΂���make_display_list()�������Ȃ���Ă��Ȃ��Ƃ�false���Ԃ����
bool is_made(MGCL::VIEWMODE viewMode=MGCL::DONTCARE){return m_vbo->is_made(viewMode);};


///���łɍ쐬�ς݂ł����Ă������I�ɍč쐬���s���B
void make_display_list(MGCL::VIEWMODE vmode=MGCL::DONTCARE){
	m_vbo->make_display_list(vmode);
};

///�`��֐�draw()�́Ais_made()==false�ł���΁A�쐬���A�\������������B
///is_made()(�`��f�[�^�쐬�ς݁j�ł���΁A���łɍ쐬���ꂽmgVBOElement�̕`����s���B
void draw(MGCL::VIEWMODE viewMode=MGCL::DONTCARE);

///draw()��mgVBOLeaf���쐬�ς݁inot null)�ł���΍쐬�������s��Ȃ����A
///redraw()�͋����I�ɍč쐬���s���`�揈���������Ȃ��B
void redraw(MGCL::VIEWMODE viewMode=MGCL::DONTCARE){m_vbo->redraw(viewMode);};
	
///�`��֐�selectionDraw()�́AObject�I���̂��߂̕\������������B
///�ʏ��draw�Ƃ̑���F///Color�Ƃ���m_bufferID��p���Asize�����ȊO��
///attributes�̏����inormal, texture, color)�����Ȃ��B
virtual void selectionDraw(MGCL::VIEWMODE viewMode=MGCL::DONTCARE);

///Obtain display list name.
///0(null) �͂���mgVBOElement��mgVBOLeaf�ł���A���O�������Ȃ����Ƃ������B
///���O��mgVBO����������
unsigned getDName()const{return 0;};

///Selection�ɐݒ肷�閼�O�����߂�B=0�̂Ƃ��A���O�̐ݒ菈�������Ȃ��B
///0(null) �͂���mgVBOElement��mgVBOLeaf�ł���A���O�������Ȃ����Ƃ������B
///���O��mgVBO����������
virtual GLuint getSelectionName()const{return m_vbo->getSelectionName();};

mgVBO* vboPointer(){return m_vbo;};

protected:

private:
	///mgVBOPointer��vbo���Q�Ƃ���̂�
	mgVBO* m_vbo;
};

/////////////////////////////////////////////////////////////////////////////
// mgVBOLeafPointer

///���łɍ쐬�ς݂�mgVBOLeaf��VBO�̃����o�[�ielement�j�Ƃ��ĕێ����邽�߂̃N���X.

///�����ȂǁAmgVBOLeaf�����łɍ쐬�ς݂̂��̂�VBO�̃����o�[�ielement�j�Ƃ��ĕ`�悷�邽�߂̃N���X.
///����mgVBOLeaf��pointer��ΏۂƂ���mgVBOLeafPointer���쐬��mgVBO�̃����o�[�Ƃ���B
///mgVBOLeafPointer��mgVBOLeaf���Q�Ƃ���̂݁B����instance�͗��p�҂̊Ǘ��ƂȂ�B
class MG_DLL_DECLR mgVBOLeafPointer:public mgVBOElement{
public:

///MGAttribedGel�p��constructor.
///mgVBOLeafPointer��vbo���Q�Ƃ���̂�
mgVBOLeafPointer(const mgVBOLeaf& leaf);

~mgVBOLeafPointer(){;};

///����mgVBOElement��null(���܂�draw/make_display_list()��������Ă��Ȃ��j���𔻒�
///mgVBOLeafPointer�͏��true���Ԃ����
bool is_made(MGCL::VIEWMODE viewMode=MGCL::DONTCARE){return true;};

///���łɍ쐬�ς݂ł����Ă������I�ɍč쐬���s���B
///mgVBOLeafPointer�͂Ȃɂ����Ȃ��B���̐����Ɋւ��Ă͗��p�҂̐ӔC�ƂȂ�
void make_display_list(MGCL::VIEWMODE vmode=MGCL::DONTCARE){;};

///�`��֐�draw()�́Ais_made()==false�ł���΁A�쐬���A�\������������B
///is_made()(�`��f�[�^�쐬�ς݁j�ł���΁A���łɍ쐬���ꂽmgVBOElement�̕`����s���B
void draw(MGCL::VIEWMODE viewMode=MGCL::DONTCARE);

///draw()��mgVBOLeaf���쐬�ς݁inot null)�ł���΍쐬�������s��Ȃ����A
///redraw()�͋����I�ɍč쐬���s���`�揈���������Ȃ��B
///mgVBOLeafPointer�ł�draw�Ɠ����B���̐����Ɋւ��Ă͗��p�҂̐ӔC�ƂȂ�
void redraw(MGCL::VIEWMODE viewMode=MGCL::DONTCARE){draw(viewMode);};
	
///�`��֐�selectionDraw()�́AObject�I���̂��߂̕\������������B
///�ʏ��draw�Ƃ̑���F///Color�Ƃ���m_bufferID��p���Asize�����ȊO��
///attributes�̏����inormal, texture, color)�����Ȃ��B
void selectionDraw(MGCL::VIEWMODE viewMode=MGCL::DONTCARE);

///Obtain display list name.
///0(null) �͂���mgVBOElement��mgVBOLeaf�ł���A���O�������Ȃ����Ƃ������B
///���O��mgVBO����������
unsigned getDName()const{return 0;};

///Selection�ɐݒ肷�閼�O�����߂�B=0�̂Ƃ��A���O�̐ݒ菈�������Ȃ��B
///0(null) �͂���mgVBOElement��mgVBOLeaf�ł���A���O�������Ȃ����Ƃ������B
///���O��mgVBO����������
GLuint getSelectionName()const{return 0;};

const mgVBOLeaf* leafPointer()const{return m_VBOLeaf;};
mgVBOLeaf* leafPointer(){return m_VBOLeaf;};

protected:

private:
	///mgVBOLeafPointer��vbo���Q�Ƃ���̂�
	mgVBOLeaf* m_VBOLeaf;
};

/** @} */ // end of DisplayHandling group
#endif // !defined(_MGVBOPOINTER__INCLUDED_)
