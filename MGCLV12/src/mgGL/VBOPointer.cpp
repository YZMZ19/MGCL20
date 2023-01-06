/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#include "StdAfx.h"
#include "mgGL/VBOPointer.h"
#include "mgGL/VBOLeaf.h"

////////////////////////  mgVBOPointer  //////////////////////////////

///mgVBOElement:
///MGGroup�Ȃǂ�member��MGGel�̕`���mgVBO��element�Ƃ��邽�߂̃N���X
///MGAttribedGel�ɑ΂���`���MGAttribedGel��unique_ptr<mgVBO>�Ƃ��ĕێ������B
///����pointer��ΏۂƂ���mgVBOPointer���쐬��mgVBO�̃����o�[�Ƃ���B
///mgVBOPointer��vbo���Q�Ƃ���̂݁B
///class mgVBOPointer:public mgVBOElement
namespace{	GLenum glErr;};

///�`��֐�draw()�́Ais_made()==false�ł���΁A�쐬���A�\������������B
///is_made()(�`��f�[�^�쐬�ς݁j�ł���΁A���łɍ쐬���ꂽmgVBOElement�̕`����s���B
void mgVBOPointer::draw(MGCL::VIEWMODE viewMode){
	if(is_display())
		m_vbo->draw(viewMode);
}
	
///�`��֐�selectionDraw()�́AObject�I���̂��߂̕\������������B
///�ʏ��draw�Ƃ̑���F///Color�Ƃ���m_bufferID��p���Asize�����ȊO��
///attributes�̏����inormal, texture, color)�����Ȃ��B
void mgVBOPointer::selectionDraw(MGCL::VIEWMODE viewMode){
	if(is_display())
		m_vbo->selectionDraw(viewMode);
}

/////////////////////////
///mgVBOLeafPointer

///MGAttribedGel�p��constructor.
///mgVBOLeafPointer��vbo���Q�Ƃ���̂�
mgVBOLeafPointer::mgVBOLeafPointer(const mgVBOLeaf& leaf){
	mgVBOLeaf* leafP=const_cast<mgVBOLeaf*>(&leaf);
	m_VBOLeaf=leafP;
}

///�`��֐�draw()�́Ais_made()==false�ł���΁A�쐬���A�\������������B
///is_made()(�`��f�[�^�쐬�ς݁j�ł���΁A���łɍ쐬���ꂽmgVBOElement�̕`����s���B
void mgVBOLeafPointer::draw(MGCL::VIEWMODE viewMode){
	m_VBOLeaf->draw(viewMode);
}
	
///�`��֐�selectionDraw()�́AObject�I���̂��߂̕\������������B
///�ʏ��draw�Ƃ̑���F///Color�Ƃ���m_bufferID��p���Asize�����ȊO��
///attributes�̏����inormal, texture, color)�����Ȃ��B
void mgVBOLeafPointer::selectionDraw(MGCL::VIEWMODE viewMode){
	if(is_display())
		m_VBOLeaf->selectionDraw(viewMode);
}
