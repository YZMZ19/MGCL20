/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#include "StdAfx.h"
#include "mgGL/VBOPointer.h"
#include "mgGL/VBOLeaf.h"

////////////////////////  mgVBOPointer  //////////////////////////////

///mgVBOElement:
///MGGroupなどでmemberのMGGelの描画をmgVBOのelementとするためのクラス
///MGAttribedGelに対する描画はMGAttribedGelにunique_ptr<mgVBO>として保持される。
///このpointerを対象としてmgVBOPointerを作成しmgVBOのメンバーとする。
///mgVBOPointerはvboを参照するのみ。
///class mgVBOPointer:public mgVBOElement
namespace{	GLenum glErr;};

///描画関数draw()は、is_made()==falseであれば、作成し、表示処理をする。
///is_made()(描画データ作成済み）であれば、すでに作成されたmgVBOElementの描画を行う。
void mgVBOPointer::draw(MGCL::VIEWMODE viewMode){
	if(is_display())
		m_vbo->draw(viewMode);
}
	
///描画関数selectionDraw()は、Object選択のための表示処理をする。
///通常のdrawとの相違：///Colorとしてm_bufferIDを用い、size処理以外の
///attributesの処理（normal, texture, color)をしない。
void mgVBOPointer::selectionDraw(MGCL::VIEWMODE viewMode){
	if(is_display())
		m_vbo->selectionDraw(viewMode);
}

/////////////////////////
///mgVBOLeafPointer

///MGAttribedGel用のconstructor.
///mgVBOLeafPointerはvboを参照するのみ
mgVBOLeafPointer::mgVBOLeafPointer(const mgVBOLeaf& leaf){
	mgVBOLeaf* leafP=const_cast<mgVBOLeaf*>(&leaf);
	m_VBOLeaf=leafP;
}

///描画関数draw()は、is_made()==falseであれば、作成し、表示処理をする。
///is_made()(描画データ作成済み）であれば、すでに作成されたmgVBOElementの描画を行う。
void mgVBOLeafPointer::draw(MGCL::VIEWMODE viewMode){
	m_VBOLeaf->draw(viewMode);
}
	
///描画関数selectionDraw()は、Object選択のための表示処理をする。
///通常のdrawとの相違：///Colorとしてm_bufferIDを用い、size処理以外の
///attributesの処理（normal, texture, color)をしない。
void mgVBOLeafPointer::selectionDraw(MGCL::VIEWMODE viewMode){
	if(is_display())
		m_VBOLeaf->selectionDraw(viewMode);
}
