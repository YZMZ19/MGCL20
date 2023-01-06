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

///MGGroupなどでmemberのMGGelの描画をmgVBOのelementとするためのクラス.

///MGAttribedGelに対する描画はMGAttribedGelにunique_ptr<mgVBO>として保持される。
///このpointerを対象としてmgVBOPointerを作成しmgVBOのメンバーとする。
///mgVBOPointerはvboを参照するのみ。
class MG_DLL_DECLR mgVBOPointer:public mgVBOElement{
public:

///MGAttribedGel用のconstructor.
///mgVBOPointerはvboを参照するのみ
mgVBOPointer(mgVBO& vbo):m_vbo(&vbo){;};

~mgVBOPointer(){;};

///このmgVBOElementがnull(いまだdraw/make_display_list()処理されていない）かを判定
///mgVBOPointerはm_vboに対するmake_display_list()処理がなされていないときfalseが返される
bool is_made(MGCL::VIEWMODE viewMode=MGCL::DONTCARE){return m_vbo->is_made(viewMode);};


///すでに作成済みであっても強制的に再作成を行う。
void make_display_list(MGCL::VIEWMODE vmode=MGCL::DONTCARE){
	m_vbo->make_display_list(vmode);
};

///描画関数draw()は、is_made()==falseであれば、作成し、表示処理をする。
///is_made()(描画データ作成済み）であれば、すでに作成されたmgVBOElementの描画を行う。
void draw(MGCL::VIEWMODE viewMode=MGCL::DONTCARE);

///draw()はmgVBOLeafが作成済み（not null)であれば作成処理を行わないが、
///redraw()は強制的に再作成を行い描画処理をおこなう。
void redraw(MGCL::VIEWMODE viewMode=MGCL::DONTCARE){m_vbo->redraw(viewMode);};
	
///描画関数selectionDraw()は、Object選択のための表示処理をする。
///通常のdrawとの相違：///Colorとしてm_bufferIDを用い、size処理以外の
///attributesの処理（normal, texture, color)をしない。
virtual void selectionDraw(MGCL::VIEWMODE viewMode=MGCL::DONTCARE);

///Obtain display list name.
///0(null) はこのmgVBOElementはmgVBOLeafであり、名前をもたないことを示す。
///名前はmgVBOだけが持つ
unsigned getDName()const{return 0;};

///Selectionに設定する名前を求める。=0のとき、名前の設定処理をしない。
///0(null) はこのmgVBOElementはmgVBOLeafであり、名前をもたないことを示す。
///名前はmgVBOだけが持つ
virtual GLuint getSelectionName()const{return m_vbo->getSelectionName();};

mgVBO* vboPointer(){return m_vbo;};

protected:

private:
	///mgVBOPointerはvboを参照するのみ
	mgVBO* m_vbo;
};

/////////////////////////////////////////////////////////////////////////////
// mgVBOLeafPointer

///すでに作成済みのmgVBOLeafをVBOのメンバー（element）として保持するためのクラス.

///文字など、mgVBOLeafをすでに作成済みのものをVBOのメンバー（element）として描画するためのクラス.
///このmgVBOLeafのpointerを対象としてmgVBOLeafPointerを作成しmgVBOのメンバーとする。
///mgVBOLeafPointerはmgVBOLeafを参照するのみ。そのinstanceは利用者の管理となる。
class MG_DLL_DECLR mgVBOLeafPointer:public mgVBOElement{
public:

///MGAttribedGel用のconstructor.
///mgVBOLeafPointerはvboを参照するのみ
mgVBOLeafPointer(const mgVBOLeaf& leaf);

~mgVBOLeafPointer(){;};

///このmgVBOElementがnull(いまだdraw/make_display_list()処理されていない）かを判定
///mgVBOLeafPointerは常にtrueが返される
bool is_made(MGCL::VIEWMODE viewMode=MGCL::DONTCARE){return true;};

///すでに作成済みであっても強制的に再作成を行う。
///mgVBOLeafPointerはなにもしない。その生成に関しては利用者の責任となる
void make_display_list(MGCL::VIEWMODE vmode=MGCL::DONTCARE){;};

///描画関数draw()は、is_made()==falseであれば、作成し、表示処理をする。
///is_made()(描画データ作成済み）であれば、すでに作成されたmgVBOElementの描画を行う。
void draw(MGCL::VIEWMODE viewMode=MGCL::DONTCARE);

///draw()はmgVBOLeafが作成済み（not null)であれば作成処理を行わないが、
///redraw()は強制的に再作成を行い描画処理をおこなう。
///mgVBOLeafPointerではdrawと同じ。その生成に関しては利用者の責任となる
void redraw(MGCL::VIEWMODE viewMode=MGCL::DONTCARE){draw(viewMode);};
	
///描画関数selectionDraw()は、Object選択のための表示処理をする。
///通常のdrawとの相違：///Colorとしてm_bufferIDを用い、size処理以外の
///attributesの処理（normal, texture, color)をしない。
void selectionDraw(MGCL::VIEWMODE viewMode=MGCL::DONTCARE);

///Obtain display list name.
///0(null) はこのmgVBOElementはmgVBOLeafであり、名前をもたないことを示す。
///名前はmgVBOだけが持つ
unsigned getDName()const{return 0;};

///Selectionに設定する名前を求める。=0のとき、名前の設定処理をしない。
///0(null) はこのmgVBOElementはmgVBOLeafであり、名前をもたないことを示す。
///名前はmgVBOだけが持つ
GLuint getSelectionName()const{return 0;};

const mgVBOLeaf* leafPointer()const{return m_VBOLeaf;};
mgVBOLeaf* leafPointer(){return m_VBOLeaf;};

protected:

private:
	///mgVBOLeafPointerはvboを参照するのみ
	mgVBOLeaf* m_VBOLeaf;
};

/** @} */ // end of DisplayHandling group
#endif // !defined(_MGVBOPOINTER__INCLUDED_)
