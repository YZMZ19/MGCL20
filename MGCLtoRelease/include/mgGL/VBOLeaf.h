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

///mgVBOに対して描画データ作成後の情報を保持するためのクラス.

///mgVBOはmgVBOLeafとmgVBO自身をelementとして保持する。
///mgVBOをelementとするのはMGGroupまたはMGShellのメンバーデータの表示情報を
///保持するため
class MG_DLL_DECLR mgVBOLeaf:public mgVBOElement{

///Maximum line width.

///Used to judge the polygon mode of the polygon if this mgVBOLeaf is a triangle.
///When polygon mode = GL_LINE, m_size=the line width.
///When polygon mode = GL_POINT, the point size=m_size-m_PolygonModePointSizeBase.
///That is, if m_size>m_PolygonModePointSizeBase, its polygon mode is GL_POINT and
///the point size=m_size-m_PolygonModePointSizeBase.
static const GLfloat m_PolygonModePointSizeBase;

private:

	mgGLSL::DrawType m_drawType;/// 描画モード mgGLSL::DrawTypeの値が入る

	GLuint m_bufferID;///<glGenBuffersでで作成されたBuffer id。未作成のときは=0.
				///m_bufferIDはm_vertexArrayIDにBindされる。

	///頂点ごとの(Staticでない）値が指定されているか否かを示す
	bool m_ColorSpecified:1, m_NormalSpecified:1, m_TextureSpecified:1;

	mutable GLuint m_vertexArrayID;///<glGenVertexArraysで作成されたVertexArray.未作成のときは=0.

	///(m_primitiveMode,m_count)は簡易表示(VIEWMODE==WIRE)または強調表示(VIEWMODE==HIGHLIGHT)に利用される。
	///(m_primitiveMode,m_count)はglDrawArraysにm_bufferIDとfirst=0で使用される。
	GLenum m_primitiveMode;///<glDrawArraysに入力とするprimitive種類
	unsigned m_count;///<glDrawArraysに入力とするVertexの数

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

	///light modeはm_elementsShadeに対してのみ有効。m_elementsに対しては常にlightはオフ
	int m_lightMode;///<  <0: undefined, =0:Light is disabled, >0:Light is enabled.

friend class mgVBO;

public:

mgVBOLeaf();
mgVBOLeaf(GLfloat size,MGColor& color);
mgVBOLeaf(const mgVBOLeaf& vbol);
mgVBOLeaf(const mgVBOLeafBuilder& builder);
~mgVBOLeaf();

///このmgVBOElementがnull(いまだdraw/make_display_list()処理されていない）かを判定

///mgVBOはm_gelに対するmake_display_list()処理がなされていないときfalseが返される
///mgVBOLeafはmake_display_list()されたときのみ存在するため常にtrueが返る。
bool is_made(MGCL::VIEWMODE viewMode=MGCL::DONTCARE){return true;};

///Remake the display list.

///すでに作成済みであっても強制的に再作成を行う。
///mgVBOLeafは作成済みの場合のみ存在するため何もしない。
void make_display_list(MGCL::VIEWMODE vmode=MGCL::DONTCARE){;};

///Draw.

///描画関数draw()は、is_made()であれば、作成し、表示処理をする。
///!is_made()(描画データ作成済み）であれば、すでに作成されたmgVBOElementの描画を行う。
///mgVBOLeafでは常に作成すみであり表示処理をおこなう。
void draw(MGCL::VIEWMODE viewMode=MGCL::DONTCARE);

///Redraw.

///draw()はmgVBOLeafが作成済み（not null)であれば作成処理を行わないが、
///redraw()は強制的に再作成を行い描画処理をおこなう。
///mgVBOLeafでは常に作成すみであり表示処理をおこなう。
void redraw(MGCL::VIEWMODE viewMode=MGCL::DONTCARE){draw(viewMode);};

///Draw for the Selection.

///描画関数selectionDraw()は、Object選択のための表示処理をする。
///通常のdrawとの相違：///Colorとしてm_bufferIDを用い、size処理以外の
///attributesの処理（normal, texture, color)をしない。
virtual void selectionDraw(MGCL::VIEWMODE viewMode=MGCL::DONTCARE);

///DrawTypeをセットする
void setDrawType(mgGLSL::DrawType drawType){m_drawType=drawType;};
mgGLSL::DrawType getDrawType()const{return m_drawType;};

///Textureをsetする。
void setTexture(mgTexture* texture);

///Textureをgetする。
const mgTexture* getTexture()const;
mgTexture* getTexture();

///Test if a texture is binded.
bool texture_is_bound()const;

///作成済のmgVBOLeafのIDの位置のVertex dataをPのデータに置き換える。
///0<= ID <getVerticesNumber();
void updateVertex(unsigned ID, const MGPosition& P);

///作成済のmgVBOLeafのstartIDの位置からnum個のVertex dataをPsのデータに
///置き換える。 0<= startID , startID+num<=getVerticesNumber().
///Here num=Ps.size();
void updateVertices(unsigned startID, const std::vector<MGPosition>& Ps);

///作成済のmgVBOLeafのstartIDの位置からnum個のVertex dataをPsのデータに
///置き換える。 0<= startID , startID+num<=getVerticesNumber().
///Here num=Ps.length();
void updateVertices(unsigned startID, const MGBPointSeq& Ps);

///作成済のmgVBOLeafのstartIDの位置からnumVertices個の頂点dataをareaのデータに
///置き換える。 0<= startID , startID+numVertices<=getVerticesNumber().
///areaの配列の長さ＝numVertices*3となる
void updateVertices(unsigned startID, unsigned numVertices, const float* area);

///Vertexの数を求める
unsigned getVerticesNumber()const{return m_count;};

///Obtain display list name.

///0(null) はこのmgVBOElementはmgVBOLeafであり、名前をもたないことを示す。
///名前はmgVBOだけが持つ
unsigned getDName()const{return 0;};

///Selectionに設定する名前を求める。

///=0のとき、名前の設定処理をしない。
///0(null) はこのmgVBOElementはmgVBOLeafであり、名前をもたないことを示す。
///名前はmgVBOだけが持つ
virtual GLuint getSelectionName()const{return 0;};

///m_primitiveMode=GL_QUAD_STRIP, GL_TRIANGLES, GL_TRIANGLE_STRIP, GL_TRIANGLE_FAN
///のときにのみ有効でそのPolygonMode（GL_POINT, GL_LINE, GL_FILL)を指定する.
///setPolygonModeしていないmgVBOLeafはGL_FILLとされる。
void setPolygonMode(GLenum mode);

///Set color.
void setStaticAttribColor(const MGColor& color){
	m_color=color;
};

///Set primitive mode.
void setPrimitiveMode(GLenum primitiveMode){m_primitiveMode=primitiveMode;};

///Set light mode. mode=-1:undefined, =0:disabled, =1:enabled.
void setLightMode(int mode);

///Line stipple属性をセットする。

///When factor=0 is input, line pattern is disabled. This means the line is solid.
///When factor<0, the stipple attribute is undefined. This means the attribute
///is defined by the environment.
///When factor<=0, pattern is unnecessary.
void setLineStipple(short int factor, GLushort pattern);

///Set line width.
void setStaticAttribLineWidth(GLfloat size);///size<=0. はundefinedを示す

///Set Point size.
void setStaticAttribPointSize(GLfloat size);///size<=0. はundefinedを示す

protected:
	

private:
	void execStaticAttrib(MGCL::VIEWMODE viewMode, bool selection=false);

};

/** @} */ // end of DisplayHandling group
#endif // !defined(_MGVBOLEAF__INCLUDED_)
