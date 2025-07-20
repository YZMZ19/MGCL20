/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#pragma once

class MGColor;
class mgTL2Triangles;
#include "StdAfx.h"
#include "mg/drawParam.h"
#include "mgGL/GLAttrib.h"
#include "mgGL/Color.h"

#ifndef _CONSOLE

#include "mgGL/Appearance.h"
#include "mgGL/glslprogram.h"
#include "mgGL/StaticGLAttrib.h"
#include "mgGL/VBOElement.h"
#include "mgGL/VBOLeafBuilder.h"

class MGPosition;
class MGBox;
class MGBPointSeq;
class MGSPointSeq;
class MGCurve;
class MGSPointSeq;
class MGCell;
class MGStl;
class MGAttribedGel;
class mgVBO;
class MGComplex;
class MGImage;
class MGObject;
class MGGelPosition;

/** @file */
/** @addtogroup DisplayHandling
 *  @{
 */

/// Vertex Buffer Object Class, a class to draw as a display list in OpenGL4 

/// 基本的な描画のためのツールとOpenGL4を使用した描画機構を提供する。
///
///******利用上の手順******
///(1) MGAttribedGel用のmgVBOは、mgVBO(const MGAttribedGel&)のconstructorを利用して構築し、
///   必要に応じてdrawXXXまたはBegin-EndによりmgVBOに追加のmgVBOLeafを作成しm_gelの描画データを作成する。
///　 非MGAttribedGelのmgVBOはdrawXXX, begin()-end()を利用して描画データ作成処理(make_display_list())を行う。
///   作成が終了したらsetDirty(false)により、作成完了を宣言する。
///   mgVBOはsetDirty(false)により描画データ作成完了宣言されていない場合にdraw()が呼ばれたとき、
///   必ずmake_dipslay_list()を呼び出す。
///(2) draw()により描画処理をする。
///   MGAttribedGel用のmgVBOに対してはmake_display_list()をせず、いきなりdraw()を呼び出してもよい。
///   setDirty(false)が発信されていない場合、mgVBOはmake_display_list()を発信し、描画データ
///   作成処理を行ったのち、drawを行う。
///
///   draw()/make_display_list()を呼ぶ前には種類の異なる複数の図形の作成を行っても
///   よい(drawXXX, begin()-end())が、draw()/make_display_list()までに作成した図形は
///   draw()/make_display_list()によりひとまとめの扱いとなる。
///   描画関数draw()はメンバーデータの各mgVBOElementに対してdraw()を発信する。
///
/// Begin()-Bnd()で囲まれた場合、必ずひとつのmgVBOLeafが作成される。ひとつであることは
/// 保障される。すべてのdrawXXXは、一つまたは一つ以上のmgVBOLeafをdraw()のために
/// 作成する（make_display_list）だけで実際の描画は行わない。
///
///mgVBOはmgVBOLeafとmgVBO自身(VBOPointer)をelementとして保持する。
///mgVBOをelementとするのはMGGroupまたはMGShellのメンバーデータの表示情報を
///保持するため。これ以外のMGObjectはmgVBOLeaf(一般には複数）をelementとして保持する
class MG_DLL_DECLR mgVBO:public mgVBOElement{

friend class MGDNameControl;
friend class MGAttribedGel;
friend class mgVBOPointer;
friend class MGOpenGLView;

public:

///Define iterator.
typedef std::vector<UniqueVBOElement> container_type;
typedef container_type::iterator iterator;
typedef container_type::const_iterator const_iterator;

protected:

	unsigned m_dname;///<name to register in MGDNameControl.
	MGAttribedGel* m_gel;///<MGAttribedGel of this mgVBO.
		///<When m_gel=0, this mgVBO is not for an MGAttribedGel but a temporary or non gel's VBO.

	std::vector<UniqueVBOElement>* m_target_elements;//=&m_elements, or &m_elementsShade.

	mutable bool m_elementsDirty:1;///<  =true if dirty and need to remake.
	std::vector<UniqueVBOElement> m_elements;///<VBO data for drawWire or highlight. Shading data
		///<is stored in m_elementsShade if exist.
		///<When m_elements.size()=0, this means this VBO's display data is not made yet.
		///<m_gel's display list will be made to draw elements by calling m_gel->make_display_list().
		///<When m_elements.size()!=0, the display lists are already made and m_elements includes
		///<the elements.

	mutable bool m_elementsShadeDirty:1;///<   =true if dirty and need to remake.
	std::vector<UniqueVBOElement> m_elementsShade;///< VBO data for shading. Ths highlight data
		///<is stored in m_elements.
		///<When m_elementsShade.size()=0 and m_gel->manifold_dimension()>=2, this means
		///<m_gel's shading data is not made yet.

	MGColor m_colorStatic;///<Color specified by setStaticAttribColor().
		///<The color of the following Begin()(mgVBOLeaf generated) is set to this color.

	GLfloat m_lineWidthStatic;///<The line width specified by setStaticAttribLineWidth.
		///<The size of the following Begin()(mgVBOLeaf generated) is set to this size.
		///<When m_lineWidthStatic<=0., the line width is undefined.

	GLfloat m_pointSizeStatic;///<The point size specified by setStaticAttribPointSize.
		///<The size of the following Begin()(mgVBOLeaf generated) is set to this size.
		///<When m_pointSizeStatic<=0., the point size is undefined.

	short int m_stippleFactor;///<Line stipple factor.
		///<If m_stippleFactor=0, line Stipple is disabled.
		///<   m_stippleFactor<0, line stipple is undefined.
	GLushort m_LineStipplePattern;///< Indicates the pattern of stipple.

	///light modeはm_elementsShadeに対してのみ有効。m_elementsに対しては常にlightはオフ
	int m_lightMode;///< <0: undefined, =0:Light is disabled, >0:Light is enabled.

	/// CoorinateType
	mgGLSL::CoordinateType m_coordinateType;

	///Buffer to store data in Begin() to End().
	std::unique_ptr<mgVBOLeafBuilder> m_builder;

public:

///非MGAttribedGel用のconstructor.
mgVBO(mgGLSL::CoordinateType coordType= mgGLSL::CoordinateType::World);

///MGAttribedGel用のconstructor.
mgVBO(const MGAttribedGel& gel);

///Copy constructor.
mgVBO(const mgVBO& vbo);

///Assignment.
mgVBO& operator=(const mgVBO& vbo);

virtual ~mgVBO();

///VBO dataを初期化する。

///(1) display modeを表示にする
///(2) viewModeに従い、clearElements()をする
///(3) clearStaticAttributes()発信
///(4) gel用のvboであれば、drawAttrib()発信
void initializeVBO(MGCL::VIEWMODE viewMode=MGCL::DONTCARE);

///描画関数draw()は、!is_made()(dirty)であれば、make_display_list()を呼び表示データ作成処理をする。

///is_made()(not dirty, 描画データ作成済み）であれば、make_display_list()を呼ばず、
///すでに作成されたmgVBOElementの描画を行う。
///MGAttribedGelはすべて専用のVBOを所有し、draw()の呼び出しにより自動的に描画データ作成される。
///MGAttribedGel以外は上記のmake_display_list()をoverrideし自分の描画を行う。
virtual void draw(MGCL::VIEWMODE viewMode=MGCL::DONTCARE);

///draw()はmgVBOLeafが作成済み（not null)であれば作成処理を行わないが、
///redraw()は強制的に再作成を行い描画処理をおこなう。
virtual void redraw(MGCL::VIEWMODE viewMode=MGCL::DONTCARE);

///描画関数selectionDraw()は、Object選択のための表示処理をする。

///通常のdrawとの相違：///Colorとしてm_bufferIDを用い、size処理以外の
///attributesの処理（normal, texture, color)をしない。
virtual void selectionDraw(MGCL::VIEWMODE viewMode=MGCL::DONTCARE);

/// <summary>
/// Given top Group(File MGGroup) mgVBO, build MGGelPosition of this,
/// This must be built by MGOpenGLView.
/// </summary>
void buildGelPosition(mgVBO* vboGroup, MGGelPosition& gelp);

///highlight属性で表示する
void highlight();

///Clear all the data.
virtual void clearElements(MGCL::DRAW_TARGET target= MGCL::WIRE_AND_SHADING);

///Static Attributes をすべてdefaultにもどす(display/noDisplayは対象外）
virtual void clearStaticAttributes();

///Get the reference of the last element.
UniqueVBOElement& back()const{return m_target_elements->back();};

///Get the reference of the fast element.
UniqueVBOElement& front()const{return m_target_elements->front();};

///Get the iterator of the fast element.
iterator begin_element(){return m_target_elements->begin();};

///Get the end() iterator.
iterator end_element(){return m_target_elements->end();};

///Get the begin() iterator.
const_iterator begin_element()const{return m_target_elements->begin();};

///Get the end() iterator.
const_iterator end_element()const{return m_target_elements->end();};

///Get the number of elements.
int elementNumber()const{return (int)m_target_elements->size();};

///Pop back the element.
void pop_back_element(){m_target_elements->pop_back();};

///Push back an element.
void push_back_element(mgVBOElement* elm){m_target_elements->emplace_back(elm);};

///Obtain display list name.
///0(null) はこのmgVBOElementはmgVBOLeafであり、名前をもたないことを示す。
///名前はmgVBOだけが持つ
unsigned getDName()const{return m_dname;};

///Get the target gle of this.
MGAttribedGel* gel()const{return m_gel;};

///Selectionに設定する名前を求める。

///=0のとき、名前の設定処理をしない。
///0(null) はこのmgVBOElementはmgVBOLeafであり、名前をもたないことを示す。
///名前はmgVBOだけが持つ
virtual GLuint getSelectionName()const;

///このmgVBOElementがnull(いまだdraw/make_display_list()処理されていない）かを判定
///mgVBOはm_gelに対するmake_display_list()処理がなされていないときfalseが返される
virtual bool is_made(MGCL::VIEWMODE viewMode=MGCL::DONTCARE);

///現在Begin()されているかどうかを調べる。

///Begin()され、End()の前であればtrueを返す
bool is_InBegin();

///typeのmgVBOLeafをひとつ(Begin-Endでひとつ）作成して追加する。

///targetはBegin()-End()で作成されるmgVBOLeafをどのelementに入れるかを示す。
///target=WIRE:Wire(highlight用のelement(m_elements)
///target=SHADING:Shading用のelement(m_elementsShade).
///typeは次のものを可とする(GL_xxxxのxxxxを示す）
///POINTS,LINES,LINE_STRIP,LINE_LOOP,TRIANGLE_FAN,TRIANGLE_STRIP,QUAD_STRIP.
///Begin()--End()でひとつのmgVBOLeafとされることが保証される。
///ここで作成されるmgVBOLeafのstaticAttributesは現在このmgVBOに設定されている
///属性とされる。
virtual void Begin(GLenum type, MGCL::DRAW_TARGET target= MGCL::WIRE);

///Begin()で始めたVBOをひとつのmgVBOLeafとして作成、追加する。

///polygonModeはBegin()のtypeがGL_QUAD_STRIP, GL_TRIANGLES, GL_TRIANGLE_STRIP, GL_TRIANGLE_FAN
///の際にのみ有効で、そのPolygonMode（GL_POINT, GL_LINE, GL_FILL)を指定する。
///作成されたmgVBOElement*が返されるがこのmgVBOElementはmgVBOが所有している。
///nullが返却されたとき、作成に失敗。
virtual mgVBOLeaf* End(GLenum polygonMode=GL_FILL);

///Set dirty flag(s) of m_elements or m_elementsShade.

///is_dirty is true if dirty and need to remake.
void setDirty(bool is_dirty=true, MGCL::DRAW_TARGET target= MGCL::WIRE_AND_SHADING)const;

///Function IDをセットする.

///setFunctionID()はBegin()-End()の間でのみ有効。
//void setFunctionID(int functionID);
void setDrawType(mgGLSL::DrawType drawType);

///Textureをsetする。

///textureは参照するのみ
///setTexture()はBegin()-End()の間でのみ有効。
void setTexture(mgTexture* texture);

///頂点ごとの色を指定する
///Begin()-End()の間でのみ有効。
void Color(const MGColor& colr);
void Color3fv(const float colr[3]);
void Color3dv(const double colr[3]);
void Color4fv(const float colr[4]);
void Color4ubv(const unsigned char rgba[4]);

///頂点ごとのnormalを指定する
///Begin()-End()の間でのみ有効。
void Normal(const MGVector& norml);
void Normal(float x, float y, float z);
void Normal3d(double x, double y, double z);
void Normal3fv(const float norml[3]);
void Normal3dv(const double norml[3]);

///頂点の座標値を指定する
///Begin()-End()の間でのみ有効。
virtual void Vertex(const MGPosition& v);
virtual void Vertex(float x, float y, float z=0.0f);
virtual void Vertex3d(double x, double y, double z=0.0);
virtual void Vertex2fv(const float v[2]);
virtual void Vertex3fv(const float v[3]);
virtual void Vertex2dv(const double v[2]);
virtual void Vertex3dv(const double v[3]);

///頂点ごとのTextureの座標値を指定する
///Begin()-End()の間でのみ有効。
void TexCoord(const MGPosition& v);
void TexCoord(float x, float y);
void TexCoord2d(double x, double y);
void TexCoord2fv(const float v[2]);
void TexCoord2dv(const double v[2]);

///Static attributeを設定する。

///begin()-end()の間であればそのmgVBOLeafに対してだけ適用され、
///begin()-end()の外であればその後に作成される(setStaticAttribXXXをうちに持たない）
///begin()-end()のmgVBOLeafすべてに適用される。
void setStaticAttribColor(const MGColor& color);//colorにはundefinedのものも許される
void setStaticAttribColor(const float color[4]);
void setStaticAttribColor(float r, float g, float b);
void setStaticAttribColor(MGColor::ColorID id);
void setStaticAttribLineWidth(GLfloat size);///size<=0. はundefinedを示す
void setStaticAttribPointSize(GLfloat size);///size<=0. はundefinedを示す

///Line stipple属性をセットする。

///When factor=0 is input, line pattern is disabled. 実線となる
///When factor<0, the stipple attribute is undefined. This means the attribute
///is defined by the environment.
///When factor<=0, pattern is unnecessary.
void setLineStipple(short int factor, GLushort pattern=0);

///Set line width.
void LineWidth(GLfloat size){setStaticAttribLineWidth(size);};

///Line patternをdisableにして実線とする
void disableLinePattern();

///Set the draw param. This is applied to all the make_display_list ofmgVBO
///after setDrawParam().
void setDrawParam(const MGDrawParam& dpara){ mgVBOElement::setDrawParam(dpara); };

///Get the point size.
GLfloat getPointSize()const;

///Get the line width.
GLfloat getLineWidth()const;

///Get the static color.
const MGColor& staticColor()const{ return m_colorStatic; };
MGColor& staticColor(){ return m_colorStatic; };

///Set light mode. mode=-1:undefined, =0:disabled, =1:enabled.
void setLightMode(int mode){ m_lightMode = mode; };
int getLightMode(){ return m_lightMode; };

//////////////////// draw functions//////////////////////////////////
//以下のdrawXXXはすべてひとつまたは一つ以上のmgVBOLeaf(またはmgVBO)をelementとして作成する。

///gelのmgVBOPointerを作成しmgVBOElementとして追加する。
///gel must be valid when draw event happens since gel is referenced at that time.
void drawGel(const MGAttribedGel& gel);

///gelのmgVBOLeafPointerを作成しmgVBOElementとして追加する。
void drawVBOLeaf(
	const mgVBOLeaf& leaf,///<leaf to add as mgVBOElement.
	MGCL::DRAW_TARGET target= MGCL::SHADING
	    ///<When target=WIRE, built elements are
		///<stored as wire mode display, else as shading mode display
);

///gelのmgVBOPointerをメンバーから外す
void deleteGel(const MGAttribedGel& gel);

///Draw a string with a color.
void drawString(const char* str, const MGPosition& P, const MGColor* colr = nullptr);
void drawString(const wchar_t* str, const MGPosition& P, const MGColor* colr = nullptr);
void drawString(const CString& str,const MGPosition& P, const MGColor* colr=nullptr);
void drawString(const std::string& str, const MGPosition& P, const MGColor* colr = nullptr);

///Screen coordinate type.
void drawStringByScreen(const char* str, const MGPosition& P, const MGColor* colr = nullptr);
void drawStringByScreen(const wchar_t* str, const MGPosition& P, const MGColor* colr = nullptr);
void drawStringByScreen(const CString& str, const MGPosition& P, const MGColor* colr = nullptr);
void drawStringByScreen(const std::string& str, const MGPosition& P, const MGColor* colr = nullptr);

//Draw an arrow symbol with implementation of OpenGL.
//data[0] is the origin of the arrow, data[1] is the top of the arrow,
//data[2], [3] are two bottoms of arrowhead.
void drawArrow(const MGPosition pos[4]);

/// Draw an object of class MGBox, by wireframe.
void drawBox(const MGBox& box);

/// Draw a control points, dotted lines shall be drawn
/// between point[i-1] and point[i], for i = 1, .., length()-1.
void drawPointSeq(
	const MGBPointSeq& bp,///<Target points.
	bool draw_points=true	///<True if points be drawn.
);
void drawPointSeq(
	const MGSPointSeq& sp,///<Target points.
	bool draw_points=true	///<True if points be drawn.
);

///Draw 3D curve in the topology's star cell world coordinates.

//Draw 3D curve in the topology's star cell world coordinates.
///obj is a boundary of the star cell, and the curves are extracted from the
///boundary of the star cell and drawn.
void drawWire_in_star(const MGLoop& obj);

///Draw curvature variation graph so-called Hige.
void drawCurvaGraph(
	const MGCurve& curve,///<The target curve.
	int density,///<Dinsity of the graph.
	bool use_radius=false,///<Indicates if curvature is used(=false) or curvature radius(=true). 
	double scaleRelative =1.,///<Scale of the graph.
	double lengthBase=-1.///base length of curvature graph. If lengthBase<=0.,
				///the length is obtained from input curve by curvatureLengthDisplay().
);

///Draw curvature variation graph so-called Hige.
void drawCurvaGraph(
	std::vector<const MGCurve*>& curves,///<The target curves.
	int density,///<Dinsity of the graph.
	bool use_radius = false,///<Indicates if curvature is used(=false) or curvature radius(=true). 
	double scaleRelative = 1.,///<Scale of the graph.
	double lengthBase = -1.///base length of curvature graph. If lengthBase<=0.,
	///the length is obtained from input curve by curvatureLengthDisplay().
);

///Draw a point using openGL functions.

///Drawing is done twice for inner and outer.
///size is the outer point size. The inner point size is outer size-2.
///If size<=0., mgVBOElement::getDefaultPointSize() is used.
///The outer point color is the VBO's color, and the inner is white.
void drawPoint(double x,double y,double z, double size=-1.);
void drawPoint(const MGPosition& pos, double size=-1.);
void drawPointInverseColor(double x,double y,double z, double size=-1.);
void drawPointInverseColor(const MGPosition& pos, double size=-1.);

///Draw a point with color.

///When innerSize<=0., only outer point is drawn.
///The outer point color is the VBO's color if not input,
///and the inner is white if not input.
void drawPointWithColor(double x,double y,double z,
	double outerSize, double innerSize,
	const MGColor* colorInner=0, const MGColor* colorOuter=0
);

///Draw a point with color.

///When innerSize<=0., only outer point is drawn.
///The outer point color is the VBO's color if not input,
///and the inner is white if not input.
void drawPointWithColor(const MGPosition& pos,
	double outerSize, double innerSize,
	const MGColor* colorInner=0, const MGColor* colorOuter=0
);

///draw points sequence ipos with 2 colors, inner and outer.

///The outer point size is obtained from the VBO's StaticAttribSize.
///The inner point size is outer size-2.
void drawPoints(
	const MGColor& boundary_color,
	const MGColor& inner_color,
	const std::vector<MGPosition>& ipos,
	double size=-1.
);

///Draw a polyline using openGL functions.
///typeは次のものを可とする(GL_xxxxのxxxxを示す）
///POINTS,LINES,LINE_STRIP,LINE_LOOP,TRIANGLE_FAN,TRIANGLE_STRIP,QUAD_STRIP
template<class... Args>
void drawTypedPolyline(GLenum type, const Args*... Ps) {
	Begin(type);
	for (const MGPosition* points[]{ Ps... }; auto & point : points) {
		const MGPosition& P = *point;
		Vertex3d(P[0], P[1], P[2]);
	}
	End();
};

///Draw a polyline using openGL functions.

/// <summary>
/// Draw open or closed polyline of any number of points:
/// e.g. drawOpenPolyline(p0, p1, p2, ...);
/// </summary>
template<class... Args>
void drawOpenPolyline(const Args*... Ps) {
	drawTypedPolyline(GL_LINE_STRIP, Ps...);
};
template<class... Args>
void drawClosedPolyline(const Args*... Ps) {
	drawTypedPolyline(GL_LINE_LOOP, Ps...);
};

///Draw a polyline using openGL functions.

///When clodes=true, 1st and last points will be connected.
void drawPolyline(const MGBPointSeq& line, bool closed=false);

///Draw a polyline using openGL functions.

///When cloded=true, 1st and last points will be connected.
void drawPolyline(const std::vector<MGPosition>& line, bool closed=false);

///Draw a polyline using openGL functions.

///When cloded=true, 1st and last points are connected.
void drawPolyline(size_t nPoints, const MGPosition line[], bool closed=false);

///Draw a line from start to end.
void drawStraight(const MGPosition& end, const MGPosition& start);

//Draw an object in its parameter space(MGDraw_in_parameter_space).

//This is valid only for Surface, Face, Loop, Edge.
void drawObjInParameterSpace(const MGObject& obj);

///Draw the rectangle of a box.
void drawRectangle(
	const MGBox& box	//Box to draw.
);

/// Renders curvatures mapping that comes into colorful image.

/// A point whose curvature is within [lower, upper], the color varies.
void drawSurfaceCurvature(
	const mgTL2Triangles& tld,
	MGCL::SURFACE_CURVATURE_KIND kind,
	double lower, double upper //minimum and maximum value of the curvatures of the kind.
);

/// Renders curvatures mapping that comes into colorful image.

/// A point whose curvature is within [lower, upper], the color varies.
void drawSurfaceCurvature(
	const std::vector<mgTL2Triangles>& tldvec,///<target triangulated data to draw.
	MGCL::SURFACE_CURVATURE_KIND kind,///<Curvature kind, see MGCL.h.
	double lower,///<Minimum value of the curvatures.
	double upper ///< Maximum value of the curvatures of the kind.
);

///OpenGL shading display of a tesselated data tris.
void drawShade(
	const mgTL2Triangles& tris,///<target triangulated data to draw.
	MGCL::DRAW_TARGET target= MGCL::SHADING,///<When target=WIRE, built elements are
			///stored as wire mode display, else as shading mode display
	GLenum polygonMode=GL_FILL//Polygon mode to draw, GLPOINT, GL_LINE, or GLFILL.
);

///OpenGL shading display of a tesselated data tris.
void drawShade(
	const std::vector<mgTL2Triangles>& trisVector,///<target triangulated data to draw.
	MGCL::DRAW_TARGET target= MGCL::SHADING,///<When target=WIRE, built elements are
			///stored as wire mode display, else as shading mode display
	GLenum polygonMode=GL_FILL//Polygon mode to draw, GLPOINT, GL_LINE, or GLFILL.
);

/// MGStlオブジェクトを描画する
void drawSTL(
	const MGStl& stl, ///< 描画するMGStlオブジェクト
	MGCL::DRAW_TARGET target= MGCL::SHADING,///<When target=WIRE, built elements are
			///<stored as wire mode display, else as shading mode display
	GLenum polygonMode=GL_FILL///<Polygon mode to draw, GLPOINT, GL_LINE, or GLFILL.
);

protected:

///Build VBO hierarchy in vbos;

///Let n=vbos.size(), then 
///vbos[i] includes vbos[i+1] as mgVBOPointer for i=0,...,n-2.
///vbos[0]=&parent, and vbos[n-1] = this mgVBO.
///Function's return value is true if found, false if not.
bool buildVBOHierarchy(
	mgVBO& parent,///<Parent mgVBO that may include this vbo as its child
		///< or the descendant.
	std::vector<mgVBO*>& vbos///<VBO hierarchy is output.
);

///Exec gel's draw action.

///mgVBOの既定の処理は次の通り：
///(1) initializeVBO
///(2) m_gelの描画データ作成のみをおこなう。
///すでに作成済みであっても強制的に再作成を行う。
///m_gel=0のときはなにもしない。
virtual void make_display_list(MGCL::VIEWMODE vmode = MGCL::DONTCARE);

///Exec color data action.
void execStaticColorAttrib(){mgGLSL::execStaticColorAttrib(m_colorStatic);};

///Exec static attributes.

/// execStaticGLAttribは固有のuniform値をセットするために
/// Overrideする場合があるので、virtualはつけておいてください。
virtual void execStaticGLAttrib();

/// Get anchor position for mgCoordinateTypeSwitcher.

/// mgVBO that has anchor position must return an adequate position pointer.
virtual const MGPosition* getAnchorPosition(){return nullptr;};

///Set display list name.
void setDlName(unsigned name){m_dname=name;};

///Set gel pointer.
void setGel(const MGAttribedGel* gel);

///Set Elements target.
void setElementTarget(MGCL::DRAW_TARGET target);

};

///Utility class to invoke mgVBO::setLightMode().

///mgLightModeSwitcher saves the current lightmode and invoke input mode.
///When mgLightModeSwitcher is destructed, the saved original mode is restored.
class MG_DLL_DECLR mgLightModeSwitcher{
public:
	mgLightModeSwitcher(mgVBO& vbo, mgGLSL::ShadeMode mode):m_vbo(vbo){
		m_orgMode = m_vbo.getLightMode();
		m_vbo.setLightMode(mode);
	};

	~mgLightModeSwitcher(){
		m_vbo.setLightMode(m_orgMode);
	};

private:
	mgVBO& m_vbo;
	int m_orgMode;
};

#else //_CONSOLE

class mgVBOLeaf;
#include "mg/Position.h"
#include "mgGL/glslprogram.h"

class MG_DLL_DECLR mgVBOElement {
public:
	static const MGDrawParam& getDrawParam() {
		static const MGDrawParam para;
		return para;
	};
};
class MG_DLL_DECLR mgVBO {
public:
	void initializeVBO(MGCL::VIEWMODE viewMode = MGCL::DONTCARE) {
		std::cout << "mgVBO::initializeVBO, viewMode="<< viewMode << std::endl;
	}
	bool is_made(MGCL::VIEWMODE viewMode = MGCL::DONTCARE) { return true; };
	void setDirty(bool is_dirty = true, MGCL::DRAW_TARGET target = MGCL::WIRE_AND_SHADING)const{
		std::cout << "mgVBO::setDirty, is_dirty="<< is_dirty
			<<", target="<<target<< std::endl;
	}
	void setStaticAttribColor(const MGColor& color) { ; };
	void setStaticAttribColor(const float color[4]) { ; };
	void LineWidth(float size) { ; };
	void Begin(GLenum type, MGCL::DRAW_TARGET target = MGCL::WIRE) {
		std::cout << "mgVBO::Begin" << std::endl;
	}
	mgVBOLeaf* End(GLenum polygonMode = GL_FILL) {
		std::cout << "mgVBO::End" << std::endl; return nullptr;
	}
	void Vertex3d(double x, double y, double z) {
		std::cout << "mgVBO::Vertex3d(" << x << ", " << y << ", " << z << ")" << std::endl;
	}
	virtual void Vertex3dv(const double v[3]) {
		std::cout << "mgVBO::Vertex3dv("<<v[0]<<", "<<v[1]<<", "<<v[2]<<")"<<std::endl;
	}
	void draw(MGCL::VIEWMODE viewMode = MGCL::DONTCARE) {
		std::cout << "mgVBO::draw, viewMode=" << viewMode<< std::endl;
	}
	void drawGel(const MGAttribedGel& gel) {
		std::cout << "mgVBO::drawGel, gel=" << &gel << std::endl;
	}
	void drawStraight(const MGPosition& end, const MGPosition& start) {
		std::cout << "mgVBO::drawStraight,end=" << end << ", start=" << start << std::endl;
	}
	void drawPoint(double x, double y, double z, double size = -1.) {
		std::cout << "mgVBO::drawPoint(" << x << ", " << y << ", " << z << ")"
			<<",size="<<size<<std::endl;
	}
	void drawShade(
		const mgTL2Triangles& tris,///<target triangulated data to draw.
		MGCL::DRAW_TARGET target = MGCL::SHADING,///<When target=WIRE, built elements are
				///stored as wire mode display, else as shading mode display
		GLenum polygonMode = GL_FILL//Polygon mode to draw, GLPOINT, GL_LINE, or GLFILL.
	){
		std::cout << "mgVBO::drawShade, tris=" << &tris << ",target="<<target
			<< ", polygonMode=" << polygonMode << std::endl;
	}

};
class MG_DLL_DECLR mgLightModeSwitcher {
public:
	mgLightModeSwitcher(mgVBO& vbo, mgGLSL::ShadeMode mode) { ; };
};

#endif //_CONSOLE


/** @} */ // end of DisplayHandling group
