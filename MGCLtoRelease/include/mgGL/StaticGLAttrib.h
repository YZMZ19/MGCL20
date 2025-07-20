/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#ifndef _MGStaticGLAttrib_HH_
#define _MGStaticGLAttrib_HH_

class MGColor;
#include "mg/MGCL.h"

/** @addtogroup GLAttrib
 *  @{
 */

//
///Define mgStaticGLAttrib Class.

///mgStaticGLAttrib defines MGColor and line width data of OpenGL.
class mgStaticGLAttrib{

public:
	mgStaticGLAttrib();
	mgStaticGLAttrib(const MGColor& color, float lineWidth);

	///Set color.
	void setColor(const MGColor& color);
	void setColor(const float color[4]);

	///Set line width.
	void setLineWidth(float lineWidth){m_lineWidth=lineWidth;};

	///Line stipple属性をセットする。

	///When factor=0 is input, line pattern is disabled. 実線となる
	///When factor<0, the stipple attribute is undefined. This means the attribute
	///is defined by the environment.
	///When factor<=0, pattern is unnecessary.
	void setLineStipple(short int factor, GLushort pattern);
	
	///Set light mode. mode=-1:undefined, =0:disabled, =1:enabled.
	void setLightMode(int mode){m_lightMode=mode;};

	///Get color.
	const float* color()const{return m_color;};
	void getColor(MGColor& color);

	///Get line width.
	float getLineWidth()const{return m_lineWidth;};

	///Get line stipple.
	void getLineStipple(short int& factor, GLushort& pattern)const;

	///Get light mode.
	int getLightMode()const{return m_lightMode;};

private:

	float m_color[4];	///<color data (r,g, b, a)
	float m_lineWidth;	///Line width.
	short int m_stippleFactor;///Line stipple factor. If m_stippleFactor=0, line Stipple is disabled,
	GLushort m_LineStipplePattern;///m_LineStipplePatternindicates the pattern.

	///light modeはm_elementsShadeに対してのみ有効。m_elementsに対しては常にlightはオフ
	int m_lightMode;/// =0:Light is disabled, >0:Light is enabled.

};

/** @} */ // end of GLAttrib group
#endif //#ifndef _MGStaticGLAttrib_HH_
