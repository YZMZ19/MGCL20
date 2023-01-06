/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#pragma once

#include "mg/MGCL.h"
#include "mg/Position.h"
#include "mgGL/VBObyAnchorPt.h"

class FTFont;
class MGColor;

/** @addtogroup DisplayHandling
 *  @{
 */

///Defines String writer class using mgVBO(OpenGL infrastructure).

///MGStringWriter makes an mgVBO(display list) to write the input string
///at the input position pos.
///Before use of MGStringWriter, Init() must be invoked once.
class MG_DLL_DECLR MGStringWriter{
public:


	///const chat* and pos as world coordinate version.
	static VBObyAnchorPt* Draw(
		const char *str,	///<String to display
		const MGPosition& pos,///<position for the string to display at.
		const MGColor* color=NULL///<Color of the string.
	);

	///const wchar_t* and pos as world coordinate version.
	static VBObyAnchorPt* Draw(
		const wchar_t *str,	///<String to display
		const MGPosition& pos,///<position for the string to display at.
		const MGColor* color=NULL///<Color of the string.
	);

	///const char* and pos as screen coordinate version.
	static VBObyAnchorPt* DrawByScreen(
		const char *str,	///<String to display
		const MGPosition& pos,///<position for the string to display at.
		const MGColor* color=NULL///<Color of the string.
	);

	///const wchar_t* and pos as screen coordinate version.
	static VBObyAnchorPt* DrawByScreen(
		const wchar_t *str,	///<String to display
		const MGPosition& pos,///<position for the string to display at.
		const MGColor* color=NULL///<Color of the string.
	);

	///Set font data.
	///The default font is "C:\\Windows\\Fonts\\arial.ttf".
	///If this is not the case, setFont must be invoked.
	static void setFont(
		const char* fontPath, ///< font file path.
		unsigned int faceSize = 12,///< the face size in points(1/72 inch).
		unsigned int resolution =96///<the resolution of the target device.
	);

	///Initialize MGStringWriter. This must be invoked once before use.
	static void Init();

private:
	static MGStringWriter* getInstance();

	MGStringWriter(void);
	virtual ~MGStringWriter(void);

	template<typename T>
	VBObyAnchorPt* DrawString(
		mgGLSL::CoordinateType type,
		const T* str,//target string to draw
		const MGPosition& pos,//Start position to draw in type coordinates
		const MGColor* color)
	{
		VBObyAnchorPt* pVBO = new VBObyAnchorPt(type);
		pVBO->setStaticAttribAnchorPoint(pos);
		if(color)
			pVBO->setStaticAttribColor(*color);

		m_pFont->Render(*pVBO, str);
		pVBO->setDirty(false);
		return pVBO;
	}

	FTFont* m_pFont;
};

/** @} */ // end of DisplayHandling group
