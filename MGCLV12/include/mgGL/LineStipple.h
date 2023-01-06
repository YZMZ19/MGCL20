/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#pragma once

#include "mgGL/GLAttrib.h"

class MGOfstream;
class MGIfstream;

//
//Define MGLineStipple Class.

/** @addtogroup GLAttrib
 *  @{
 */

///Define line font enum.
enum class LineFont {
	UndefinedFont = -1,
	Solid = 1,		///<factor=2, pattern=0xFFFF
	Dashed = 2,		///<facotr=2, pattern=0x3333
	Phantom = 3,		///<facotr=2, pattern=0x5757
	CenterLine = 4,	///<facotr=2, pattern=0x5f5f
	Dotted = 5		///<facotr=2, pattern=0x1111
};

///MGLineStipple defines line stipple patters.

///The pattern is defined as a binary data of unsigned short.
class MG_DLL_DECLR MGLineStipple:public MGGLAttrib{

public:

MGLineStipple():MGGLAttrib(static_cast<int>(mgGLMode::UNDEFINED)),m_pattern(0xffff){;};
MGLineStipple(LineFont font);
MGLineStipple(int factor,unsigned short pattern)
:MGGLAttrib(factor){m_pattern=pattern;};

///Assignment
MGLineStipple& operator=(const MGGel& gel2);
MGLineStipple& operator=(const MGLineStipple& gel2);

///comparison
bool operator<(const MGLineStipple& gel2)const;
bool operator<(const MGGel& gel2)const;

////////////Member Function////////////

///Generate a newed clone object.
MGLineStipple* clone()const;

///Invoke appropriate OpenGL fucntion to the drawing environment.
void exec()const;

///Invoke appropriate OpenGL fucntion to this attributefor the vbo.
void exec(mgVBO& vbo)const;

///Set factor and pattern.
void set(int factor,unsigned short pattern){data()=factor;m_pattern=pattern;};

///Get factor.
int get_factor()const{return data();};

///Get pattern.
unsigned short get_pattern()const{return m_pattern;};

///Get the font number
LineFont get_font_number()const;

///draw GLAttribute process.
void drawAttrib(
	mgVBO& vbo,///<The target graphic object.
	bool no_color=false	///<if true, color attribute will be neglected.
)const{exec(vbo);};

///Test if this is highlight attrib or not.
bool is_highlight_attrib()const{return true;};

///render GLAttribute process.
void render(mgVBO& vbo)const{exec(vbo);};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
void set_draw_attrib_mask(unsigned int& mask)const{set_Amask(mask,LINE_BIT);};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
void reset_draw_attrib_mask(unsigned int& mask)const{reset_Amask(mask,LINE_BIT);};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
void set_render_attrib_mask(unsigned int& mask)const{set_Amask(mask,LINE_BIT);};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
void reset_render_attrib_mask(unsigned int& mask)const{reset_Amask(mask,LINE_BIT);};

/// Return This object's typeID
long identify_type() const{return MGLINE_STIPPLE_TID;};

///Get the name of the class.
std::string whoami()const{return "LineStipple";};

///Read all member data.
void ReadMembers(MGIfstream& buf);
///Write all member data
void WriteMembers(MGOfstream& buf)const;

/// Output function.
std::ostream& toString(std::ostream&) const;

private:
	unsigned short m_pattern;///<line stipple pattern.

};

/** @} */ // end of GLAttrib group

