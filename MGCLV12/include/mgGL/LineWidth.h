/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#pragma once

#include "mgGL/GLAttrib.h"

class MGOfstream;
class MGIfstream;

//
//Define MGLineWidth Class.

/** @addtogroup GLAttrib
 *  @{
 */

///MGLineWidth defines line width of a curve.

///Line width is defined by a float data.
class MG_DLL_DECLR MGLineWidth:public MGGLAttrib{

public:

MGLineWidth(mgGLMode m= mgGLMode::UNDEFINED):MGGLAttrib(static_cast<int>(m)){;};
MGLineWidth(float width)
	:MGGLAttrib(static_cast<int>(mgGLMode::ENABLED)),m_line_width(width){;};

///Assignment
MGLineWidth& operator=(const MGGel& gel2);
MGLineWidth& operator=(const MGLineWidth& gel2);

///comparison
bool operator<(const MGLineWidth& gel2)const;
bool operator<(const MGGel& gel2)const;

////////////Member Function////////////

///Generate a newed clone object.
MGLineWidth* clone()const;

///Invoke appropriate OpenGL fucntion to the drawing environment.
void exec()const;

///Invoke appropriate OpenGL fucntion to this attribute.
void exec(mgVBO& vbo)const;

///Set width.
void set_width(float width);

///Get width.
float get_width()const{return m_line_width;};

///draw GLAttribute process.
void drawAttrib(
	mgVBO& vbo,///<The target graphic object.
	bool no_color=false	///<if true, color attribute will be neglected.
)const{exec(vbo);};

///get maximum line width
float get_maximum_width()const;

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
long identify_type() const{return MGLINE_WIDTH_TID;};

///Get the name of the class.
std::string whoami()const{return "LineWidth";};

///Read all member data.
void ReadMembers(MGIfstream& buf);

///Write all member data
void WriteMembers(MGOfstream& buf)const;

/// Output function.
std::ostream& toString(std::ostream&) const;

private:

	float m_line_width;	///<line width.

};

/** @} */ // end of GLAttrib group
