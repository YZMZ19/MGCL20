/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#pragma once

#include <iosfwd>
#include "mg/MGCL.h"
#include "mg/Attrib.h"

//
//Define MGGLAttrib Class.

class MGOfstream;
class MGIfstream;
class MGColor;
class MGLight;
class MGRenderAttr;
class MGLineStipple;
class MGLights;
class MGLineWidth;

/** @addtogroup GLAttrib
 *  @{
 */

enum class mgGLMode{
	UNDEFINED = -3,
	DISABLED = -2,
	ENABLED = 1,
};

///MGGLAttrib is an abstract class which defines the enum of undefined or disabled.

///Subclass of MGGLAttrib can use m_flag as a part of its own class attributes.
///In this case, -3 and -2 must be avoided. -3 and -2 do not appear in OpenGL attributes.
class MG_DLL_DECLR MGGLAttrib:public MGAttrib{

public:

///Defines MGGLAttrib bit positions, which are equal to OpenGL's ones.
enum ATTRIB_MASK {
	CURRENT_BIT =		0x00000001,//GL_CURRENT_BIT,
	POINT_BIT=			0x00000002,//GL_POINT_BIT,
	LINE_BIT=			0x00000004,//GL_LINE_BIT,
	POLYGON_BIT=		0x00000008,//GL_POLYGON_BIT,
	POLYGON_STIPPLE_BIT=0x00000010,//GL_POLYGON_STIPPLE_BIT,
	PIXEL_MODE_BIT=		0x00000020,//GL_PIXEL_MODE_BIT,
	LIGHTING_BIT=		0x00000040,//GL_LIGHTING_BIT,
	FOG_BIT=			0x00000080,//GL_FOG_BIT,
	DEPTH_BUFFER_BIT=	0x00000100,//GL_DEPTH_BUFFER_BIT,
	ACCUM_BUFFER_BIT=	0x00000200,//GL_ACCUM_BUFFER_BIT,
	STENCIL_BUFFER_BIT= 0x00000400,//GL_STENCIL_BUFFER_BIT,
	VIEWPORT_BIT=		0x00000800,//GL_VIEWPORT_BIT,
	TRANSFORM_BIT=		0x00001000,//GL_TRANSFORM_BIT,
	ENABLE_BIT=			0x00002000,//GL_ENABLE_BIT,
	COLOR_BUFFER_BIT=	0x00004000,//GL_COLOR_BUFFER_BIT,
	HINT_BIT=			0x00008000,//GL_HINT_BIT,
	EVAL_BIT=			0x00010000,//GL_EVAL_BIT,
	LIST_BIT=			0x00020000,//GL_LIST_BIT,
	TEXTURE_BIT=		0x00040000,//GL_TEXTURE_BIT,
	SCISSOR_BIT=		0x00080000 //GL_SCISSOR_BIT
};

MGGLAttrib(int flag= static_cast<int>(mgGLMode::UNDEFINED)):m_flag(flag){;};

///Assignment operator.
virtual MGGLAttrib& operator=(const MGGLAttrib& gel2){set_glattrib(gel2);return *this;};

////////////Member Function////////////

///Generate a newed clone object.
virtual MGGLAttrib* clone()const=0;

///Test if this is defined data or not.
bool undefined()const{return m_flag== static_cast<int>(mgGLMode::UNDEFINED);};
bool defined()const{return !undefined();};

///Test if this is enabled or not.
bool disabled()const{return m_flag== static_cast<int>(mgGLMode::DISABLED);};
bool enabled()const {return !undefined() && !disabled();}

///Set as undefined.
void set_undefined(){setFlag(mgGLMode::UNDEFINED);};

///Set as undefined.
void set_disabled(){ setFlag(mgGLMode::DISABLED);};

///retrieve the data.
int data()const{return m_flag;};
int& data(){return m_flag;};

///draw GLAttribute process.
virtual void drawAttrib(
	mgVBO& vbo,///<The target graphic object.
	bool no_color=false	///<if true, color attribute will be neglected.
)const=0;

///Test if this is highlight attrib or not.
virtual bool is_highlight_attrib()const{return false;};

///render GLAttribute process.
virtual void render(mgVBO& vbo)const=0;

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
virtual void set_draw_attrib_mask(unsigned int& mask)const=0;

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
virtual void reset_draw_attrib_mask(unsigned int& mask)const=0;

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
virtual void set_render_attrib_mask(unsigned int& mask)const=0;

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
virtual void reset_render_attrib_mask(unsigned int& mask)const=0;

/// Return This object's typeID
virtual long identify_type() const{return MGGLATTRIBUTE_TID;};

///Read all member data.
virtual void ReadMembers(MGIfstream& buf);
///Write all member data
virtual void WriteMembers(MGOfstream& buf)const;

/// Output virtual function.
virtual std::ostream& toString(std::ostream&) const;

///Compare if this and at2 are the same leaf MGGLAttrib class.
bool same_type(const MGGLAttrib& at2)const;

protected:

int m_flag;	///< =-3:undefined, will be inheritted.
			///< =-2:disabled.
			///< =other:each subclass's data(enabled).

///set m_flag.
void setFlag(mgGLMode flag) { m_flag = static_cast<int>(flag); };

///Assignment
MGGLAttrib& set_glattrib(const MGGLAttrib& gel2);

};

///Set the bit of mask.
void set_Amask(unsigned int& mask, MGGLAttrib::ATTRIB_MASK bit);

///Reset the bit of mask.
void reset_Amask(unsigned int& mask, MGGLAttrib::ATTRIB_MASK bit);

/** @} */ // end of GLAttrib group

