/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGNAME_HH_
#define _MGNAME_HH_

#include <iosfwd>
#include <string>
#include "mgGL/GLAttrib.h"
#include "mgGL/Color.h"
#include "mgGL/Appearance.h"

class MGOfstream;
class MGIfstream;

/** @addtogroup GLAttrib
 *  @{
 */

///Defines MGAttribedGel's Name data.

///MGName provides naming facility for MGAttribedGel, that is,
///MGGroup or MGObject.
///Add a name data of MGName as an MGGLAttrib data to MGAttribedGel, 
///using set_name, then search of it can be performed by giving the name data.
///MGName inherits MGGLAttrib, it does not have no effect to OpenGL attributes.
///MGName uses only MGAppearqance facilities.
///TID of MGName:	MGNAME_TID=			0x02010300L;
class MG_DLL_DECLR MGName : public MGGLAttrib{

public:

MGName();
MGName(const std::string& name);

///Test if this is the same name as name2.
bool operator ==(const MGName& name2)const{return m_name==name2.m_name;};

///Generate a newed clone object.
MGName* clone()const;

///Set name.
void setName(std::string& name);
void setName(char* name);

///Get name.
std::string& getName(){return m_name;};
const std::string& getName()const{return m_name;};

///render GLAttribute process.
void exec(mgVBO& vbo)const{;};

///draw GLAttribute process.
void drawAttrib(
	mgVBO& vbo,///<The target graphic object.
	bool no_color=false	///<if true, color attribute will be neglected.
)const{;};

///render GLAttribute process.
void render(mgVBO& vbo)const{;};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
void set_draw_attrib_mask(unsigned int& mask)const{mask=0;};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
void reset_draw_attrib_mask(unsigned int& mask)const{;};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
void set_render_attrib_mask(unsigned int& mask)const{mask=0;};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
void reset_render_attrib_mask(unsigned int& mask)const{;};

/// Return This object's typeID
long identify_type() const{return MGNAME_TID;};

///Get the name of the class.
std::string whoami()const{return "Name";};

///Read all member data.
void ReadMembers(MGIfstream& buf);
///Write all member data
void WriteMembers(MGOfstream& buf)const;

/// Output function.
std::ostream& toString(std::ostream&) const;

private:

    std::string m_name;//Name data will be stored.

};

/** @} */ // end of GLAttrib group
#endif // _MGNAME_HH_
