/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#pragma once

class MGObject;
class MGGLAttrib;
class MGColor;
class MGAppearance;
class MGName;
class mgVBO;

#include "mg/MGCL.h"
#include "mg/Gel.h"

#ifndef _CONSOLE

//
//Define MGAttribedGel Class.


/** @addtogroup GelRelated
 *  @{
 */

///Is an abstract class which provides interfaces of MGGel that has MGAppearance.

///MGAppearance has MGAttrib's that decorate MGGel. MGAttribedGel is MGGel
///that has MGAppearance.
class MG_DLL_DECLR MGAttribedGel:virtual public MGGel{
protected:
	mutable std::unique_ptr<mgVBO> m_VBO;///<display name if m_VBO!=0;

public:

////////Special member functions/////////
MGAttribedGel()=default;///void constructor.
virtual ~MGAttribedGel();
MGAttribedGel(const MGAttribedGel& rhs);//copy constructor.
MGAttribedGel& operator=(const MGAttribedGel& rhs);//copy assignment.
MGAttribedGel(MGAttribedGel&& rhs);//move constructor.
MGAttribedGel& operator=(MGAttribedGel&& rhs);//move assignment.

///Get the MGAppearance pointer in this MGAttribedGel.
///If not defined, null will be returned.
virtual MGAppearance* appearance()=0;
virtual const MGAppearance* appearance()const=0;

///copy the appearance of gel2 to this.
void copy_appearance(const MGAttribedGel& gel2);

///Obtain display list name.
mgVBO* dlist_name()const;

///Judge if the display list for vmode is made or not.
virtual bool displayList_is_made(MGCL::VIEWMODE vmode)const;

///Delete VBO and remove dlist name from MGDNameControl.
void deleteVBO()const;

///Process of draw or render attributes.
virtual void drawAttrib(
	mgVBO& vbo,///<The target graphic object.
	bool no_color=false	///if true, color attribute will be neglected.
)const;
virtual void render_attribute()const;

///Make sure that this MGAttribedGel has appearance,
/// and get the MGAppearance pointer.
virtual MGAppearance* ensure_appearance()=0;

///Obtain attribute mask for glPushAttrib().
virtual int get_draw_attrib_mask()const;
virtual int get_render_attrib_mask()const;

///Get the number of elements of m_VBO.
int getVBOElementsNumber()const;

///Get the number of shading elements of m_VBO.
int getVBOShaderElementsNumber()const;

///Make a display list of this gel.
virtual void make_display_list(MGCL::VIEWMODE vmode=MGCL::DONTCARE)const{;};

///Test if this group is attributed  as no display.
///true if attributed as no display.
virtual bool no_display()const;

///Remove the MGAppearance of this MGAttribedGel.
virtual std::unique_ptr<MGAppearance> remove_appearance()=0;

///Removed the attribute of specified type.
std::unique_ptr<MGGLAttrib> remove_GLattrib(long tid);

///Set the attribute in this list. attr must be a newed object, and the
///ownership will be transfered to this MGAppearance.
virtual void set_GLattrib(MGGLAttrib* attr);

///Set this group as display or no display group.
virtual void set_display();
virtual void set_no_display();

///Test if this is visible or not.
bool visible()const{return !no_display();};

///set the copy of appearance appr2 to this MGAttribedgel.
virtual void set_appearance(const MGAppearance& appr2)=0;//Copy.
virtual void set_appearance(MGAppearance* appr2)=0;//Ownership transfer.

///Set name to this MGAttribedGel.
///If this had a name already, the old name is replace with newName;
void set_name(const MGName& newName);

///Get the name of this MGAttribedGel.
///If this had no name, null will be returned.
///If this did have a name, the pointer is returned.
const MGName* get_name()const;

///Set color to this MGAttribedGel.
///If this had a color already, the old color is replace with newColor;
void set_color(const MGColor& newColor);

///Get the color of this MGAttribedGel.
///If this had no color, null will be returned.
///If this did have a color, the pointer is returned.
const MGColor* get_color()const;

///Set DlistName.
///vbo must be newed one, and the ownership is transfered to this MGAttribedGel.
void setDlistName(mgVBO* vbo=0)const;

///Set dirty flag(s) of this VBO(m_VBO).
void setDirty(
	bool is_dirty=true	///< =true if dirty and need to remake.
)const;

friend class MGGroup;
friend class MGObject;
friend class mgVBO;
friend class MGDNameControl;
};

#else //_CONSOLE

#include "mgGL/VBO.h"
class MG_DLL_DECLR MGAttribedGel :virtual public MGGel {
public:
	mgVBO* dlist_name()const {
		static mgVBO vbo;
		return &vbo;
	}
	virtual MGAppearance* appearance() = 0;
	virtual const MGAppearance* appearance()const = 0;
	void copy_appearance(const MGAttribedGel& gel2) { ; };
	const MGName* get_name()const { return nullptr; };
	virtual std::unique_ptr<MGAppearance> remove_appearance() = 0;
	virtual void set_appearance(const MGAppearance& appr2) = 0;//Copy.
	virtual void set_appearance(MGAppearance* appr2) = 0;//Ownership transfer.
	void setDirty(bool is_dirty = true)const { ; };
	virtual void set_display() { ; };
	virtual void set_no_display() {	;};
	virtual bool no_display()const { return true; };
	void set_GLattrib(MGGLAttrib* attr) { ; };
	virtual bool displayList_is_made(MGCL::VIEWMODE vmode)const { return true; };
};

#endif //_CONSOLE

/** @} */ // end of GelRelated group
