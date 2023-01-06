/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#include "StdAfx.h"
#include "mg/DNameControl.h"
#include "mg/AttribedGel.h"
#include "mgGL/Appearance.h"
#include "mgGL/GLAttrib.h"
#include "mgGL/Name.h"
#include "mgGL/VBO.h"

//
//Define MGAttribedGel Class.
//MGAttribedGel is an abstract class which provides function interfaces of
//MGGroup and MGObject that have MGAppearance. 
//

///copy constructor.
MGAttribedGel::MGAttribedGel(const MGAttribedGel& gel2)
:MGGel(gel2){
}

//Move constructor.
MGAttribedGel::MGAttribedGel(MGAttribedGel&& gel2)
:MGGel(std::move(gel2)), m_VBO(std::move(gel2.m_VBO)){
	mgVBO* vbo=m_VBO.get();
	if(vbo)
		vbo->setGel(this);
}

MGAttribedGel::~MGAttribedGel() { ; }

//Copy assignment.
MGAttribedGel& MGAttribedGel::operator=(const MGAttribedGel& gel2){
	if(this!=&gel2){
		MGGel::operator=(gel2);
		deleteVBO();
	}
	return *this;
}

//Move assignment.
MGAttribedGel& MGAttribedGel::operator=(MGAttribedGel&& gel2){
	MGGel::operator=(std::move(gel2));
	m_VBO=std::move(gel2.m_VBO);
	mgVBO* vbo=m_VBO.get();
	if(vbo)
		vbo->setGel(this);
	return *this;
}

//copy the appearance of gel2 to this.
void MGAttribedGel::copy_appearance(const MGAttribedGel& gel2){
	remove_appearance();
	const MGAppearance* appr2=gel2.appearance();
	if(appr2){
		set_appearance(*appr2);
	}
}

///Obtain display list name.
mgVBO* MGAttribedGel::dlist_name()const{
	mgVBO* vbo=m_VBO.get();
	if(!vbo){
		//When this vbo is null, make vbo.
		MGDNameControl& dnc=getDNameControlInstance();
		vbo=dnc.insertDlistMap(this);
		m_VBO.reset(vbo);
	}
	return vbo;
}

void MGAttribedGel::deleteVBO()const{
	if(m_VBO.get()){
		MGDNameControl& dnc=getDNameControlInstance();
		mgVBO* vbo=dnc.deleteDlistMap(this);
		m_VBO.reset(0);
		delete vbo;
	}
}

///Process of draw or render attributes.
void MGAttribedGel::drawAttrib(
	mgVBO& vbo,///<The target graphic object.
	bool no_color///<if true, color attribute will be neglected.
)const{
	const MGAppearance* app=appearance();
	if(app)
		app->drawAttrib(vbo,no_color);
}
void MGAttribedGel::render_attribute()const{
	mgVBO* vbo=dlist_name();
	const MGAppearance* app=appearance();
	if(app)
		app->render(*vbo);
}

//Obtain attribute mask for glPushAttrib().
int MGAttribedGel::get_draw_attrib_mask()const{
	int mask=0;
	const MGAppearance* app=appearance();
	if(app)
		mask=app->get_draw_attrib_mask();
	return mask;
}
int MGAttribedGel::get_render_attrib_mask()const{
	int mask=0;
	const MGAppearance* app=appearance();
	if(app)
		mask=app->get_render_attrib_mask();
	return mask;
}

//Test if this group is attributed  as no display.
//true if attributed as no display.
bool MGAttribedGel::no_display()const{
	const MGAppearance* app=appearance();
	if(app)
		return app->no_display();
	return false;
}

//Removed the attribute of specified type.
std::unique_ptr<MGGLAttrib> MGAttribedGel::remove_GLattrib(long tid){
	std::unique_ptr<MGGLAttrib> atr;
	MGAppearance* app=appearance();
	if(app){
		atr=app->release_attrib(tid);
		if(app->can_be_removed())
			remove_appearance();
	}
	return atr;
}

//Set the attribute in this list. attr must be a newed object, and the
//ownership will be transfered to this MGAppearance.
void MGAttribedGel::set_GLattrib(MGGLAttrib* attr){
	if(!attr)
		return;
	MGAppearance* app=ensure_appearance();
	app->set_attrib(attr);
}

//Set this group as display or no display group.
void MGAttribedGel::set_display(){
	MGAppearance* app=appearance();
	if(!app)
		return;
	app->set_display();
}
void MGAttribedGel::set_no_display(){
	MGAppearance* app=ensure_appearance();
	app->set_no_display();
}

///Get the name of this MGAttribedGel.
///If this had no name, null will be returned.
///If this did have a name, the pointer is returned.
const MGName* MGAttribedGel::get_name()const{
	const MGAppearance* apr=appearance();
	if(!apr)
		return 0;
	MGAppearance::const_iterator i=apr->search_by_id(MGNAME_TID);
	if(i==apr->end())
		return 0;

	const MGName* nm=static_cast<const MGName*>(i->get());
	return nm;
}

///Set name to this MGAttribedGel.
///If this had a name already, the old name is replace with newName;
void MGAttribedGel::set_name(const MGName& newName){
	MGAppearance* apr=ensure_appearance();
	MGAppearance::iterator i=apr->search_by_id(MGNAME_TID);
	if(i==apr->end()){
		apr->set_attrib(new MGName(newName));
	}else{
		MGName& oldName=static_cast<MGName&>(**i);
		oldName=newName;
	}
}

///Set color to this MGAttribedGel.
///If this had a color already, the old color is replace with newColor;
void MGAttribedGel::set_color(const MGColor& newColor){
	MGAppearance* apr=ensure_appearance();
	MGAppearance::iterator i=apr->search_by_id(MGCOLOR_TID);
	if(i==apr->end()){
		apr->set_attrib(new MGColor(newColor));
	}else{
		MGColor& oldColor=static_cast<MGColor&>(**i);
		oldColor=newColor;
	}
}

///Get the color of this MGAttribedGel.
///If this had no color, null will be returned.
///If this did have a color, the pointer is returned.
const MGColor* MGAttribedGel::get_color()const{
	const MGAppearance* apr=appearance();
	if(!apr)
		return 0;
	MGAppearance::const_iterator i=apr->search_by_id(MGCOLOR_TID);
	if(i==apr->end())
		return 0;

	const MGColor* color=static_cast<MGColor*>(i->get());
	return color;
}

///Set DlistName.
///vbo must be newed one, and the ownership is transfered to this MGAttribedGel.
void MGAttribedGel::setDlistName(mgVBO* vbo)const{
	m_VBO.reset(vbo);
}

///Set dirty flag(s) of this VBO(m_VBO).
void MGAttribedGel::setDirty(bool is_dirty)const{
	if(m_VBO.get()){
		m_VBO->setDirty(is_dirty);
	}
}
