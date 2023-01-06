/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include <bitset>
#include "mg/GelFactory.h"
#include "mg/Ofstream.h"
#include "mg/Ifstream.h"
#include "mgGL/Color.h"
#include "mgGL/LineWidth.h"
#include "mgGL/Appearance.h"
#include "mgGL/LineStipple.h"
#include "mgGL/Name.h"
#include "mgGL/VBO.h"

#include "mgGL/DirectionalLight.h"
#include "mgGL/SpotLight.h"
#include "mgGL/RenderAttr.h"
#include "mgGL/GLSLProgram.h"
#include "mgGL/Lights.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//Set or reset the bit of mask.
void set_Amask(unsigned int& mask, MGGLAttrib::ATTRIB_MASK bit){mask|=bit;}
void reset_Amask(unsigned int& mask, MGGLAttrib::ATTRIB_MASK bit){mask&=~bit;}

//////////////  MGGLAttrib  ////////////////

//Assignment
MGGLAttrib& MGGLAttrib::set_glattrib(const MGGLAttrib& gel2){
	if(this==&gel2)
		return *this;

	MGAttrib::operator=(gel2);
	m_flag=gel2.m_flag;
	return *this;
}

//Read all member data.
//Write all member data
void MGGLAttrib::WriteMembers(MGOfstream& buf)const{
	buf << m_flag;
}
void MGGLAttrib::ReadMembers(MGIfstream& buf){
	buf >> m_flag;
}

std::ostream& MGGLAttrib::toString(
	std::ostream& ostrm
)const{
	ostrm<<this;
	if(undefined()) ostrm<<":UNDEFINED";
	else if(disabled()){
		ostrm<<":DISABLED";
	}
	return ostrm;
}

//Compare if this and at2 are the same leaf MGGLAttrib class.
bool MGGLAttrib::same_type(const MGGLAttrib& at2)const{
	return identify_type()==at2.identify_type();
}

//////////////  MGLineWidth  ////////////////

//assignment
MGLineWidth& MGLineWidth::operator=(const MGLineWidth& gel2){
	if(this==&gel2)
		return *this;

	MGGLAttrib::operator=(gel2);
	m_line_width=gel2.m_line_width;
	return *this;
}
MGLineWidth& MGLineWidth::operator=(const MGGel& gel2){
	const MGLineWidth* gel2_is_this=dynamic_cast<const MGLineWidth*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGLineWidth::operator<(const MGLineWidth& gel2)const{
	if(m_flag==gel2.m_flag)
		return m_line_width<gel2.m_line_width;
	return m_flag<gel2.m_flag;
}
bool MGLineWidth::operator<(const MGGel& gel2)const{
	const MGLineWidth* gel2_is_this=dynamic_cast<const MGLineWidth*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

MGLineWidth* MGLineWidth::clone()const{
	return new MGLineWidth(*this);
}

//get maximum line width
float MGLineWidth::get_maximum_width()const{
	float wd[2] = { 0.,0. };
#ifndef _CONSOLE
	glGetFloatv(GL_LINE_WIDTH_RANGE,wd);
#endif
	return wd[1];
}

// Output function.
std::ostream& MGLineWidth::toString(std::ostream& ostrm) const{
	ostrm<<"LineWidth="; MGGLAttrib::toString(ostrm);
	if(undefined())
		return ostrm;
	ostrm<<"="<<m_line_width;
	return ostrm;
}

void MGLineWidth::set_width(float width){
	m_line_width=width;
	setFlag(mgGLMode::ENABLED);
}

///環境属性として属性をセットする
void MGLineWidth::exec()const{
	if(undefined())
		return;

	mgGLSL::execStaticLineWidth(get_width());//glLineWidth(get_width());
}

///vboに対して属性をセットする
void MGLineWidth::exec(mgVBO& vbo)const{
	if(undefined())
		return;
	vbo.LineWidth(get_width());//glLineWidth(get_width());
}

void MGLineWidth::WriteMembers(MGOfstream& buf)const{
	MGGLAttrib::WriteMembers(buf);
	if(undefined())
		return;
	buf<<m_line_width;
}
void MGLineWidth::ReadMembers(MGIfstream& buf){
	MGGLAttrib::ReadMembers(buf);
	if(undefined())
		return;
	buf>>m_line_width;
}

//////////////////MGLineStipple//////////////////

//assignment
MGLineStipple& MGLineStipple::operator=(const MGLineStipple& gel2){
	if(this==&gel2)
		return *this;

	MGGLAttrib::operator=(gel2);
	m_pattern=gel2.m_pattern;
	return *this;
}
MGLineStipple& MGLineStipple::operator=(const MGGel& gel2){
	const MGLineStipple* gel2_is_this=dynamic_cast<const MGLineStipple*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGLineStipple::operator<(const MGLineStipple& gel2)const{
	if(m_flag==gel2.m_flag)
		return m_pattern<gel2.m_pattern;
	return m_flag<gel2.m_flag;
}
bool MGLineStipple::operator<(const MGGel& gel2)const{
	const MGLineStipple* gel2_is_this=dynamic_cast<const MGLineStipple*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

MGLineStipple::MGLineStipple(LineFont font):MGGLAttrib(2){
	switch(font){
	//case LineFont::Solid: m_pattern=0xffff; break;
	case LineFont::Dashed: m_pattern=0x3333; break;
	case LineFont::Phantom: m_pattern=0x5757; break;
	case LineFont::CenterLine: m_pattern=0x5f5f; break;
	case LineFont::Dotted: m_pattern=0x1111; break;
	default: m_pattern = 0xffff;
	}
}

MGLineStipple* MGLineStipple::clone()const{
	return new MGLineStipple(*this);
}

///環境属性として属性をセットする
void MGLineStipple::exec()const{
	if(undefined())
		return;
	if(disabled()){
		mgGLSL::execStaticLineStipple(0,0);//glDisable(GL_LINE_STIPPLE);
	}else{
		mgGLSL::execStaticLineStipple((short)get_factor(),get_pattern());
	}
}

///vboに対して属性をsetする
void MGLineStipple::exec(mgVBO& vbo)const{
	if(undefined())
		return;
#ifndef _CONSOLE
	vbo.setLineStipple((short)get_factor(),get_pattern());
#endif
}

//Get the font number
LineFont MGLineStipple::get_font_number()const{
	int factr=get_factor();
	if(factr!=2)
		return LineFont::UndefinedFont;

	unsigned short font=get_pattern();
	switch(font){
	case 0xFFFF: return LineFont::Solid;
	case 0x3333: return LineFont::Dashed;
	case 0x5757: return LineFont::Phantom;
	case 0x5f5f: return LineFont::CenterLine;
	case 0x1111: return LineFont::Dotted;
	}

	return LineFont::UndefinedFont;
}

std::ostream& MGLineStipple::toString(std::ostream& ostrm) const{
	ostrm<<"LineStipple="; MGGLAttrib::toString(ostrm);
	if(undefined()) return ostrm;
	ostrm<<",pattern="<<m_pattern;
	return ostrm;
}

void MGLineStipple::WriteMembers(MGOfstream& buf)const{
	MGGLAttrib::WriteMembers(buf);
	if(undefined()) return;
	buf<<m_pattern;
}
void MGLineStipple::ReadMembers(MGIfstream& buf){
	MGGLAttrib::ReadMembers(buf);
	if(undefined()) return;
	buf>>m_pattern;
}

AUTO_GEL_REGISTER(MGColor, MGCOLOR_TID);
AUTO_GEL_REGISTER(MGLineWidth, MGLINE_WIDTH_TID);
AUTO_GEL_REGISTER(MGLineStipple, MGLINE_STIPPLE_TID);
AUTO_GEL_REGISTER(MGName, MGNAME_TID);
#ifndef _CONSOLE
AUTO_GEL_REGISTER(MGLight, MGLIGHT_TID);
AUTO_GEL_REGISTER(MGDirectionalLight, MGDIRECTIONAL_LIGHT_TID);
AUTO_GEL_REGISTER(MGPointLight, MGPOINT_LIGHT_TID);
AUTO_GEL_REGISTER(MGSpotLight, MGSPOT_LIGHT_TID);
AUTO_GEL_REGISTER(MGLights, MGLIGHTS_TID);
AUTO_GEL_REGISTER(MGRenderAttr, MGRENDER_ATTR_TID);
#endif