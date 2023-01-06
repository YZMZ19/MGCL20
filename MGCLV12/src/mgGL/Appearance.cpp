/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#include "StdAfx.h"
#include "mg/Attrib.h"
#include "mg/Ofstream.h"
#include "mg/Ifstream.h"
#include "mgGL/Color.h"
#include "mgGL/LineWidth.h"
#include "mgGL/Appearance.h"
#include "mgGL/LineStipple.h"
#include "mgGL/VBO.h"
#include "mgGL/Light.h"
#include "mgGL/Lights.h"

//
//Define MGAppearance Class.
//MGAppearance is a class to contain MGGLAttrib objects.
//MGAppearance acts just like as std::auto_ptr.
//That is, MGAppearance holds newed object pointers of MGGLAttrib,
//and when copy constructor or assignment operator is invoked,
//the pointer ownership is transfered to the new MGAppearance object.

MGAppearance& MGAppearance::operator=(MGAppearance&& gel2){
	m_glattribs=std::move(gel2.m_glattribs);
	m_no_display=gel2.m_no_display;
	return *this;
}

bool MGAppearance::operator<(const MGAppearance& gel2)const{
	size_t n1=size(), n2=gel2.size();
	if(n1==n2){
		return size_t(this)<size_t(&gel2);
	}else
		return n1<n2;
}
bool MGAppearance::operator<(const MGGel& gel2)const{
	const MGAppearance* gel2_is_this=dynamic_cast<const MGAppearance*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

//Test if this MGAppearance can be removed or not.
bool MGAppearance::can_be_removed()const{
	return (size()==0 && !no_display());
}

//Generate copied gel of this gel.
//Returned is a newed object. User must delete the object.
MGAppearance* MGAppearance::clone()const{
	MGAppearance* appr2=new MGAppearance;
	const_iterator i=begin(), ie=end();	
	for(; i!=ie; i++){
		MGAttrib* attri=static_cast<MGAttrib*>((**i).clone());
		appr2->push_back(attri);
	}
	appr2->m_no_display=m_no_display;
	return appr2;
}

std::ostream& MGAppearance::toString(
	std::ostream& ostrm
)const{
	ostrm<<"MGAppearance="<<this<<",no_display="<<m_no_display;
	ostrm<<", number of attribs="<<size()<<std::endl;
	const_iterator i=begin(), ie=end();	
	for(int j=0; i!=ie; i++, j++){
		ostrm<<"  attr["<<j<<"]::"<<(**i)<<std::endl;
	}
	return ostrm;
}

//////////Member Function//////////

//Release the specified attribute.
//Function's return value is the MGGLAttrib* that is released.
std::unique_ptr<MGGLAttrib> MGAppearance::release_attrib(long tid){
	iterator i=begin(), ie=end();
	for(;i!=ie; i++){
		if((*i)->identify_type()==tid){
			std::unique_ptr<MGGLAttrib> gla
				=std::unique_ptr<MGGLAttrib>(static_cast<MGGLAttrib*>(i->release()));
			m_glattribs.erase(i);
			return gla;
		}
	}
	return 0;
}

//Set the attribute in this list. attr must be a newed object, and the
//ownership will be transfered to this MGAppearance.
void MGAppearance::set_attrib(MGGLAttrib* attr){
	MGGLAttrib* olda=set_attrib_with_old(attr);
	if(olda)
		delete olda;
}
void MGAppearance::set_attrib(UniqueGLAttribVec& attrs){
	int n=(int)attrs.size();
	for(int i=0; i<n; i++){
		MGGLAttrib* attr=attrs[i].release();
		MGGLAttrib* olda=set_attrib_with_old(attr);
		if(olda)
			delete olda;
	}
}

//Set the attribute in this list. attr must be a newed object, and the
//ownership will be transfered to this MGAppearance.
//When the appearance held an attribute, the old one will be returned
//as the function's return value. Users must delete it.
MGGLAttrib* MGAppearance::set_attrib_with_old(MGGLAttrib* attr){
	if(!attr)
		return 0;

	iterator i=search(attr);
	if(i!=end()){//If found.
		MGGLAttrib* gla=static_cast<MGGLAttrib*>(i->release());
		iterator j=m_glattribs.erase(i); insert(j,attr);
		return gla;
	}else{//If not found.
		push_back(attr);
		return 0;
	}
}

//Functional object for find_if.
class MGAppearanceSearch{
public:
	MGAppearanceSearch(const MGGLAttrib* atr):m_attr(atr){;};
	bool operator()(const std::unique_ptr<MGGel>& atr2){
		return m_attr->same_type(*(static_cast<const MGGLAttrib*>(atr2.get())));
	};
	const MGGLAttrib* m_attr;
};

//Search the same MGGLAttrib leaf class object in this list.
//If not found, end() will be returned.
MGAppearance::iterator MGAppearance::search(const MGGLAttrib* atr){
	return std::find_if(begin(), end(), MGAppearanceSearch(atr));
}
MGAppearance::const_iterator MGAppearance::search(const MGGLAttrib* atr)const{
	return std::find_if(begin(), end(), MGAppearanceSearch(atr));
}
MGAppearance::iterator MGAppearance::search_by_id(MGGEL_TID tid){
	MGAppearance::iterator i=begin(), ie=end();
	for(;i!=ie; i++){
		const MGGLAttrib* gla=static_cast<const MGGLAttrib*>(i->get());
		if(gla->identify_type()==tid) return i;
	}
	return ie;
}
MGAppearance::const_iterator MGAppearance::search_by_id(MGGEL_TID tid)const{
	MGAppearance::const_iterator i=begin(), ie=end();
	for(;i!=ie; i++){
		const MGGLAttrib* gla=static_cast<const MGGLAttrib*>(i->get());
		if(gla->identify_type()==tid) return i;
	}
	return ie;
}
//Turn on the appropriate mask bit for this attribute. See glPushAttrib().
int MGAppearance::get_draw_attrib_mask()const{
	unsigned int mask=0;
	MGAppearance::const_iterator i=begin(), ie=end();
	for(;i!=ie; i++){
		const MGGLAttrib* gla=static_cast<const MGGLAttrib*>(i->get());
		gla->set_draw_attrib_mask(mask);
	}
	return mask;
}

//Turn on the appropriate mask bit for this attribute. See glPushAttrib().
int MGAppearance::get_render_attrib_mask()const{
	unsigned int mask=0;
	MGAppearance::const_iterator i=begin(), ie=end();
	for(;i!=ie; i++){
		const MGGLAttrib* gla=static_cast<const MGGLAttrib*>(i->get());
		gla->set_render_attrib_mask(mask);
	}
	return mask;
}

//メンバデータを書き込む関数
void MGAppearance::WriteMembers(MGOfstream& buf)const{
	m_glattribs.WriteMembers(buf);
	int no_disp=1;
	if(!no_display())
		no_disp=0;
	buf<<no_disp;
}

//メンバデータを読み出す関数
void MGAppearance::ReadMembers(MGIfstream& buf){
	m_glattribs.ReadMembers(buf);
	int no_disp;
	buf>>no_disp;
	m_no_display=false; if(no_disp) m_no_display=true;
}

///draw GLAttributes process.
void MGAppearance::drawAttrib(
	mgVBO& vbo,///<The target graphic object.
	bool no_color	//if true, color attribute will be neglected.
)const {
	if (no_display())
		return;

	const_iterator i = begin(), ie = end();
	for (; i != ie; i++) {
		const MGGLAttrib* gla = dynamic_cast<const MGGLAttrib*>(i->get());
		if (gla) {
			if (no_color && gla->is_highlight_attrib())
				continue;
			gla->drawAttrib(vbo);
		}
	}
}

//render GLAttributes process.
void MGAppearance::render(mgVBO& vbo)const {
	if (no_display())
		return;

	const_iterator i = begin(), ie = end();
	for (; i != ie; i++) {
		const MGGLAttrib* gla = dynamic_cast<const MGGLAttrib*>(i->get());
		if (gla)
			gla->render(vbo);
	}
}

#ifndef _CONSOLE

//Set the material. When rs=FRONT_AND_BACK and different material for the back side
//is used, set_back_material must be invoked after invoking set_material.
//Else the same material will be appllied for the both sides.
void MGAppearance::set_material(
	MGRenderAttr::RENDERSIDE rs,
	const float ambient[3],
	const float diffuse[3],
	const float specular[3],
	const float emission[3],
	float shininess,
	float transparency
){
	MGRenderAttr* ra;
	iterator i=search_by_id(MGRENDER_ATTR_TID);
	if(i==end()){
		ra=new MGRenderAttr();
		push_back(ra);
	}else ra=static_cast<MGRenderAttr*>(i->get());
	ra->set_material(rs,ambient,diffuse,specular,emission,shininess,transparency);
}

//Set the back side material. Invoking set_back_material means two sided material
//and setting different material to the back side.
//Before use of set_back_material, set_material must be invoked first.
//set_back_material will set two sided material.
void MGAppearance::set_back_material(
	const float ambient[3],
	const float diffuse[3],
	const float specular[3],
	const float emission[3],
	float shininess,
	float transparency
){
	MGRenderAttr* ra;
	iterator i=search_by_id(MGRENDER_ATTR_TID);
	if(i==end()){
		ra=new MGRenderAttr();
		push_back(ra);
	}else ra=static_cast<MGRenderAttr*>(i->get());
	ra->set_back_material(
		ambient,diffuse,specular,emission,shininess,transparency);
}

//Add a light. light must be a newed object, and the ownership will be
//transfered to this object.
//Function's return value is the number of lights after added.
int MGAppearance::add_light(MGLight* light) {
	iterator i = search_by_id(MGLIGHTS_TID);
	MGLights* lights;
	if (i == end()) {
		lights = new MGLights;
		push_back(lights);
	}
	else
		lights = static_cast<MGLights*>(i->get());
	lights->push_back(light);
	return lights->size();
}

#endif //_CONSOLE

void MGAppearance::set_color(const MGColor& color) {
	MGColor* colr = new MGColor(color);
	set_attrib(colr);
}
void MGAppearance::set_color(const float color[4]) {
	MGColor* colr = new MGColor(color[0], color[1], color[2], color[3]);
	set_attrib(colr);
}
void MGAppearance::set_color(float red, float green, float blue, float alpha) {
	MGColor* colr = new MGColor(red, green, blue, alpha);
	set_attrib(colr);
}

void MGAppearance::setLineWidth(float width) {
	MGLineWidth* lw = new MGLineWidth(width);
	set_attrib(lw);
}

///Line stipple属性をセットする。
///When factor=0 is input, line pattern is disabled. This means the line is solid.
///When factor<0, the stipple attribute is undefined. This means the attribute
///is defined by the environment.
///When factor<=0, pattern is unnecessary.
void MGAppearance::setLineStipple(short int factor,unsigned short pattern){
	MGLineStipple* ls=new MGLineStipple(factor,pattern);
	set_attrib(ls);
}

#ifdef _CONSOLE

int MGAppearance::add_light(MGLight* light) {
	std::cout << "MGAppearance::add_light, light=" << light << std::endl;
	return 0;
}

namespace mgGLSL {
	void execStaticColorAttrib(const MGColor& color) {
		std::cout << "execStaticColorAttrib, color="<<&color<<std::endl;
	}
	void execStaticGLAttrib(const mgStaticGLAttrib& attrib) {
		std::cout << "execStaticGLAttrib, attrib=" << &attrib << std::endl;
	}
	void execStaticColorAttrib(const float color[4]) {
		std::cout << "execStaticColorAttrib, color=" << color << std::endl;
	}
	void execStaticLineWidth(float lineWidth) {
		std::cout << "execStaticLineWidth, lineWidth=" << lineWidth << std::endl;
	}
	void execStaticLineStipple(short int factor, GLuint pattern) {
		std::cout << "execStaticLineStipple, factor="
			<< factor<<", pattern="<<pattern << std::endl;
	}
	void execLightMode(int mode) {
		std::cout << "execLightMode, mode=" << mode << std::endl;
	}
};

#endif // _CONSOLE
