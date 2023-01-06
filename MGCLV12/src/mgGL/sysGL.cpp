/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
// mgSysGLList.cpp : mgSysGLList クラスのimplementation。
#include "StdAfx.h"
#include "mg/DNameControl.h"
#include "mgGL/OpenGLView.h"
#include "mgGL/SysGL.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif


mgSysGL::mgSysGL():m_fucntion_id(0),m_gel(0){;}

mgSysGL::mgSysGL(int fucntion_code,const MGGel* object_id)
:m_fucntion_id(fucntion_code),m_gel(const_cast<MGGel*>(object_id)){
}

//Copy constructor, replacing gel_old to gel_new.
mgSysGL::mgSysGL(
	mgSysGL& glold,
	const MGGel* gel_old,
	const MGGel* gel_new
):mgVBO(glold),m_fucntion_id(glold.m_fucntion_id),m_gel(glold.m_gel){
	replace(gel_old,gel_new);	
}

///Copy constructor.
mgSysGL::mgSysGL(const mgSysGL& sp)
:mgVBO(sp),m_fucntion_id(sp.m_fucntion_id),m_gel(sp.m_gel){
}

///Assignment
mgSysGL& mgSysGL::operator=(
	const mgSysGL& sp
){
	mgVBO::operator=(sp);
	m_fucntion_id=sp.m_fucntion_id;
	m_gel=sp.m_gel;
	return *this;
}

mgSysGL::~mgSysGL(){
}

//Construct new object by copying to newed area.
//User must delete this copied object by "delete".
mgSysGL* mgSysGL::clone()const{
	return new mgSysGL(*this);
}

///Test if this mgSysGL includes gel(return true) or not.
///The default includes tests if the input gel is m_gel of this
///member data.
bool mgSysGL::includes(const MGGel* gel)const{
	if(!gel || !m_gel)
		return false;
	return gel==m_gel;
}

///(1) initializeVBO
///(2)drawSysGLで描画データを作成しなおす。
void mgSysGL::make_display_list(MGCL::VIEWMODE vmode){
	initializeVBO(vmode);
	drawSysGL();
	setDirty(false);
}

//Make system display list in glv.
//This must be a newed object and the ownership will be transfered to
//glv(glv.m_sysgllist).
void mgSysGL::makeSysGLDisplayList(
	MGOpenGLView& glv
){
	initializeVBO();
	drawSysGL();
	glv.push_back_to_sysgl(this);
	setDirty(false);
}

//replace gel_old to gel_new.
//If gel_old is not included in this, do nothing.
void mgSysGL::replace(
	const MGGel* gel_old, const MGGel* gel_new
){
	if(m_gel==gel_old)
		m_gel=const_cast<MGGel*>(gel_new);
}

// Output virtual function.
//Output to stream file:メンバデータを標準出力に出力する。
std::ostream& mgSysGL::toString(std::ostream& ostrm) const{
//	ostrm.setf(ios::scientific,ios::floatfield);
//	ostrm.precision(10);
	ostrm<<"mgSysGL::"<<this;
	ostrm<<",m_fucntion_id="<<m_fucntion_id
		<<",m_gel="<<m_gel;
	return ostrm;
}

//////////// mgSysGL output ////////////
std::ostream& operator<<(std::ostream& outp, const mgSysGL& sysgl){
	sysgl.toString(outp);
	return outp;
}