/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Ofstream.h"
#include "mg/Ifstream.h"
#include "mgGL/Light.h"
#include "mgGL/Lights.h"
#include "mgGL/GLSLProgram.h"
#include "mgGL/VBO.h"
#include "mgGL/Appearance.h"


#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

////////////////////////////////////////////////////////////////////

MGLights::MGLights(const MGLights & rhs):MGGLAttrib(rhs),m_lights(rhs.m_lights.size()){
	size_t i=0;
	for(auto& j:rhs.m_lights){
		m_lights[i++].reset(j->clone());
	}
}
//assignment
MGLights& MGLights::operator=(const MGLights& gel2){
	if(this==&gel2)
		return *this;

	MGGLAttrib::operator=(gel2);
	m_lights.clear();
	int n=(int)gel2.m_lights.size();
	for(int i=0; i<n; i++){
		m_lights.emplace_back(gel2.m_lights[i]->clone());
	}
	return *this;
}
MGLights& MGLights::operator=(const MGGel& gel2){
	const MGLights* gel2_is_this=dynamic_cast<const MGLights*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGLights::operator<(const MGLights& gel2)const{
	return m_lights.size()<gel2.m_lights.size();
}
bool MGLights::operator<(const MGGel& gel2)const{
	const MGLights* gel2_is_this=dynamic_cast<const MGLights*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

MGLights::MGLights():MGGLAttrib(1){
	//m_lights.push_back(new MGLight);
}

MGLights* MGLights::clone()const{
	MGLights* lights=new MGLights;
	const_iterator i=begin(), ie=end();	
	for(; i!=ie; i++){
		lights->push_back(static_cast<MGLight*>((**i).clone()));
	}
	return lights;
}

void MGLights::exec()const{

	mgGLSLProgram* glsl=mgGLSLProgram::getCurrentGLSLProgram();
	glsl->setUniform(mgGLSLProgram::ShaderMode, mgGLSL::NoShading);
	if(undefined())
		return;

	for(GLint idx = 0; idx<10;++idx){
		MGLight::ResetLight(idx);
	}

	int n=size();
	if(!n)
		return;

	for(int i=0; i<n; i++)
		m_lights[i]->exec();
//	glEnable(GL_LIGHTING);// ライティング開始
	glsl->setUniform(mgGLSLProgram::ShaderMode, mgGLSL::Shading);
	std::cerr << *this;
}

// Output function.
std::ostream& MGLights::toString(std::ostream& ostrm) const{
	ostrm<<"Lgihts=";MGGLAttrib::toString(ostrm);
	int n=size();
	ostrm<<",number of lights="<<n<<std::endl;
	for(int i=0; i<n; i++) m_lights[i]->toString(ostrm);
	return ostrm;
}

//add one light.
//Function's return value is the numbe of lights defined.
int MGLights::push_back(MGLight* light){
	int lnum=(int)m_lights.size();
	light->set_light_number(lnum);
	m_lights.emplace_back(light);
	return (int)m_lights.size();
}

//Read all member data.
void MGLights::ReadMembers(MGIfstream& buf){
	MGGLAttrib::ReadMembers(buf);
	int n=size();
	buf>>n;
	for(int i=0; i<n; i++)
		push_back(static_cast<MGLight*>(buf.ReadPointer()));
}
//Write all member data
void MGLights::WriteMembers(MGOfstream& buf)const{
	MGGLAttrib::WriteMembers(buf);
	int n=size();
	buf<<n;
	const_iterator i=begin(), ie=end();	
	for(; i!=ie; i++)
		buf.WritePointer(i->get());
}
