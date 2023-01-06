/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Ofstream.h"
#include "mg/Ifstream.h"
#include "mgGL/Name.h"
#include "mgGL/VBO.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

MGName::MGName():MGGLAttrib(static_cast<int>(mgGLMode::UNDEFINED)){;};
MGName::MGName(const std::string& name):MGGLAttrib(1),
m_name(name){;}
	
MGName* MGName::clone()const{
	return new MGName(*this);
}

void MGName::setName(std::string& name){
	m_flag=1;
	m_name=name;
}
void MGName::setName(char* name){
	m_flag=1;
	m_name=name;
}

// Output function.
std::ostream& MGName::toString(std::ostream& ostrm) const{
	ostrm<<"Name="; MGGLAttrib::toString(ostrm);
	if(enabled()){
		ostrm<<"="<<m_name;
	}
	return ostrm;
}

//Write all member data
void MGName::WriteMembers(MGOfstream& buf)const{
	MGGLAttrib::WriteMembers(buf);
	if(undefined()) return;

	int nchar=(int)m_name.size();
	buf<<nchar;
	buf.writenChar(m_name.data(),nchar);
}
//Read all member data.
void MGName::ReadMembers(MGIfstream& buf){
	MGGLAttrib::ReadMembers(buf);
	if(undefined()) return;

	int nchar;
	buf>>nchar;
	char* c=new char[nchar+1];
	buf.readnChar(c,nchar);
	c[nchar]=0;
	m_name=std::string(c);
	delete[] c;
}
