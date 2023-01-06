/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Ofstream.h"
#include "mg/Ifstream.h"
#include "mg/DNameControl.h"
#include "mgGL/OpenGLView.h"
#include "mgGL/VBO.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

MGDNameControl::MGDNameControl():m_nextName(OpenGLStartDisplayName()){;}

///Get the singleton instance of the MGDNameControl.
MGDNameControl& getDNameControlInstance(){
	static MGDNameControl m_displayName;
	return m_displayName;
}

///Get mgVBO* from the map.
///Null will be returned when dlname is not registered.
mgVBO* MGDNameControl::VBO_from_dlistName(unsigned dlname){
	Dlist2VBOMap& d2g=m_dlist2VBOMap;
	Dlist2VBOMap::iterator i, iend=d2g.end();
	i=d2g.find(dlname);
	if(i!=iend){
		mgVBO* vbo=i->second;
		return vbo;
	}
	return 0;
}

MGAttribedGel* MGDNameControl::Gelpointer_from_dlistName(unsigned dlname){
	Dlist2VBOMap& d2g=m_dlist2VBOMap;
	Dlist2VBOMap::iterator i, iend=d2g.end();
	i=d2g.find(dlname);
	if(i!=iend){
		mgVBO* vbo=i->second;
		return vbo->gel();
	}
	return 0;
}

///Insert vbo to dlist map.
void MGDNameControl::insertVBO2DlistMap(unsigned dlistName,mgVBO* vbo){
	m_dlist2VBOMap.insert(Dlist2VBOMapPair(dlistName,vbo));
}

///Insert gel to dlist map. Function's return value is the mgVBO generated.
mgVBO* MGDNameControl::insertDlistMap(const MGAttribedGel* gel){
	MGAttribedGel* gel2=const_cast<MGAttribedGel*>(gel);
	mgVBO* vbo=getNewName(gel);
	return vbo;
}

mgVBO* MGDNameControl::deleteDlistMap(unsigned dlname){
	Dlist2VBOMap& d2g=m_dlist2VBOMap;
	Dlist2VBOMap::iterator i, iend=d2g.end();
	mgVBO* vbo=0;
	i=d2g.find(dlname);
	if(i!=iend){
		vbo=i->second;
		MGAttribedGel* gel=vbo->gel();
		if(gel){
			gel->m_VBO.release();
		}
		d2g.erase(i);
	}
	return vbo;
}
mgVBO* MGDNameControl::deleteDlistMap(const MGAttribedGel* gel){
	const mgVBO* vbo=gel->m_VBO.get();
	if(!vbo)
		return 0;
	unsigned nm=vbo->getDName();
	return deleteDlistMap(nm);
}

///Get the next available name inputing vbo
///The name will be returned.
unsigned MGDNameControl::getNewNamebyVBO(mgVBO* vbo){
	unsigned name=m_nextName;
	m_nextName++;
	insertVBO2DlistMap(name,vbo);
	return name;
}

///Get the next available name.
mgVBO* MGDNameControl::getNewName(const MGAttribedGel* gel){
	mgVBO* vbo;
	if(gel==0){
		vbo=new mgVBO();
	}else{
		vbo=new mgVBO(*gel);
	}
	unsigned name=getNewNamebyVBO(vbo);
	vbo->setDlName(name);
	return vbo;
}
