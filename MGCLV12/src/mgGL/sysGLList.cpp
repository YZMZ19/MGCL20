/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
// mgSysGLList.cpp : mgSysGLList クラスのimplementation。
#include "StdAfx.h"
#include "mg/Curve.h"
#include "mgGL/SysGLList.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

///////////////Member fucntions/////////////////

//Clear the list.
void mgSysGLList::clear(){m_sysgls.clear();}

//Delete all the display lists that have the fucntion_code fc.
//function's retrurn value is true if any one is deleted.
bool mgSysGLList::delete_lists_by_function_code(int fc){
	bool deleted=false;
	iterator i=m_sysgls.begin(), in, ie=m_sysgls.end();
	in=i;
	while(i!=ie){
		in++;
		mgSysGL& sysgl=**i;
		if(sysgl.function_code()==fc){
			deleted=true;
			m_sysgls.erase(i);
		}
		i=in;
	}
	return deleted;
}

//Delete all the display lists that have the object_id oi.
//function's retrurn value is true if any one is deleted.
bool mgSysGLList::delete_lists_by_object_id(
	const MGGel* oi,
	std::vector<UniqueSysGL>& functions	//Deleted mgSysGL pointer  appended.
){
	bool deleted=false;
	iterator i=m_sysgls.begin(), in, ie=m_sysgls.end();
	in=i;
	while(i!=ie){
		in++;
		UniqueSysGL& sysgl=*i;
		if(sysgl->includes(oi)){
			functions.emplace_back(std::move(sysgl)); 
			deleted=true;
			m_sysgls.erase(i);
		}
		i=in;
	}

	return deleted;
}

//Delete the display list that have the fucntion_code fc and the object id gel.
//function's retrurn value is true if any one is deleted.
bool mgSysGLList::delete_lists_by_function_object_code(
	int fc, const MGGel* gel
){
	bool deleted=false;
	iterator i=m_sysgls.begin(), in, ie=m_sysgls.end();
	in=i;
	while(i!=ie){
		in++;
		UniqueSysGL& sysgl = *i;
		if(sysgl->function_code()==fc && sysgl->includes(gel)){
			deleted=true;
			m_sysgls.erase(i);
		}
		i=in;
	}
	return deleted;
}

//Draw all the objects by calling glCallList in this list.
void mgSysGLList::draw_list(MGCL::VIEWMODE vMode){
	for (auto& sg : m_sysgls)
		sg->draw(vMode);
}

//ファンクションコードに対応するSysGLを取得する
mgSysGL* mgSysGLList::getSysGLByFunctionCode(int fc){
	iterator i = m_sysgls.begin(), ie = m_sysgls.end();
	for(;i != ie; ++i){
		UniqueSysGL& sysgl = *i;
		if(sysgl->function_code()==fc){
			return sysgl.get();
		}
	}
	return 0;
}

//Test if this list includes the fucntion code fc's SysGL or not.
bool mgSysGLList::includes(int fc)const{
	const_iterator i=m_sysgls.begin(), ie=m_sysgls.end();
	for(;i!=ie; i++){
		if((**i).function_code()==fc)
			return true;
	}
	return false;
}

mgSysGL* mgSysGLList::push_back(int fc,const MGGel* oi){
	mgSysGL* sgl=new mgSysGL(fc,oi);
	m_sysgls.emplace_back(sgl);
	return sgl;
}
mgSysGL* mgSysGLList::push_front(int fc,const MGGel* oi){
	mgSysGL* sgl=new mgSysGL(fc,oi);
	m_sysgls.emplace_front(sgl);
	return sgl;
}

//sysgl must be a newed object, and the ownership will be 
//transfered to this.
void mgSysGLList::push_back(mgSysGL*&& sysgl){
	for(auto& sg:m_sysgls)
		if(sg.get() == sysgl)
			return;
	m_sysgls.emplace_back(std::move(sysgl));
}
