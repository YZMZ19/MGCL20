/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGDNAMECONTROL_HH_
#define _MGDNAMECONTROL_HH_


#include <map>
#include "MGCL.h"
#include "mgGL/VBO.h"
class MGAttribedGel;

/// @cond

///MGDNameControlはOpenGLに利用するDisplay Listの名前を制御します。
///OpenGLのglGenListsの代替となります。
///MGDNameControlはいわゆるsingletonでひとつのプログラムにひとつの
///オブジェクトのみ存在します。
///最初の値はOpenGLStartDisplayName()
///(DISPNAME_BASE_ID)からとなります。
///

class MG_DLL_DECLR MGDNameControl{

public:
	
typedef std::map<unsigned, mgVBO*> Dlist2VBOMap;
typedef Dlist2VBOMap::value_type Dlist2VBOMapPair;

///Get the singleton instance of the MGDNameControl.
MG_DLL_DECLR friend MGDNameControl& getDNameControlInstance();

MGDNameControl();

///Get MGAttribedGel* from the map.
///Null will be returned when dlname is not registered.
MGAttribedGel* Gelpointer_from_dlistName(unsigned dlname);

///Get mgVBO* from the map.
///Null will be returned when dlname is not registered.
mgVBO* VBO_from_dlistName(unsigned dlname);

///Get the next available name.
mgVBO* getNewName(const MGAttribedGel* gel=0);

///Get the next available name inputing vbo
///The name will be returned.
unsigned getNewNamebyVBO(mgVBO* vbo);

///Insert gel to dlist map. Function's return value is the mgVBO generated.
mgVBO* insertDlistMap(const MGAttribedGel* gel);

///Delete gel's vbo and the map of (dlname, vbo).
///VBO of this name will be returned(must be deleted).
mgVBO* deleteDlistMap(unsigned dlname);

///Delete dlname's vbo and the map of (dlname, vbo).
///VBO of this gel will be returned(must be deleted).
mgVBO* deleteDlistMap(const MGAttribedGel* gel);

private:

unsigned m_nextName;///<次に払いだされる名前（unsigned integer)。
Dlist2VBOMap m_dlist2VBOMap;///<Map to obtain MGAttribedGel* from the dlistName.

///Insert vbo to dlist map.
void insertVBO2DlistMap(unsigned dlistName,mgVBO* vbo);

};

///Get the singleton instance of the MGDNameControl.
MG_DLL_DECLR MGDNameControl& getDNameControlInstance();

/// @endcond
#endif // _MGDNAMECONTROL_HH_
