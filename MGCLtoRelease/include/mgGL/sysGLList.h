/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGSYSGLList_HH_
#define _MGSYSGLList_HH_

#include <vector>
#include <utility>
#include <bitset>
#include <list>
#include "mg/MGCL.h"
#include "mgGL/sysGL.h"

class MGGel;
class MGCurve;

/** @addtogroup DisplayHandling
 *  @{
 */

///Defines a list of mgSysGL.

///mgSysGLList is a class to constrol system generated display list.
///System generated display list means the display lists not generated
///from the document, m_display_list of the document.
///mgSysGLList returns a display list name(mgVBO*) when an mgSysGL(fucntion code, object id)
///is added(push_back, or push_front). Using the display list name, mgSysGLList invokes
///mgVBO->draw().
///Then mgSysGLList will manage to delete the display list,
/// giving the fucntion code or the object id.
class MG_DLL_DECLR mgSysGLList{

public:

typedef std::list<mgSysGL*> container_type;

typedef container_type::iterator iterator;
typedef container_type::const_iterator const_iterator;

///////////////Constructor/////////////

mgSysGLList();

///Copy constructor.
mgSysGLList(const mgSysGLList& list2);

////////////Destructor////////////////
~mgSysGLList();

////////////////Operator overload//////////////////

///Assignment operator.
mgSysGLList& operator=(const mgSysGLList& list2);

/// オペレーション

///const_iterator begin()const{return m_sysgls.begin();};
///iterator begin(){return m_sysgls.begin();};

///Clear the list. Delete all the display lists.
void clear();

///Delete all the display lists that have the fucntion_code fc.
///function's retrurn value is true if any one is deleted.
bool delete_lists_by_function_code(int fc);

///Delete all the display lists that have the object_id oi.
///function's retrurn value is true if any one is deleted.
bool delete_lists_by_object_id(
	const MGGel* oi,///<target gel to delete.
	std::vector<UniqueSysGL>& functions	///<mgSysGL pointer will be appended.
);

///Delete the display list that have the fucntion_code fc and the object id gel.
///function's retrurn value is true if any one is deleted.
bool delete_lists_by_function_object_code(int fc, const MGGel* gel);

///Draw all the objects by calling glCallList in this list.
void draw_list(MGCL::VIEWMODE vMode);

//ファンクションコードに対応するSysGLを取得する
mgSysGL* getSysGLByFunctionCode(int fc);

///Test if this list includes the fucntion code fc's SysGL or not.
bool includes(int fc)const;

///const_iterator end()const{return m_sysgls.end();};
///iterator end(){return m_sysgls.end();};

///Push back (function_code, object id) to the end or the beginning of the list.
///Function's return value is the mgSysGL pointer pushed.
mgSysGL* push_back(int fc, const MGGel* oi);
mgSysGL* push_front(int fc, const MGGel* oi);

///sysgl must be a newed object, and the ownership will be 
///transfered to this.
void push_back(mgSysGL* sysgl);

///Get the size of this list.
int size()const{return int(m_sysgls.size());};

private:

///In m_sysgls mgSysGL pointers are stored.
container_type m_sysgls;			///<list of mgSysGL.

};

/** @} */ // end of DisplayHandling group
#endif
