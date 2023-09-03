/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGSYSGL_HH_
#define _MGSYSGL_HH_

#include <iosfwd>
#include "mg/MGCL.h"
#include "mgGL/VBO.h"
#include "mgGL/Appearance.h"

class MGGel;
class MGOpenGLView;

/** @file */
/** @addtogroup DisplayHandling
 *  @{
 */

///mgSysGL is a class to provide a facility to draw temporal pictures.

///MGOpenGLView holds a list of mgSysGL and draws the pictures by invoking
///display list drawer.
///As long as the codes are unique, function codes can be
///any numbers. Usually the id is a command id.
class MG_DLL_DECLR mgSysGL: public mgVBO{

public:

///Output to stream.
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream& outp, const mgSysGL& sysgl);

///////////////Constructor/////////////

mgSysGL();
mgSysGL(int fucntion_code,const MGGel* object_id);

///Copy constructor.
mgSysGL(const mgSysGL& sp);

///Copy constructor, replacing gel_old to gel_new.
mgSysGL(mgSysGL& glold,const MGGel* gel_old, const MGGel* gel_new);

virtual ~mgSysGL();

////////////////Operator overload//////////////////

///Assignment
mgSysGL& operator=(const mgSysGL& sgl);

/////////////

///Construct new object by copying to newed area.
///User must delete this copied object by "delete".
virtual mgSysGL* clone()const;

///Draw this Sysgl.
///This drawSysGL is used to drawSysGL the pictures for Undo(, Redo) operations.
virtual void drawSysGL(){;};

///(1) initializeVBO
///(2)drawSysGLで描画データを作成しなおす。
virtual void make_display_list(MGCL::VIEWMODE vmode=MGCL::DONTCARE);

///Make system display list in glv  by invoking the above drawSysGL() function.
virtual void makeSysGLDisplayList(MGOpenGLView& glv);

///Get the i-th ellement's object_id or fucntion_code.
int function_code()const{return m_fucntion_id;};

///Test if this mgSysGL includes gel(return true) or not.
///The default includes tests if the input gel is m_gel of this
///member data.
virtual bool includes(const MGGel* gel)const;

/// Output virtual function.
///Output to stream file:メンバデータを標準出力に出力する。
virtual std::ostream& toString(std::ostream& ostrm) const;

///replace gel_old to gel_new.
///If gel_old is not included in this, do nothing.
virtual void replace(const MGGel* gel_old, const MGGel* gel_new);

///Set function code.
void set_function_code(int fc){m_fucntion_id=fc;};

protected:

///Get gel.
const MGGel* object_id()const{return m_gel;};

///Set gel.
void set_object_id(const MGGel* oi){m_gel=oi;};

private:

	int m_fucntion_id;///<fucntion code.
	const MGGel* m_gel;	///<Object id. When more than 1 objects are concerned, 
						///<mgSysGL class will be inheritted and the subclass
						///<will retain the objcets.
};

/** @} */ // end of DisplayHandling group
#endif
