/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#if !defined(AFX_SYS2GL_H__1EC223A1_C0D3_4AA7_80BD_89392C970098__INCLUDED_)
#define AFX_SYS2GL_H__1EC223A1_C0D3_4AA7_80BD_89392C970098__INCLUDED_

#include "mgGL/SysGL.h"

/** @addtogroup DisplayHandling
 *  @{
 */

///Defines two mgSysGL objects.
class MG_DLL_DECLR mgSys2GL: public mgSysGL{

public:

mgSys2GL(
	int function_code,///<Function code.
	const MGGel* gel1,///<1st gel.
	const MGGel* gel2 ///<2nd gel.
);

///Get gel1.
const MGGel* gel1()const{return object_id();};

///Get gel2.
const MGGel* gel2()const{return m_gel2;};

///Test if this mgSysGL includes gel(return true) or not.
virtual bool includes(const MGGel* gel)const;

///replace gel_old to gel_new.
///If gel_old is not included in this, do nothing.
virtual void replace(
	const MGGel* gel_old,	///<gel_old must be a MGCurve.
	const MGGel* gel_new	///<gel_new must be a MGFSurface.
);

///Output to stream file.
virtual std::ostream& toString(std::ostream& ostrm) const;

private:
	const MGGel* m_gel2;///<MGGel. 1st gel is stored in mgSysGL's m_gel.

};

/** @} */ // end of DisplayHandling group

#endif // !defined(AFX_SYS2GL_H__1EC223A1_C0D3_4AA7_80BD_89392C970098__INCLUDED_)
