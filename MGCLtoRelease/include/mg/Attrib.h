/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGAttrib_HH_
#define _MGAttrib_HH_

#include <vector>
#include <memory>
#include "mg/Gel.h"

//
//Define MGAttrib Class.

class MGIfstream;
class MGOfstream;
class MGAttrib;
class MGGLAttrib;
class MGContext;

/** @addtogroup GelRelated
 *  @{
 */


///MGAttrib is an abstract class that defines attribute elements of MGGel.

///Currently main subclasses of MGAttrib are MGGLAttrib, MGAppearance, and MGContex.
class MG_DLL_DECLR MGAttrib: public MGGel{

public:

////////Special member functions/////////
MGAttrib()=default;
virtual ~MGAttrib()=default;
MGAttrib(const MGAttrib&)=default;///Copy constructor.
MGAttrib& operator=(const MGAttrib&)=default;
MGAttrib(MGAttrib&&)=default;		///Move constructor.
MGAttrib& operator= (MGAttrib&&)=default;///Move assignment.

///draw attribute data.
virtual void drawAttrib(
	mgVBO& vbo,///<The target graphic object.
	bool no_color=false	///<if true, color attribute will be neglected.
)const=0;

/// Return This object's typeID
///Sub class of MGAttrib must have the id as 0x02nnnnnnL.
virtual long identify_type() const{return MGATTRIB_TID;};

friend class MGIfstream;
friend class MGOfstream;

};

/** @} */ // end of GelRelated group
#endif
