/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGGel_HH_
#define _MGGel_HH_

#include "StdAfx.h"
#include "mg/MGCL.h"
#include "mg/types.h"

class MGDrawParam;
class MGAbstractGels;
class MGAttribedGel;
class MGIfstream;
class MGOfstream;
class MGFace;
class MGShell;
class MGOpenGLView;
class MGIgesOfstream;
class mgSysGL;
class mgVBO;

/** @defgroup GelRelated Gel Related class
 *  MGGel is top abstract class for MGObject, MGGroup, and MGGLAttrib.
 *  Interface to store data in MGGroup.
 *  @{
 */

///MGGel is an abstract class which represents a group element.

///Gel is the abbreviation of group element, is designed to store
///in MGGroup as an element.
///Subclasses of MGGel are:
///(1) MGAttribedGel(whose sub are MGObject, MGGroup), or (2) MGAttrib.
///MGGel provides functions of serialization of objects.
///All the objects of MGGel subclasses can be serialized using
///MGGroup::make_file(), and MGGroup constructor.
class MG_DLL_DECLR MGGel{

public:

///Virtual Destructor
virtual ~MGGel()=default;

///Assignment.

///When the leaf objects of this and gel2 are not equal, this assignment
///does nothing.
virtual MGGel& operator=(const MGGel& gel2){return *this;};

///Comparison.
//gel2 must be the same class as this.
virtual bool equal_test(const MGGel& gel2)const {
	return false;
};
//gel2 must be the same class as this.
virtual std::partial_ordering ordering_test(const MGGel& gel2)const {
	return std::partial_ordering::unordered;
};

std::strong_ordering  typeCompare(const MGGel& gel2)const;

/// Output virtual function.
virtual std::ostream& toString(std::ostream&) const=0;

///Output the content as std::string.
///The output string is the same as std::cout<<MGGel.
std::string string_content()const;

///IGES output function
///(Default function is no operation to output)
///Function's return value is the directory entry id created.
virtual int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const{return 0;};

///Generate copied gel of this gel.
///Returned is a newed object. User must delete the object.
virtual MGGel* clone()const=0;

///Get manifold dimension.
///MGGroup returns right one, MGGroup return 2, and others return -1.
virtual int manifold_dimension() const{return -1;};

/// Return This object's typeID
virtual long identify_type() const = 0;

///Determine if this is one of the input types or not.
///Function's return value is true if this is one of the input types.
bool type_is(const MGAbstractGels& types)const;

///Get the class name.
virtual std::string whoami()const=0;

protected:

///Read all member data.
virtual void ReadMembers(MGIfstream& buf)=0;

///Write all member data
virtual void WriteMembers(MGOfstream& buf)const=0;

private:

///string stream function
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream& ostrm, const MGGel& gel);
#ifdef FALSE__UNICODE
MG_DLL_DECLR friend std::wostream& operator<< (std::wostream& ostrm, const MGGel& gel);
#endif

friend class MGIfstream;
friend class MGOfstream;

};

long abstractGelId(long typeId);

///Construct a null newed MGGel from the type id TID.
MG_DLL_DECLR MGGel* MGNullGel(long TID);

/** @} */ // end of GelRelated group
#endif
