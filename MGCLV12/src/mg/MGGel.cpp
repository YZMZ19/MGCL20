/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Gel.h"
#include "mg/GelFactory.h"
#include "mg/AbstractGels.h"
#include "mg/Attrib.h"
#include "mg/Ifstream.h"
#include "mg/Ofstream.h"
#include "mg/Group.h"
#include "mgGL/VBO.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

// MGGel
// Implementation of MGGel.
//MGGel is an abstract class which represents a group element.
//Gel is the abbreviation of group element.
//Subclasses of MGGel are:
//(1) MGAttribedGel, or (2) MGAttrib.
//MGGel provides functions of serialization of objects.
//All the objects of MGGel subclasses can be serialized using
//MGGroup::make_file(), and MGGroup constructor.

long abstractGelId(long typeId){
	long id = typeId; id &= 0xff000000;
	switch(id){
		case MGOBJECT_TID:	return id;
		case MGGROUP_TID:	return id;
		case MGATTRIB_TID:	return id;
	}
	return -1;
}

//Construct a null newed MGGel from the type id TID.
MGGel* MGNullGel(long TID){
	MGGelFactoryRegistry* reg = MGGelFactoryRegistry::get_instance();
	return reg->create_gel(TID);
}

bool MGGel::operator<(const MGGel& gel2)const{
	return identify_type() < gel2.identify_type();
}

//Determine if this is one of the input types or not.
//Function's return value is true if this is one of the input types.
bool MGGel::type_is(const MGAbstractGels& types)const{
	MGAbstractGels::const_iterator i=types.begin(), ie=types.end();
	for(; i!=ie; i++){
		MGGEL_KIND agel=(*i).first;
		long tid1=identify_type()&agel;
		long tid2=(*i).second&agel;
		if(tid1==tid2)
			return true;
	}
	return false;
}

//Output the content as std::string.
//The output string is the same as std::cout<<MGGel.
std::string MGGel::string_content()const{
	std::ostringstream s;
#ifdef WIN32
	//added by Tomoko.
	//Global Localeをセットすること。
	s.imbue(std::locale::empty());
#endif
	toString(s);
	return s.str();
}

//////////// MGGel ////////////
ostream& operator<< (ostream& ostrm, const MGGel& gel){
	return gel.toString(ostrm);
}

#ifdef FALSE__UNICODE
std::wostream& operator<< (std::wostream& ostrm, const MGGel& gel){
	std::ostringstream buf;
	gel.toString(buf);
	buf.flush();

	std::string tmp = buf.str();
	return ostrm<<CA2T(tmp.c_str());
}
#endif